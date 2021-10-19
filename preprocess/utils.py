from cffi import FFI
import logging
import os
import sys
import re
import subprocess
import shlex
logging.getLogger().setLevel(logging.INFO)
from shared.utils import file_path_from, subprocess_popen

use_mpmath = True
try:
    import mpmath as math
except ImportError:
    import math
    use_mpmath = False

LOG_10 = 2.3025
LOG_2 = 0.3010

# compress intermediate gvcf row using lz4, issue:https://github.com/HKU-BAL/Clair3/issues/48
lz4_path = subprocess.run("which lz4", stdout=subprocess.PIPE, shell=True).stdout.decode().rstrip()
COMPRESS_GVCF = True if lz4_path != "" else False
LZ4_COMPRESS = "lz4 -c"
LZ4_DECOMPRESS = "lz4 -fdc"
GVCF_SUFFIX = ".tmp.gvcf"

class compressReaderWriter(object):
    def __init__(self, input_path=None, output_path=None, compress=False):
        self.input_path = input_path
        self.output_path = output_path
        self.compress = compress
        self.read_proc = None
        self.reader = None

        self.writer = None
        self.write_proc = None
        self.write_fpo = None

    def read_input(self):
        if self.compress:
            self.read_proc = subprocess_popen(shlex.split("{} {}".format(LZ4_DECOMPRESS, self.input_path)), stderr=subprocess.DEVNULL)
            a = subprocess_popen(shlex.split("{} {}".format(LZ4_DECOMPRESS, self.input_path)), stderr=subprocess.DEVNULL)
            streamdata = a.communicate()[0]
            rc = a.returncode
            self.reader = self.read_proc.stdout
        else:
            self.reader = open(self.input_path, 'r')
        return self.reader

    def close_reader(self):
        if self.compress:
            self.read_proc.stdout.close()
            self.read_proc.wait()
        else:
            self.reader.close()

    def write_output(self):
        if self.compress:
            self.write_fpo = open(self.output_path, 'w')
            self.write_proc = subprocess_popen(shlex.split(LZ4_COMPRESS), stdin=subprocess.PIPE, stdout=self.write_fpo, stderr=subprocess.DEVNULL)
            self.writer = self.write_proc.stdin

        else:
            self.writer = open(self.output_path, 'w')
        return self.writer

    def close_writer(self):
        if self.compress:
            self.write_proc.stdin.close()
            self.write_proc.wait()
            self.write_fpo.close()
        else:
            self.writer.close()


class gvcfGenerator(object):

    def __init__(self, ref_path, samtools='samtools'):

        self.reference_file_path = ref_path
        self.samtools = samtools
        pass

    def readCalls(self, callPath, callType='variant', ctgName=None, ctgStart=None, ctgEnd=None, add_header=False, writer=None):

        CR = compressReaderWriter(input_path=callPath, compress=COMPRESS_GVCF)
        reader = CR.read_input()
        need_write_header = True
        header = []
        for line in reader:
            if (line.startswith('#')):
                if add_header and line not in header:
                        header.append(line)

                continue
            if add_header and len(header) and need_write_header:
                print(''.join(header).rstrip(), file=writer)
                need_write_header = False
            if (callType == 'non-variant'):
                cur_non_variant_start = int(line.strip('\n').split('\t')[1])
                cur_non_variant_end = int(re.search(r'.*END=(.*)\tGT.*', line).group(1))
                cur_non_variant_chr = line.strip('\n').split('\t')[0]
                if ((ctgName and cur_non_variant_chr == ctgName) or (not ctgName)):
                    if ((ctgStart and cur_non_variant_start >= ctgStart) or (not ctgStart)):
                        if ((ctgEnd and cur_non_variant_end <= ctgEnd) or (not ctgEnd)):
                            yield line.strip('\n'), cur_non_variant_start, cur_non_variant_end, 'original'
            else:
                # for variant calls, return "pos"
                # DEL and INS should be considered here
                tmp = line.strip('\n').split('\t')
                ref = tmp[3]
                alt = tmp[4]
                n_alt = len(alt.split(','))
                cur_variant_start = int(line.strip('\n').split('\t')[1])
                cur_variant_end = cur_variant_start - 1 + len(ref)
                is_reference_call = (alt == '.') or (ref == alt)
                if not is_reference_call:
                # assuming AD is at the columns [-3], add 0 to AD for gVCF
                    ori_info = tmp[-1].split(':')
                    ori_info[-3] += ',0'
                    tmp[-1] = ':'.join(ori_info)

                    # assumeing PL is at the last column
                    # add <NON_REF> to variant calls
                    tmp[4] = tmp[4] + ',<NON_REF>'
                    if (n_alt == 1):

                        tmp[-1] = tmp[-1] + ',990,990,990'

                    elif (n_alt == 2):
                        tmp[-1] = tmp[-1] + ',990,990,990,990'
                else:
                    # skip reference calls
                    continue
                new_line = '\t'.join(tmp)

                cur_variant_chr = tmp[0]


                if ((ctgName and cur_variant_chr == ctgName) or (not ctgName)):
                    if ((ctgStart and cur_variant_start >= ctgStart) or (not ctgStart)):
                        if ((ctgEnd and cur_variant_end <= ctgEnd) or (not ctgEnd)):
                            yield new_line, cur_variant_start, cur_variant_end

        CR.close_reader()

    def readReferenceBaseAtPos(self, pos):

        cmd = self.samtools + ' faidx ' + self.reference_file_path + ' ' + pos

        reader = os.popen(cmd)
        for line in reader:
            if (line.startswith('>')):
                continue
            else:
                ref_base = line.strip('\n').upper()
                reader.close()
                return ref_base
        
    def _writeRightBlock(self, block_new_start, curNonVarEnd, curNonVarCall, save_writer):

        pos_cmd = str(curNonVarCall.split('\t')[0]) + ':' + str(block_new_start) + '-' + str(block_new_start)
        new_ref = self.readReferenceBaseAtPos(pos_cmd)
        tmp = curNonVarCall.split('\t')
        tmp[1] = str(block_new_start)
        tmp[3] = str(new_ref)
        print('\t'.join(tmp), file=save_writer)

    def _writeLeftBlock(self, end_pos, curNonVarCall, save_writer):

        new_left_block = re.sub("END=[0-9]*\t", "END=" + str(end_pos) + '\t', curNonVarCall)
        print(new_left_block, file=save_writer)
        pass

    def writeNonVarBlock(self, start, end, pos_flag, curNonVarCall, save_writer):

        if (pos_flag == 'left'):
            self._writeLeftBlock(end, curNonVarCall, save_writer)
        elif (pos_flag == 'right'):
            self._writeRightBlock(start, end, curNonVarCall, save_writer)
        else:
            print(curNonVarCall, file=save_writer)
    def mergeCalls(self, variantCallPath, nonVarCallPath, savePath, sampleName, ctgName=None, ctgStart=None,
                   ctgEnd=None):

        '''
        merge calls between variant and non-variant
        '''

        varCallStop = False
        nonVarCallStop = False

        #output writer
        CW = compressReaderWriter(output_path=savePath, compress=COMPRESS_GVCF)
        save_writer = CW.write_output()

        varCallGenerator = self.readCalls(variantCallPath, 'variant', ctgName, ctgStart, ctgEnd)
        nonVarCallGenerator = self.readCalls(nonVarCallPath, 'non-variant', ctgName, ctgStart, ctgEnd, add_header=True, writer=save_writer)
        hasVar = True
        # in case of empty file
        try:
            curVarCall, curVarStart, curVarEnd = next(varCallGenerator)
        except StopIteration:
            varCallStop = True
            hasVar = False
        try:
            curNonVarCall, curNonVarStart, curNonVarEnd, curNonVarPos = next(nonVarCallGenerator)
        except StopIteration:
            nonVarCallStop = True

        while True and (not varCallStop) and (not nonVarCallStop):
            if (curNonVarEnd < curVarStart):

                '''
                |____|   {____}
                 nonVar    Var
                nonVar region is on the left, no overlapped region
                '''
                # print(curNonVarCall,file=save_writer)
                self.writeNonVarBlock(curNonVarStart, curNonVarEnd, curNonVarPos, curNonVarCall, save_writer)
                # move non variant calls to the next
                try:
                    curNonVarCall, curNonVarStart, curNonVarEnd, curNonVarPos = next(nonVarCallGenerator)
                except StopIteration:
                    nonVarCallStop = True
                    break
            elif (curVarEnd < curNonVarStart):
                '''
                {____}     |____|
                 Var         nonVar
                var region is on the left, no overlapped region
                '''
                # print("{_____} |_____|")
                print(curVarCall, file=save_writer)
                try:
                    curVarCall, curVarStart, curVarEnd = next(varCallGenerator)
                except StopIteration:
                    varCallStop = True
                    break

            elif (curVarStart <= curNonVarStart and curVarEnd >= curNonVarStart):
                '''
                {____|____}___|
                or
                {____|_______|____}
                the left point of nonvar block can be included
                var region is on the left, has overlapped region
                '''
                # write the current variant Call
                print(curVarCall, file=save_writer)
                block_new_start = curVarEnd + 1
                try:
                    curVarCall, curVarStart, curVarEnd = next(varCallGenerator)
                except StopIteration:
                    varCallStop = True
                    break

                while (block_new_start > curNonVarEnd):
                    # skip the non-variant block within the current variant block

                    try:

                        curNonVarCall, curNonVarStart, curNonVarEnd, curNonVarPos = next(nonVarCallGenerator)
                    except StopIteration:
                        nonVarCallStop = True
                        break

                if (nonVarCallStop):
                    break

                # check if the start of the current non-variant block
                if ((block_new_start - 1) >= curNonVarStart):
                    # there is overlap between variants and non-variant block
                    # just write the right part of the non-variant block
                    curNonVarStart = block_new_start
                    curNonVarPos = 'right'

            elif (curVarStart > curNonVarStart):
                '''
                |_{__________________}__|
                or 
                |__{______________|____}
                '''
                # var call is within the non-var block
                # split the non-var block
                non_var_block_left_end = curVarStart - 1
                if (non_var_block_left_end >= curNonVarStart):
                    self._writeLeftBlock(non_var_block_left_end, curNonVarCall, save_writer)
                # print out variant call
                print(curVarCall, file=save_writer)
                # take care here, whether write the left right non-variant block,
                # it dependes on the position of the next variant calls
                non_var_block_right_start = curVarEnd + 1

                try:
                    curVarCall, curVarStart, curVarEnd = next(varCallGenerator)
                except StopIteration:
                    varCallStop = True
                    break

                # still has the right left block
                if (non_var_block_right_start <= curNonVarEnd):
                    curNonVarStart = non_var_block_right_start
                    curNonVarPos = 'right'
                else:
                    # get the next non-variant block,skip the non-var block that is within the variant
                    while True:
                        try:
                            curNonVarCall, curNonVarStart, curNonVarEnd, curNonVarPos = next(nonVarCallGenerator)
                        except StopIteration:
                            nonVarCallStop = True
                            break
                        if (non_var_block_right_start <= curNonVarEnd):
                            break

                    if (nonVarCallStop):
                        break

                    curNonVarStart = non_var_block_right_start
                    curNonVarPos = 'right'

            else:
                print("[ERROR] CurVarStart", curVarStart, 'curVarEnd', curVarEnd, 'curNonVarStart', curNonVarStart,
                      'curNonVarEnd', curNonVarEnd)

        # printout the remain content
        if (not varCallStop):
            # print out the left

            print(curVarCall, file=save_writer)
            for curVarCall, curVarStart, curVarEnd in varCallGenerator:
                print(curVarCall, file=save_writer)
        if (not nonVarCallStop):
            if (hasVar and curNonVarEnd > curVarEnd):
                self.writeNonVarBlock(curVarEnd + 1, curNonVarEnd, curNonVarPos, curNonVarCall, save_writer)
            for curNonVarCall, curNonVarStart, curNonVarEnd, curNonVarPos in nonVarCallGenerator:
                print(curNonVarCall, file=save_writer)

        CW.close_writer()


class variantInfoCalculator(object):

    def __init__(self, gvcfWritePath, ref_path, p_err, gq_bin_size, ctgName, bp_resolution=False, sample_name='None', mode='L'):

        # default p_error is 0.001, while it could be set by the users' option
        self.p_error = p_err
        self.LOG_10 = LOG_10
        self.logp = math.log(self.p_error) / self.LOG_10
        self.log1p = math.log1p(-self.p_error) / self.LOG_10
        self.LOG_2 = LOG_2
        # need to check with the clair3 settings
        #self.max_gq = 255
        self.max_gq = 50
        self.variantMath = mathcalculator()
        self.constant_log10_probs = self.variantMath.normalize_log10_prob([-1.0, -1.0, -1.0])
        self.gq_bin_size = gq_bin_size
        self.CW = None
        # set by the users
        if (gvcfWritePath != "PIPE"):
            if (not os.path.exists(gvcfWritePath)):
                os.mkdir(gvcfWritePath)

            self.CW = compressReaderWriter(output_path=os.path.join(gvcfWritePath, sample_name + GVCF_SUFFIX), compress=COMPRESS_GVCF)
            self.vcf_writer = self.CW.write_output()
        else:
            self.vcf_writer = sys.stdout
        self.writePath = gvcfWritePath
        self.sampleName = sample_name.split('.')[0]
        self.bp_resolution = bp_resolution
        self.reference_file_path = ref_path

        if (mode == 'L'):
            # dictionary to store constant log values for speeding up
            self.normalized_prob_pool = {}

            self.current_block = []
            self._print_vcf_header()
            self.cur_gq_bin_index = None
            self.cur_gt = None
            self.cur_min_DP = None
            self.cur_max_DP = None
            self.cur_chr = None
            self.cur_raw_gq = None
        pass
    def write_empty_pileup(self,ctgName,ctgStart,ctgEnd):
        
        non_variant_info = {"validPL": False, "gq": 1, "binned_gq": 1, "pl": [0,0,0],
                            "chr": ctgName, 'pos': max(1,ctgStart), 'ref': 'N',
                            "gt": './.', 'min_dp': 0, 'END': ctgEnd}
        self.write_to_gvcf(non_variant_info)
    def make_gvcf_online(self, variant_summary, push_current=False):

        '''
        
        make gvcf while reading from pileup
        '''

        if (push_current):
            if (len(self.current_block) > 0):
                self.write_to_gvcf_batch(self.current_block, self.cur_min_DP, self.cur_raw_gq)
                self.current_block = []
                self.cur_gq_bin_index = None
                self.cur_gt = None
                self.cur_min_DP = None
                self.cur_max_DP = None
                self.cur_chr = None
                self.cur_raw_gq = None
            return

        cur_item = self.reference_likelihood(variant_summary)
        _gq_bin = cur_item['binned_gq']
        _gt = cur_item["gt"]
        _DP = cur_item["min_dp"]
        _chr = cur_item['chr']
        _raw_gq = cur_item['gq']
        _cur_ref = variant_summary['ref']

        if (self.cur_gq_bin_index == None):
            self.current_block, self.cur_gq_bin_index, self.cur_gt, self.cur_min_DP, self.cur_max_DP, self.cur_chr, self.cur_raw_gq = (
                [cur_item], _gq_bin, _gt, _DP, _DP, _chr, _raw_gq)
            self.cur_ref = _cur_ref
 
        elif (_gq_bin != self.cur_gq_bin_index):
            self.write_to_gvcf_batch(self.current_block, self.cur_min_DP, self.cur_raw_gq)
            self.current_block, self.cur_gq_bin_index, self.cur_gt, self.cur_min_DP, self.cur_max_DP, self.cur_chr, self.cur_raw_gq = (
                [cur_item], _gq_bin, _gt, _DP, _DP, _chr, _raw_gq)
            self.cur_ref = _cur_ref

        elif (_gt != self.cur_gt):
            self.write_to_gvcf_batch(self.current_block, self.cur_min_DP, self.cur_raw_gq)
            self.current_block, self.cur_gq_bin_index, self.cur_gt, self.cur_min_DP, self.cur_max_DP, self.cur_chr, self.cur_raw_gq = (
                [cur_item], _gq_bin, _gt, _DP, _DP, _chr, _raw_gq)
            self.cur_ref = _cur_ref

        elif (_chr != self.cur_chr):
            self.write_to_gvcf_batch(self.current_block, self.cur_min_DP, self.cur_raw_gq)
            self.current_block, self.cur_gq_bin_index, self.cur_gt, self.cur_min_DP, self.cur_max_DP, self.cur_chr, self.cur_raw_gq = (
                [cur_item], _gq_bin, _gt, _DP, _DP, _chr, _raw_gq)
            self.cur_ref = _cur_ref
        elif( (_cur_ref != self.cur_ref) and ((_cur_ref=='N') or (self.cur_ref=='N'))):
            self.write_to_gvcf_batch(self.current_block, self.cur_min_DP, self.cur_raw_gq)
            self.current_block, self.cur_gq_bin_index, self.cur_gt, self.cur_min_DP, self.cur_max_DP, self.cur_chr, self.cur_raw_gq = (
                [cur_item], _gq_bin, _gt, _DP, _DP, _chr, _raw_gq)
            self.cur_ref = _cur_ref
        else:
            '''         
            # do not consider DP 
            if(_DP < self.cur_min_DP):
                self.cur_min_DP = _DP
            if(_raw_gq < self.cur_raw_gq):
                self.cur_raw_gq = _raw_gq
            self.current_block.append(cur_item)
            '''
            
            if (_DP < self.cur_min_DP):
                tmp_cur_min_DP = _DP
                #if (self.cur_max_DP > math.ceil((tmp_cur_min_DP + min(3, tmp_cur_min_DP * 0.3)))):
                if (self.cur_max_DP > math.ceil(tmp_cur_min_DP + tmp_cur_min_DP * 0.3)):
                    self.write_to_gvcf_batch(self.current_block, self.cur_min_DP, self.cur_raw_gq)
                    self.current_block, self.cur_gq_bin_index, self.cur_gt, self.cur_min_DP, self.cur_max_DP, self.cur_chr, self.cur_raw_gq = (
                        [cur_item], _gq_bin, _gt, _DP, _DP, _chr, _raw_gq)
                else:
                    self.cur_min_DP = tmp_cur_min_DP
                    if (_raw_gq < self.cur_raw_gq):
                        self.cur_raw_gq = _raw_gq
                    self.current_block.append(cur_item)
            elif (_DP > self.cur_max_DP):
                #if (_DP <= math.ceil(self.cur_min_DP + min(3, self.cur_min_DP * 0.3))):
                if (_DP <= math.ceil(self.cur_min_DP + self.cur_min_DP * 0.3)):
                    self.cur_max_DP = _DP
                    if (_raw_gq < self.cur_raw_gq):
                        self.cur_raw_gq = _raw_gq
                    self.current_block.append(cur_item)
                else:
                    self.write_to_gvcf_batch(self.current_block, self.cur_min_DP, self.cur_raw_gq)
                    self.current_block, self.cur_gq_bin_index, self.cur_gt, self.cur_min_DP, self.cur_max_DP, self.cur_chr, self.cur_raw_gq = (
                        [cur_item], _gq_bin, _gt, _DP, _DP, _chr, _raw_gq)
            else:
                if (_raw_gq < self.cur_raw_gq):
                    self.cur_raw_gq = _raw_gq
                self.current_block.append(cur_item)
            
    def reference_likelihood(self, variant_summary):

        '''
        for non-variant sites, this function is calculate the GQ,QUAL,PL,etc for the genotype 0/0 or ./.
        '''

        n_ref = variant_summary["n_ref"]
        n_total = variant_summary['n_total']
      
        
        validPL, gq, binned_gq, log10_probs = self._cal_reference_likelihood(n_ref, n_total)
        if (validPL):
            gt = '0/0'
        else:
            gt = './.'
        _tmp_phred_probs = [-10 * x for x in log10_probs]
        min_phred_probs = min(_tmp_phred_probs)

        phred_probs = [int(x - min_phred_probs) for x in _tmp_phred_probs]

        if(variant_summary['ref'] not in ['A','T','C','G']):
            tmp_ref = 'N'
            gq = 1
            binned_gq = 1
            phred_probs = [0,0,0]
        else:
            tmp_ref = variant_summary['ref']
        non_variant_info = {"validPL": validPL, "gq": gq, "binned_gq": binned_gq, "pl": phred_probs,
                            "chr": variant_summary['chr'], 'pos': variant_summary['pos'], 'ref': tmp_ref,
                            "gt": gt, 'min_dp': variant_summary['n_total'], 'END': variant_summary['pos']}

        return non_variant_info
        pass

    def _cal_reference_likelihood(self, n_ref, n_total):

        '''
        calculate the phred genotype likelihood for a single non-variant site.
        n_ref: number of referece bases
        n_total: number of all bases by ignoring Ns
        P(hom_ref) =  (1-prr)^n_ref*prr^(n_total-n_ref)
        P(Het_alt) = (1/2)^n_total
        P(hom_alt) = prr^n_ref*(1-prr)^(n_total-n_ref)
        return flag of validPL, raw GQ, binned GQ, PLs
        '''

        validPL = True
        if (n_total == 0):
            # when the coverage is 0
            log10_probs = self.constant_log10_probs
            
            pass
        else:

            n_alts = n_total - n_ref
            
            log10_p_ref = n_ref * self.log1p + n_alts * self.logp
            
            log10_p_het = -n_total * self.LOG_2
            log10_p_hom_alt = n_ref * self.logp + n_alts * self.log1p
            
            # normalization
             
            log10_probs = self.variantMath.normalize_log10_prob([log10_p_ref, log10_p_het, log10_p_hom_alt])
            
        
        
        gq = self.variantMath.log10p_to_phred(log10_probs[0])
        
        gq = int(min(int(gq), self.max_gq))
        if (gq >= 1):
            binned_index = (gq - 1) // self.gq_bin_size
            binned_gq = binned_index * self.gq_bin_size + 1
        else:
            binned_gq = 0

        
        validPL = log10_probs[0] == max(log10_probs)
        return validPL, gq, binned_gq, log10_probs

    def _print_vcf_header(self):

        from textwrap import dedent
        print(dedent("""\
            ##fileformat=VCFv4.2
            ##FILTER=<ID=PASS,Description="All filters passed">
            ##FILTER=<ID=LowQual,Description="Low quality variant">
            ##FILTER=<ID=RefCall,Description="Reference call">
            ##INFO=<ID=P,Number=0,Type=Flag,Description="Result from pileup calling">
            ##INFO=<ID=F,Number=0,Type=Flag,Description="Result from full-alignment calling">
            ##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
            ##INFO=<ID=END,Number=1,Type=Integer,Description="End position (for use with symbolic alleles)">
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
            ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">
            ##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
            ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">
            ##FORMAT=<ID=AF,Number=1,Type=Float,Description="Estimated allele frequency in the range of [0,1]">"""), file
        =self.vcf_writer)
        if self.reference_file_path is not None:
            reference_index_file_path = file_path_from(self.reference_file_path, suffix=".fai", exit_on_not_found=True, sep='.')
            with open(reference_index_file_path, "r") as fai_fp:
                for row in fai_fp:
                    columns = row.strip().split("\t")
                    contig_name, contig_size = columns[0], columns[1]
                    print("##contig=<ID=%s,length=%s>" % (contig_name, contig_size), file=self.vcf_writer)

        print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s' % (self.sampleName), file=self.vcf_writer)

        pass

   
    def write_to_gvcf_batch(self, block, block_min_dp, block_min_raw_gq):
         
        if ((self.bp_resolution or block[0]['gt'] == "./.") and block[0]['ref']!='N'):
            # write it to VCF
            for item in block:
                self.write_to_gvcf(item)
        
        else:
            start_pos = block[0]['pos']
            end_pos = block[-1]['pos']
            first_PL = block[0]['pl']
            first_gq = block[0]['gq']
            first_binned_gq = block[0]['binned_gq']
            first_gt = block[0]['gt']
            first_ref = block[0]['ref']
            first_chr = block[0]['chr']
            if(first_ref=='N'):
                # special case for reference N
                non_variant_info = {"gq": 1, "binned_gq": 1, "pl": [0,0,0], "chr": first_chr,
                                'pos': start_pos, 'ref': first_ref, "gt": "./.", 'min_dp': block_min_dp,
                                'END': end_pos}
 
            # write binned_gq
            # non_variant_info = { "gq":first_gq, "binned_gq":first_binned_gq, "pl":first_PL,"chr":first_chr,'pos':start_pos,'ref':first_ref,"gt":first_gt,'min_dp':block_min_dp,'END':end_pos}
            # write min raw gq
            else:
                non_variant_info = {"gq": first_gq, "binned_gq": block_min_raw_gq, "pl": first_PL, "chr": first_chr,
                                'pos': start_pos, 'ref': first_ref, "gt": first_gt, 'min_dp': block_min_dp,
                                'END': end_pos}
            self.write_to_gvcf(non_variant_info)

    def write_to_gvcf(self, variant_info):

        '''
        write a temporary file gvcf. This file is needed to be merged with model variant calls.
        '''

        _tmpLine = str(variant_info["chr"]) + '\t' + str(variant_info["pos"]) + "\t.\t" + variant_info[
            'ref'] + '\t<NON_REF>\t0\t.\tEND=' + str(variant_info['END']) + '\tGT:GQ:MIN_DP:PL\t' + variant_info[
                       'gt'] + ':' + str(variant_info['binned_gq']) + ':' + str(variant_info['min_dp']) + ':' + str(
            variant_info['pl'][0]) + ',' + str(variant_info['pl'][1]) + ',' + str(variant_info['pl'][2])
        print(_tmpLine, file=self.vcf_writer)


    def close_vcf_writer(self):
        self.CW.close_writer()

class mathcalculator(object):


    def __init__(self,speedUp=True):

        self.LOG_10 = LOG_10
        self.maxPhredScore = 255
        self.speedUp = speedUp
        if(speedUp):
            try:
                self._creatCFFIFunc()
            except:
                self.speedUp = False
                pass
            

    
    def _creatCFFIFunc(self):
        self.ffi = FFI()
        self.ffi.cdef("""
                double log10p_to_phred(double log10p);
                double log10sumexp(double log10_array[],int n_array);
                double getMyMaxItem(double list[],int n_list);
                """)
        self.lib = self.ffi.verify("""
                        #include<math.h>
                        #include<stdio.h>
                        double LOG_10 = 2.3025;
                        double log10p_to_phred(double log10p){
                            double ptrue;
                            ptrue = pow(10,log10p);
                            if(ptrue==1){
                                
                                return 50;
                            }
                            return -10*(log(1.0-ptrue)/LOG_10);
                        }
                       double f_log10(double myInput){
                            double res;
                            res=log(myInput)/LOG_10;
                            return res;
                       }
                       double getMyMaxItem(double list[],int n_list){
                           double curMax;
                           int i;
                           curMax = list[0];
                           for(i=1;i<=n_list;i++){
                               if(list[i]>curMax){
                                   curMax= list[i];
                               }
                            }
                           return curMax;
                       } 
                       double log10sumexp(double log10_array[],int n_array){
                           double m, mySum,tmp;
                           int i;
        
                           m = getMyMaxItem(log10_array,n_array);
                           mySum = 0.0;
                           for(i=0;i<n_array;i++){
                               tmp = pow(10,log10_array[i]-m);
                               mySum += tmp;
                           }
                           return m+log(mySum)/LOG_10;
                       }
                       """
                       )

          
    def log10p_to_phred(self,log10p):

        '''
        log10(p) to -10log10(1-p)
        '''
        
        if(self.speedUp):
            return round(self.lib.log10p_to_phred(log10p),6)

        ptrue = math.power(10, log10p) if use_mpmath else math.pow(10, log10p)
        if ptrue == 1:
            return 50
        return round(-10*(math.log(1-ptrue)/self.LOG_10),6)
        
    def log10sumexp(self,log10_array):

        '''
        return value approxiameting to log10(sum(10^log10_probs))
        '''

        if(self.speedUp):
            
            n_array = self.ffi.cast("int", len(log10_array))
            log10_array_c = self.ffi.new("double []",log10_array)
            return self.lib.log10sumexp(log10_array_c,n_array) 
            
        m = max(log10_array)
        return m + math.log10(sum(pow(10.0, x - m) for x in log10_array))
   

    def normalize_log10_prob(self,log10_prob):

        '''
        normalize the log10 probs that make sum of the probabilities close to 1
    
        '''
        # keep 6 digits
        lse = round(self.log10sumexp(log10_prob),6)
         
        normalized_log10_probs = []
         
        
        for x in log10_prob:
            normalized_log10_probs.append(min(x-lse,0))
   
        return normalized_log10_probs