"""Pangenome haplotype query and tensor construction for Clair3.

Loads a pangenome wave VCF (e.g., HPRC v2.0) and constructs additional
haplotype rows to append to the full-alignment tensor. Each row represents
one pangenome haplotype's allele at the candidate window positions.

Haplotype selection strategy: for each candidate window, prioritize haplotypes
that carry ALT alleles within the window (most informative), then fill remaining
slots with REF-only haplotypes for diversity context.
"""

import logging
import numpy as np
from collections import defaultdict

logging.basicConfig(format='%(message)s', level=logging.INFO)

# Same base encoding as Clair3 FA tensors (from clair3_full_alignment_dwell.h)
BASE_ENCODE = {'A': 100, 'C': 25, 'G': 75, 'T': -50, 'N': 0}
PANGENOME_MARKER = 25  # Channel 2 value to distinguish pangenome rows from reads


class PangenomeDB:
    """Load and query a pangenome wave VCF for haplotype information.

    For each candidate position's 33bp window, selects the most informative
    haplotypes: those carrying ALT alleles in the window are prioritized,
    then REF-only haplotypes fill remaining slots.
    """

    def __init__(self, vcf_path, ref_fasta_path=None, num_haplotypes=32):
        """
        Args:
            vcf_path: Path to tabix-indexed pangenome VCF
            ref_fasta_path: Path to reference FASTA (for base lookups)
            num_haplotypes: Number of haplotype rows to include in tensor
        """
        import pysam
        self.vcf = pysam.VariantFile(vcf_path)
        self.ref_fasta = pysam.FastaFile(ref_fasta_path) if ref_fasta_path else None
        self.num_haplotypes = num_haplotypes
        self.sample_names = list(self.vcf.header.samples)
        self.num_samples = len(self.sample_names)
        # Total available haplotypes = samples * 2
        self.total_haplotypes = self.num_samples * 2
        logging.info("[INFO] PangenomeDB: loaded %d samples (%d haplotypes), will select %d per site",
                     self.num_samples, self.total_haplotypes, num_haplotypes)

        # Cache for variant lookups
        self._variant_cache = {}
        self._cache_chrom = None
        self._cache_start = None
        self._cache_end = None

    def _ensure_cache(self, chrom, start, end):
        """Load variants for a region into cache. Positions are 1-based."""
        cache_margin = 100000
        if (self._cache_chrom == chrom and
            self._cache_start is not None and
            start >= self._cache_start and end <= self._cache_end):
            return

        # Convert 1-based input to 0-based for pysam.fetch()
        cache_start = max(0, start - 1 - cache_margin)
        cache_end = end + cache_margin
        self._variant_cache = {}
        try:
            for rec in self.vcf.fetch(chrom, cache_start, cache_end):
                pos = rec.pos  # 1-based
                if pos not in self._variant_cache:
                    self._variant_cache[pos] = rec
        except Exception:
            pass
        self._cache_chrom = chrom
        self._cache_start = max(1, start - cache_margin)
        self._cache_end = end + cache_margin

    def get_ref_base(self, chrom, pos_1based):
        """Get reference base at 1-based position."""
        if self.ref_fasta:
            try:
                return self.ref_fasta.fetch(chrom, pos_1based - 1, pos_1based).upper()
            except Exception:
                return 'N'
        return 'N'

    def _select_haplotypes_for_window(self, chrom, start_pos, end_pos):
        """Select the most informative haplotypes for a 33bp window.

        Strategy:
        1. Find all variants in the window
        2. For each haplotype (sample_idx * 2 + allele_idx), count how many
           ALT alleles it carries in this window
        3. Rank by ALT count descending
        4. Take top num_haplotypes

        Returns:
            List of (sample_name, allele_idx) tuples, length = num_haplotypes
        """
        # Count ALT alleles per haplotype in this window
        # haplotype_id = (sample_name, allele_idx)  allele_idx in {0, 1}
        alt_counts = defaultdict(int)

        for pos in range(start_pos, end_pos):
            rec = self._variant_cache.get(pos)
            if rec is None:
                continue
            for sample_name in self.sample_names:
                try:
                    sample = rec.samples[sample_name]
                    alleles = sample.alleles
                    if alleles is None:
                        continue
                    for ai in range(min(2, len(alleles))):
                        a = alleles[ai]
                        if a is not None and a != rec.ref:
                            alt_counts[(sample_name, ai)] += 1
                except Exception:
                    continue

        # Build sorted list: ALT-carrying haplotypes first (by count desc),
        # then REF-only haplotypes
        alt_haplotypes = sorted(alt_counts.keys(), key=lambda k: -alt_counts[k])
        selected = []
        selected_set = set()

        # Phase 1: add haplotypes with ALT alleles
        for hap_id in alt_haplotypes:
            if len(selected) >= self.num_haplotypes:
                break
            selected.append(hap_id)
            selected_set.add(hap_id)

        # Phase 2: fill remaining with REF-only haplotypes
        if len(selected) < self.num_haplotypes:
            for sample_name in self.sample_names:
                for ai in range(2):
                    hap_id = (sample_name, ai)
                    if hap_id not in selected_set:
                        selected.append(hap_id)
                        selected_set.add(hap_id)
                        if len(selected) >= self.num_haplotypes:
                            break
                if len(selected) >= self.num_haplotypes:
                    break

        return selected

    def build_haplotype_tensor(self, chrom, center_pos_1based, no_of_positions, channel_size,
                                flanking_base_num, normalize_num=100):
        """Build pangenome haplotype tensor for one candidate position.

        Selects the most informative haplotypes for this specific window,
        prioritizing those with ALT alleles.

        Args:
            chrom: Chromosome name
            center_pos_1based: 1-based center position of the candidate
            no_of_positions: Window width (33)
            channel_size: Number of channels (8 or 9)
            flanking_base_num: Number of flanking bases (16)
            normalize_num: Normalization constant (100)

        Returns:
            tensor: np.ndarray [num_haplotypes, no_of_positions, channel_size] int8
        """
        N_hap = self.num_haplotypes
        tensor = np.zeros((N_hap, no_of_positions, channel_size), dtype=np.int8)

        start_pos = center_pos_1based - flanking_base_num
        end_pos = center_pos_1based + flanking_base_num + 1

        self._ensure_cache(chrom, start_pos, end_pos)

        # Select best haplotypes for this window
        selected_haplotypes = self._select_haplotypes_for_window(chrom, start_pos, end_pos)

        for offset in range(no_of_positions):
            pos = start_pos + offset
            ref_base = self.get_ref_base(chrom, pos)
            ref_val = BASE_ENCODE.get(ref_base, 0)

            variant_rec = self._variant_cache.get(pos)

            for hap_idx, (sample_name, allele_idx) in enumerate(selected_haplotypes):
                if hap_idx >= N_hap:
                    break

                hap_base_val = ref_val
                hap_af = 0
                insert_val = 0

                if variant_rec is not None:
                    try:
                        sample = variant_rec.samples[sample_name]
                        alleles = sample.alleles
                        if alleles and len(alleles) > allele_idx and alleles[allele_idx] is not None:
                            allele = alleles[allele_idx]
                            if len(allele) == 1:
                                hap_base_val = BASE_ENCODE.get(allele.upper(), 0)
                            elif len(allele) > 1:
                                hap_base_val = BASE_ENCODE.get(allele[0].upper(), 0)
                                insert_val = BASE_ENCODE.get(allele[1].upper(), 0)
                            elif len(allele) == 0:
                                hap_base_val = 0

                        af_values = variant_rec.info.get('AF')
                        if af_values:
                            hap_af = int(min(float(af_values[0]), 1.0) * normalize_num)
                    except Exception:
                        pass

                tensor[hap_idx, offset, 0] = ref_val
                tensor[hap_idx, offset, 1] = hap_base_val
                tensor[hap_idx, offset, 2] = PANGENOME_MARKER
                tensor[hap_idx, offset, 3] = 100
                tensor[hap_idx, offset, 4] = 100
                tensor[hap_idx, offset, 5] = hap_af
                tensor[hap_idx, offset, 6] = insert_val
                tensor[hap_idx, offset, 7] = 0

        return tensor

    def build_batch_haplotype_tensor(self, chrom, candidate_positions_1based,
                                      no_of_positions, channel_size,
                                      flanking_base_num, normalize_num=100):
        """Build pangenome tensors for a batch of candidate positions.

        Args:
            chrom: Chromosome name
            candidate_positions_1based: List of 1-based center positions
            no_of_positions: Window width (33)
            channel_size: Number of channels (8 or 9)
            flanking_base_num: Flanking bases (16)

        Returns:
            tensor: np.ndarray [N_candidates, num_haplotypes, no_of_positions, channel_size] int8
        """
        N = len(candidate_positions_1based)
        N_hap = self.num_haplotypes
        batch_tensor = np.zeros((N, N_hap, no_of_positions, channel_size), dtype=np.int8)

        if N > 0:
            min_pos = min(candidate_positions_1based) - flanking_base_num
            max_pos = max(candidate_positions_1based) + flanking_base_num + 1
            self._ensure_cache(chrom, min_pos, max_pos)

        for i, center_pos in enumerate(candidate_positions_1based):
            batch_tensor[i] = self.build_haplotype_tensor(
                chrom, center_pos, no_of_positions, channel_size,
                flanking_base_num, normalize_num
            )

        return batch_tensor

    def close(self):
        """Release resources."""
        if self.vcf:
            self.vcf.close()
        if self.ref_fasta:
            self.ref_fasta.close()
