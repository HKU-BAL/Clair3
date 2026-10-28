import logging
import multiprocessing
import os
import sys
from concurrent.futures import ProcessPoolExecutor


def load_optimized_core():
    path = os.environ.get("CLAIR3_FAST_CORE")
    if not path:
        return
    try:
        import ctypes
        import importlib.machinery
        import importlib.metadata
        import importlib.util
        import types

        if importlib.metadata.version("whatshap") != "2.8":
            return
        package_spec = importlib.util.find_spec("whatshap")
        package = types.ModuleType("whatshap")
        package.__path__ = list(package_spec.submodule_search_locations)
        package.__version__ = "2.8"
        sys.modules["whatshap"] = package
        ctypes.CDLL(path, mode=ctypes.RTLD_GLOBAL)
        spec = importlib.util.spec_from_file_location("whatshap.core", path)
        module = importlib.util.module_from_spec(spec)
        sys.modules["whatshap.core"] = module
        spec.loader.exec_module(module)
    except Exception as error:
        sys.modules.pop("whatshap.core", None)
        sys.modules.pop("whatshap", None)
        logging.getLogger(__name__).warning(
            "Optimized WhatsHap core unavailable; using the installed implementation: %s",
            error,
        )


load_optimized_core()

from whatshap import __version__
from whatshap.__main__ import main
from whatshap.cli import open_readset_reader, phase
from whatshap.core import NumericSampleIds, Pedigree, PedigreeDPTable, Read, ReadSet
from whatshap.graph import ComponentFinder

logger = logging.getLogger(__name__)
_original_dp_table = PedigreeDPTable
_original_input_reader = phase.PhasedInputReader
_sample_ids = {}
_context = None
_tasks = None
_read_input = None
_read_jobs = None


def solve(task_id):
    positions, reads, samples, genotypes, likelihoods = _tasks[task_id]
    numeric_sample_ids = NumericSampleIds()
    pedigree = Pedigree(numeric_sample_ids)
    for name, sample_genotypes, sample_likelihoods in zip(samples, genotypes, likelihoods):
        numeric_sample_ids[name]
        pedigree.add_individual(name, sample_genotypes, sample_likelihoods)
    table = _original_dp_table(
        reads, [10**6] * len(positions), pedigree, True, positions
    )
    superreads, _ = table.get_super_reads()
    serialized = []
    for readset in superreads:
        output = []
        for read in readset:
            variants = [(v.position, v.allele, v.quality) for v in read]
            output.append((read.name, read.sample_id, variants))
        serialized.append(output)
    cost = table.get_optimal_cost()
    logger.info("Solved phase block %d with %d variants, MEC %d", task_id, len(positions), cost)
    return serialized, cost


class ParallelDPTable:
    def __init__(self, reads, recombcost, pedigree, distrust, positions):
        self.table = None
        if not distrust or len(pedigree) != 1 or _workers() == 1:
            self.table = _original_dp_table(reads, recombcost, pedigree, distrust, positions)
            return
        finder = ComponentFinder(positions)
        for read in reads:
            read_positions = [variant.position for variant in read]
            for position in read_positions[1:]:
                finder.merge(read_positions[0], position)
        groups = {}
        for position in positions:
            groups.setdefault(finder.find(position), []).append(position)
        if len(groups) < 2:
            self.table = _original_dp_table(reads, recombcost, pedigree, distrust, positions)
            return
        self.reads = reads
        self.cost = 0
        component_reads = {component: ReadSet() for component in groups}
        for read in reads:
            component_reads[finder.find(read[0].position)].add(read)
        tasks = []
        indices = {position: index for index, position in enumerate(positions)}
        sample_name = _sample_ids[id(pedigree)].inverse_mapping()[0]
        for block in groups.values():
            readset = component_reads[finder.find(block[0])]
            genotypes = [pedigree.genotype(sample_name, indices[p]) for p in block]
            likelihoods = [
                pedigree.genotype_likelihoods(sample_name, indices[p]) for p in block
            ]
            tasks.append((block, readset, [sample_name], [genotypes], [likelihoods]))
        tasks.sort(key=lambda task: len(task[0]), reverse=True)
        global _tasks
        _tasks = tasks
        workers = min(len(tasks), _workers())
        logger.info("Solving %d independent phase blocks with %d workers", len(tasks), workers)
        merged = {}
        with ProcessPoolExecutor(max_workers=workers, mp_context=_context) as pool:
            results = pool.map(solve, range(len(tasks)), chunksize=1)
            for task_id, (serialized, cost) in enumerate(results):
                _tasks[task_id] = None
                self.cost += cost
                for readset in serialized:
                    for name, sample, variants in readset:
                        merged.setdefault((sample, name), []).extend(variants)
        _tasks = None
        combined = ReadSet()
        for (sample, name), variants in merged.items():
            variants.sort()
            read = Read(name, -1, -1, sample)
            for position, allele, quality in variants:
                read.add_variant(position, allele, quality)
            combined.add(read)
        self.superreads = [combined]

    def get_super_reads(self):
        if self.table is not None:
            return self.table.get_super_reads()
        return self.superreads, [0] * len(self.reads.get_positions())

    def get_optimal_cost(self):
        if self.table is not None:
            return self.table.get_optimal_cost()
        return self.cost

    def get_optimal_partitioning(self):
        if self.table is not None:
            return self.table.get_optimal_partitioning()
        raise RuntimeError("parallel phasing does not support --output-read-list")


def _workers():
    return max(1, int(os.environ.get("CLAIR3_PHASE_THREADS", "1")))


def read_alleles(job_id):
    source = _read_input._readset_reader
    settings = {
        name[1:]: getattr(source, name)
        for name in (
            "_mapq_threshold", "_overhang", "_use_affine", "_gap_start", "_gap_extend",
            "_default_mismatch", "_duplicates", "_use_kmerald", "_kmeralign_costs_path",
            "_kmer_size", "_kmerald_gappenalty", "_kmerald_window", "_use_supplementary",
            "_supplementary_distance_threshold", "_allow_supplementary_only_read_groups",
        )
    }
    settings["affine"] = settings.pop("use_affine")
    reader = open_readset_reader(
        source._paths, _read_input._fasta.filename, source._numeric_sample_ids, **settings
    )
    chromosome, variants, sample, reference = _read_input._parallel_read_input
    output, total, covered = [], 0, set()
    logging.getLogger("whatshap.variants").setLevel(logging.WARNING)
    with reader:
        for start, end, first, last in _read_jobs[job_id]:
            reads = reader.read(
                chromosome, variants[first:last], sample, reference, [(start, end)], None
            )
            total += len(reads)
            covered.update(reads.get_positions())
            for read in reads:
                if len(read) >= 2:
                    read.sort()
                    output.append(read.__getstate__())
    logger.info("Read BAM job %d: %d reads covering %d variants", job_id, total, len(covered))
    return output, total, covered


def parallel_read(self, chromosome, variants, sample, **kwargs):
    if (
        kwargs.get("regions") is not None
        or kwargs.get("restricted_genotypes") is not None
        or not kwargs.get("read_vcf", True)
        or self._vcfs
        or not self._ignore_read_groups
        or len(self._bam_paths) != 1
        or self._readset_reader._use_supplementary
    ):
        return self._parallel_serial_read(chromosome, variants, sample, **kwargs)
    global _read_input, _read_jobs
    _read_input = self
    self._parallel_read_input = (
        chromosome, variants, None, self._fasta[chromosome][:]
    )
    with open(os.environ["CLAIR3_PHASE_REGIONS"]) as stream:
        rows = [tuple(map(int, line.split("\t"))) for line in stream if line.strip()]
    workers = _workers()
    target = (len(variants) + workers - 1) // workers
    jobs, current, count = [], [], 0
    for row in rows:
        if current and count + row[3] - row[2] > target and len(jobs) < workers - 1:
            jobs.append(current)
            current, count = [], 0
        current.append(row)
        count += row[3] - row[2]
    if current:
        jobs.append(current)
    if not jobs:
        return ReadSet(), set()
    _read_jobs = jobs
    logger.info("Reading %d independent BAM regions with %d workers", len(rows), len(jobs))
    readset = ReadSet()
    total, covered = 0, set()
    with ProcessPoolExecutor(max_workers=len(jobs), mp_context=_context) as pool:
        for reads, count, positions in pool.map(read_alleles, range(len(jobs))):
            total += count
            covered.update(positions)
            for state in reads:
                read = Read()
                read.__setstate__(state)
                readset.add(read)
    logger.info("Found %d reads covering %d variants", total, len(covered))
    readset.sort()
    _read_input = None
    _read_jobs = None
    return readset, set()


def make_parallel_input(*args, **kwargs):
    reader = _original_input_reader(*args, **kwargs)
    reader._parallel_serial_read = reader.read
    reader.read = parallel_read.__get__(reader)
    return reader


def make_pedigree(numeric_sample_ids):
    pedigree = Pedigree(numeric_sample_ids)
    _sample_ids[id(pedigree)] = numeric_sample_ids
    return pedigree


def run():
    global _context
    if (
        sys.argv[1:2] == ["phase"]
        and tuple(int(part) for part in __version__.split(".")[:2]) >= (2, 8)
        and "fork" in multiprocessing.get_all_start_methods()
        and "--output-read-list" not in sys.argv
        and "--recombination-list" not in sys.argv
        and "--ped" not in sys.argv
    ):
        _context = multiprocessing.get_context("fork")
        phase.Pedigree = make_pedigree
        phase.PedigreeDPTable = ParallelDPTable
        if os.environ.get("CLAIR3_PHASE_REGIONS"):
            phase.PhasedInputReader = make_parallel_input
    main()


if __name__ == "__main__":
    run()
