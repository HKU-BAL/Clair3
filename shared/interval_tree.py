import shlex
import sys
from shared.intervaltree.intervaltree import IntervalTree

from shared.utils import subprocess_popen


def bed_tree_from(bed_file_path, expand_region=None, contig_name=None, bed_ctg_start=None, bed_ctg_end=None,
                  return_bed_region=False, padding=None):
    """
    0-based interval tree [start, end)
    """

    tree = {}
    if bed_file_path is None:
        if return_bed_region:
            return tree, None, None
        return tree

    bed_start, bed_end = float('inf'), 0
    unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (bed_file_path)))
    for row_id, row in enumerate(unzip_process.stdout):
        if row[0] == '#':
            continue
        columns = row.strip().split()

        ctg_name = columns[0]
        if contig_name != None and ctg_name != contig_name:
            continue
        if ctg_name not in tree:
            tree[ctg_name] = IntervalTree()

        ctg_start, ctg_end = int(columns[1]), int(columns[2])

        if ctg_end < ctg_start or ctg_start < 0 or ctg_end < 0:
            sys.exit("[ERROR] Invalid bed input in {}-th row {} {} {}".format(row_id+1, ctg_name, ctg_start, ctg_end))

        if bed_ctg_start and bed_ctg_end:
            if ctg_end < bed_ctg_start or ctg_start > bed_ctg_end:
                continue
        if padding:
            ctg_start += padding
            ctg_end -= padding
        bed_start = min(ctg_start, bed_start)
        bed_end = max(ctg_end, bed_end)
        if ctg_start == ctg_end:
            ctg_end += 1

        tree[ctg_name].addi(ctg_start, ctg_end)

    unzip_process.stdout.close()
    unzip_process.wait()
    if return_bed_region:
        return tree, bed_start, bed_end
    return tree


def is_region_in(tree, contig_name, region_start=None, region_end=None):
    if not tree or (contig_name is None) or (contig_name not in tree):
        return False

    interval_tree = tree[contig_name]
    return len(
        interval_tree.at(region_start)
        if region_end is None else
        interval_tree.overlap(begin=region_start, end=region_end)
    ) > 0
