import itertools
import os
import platform
from subprocess import run
from cffi import FFI

samver = "1.15.1"
file_directory = os.path.dirname(os.path.realpath(__file__))
htslib_dir = os.path.join(file_directory, 'samtools-{}'.format(samver), 'htslib-{}'.format(samver))

libraries = ['m', 'z', 'lzma', 'bz2', 'pthread', 'curl', 'crypto']
extra_link_args = []
library_dirs = []
if os.path.isdir(htslib_dir):
    library_dirs.append(htslib_dir)
src_dir = os.path.join(file_directory, 'src')

extra_compile_args = ['-std=c99', '-O3']
if platform.machine() in {"aarch64", "arm64"}:
    if platform.system() == "Darwin":
        pass
    else:
        extra_compile_args.append("-march=armv8-a+simd")
else:
    extra_compile_args.append("-mtune=haswell")
    libraries.append('deflate')
    try:
        conda_path = os.environ['CONDA_PREFIX']
        extra_link_args = ['-Wl,-rpath={}/lib'.format(conda_path)]
        library_dirs.append(os.path.join(conda_path, 'lib'))
    except KeyError:
        print("[WARNING] Conda prefix not found, please activate clair3 conda environment first!\n")

htslib_static = os.path.join(htslib_dir, 'libhts.a')
extra_objects = []
if os.path.exists(htslib_static):
    extra_objects = [htslib_static]
else:
    libraries.append('hts')

ffibuilder = FFI()
include_dirs = [src_dir]
if os.path.isdir(htslib_dir):
    include_dirs.append(htslib_dir)

ffibuilder.set_source("libclair3",
    r"""
    #include "kvec.h"
    #include "khash.h"
    #include "levenshtein.h"
    #include "medaka_bamiter.h"
    #include "medaka_common.h"
    #include "medaka_khcounter.h"
    #include "clair3_pileup.h"
    #include "clair3_full_alignment_dwell.h"
    """,
    libraries=libraries,
    library_dirs=library_dirs,
    include_dirs=include_dirs,
    sources=[
        os.path.join(src_dir, x) for x in (
            'levenshtein.c',
            'medaka_bamiter.c',
            'medaka_common.c',
            'medaka_khcounter.c',
            'clair3_pileup.c',
            'clair3_full_alignment_dwell.c')],
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
    extra_objects=extra_objects
)

cdef = [
    "typedef struct { ...; } bam_fset;"
    "bam_fset* create_bam_fset(char* fname, char* fasta_path);"
    "void destroy_bam_fset(bam_fset* fset);"
]
for header in ('clair3_pileup.h', 'clair3_full_alignment_dwell.h'):
    with open(os.path.join(src_dir, header), 'r') as fh:
        # remove directives
        lines = ''.join(x for x in fh.readlines() if not x.startswith('#'))
        cdef.append(lines)

ffibuilder.cdef('\n\n'.join(cdef))


if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
    #check if the libclair3.so file is already in the file_directory
    if os.path.exists(os.path.join(file_directory, "libclair3.so")):
        os.remove(os.path.join(file_directory, "libclair3.so"))
    run("cp {}/libclair3*.so {}/libclair3.so".format(file_directory, file_directory), shell=True)

