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
library_dirs = [htslib_dir]
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
    except:
        print("[WARNING] Conda prefix not found, please activate clair3 conda environment first!")

ffibuilder = FFI()
ffibuilder.set_source("libclair3",
    r"""
    #include "kvec.h"
    #include "khash.h"
    #include "levenshtein.h"
    #include "medaka_bamiter.h"
    #include "medaka_common.h"
    #include "medaka_khcounter.h"
    #include "clair3_pileup.h"
    #include "clair3_full_alignment.h"
    """,
    libraries=libraries,
    library_dirs=library_dirs,
    include_dirs=[src_dir, htslib_dir],
    sources=[
        os.path.join(src_dir, x) for x in (
            'levenshtein.c',
            'medaka_bamiter.c',
            'medaka_common.c',
            'medaka_khcounter.c',
            'clair3_pileup.c',
            'clair3_full_alignment.c')],
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
    extra_objects=[os.path.join(htslib_dir, 'libhts.a')]
)

cdef = [
    "typedef struct { ...; } bam_fset;"
    "bam_fset* create_bam_fset(char* fname);"
    "void destroy_bam_fset(bam_fset* fset);"
]
for header in ('clair3_pileup.h', 'clair3_full_alignment.h'):
    with open(os.path.join(src_dir, header), 'r') as fh:
        # remove directives
        lines = ''.join(x for x in fh.readlines() if not x.startswith('#'))
        cdef.append(lines)

ffibuilder.cdef('\n\n'.join(cdef))


if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
    run("cp {}/libclair3*.so {}/libclair3.so".format(file_directory, file_directory), shell=True)

