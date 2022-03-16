import itertools
import os
import platform
from subprocess import run
from cffi import FFI


samver = "1.10"
longphase_version = "1.0"
file_directory = os.path.dirname(os.path.realpath(__file__))
def compile_samtools_package():
    # just a simple way to compile samtools htslib
    if not os.path.exists(os.path.join(file_directory, 'libhts.a')):
        samtools_source = "samtools-{}.tar.bz2 https://github.com/samtools/samtools/releases/download/{}/samtools-{}.tar.bz2".format(samver, samver, samver)
        run("curl -L -o {}".format(samtools_source), shell=True)
        run("tar -xjf samtools-{}.tar.bz2".format(samver), shell=True)
        run("rm samtools-{}.tar.bz2".format(samver), shell=True)
        run("cd samtools-{} && autoheader && autoconf -Wno-syntax && CFLAGS='-fpic -O3' ./configure && make".format(samver), shell=True)
        run("cp samtools-{}/htslib-{}/libhts.a {}".format(samver, samver, file_directory), shell=True)


def compile_longphase_package():
    if not os.path.exists(os.path.join(file_directory, 'longphase')):
        longphase_source = "https://github.com/twolinin/longphase/archive/refs/tags/v{}.tar.gz".format(longphase_version)
        run("wget {}".format(longphase_source), shell=True)
        run("tar -zxvf v{}.tar.gz".format(longphase_version), shell=True)
        run("rm v{}.tar.gz".format(longphase_version), shell=True)
        run("cd longphase-{} && autoreconf -i && ./configure && make -j4".format(longphase_version), shell=True)
        run("mv longphase-{}/longphase {}".format(longphase_version, file_directory), shell=True)
        run("rm -r longphase-{}".format(longphase_version), shell=True)

def clean_samtools_package():
    # after ffi building, clean the samtools htslib source
    if os.path.exists(os.path.join(file_directory, 'libhts.a')):
        run("rm -r samtools-{}".format(samver), shell=True)

htslib_dir=os.path.join(file_directory, 'samtools-{}'.format(samver), 'htslib-{}'.format(samver))

libraries=['m', 'z', 'lzma', 'bz2', 'pthread', 'curl', 'crypto']
library_dirs=[htslib_dir]
src_dir=os.path.join(file_directory, 'src')

extra_compile_args = ['-std=c99', '-O3']
if platform.machine() in {"aarch64", "arm64"}:
    if platform.system() == "Darwin":
        pass
    else:
        extra_compile_args.append("-march=armv8-a+simd")
else:
    extra_compile_args.append("-mtune=haswell")

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
    extra_objects=['libhts.a']
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
    compile_samtools_package()
    compile_longphase_package()
    ffibuilder.compile(verbose=True)
    clean_samtools_package()
