OS := $(shell uname)
ARCH := $(shell arch)

PYTHON ?= python3

all : libhts.a longphase libclair3.so models
clean : clean_htslib clean_longphase clean_libclair3

SAMVER	=	1.15.1
PREFIX	?=	${CONDA_PREFIX}
LDFLAGS	=	-L ${PREFIX}/lib
CFLAGS	= -fpic -std=c99 -O3 -I ${PREFIX}/include -L ${PREFIX}/lib
CPPFLAGS	=	-std=c++11 -Wall -O3 -I ${PREFIX}/include -L ${PREFIX}/lib -Wl,-rpath=${PREFIX}/lib

ifeq ($(OS),Darwin)
    # Mac settings
    LPVER := 1.5
    CC_PATH ?= clang
    GCC ?= $(CC_PATH)
else
    # Linux settings
    LPVER := 1.7.3
    GCC ?= gcc
    GXX ?= g++
endif


samtools-$(SAMVER)/Makefile:
	curl -L -o samtools-${SAMVER}.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMVER}/samtools-${SAMVER}.tar.bz2
	tar -xjf samtools-${SAMVER}.tar.bz2
	rm samtools-${SAMVER}.tar.bz2

libhts.a: samtools-$(SAMVER)/Makefile
	# this is required only to add in -fpic so we can build python module
	@echo "\x1b[1;33mMaking $(@F)\x1b[0m"
	cd samtools-${SAMVER}/htslib-${SAMVER}; \
	if [ "$(OS)" = "Darwin" ]; then \
		autoheader && autoconf -Wno-syntax && ./configure; \
	else \
		CFLAGS="${CFLAGS}" LDFLAGS="${LDFLAGS}" ./configure; \
	fi; \
	make CFLAGS="${CFLAGS}" LDFLAGS="${LDFLAGS}"
	cp samtools-${SAMVER}/htslib-${SAMVER}/$@ $@


longphase:
	if [ "$(OS)" = "Darwin" ]; then \
		curl -L -o v${LPVER}.tar.gz https://github.com/twolinin/longphase/archive/refs/tags/v${LPVER}.tar.gz; \
		tar -zxf v${LPVER}.tar.gz; \
		cd longphase-${LPVER} && export CC=${CC_PATH} && autoreconf -i && ./configure && make -j; \
		cd .. && rm v${LPVER}.tar.gz; \
		cp longphase-${LPVER}/longphase $@; \
	else \
		curl -L -o longphase-${LPVER}.tar.xz https://github.com/twolinin/longphase/releases/download/v${LPVER}/longphase_linux-x64.tar.xz; \
		tar -xJf longphase-${LPVER}.tar.xz; \
		mv longphase_linux-x64 $@; \
		rm longphase-${LPVER}.tar.xz; \
	fi

models:
	if [ "$(OS)" = "Darwin" ]; then \
        curl -L -o clair3_models.tar.gz http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz; \
        mkdir -p ${PREFIX}/bin/models && tar -zxvf clair3_models.tar.gz -C ${PREFIX}/bin/models; \
        rm clair3_models.tar.gz; \
    fi

libclair3.so: samtools-${SAMVER}/htslib-${SAMVER} libhts.a
	${PYTHON} build.py

src/%.o: src/%.c
	$(GCC) -g -Isrc -Isamtools-${SAMVER}/htslib-${SAMVER} -c -pthread -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $^ -o $@

fa_test: src/fa_test.o src/levenshtein.o src/medaka_bamiter.o src/medaka_common.o src/medaka_khcounter.o src/clair3_pileup.o src/clair3_full_alignment.o src/csv.o
	gcc -g -Isrc -L ./ -pthread -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $^ \
		-lhts -llzma -lcurl -lcrypto -lz -lm -lbz2 -o $@


.PHONY: clean_htslib
clean_htslib:
	cd samtools-${SAMVER} && make clean || exit 0
	cd samtools-${SAMVER}/htslib-${SAMVER} && make clean || exit 0
	rm libhts.a

.PHONY: clean_longphase
clean_longphase:
	rm longphase
	if [ "$(OS)" = "Darwin" ]; then \
		rm -rf $(LONGPHASE_DIR); \
	fi

.PHONY: clean_libclair3
clean_libclair3:
	rm libclair3.*
