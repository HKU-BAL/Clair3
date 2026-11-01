OS := $(shell uname)
ARCH := $(shell arch)
CARGO ?= cargo
PREFIX ?= /usr/local
LIBTORCH ?=
LONGPHASE_CC ?= $(CC)
LONGPHASE_CXX ?= $(CXX)
PYTHON ?= python3

.PHONY: all clair3 longphase python_tools models install clean
all: clair3 longphase python_tools

python_tools:
	$(PYTHON) -m pip install --no-cache-dir whatshap==2.8
	@if [ "$(OS)" = "Darwin" ] && [ -z "$(LIBTORCH)" ]; then \
		$(PYTHON) -m pip install --no-cache-dir torch==2.13.0; \
	fi

clair3: python_tools
	@if [ -n "$(LIBTORCH)" ]; then \
		LIBTORCH="$(LIBTORCH)" $(CARGO) build --release --locked; \
	elif [ "$(OS)" = "Darwin" ]; then \
		LIBTORCH_USE_PYTORCH=1 PYTHON="$(PYTHON)" $(CARGO) build --release --locked; \
	else \
		$(CARGO) build --release --locked; \
	fi
	mkdir -p $(PREFIX)/bin $(PREFIX)/lib/clair3
	cp target/release/clair3 $(PREFIX)/bin/clair3
	@if [ -n "$(LIBTORCH)" ]; then \
		src="$(LIBTORCH)/lib"; \
	elif [ "$(OS)" = "Darwin" ]; then \
		src=$$($(PYTHON) -c 'from torch.utils import cpp_extension; print(cpp_extension.library_paths()[0])'); \
	else \
		src=$$(find target/release/build -path '*/out/libtorch/libtorch/lib' -type d | sort | tail -n 1); \
	fi; \
	if [ -z "$$src" ]; then \
		echo "LibTorch library directory not found" >&2; \
		exit 1; \
	fi; \
	for file in "$$src"/*; do \
		if [ -d "$$file" ] || { [ ! -f "$$file" ] && [ ! -L "$$file" ]; }; then \
			continue; \
		fi; \
		case "$${file##*/}" in \
			libtorch_python*|libshm_windows*) continue ;; \
		esac; \
		cp -f "$$file" $(PREFIX)/lib/clair3/; \
	done

longphase:
	mkdir -p $(PREFIX)/bin
	@if [ "$(OS)" = "Darwin" ]; then \
		curl -fLsS -o longphase.tar.gz https://github.com/twolinin/longphase/archive/refs/tags/v1.7.3.tar.gz && \
		tar -xzf longphase.tar.gz && \
		cd longphase-1.7.3 && autoreconf -i && \
		CC="$(LONGPHASE_CC)" CXX="$(LONGPHASE_CXX)" ./configure && \
		$(MAKE) -j2 CC="$(LONGPHASE_CC)" CXX="$(LONGPHASE_CXX)" && \
		cp longphase $(PREFIX)/bin/longphase; \
	else \
		curl -fLsS -o longphase.tar.xz https://github.com/twolinin/longphase/releases/download/v1.7.3/longphase_linux-x64.tar.xz && \
		tar -xJf longphase.tar.xz && \
		cp longphase_linux-x64 $(PREFIX)/bin/longphase; \
	fi
	chmod +x $(PREFIX)/bin/longphase

models:
	mkdir -p $(PREFIX)/bin/models /tmp/clair3-models
	wget -r -np -nH --cut-dirs=2 -R "index.html*" -P /tmp/clair3-models https://www.bio8.cs.hku.hk/clair3/clair3_models_pytorch/
	cp -a /tmp/clair3-models/. $(PREFIX)/bin/models/
	rm -rf /tmp/clair3-models

install: all
	@echo "Clair3 installed in $(PREFIX)/bin"

clean:
	$(CARGO) clean --release
