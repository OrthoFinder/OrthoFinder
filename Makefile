.ONESHELL:
SHELL := /bin/bash
.SHELLFLAGS := -e -o pipefail -c

.DEFAULT_GOAL := install

ifndef MAKE_VERSION
$(error This Makefile requires GNU Make. Please use gmake instead of make.)
endif

# ----------------------------------------------------------------- usage ----
# make install_tools                 the programs OrthoFinder uses by default (DIAMOND, MCL, FAMSA,
#                                    FastTree, ASTRAL-Pro, FastME, MAFFT), in ~/local/bin
# make install_all_tools             also the optional ones (IQ-TREE, RAxML, RAxML-NG, MUSCLE,
#                                    MMseqs2, BLAST+)
# make install_<program>             one of them, e.g. install_diamond
# make update_package_bin            replace the programs in src/orthofinder/bin, the ones
#                                    OrthoFinder runs (same as: make install_tools PACKAGE_BIN=true)
# make show_versions                 the versions of the programs in the install directory
# make conda_install                 OrthoFinder and the programs from bioconda, in a conda environment
#
# Where the programs are installed:
#   (default)          ~/local/bin, added to PATH in the shell's startup file; a program already
#                      on PATH is skipped (FORCE=true installs it again)
#   SYSTEM_WIDE=true   /usr/local/bin (with sudo)
#   BIN_DIR=<dir>      any directory; a program already there is skipped (FORCE=true replaces it)
#   PACKAGE_BIN=true   src/orthofinder/bin: the programs there are replaced (git keeps the old ones)
# Other options: QUIET=false shows the output of the builds; JOBS=<n> parallel build jobs;
# PORTABLE=true builds binaries that run on other machines (no CPU-specific code, static where
# possible; the default with PACKAGE_BIN=true); <PROGRAM>_VERSION=<version> another version.

QUIET ?= true
FORCE ?= false
JOBS ?= 4
CONDA_PYTHON_VERSION ?= 3.10
PYTHON_VERSION ?= python3
ORTHOFINDER_ENV_DEFAULT := of3_env
ENV_NAME ?= $(ORTHOFINDER_ENV_DEFAULT)

# ----------------------------------------------- versions of the programs ----
# The latest stable releases (checked 2026-10-11).
# OrthoFinder: empty for the latest release on Bioconda (ORTHOFINDER_VERSION=3.2 for another)
ORTHOFINDER_VERSION ?=
DIAMOND_VERSION ?= 2.2.8
MCL_VERSION ?= 22-282
CIMFOMFA_VERSION ?= 22-273
FAMSA_VERSION ?= 2.5.2
FASTME_VERSION ?= 2.1.6.3
FASTTREE_VERSION ?= 2.2.0
ASTER_VERSION ?= 1.25
MAFFT_VERSION ?= 7.526
RAXML_VERSION ?= 8.2.13
RAXMLNG_VERSION ?= 2.0.3
MUSCLE_VERSION ?= 5.3
IQTREE_VERSION ?= 3.1.4
BLAST_VERSION ?= 2.17.0
MMSEQS_VERSION ?= 18-8cc5c

# bioconda's versions, where its latest differs or it writes the version differently
CONDA_MCL_VERSION ?= 22.282
CONDA_FAMSA_VERSION ?= 2.4.1
CONDA_MAFFT_VERSION ?= 7.525
CONDA_MMSEQS_VERSION ?= 18.8cc5c

# FastTree is built single precision (as the FastTree that OrthoFinder has been
# shipping): FastTree 2.2 is double precision by default, which is ~1.5x slower
FASTTREE_CFLAGS ?= -O3 -finline-functions -funroll-loops -DUSE_SINGLE

# The BLAST+ programs installed (the release has ~30)
BLAST_PROGRAMS ?= blastp blastn makeblastdb

# RAxML: the binary installed (the names of the bioconda package), the makefile that
# builds it, and the name that makefile gives it (the AVX ones build "...-AVX")
RAXML_BINARY_VERSION ?= raxmlHPC
RAXML_MAKEFILE_MAP = \
	raxmlHPC:Makefile.gcc \
	raxmlHPC-SSE3:Makefile.SSE3.gcc \
	raxmlHPC-AVX2:Makefile.AVX.gcc \
	raxmlHPC-PTHREADS:Makefile.PTHREADS.gcc \
	raxmlHPC-PTHREADS-SSE3:Makefile.SSE3.PTHREADS.gcc \
	raxmlHPC-PTHREADS-AVX2:Makefile.AVX.PTHREADS.gcc
RAXML_MAKEFILE := $(or $(word 2, $(subst :, ,$(filter $(RAXML_BINARY_VERSION):%,$(RAXML_MAKEFILE_MAP)))),Makefile.gcc)
RAXML_BUILT_NAME := $(patsubst %-AVX2,%-AVX,$(RAXML_BINARY_VERSION))

# ------------------------------------------------------------ the sources ----
ORTHOFINDER_REPO := https://github.com/OrthoFinder/OrthoFinder.git
DIAMOND_REPO := https://github.com/bbuchfink/diamond.git
DIAMOND_RELEASE := https://github.com/bbuchfink/diamond/releases/download/v$(DIAMOND_VERSION)
MCL_SOURCE := https://micans.org/mcl/src
FAMSA_REPO := https://github.com/refresh-bio/FAMSA.git
FAMSA_RELEASE := https://github.com/refresh-bio/FAMSA/releases/download/v$(FAMSA_VERSION)
FASTME_REPO := https://gite.lirmm.fr/atgc/FastME.git
FASTTREE_SOURCE := https://raw.githubusercontent.com/morgannprice/fasttree/v$(FASTTREE_VERSION)/FastTree.c
ASTER_REPO := https://github.com/chaoszhang/ASTER.git
MAFFT_REPO := https://gitlab.com/sysimm/mafft.git
RAXML_REPO := https://github.com/stamatak/standard-RAxML.git
RAXMLNG_REPO := https://github.com/amkozlov/raxml-ng.git
RAXMLNG_RELEASE := https://github.com/amkozlov/raxml-ng/releases/download/$(RAXMLNG_VERSION)
MUSCLE_RELEASE := https://github.com/rcedgar/muscle/releases/download/v$(MUSCLE_VERSION)
IQTREE_RELEASE := https://github.com/iqtree/iqtree3/releases/download/v$(IQTREE_VERSION)
BLAST_RELEASE := https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$(BLAST_VERSION)
MMSEQS_RELEASE := https://github.com/soedinglab/MMseqs2/releases/download/$(MMSEQS_VERSION)

# -------------------------------------------------------- where they go ----
SYSTEM_WIDE ?= false
PACKAGE_BIN ?= false
HOME_DIR := $(if $(HOME),$(HOME),$(shell echo ~))
USER_INSTALL_DIR ?= $(HOME_DIR)/local/bin
SYSTEM_INSTALL_DIR := /usr/local/bin
PACKAGE_BIN_DIR := $(CURDIR)/src/orthofinder/bin
ORTHOFINDER_DIR := $(USER_INSTALL_DIR)

ifeq ($(PACKAGE_BIN),true)
    BINARY_INSTALL_DIR := $(PACKAGE_BIN_DIR)
    REPLACE := true
    PORTABLE ?= true
else ifneq ($(BIN_DIR),)
    BINARY_INSTALL_DIR := $(abspath $(BIN_DIR))
else ifeq ($(SYSTEM_WIDE),true)
    BINARY_INSTALL_DIR := $(SYSTEM_INSTALL_DIR)
    SUDO_PREFIX := $(if $(shell [ "$(USER)" = "root" ] || echo 1),sudo)
    CHECK_PATH := true
else
    BINARY_INSTALL_DIR := $(USER_INSTALL_DIR)
    CHECK_PATH := true
    ADD_TO_PATH := true
endif
PORTABLE ?= false
REPLACE ?= false
CHECK_PATH ?= false

# Compiler flags of the builds: portable binaries have no CPU-specific code and,
# on Linux, are statically linked (they then run on other, older systems)
ifeq ($(PORTABLE),true)
    NATIVE_FLAGS :=
    STATIC_FLAGS := $(if $(filter Linux,$(shell uname -s)),-static)
else
    NATIVE_FLAGS := -march=native
    STATIC_FLAGS :=
endif

PYTHON := ./$(ENV_NAME)/bin/python3
PIP := ./$(ENV_NAME)/bin/pip
VENV_BIN := ./$(ENV_NAME)/bin

# --------------------------------------------- the shell of the recipes ----
# Functions for the recipes that install a program (each recipe is one shell)
define SH_FUNCTIONS
run() { if [ "$(QUIET)" = "true" ]; then "$$@" > /dev/null 2>&1; else "$$@"; fi; }
fail() { echo "Error: $$*" >&2; exit 1; }
download() { echo "Downloading $$1"; run curl -fL --retry 3 -o "$$2" "$$1" || fail "could not download $$1"; }
clone() { echo "Getting $$1 at $$2"; run git clone --depth 1 --branch "$$2" --recursive "$$1" "$$3" || fail "could not get $$1 at $$2"; }
put() { $(SUDO_PREFIX) install -m 755 "$$1" "$(BINARY_INSTALL_DIR)/$$2" || fail "could not install $$2 in $(BINARY_INSTALL_DIR)"; echo "Installed $(BINARY_INSTALL_DIR)/$$2"; }
wanted() {
if [ "$(FORCE)" = "true" ] || [ "$(REPLACE)" = "true" ]; then return 0; fi
if [ "$(CHECK_PATH)" = "true" ]; then
if command -v "$$2" > /dev/null 2>&1; then echo "$$2 is already installed: $$(command -v $$2). Skipping it (FORCE=true installs it again)."; return 1; fi
elif [ -e "$(BINARY_INSTALL_DIR)/$$1" ]; then echo "$(BINARY_INSTALL_DIR)/$$1 is already there. Skipping it (FORCE=true replaces it)."; return 1; fi
return 0
}
platform="$$(uname -s)-$$(uname -m)"
tmp="$$(mktemp -d)"
trap 'rm -rf "$$tmp"' EXIT
endef

# The directory the programs are installed in (and, for ~/local/bin, PATH)
bin_dir:
	@$(SUDO_PREFIX) mkdir -p "$(BINARY_INSTALL_DIR)"
	if [ "$(ADD_TO_PATH)" != "true" ]; then exit 0; fi
	case "$$(basename "$${SHELL}")" in
	zsh) rc="$${HOME}/.zshrc";;
	bash) if [ "$$(uname)" = "Darwin" ]; then rc="$${HOME}/.bash_profile"; else rc="$${HOME}/.bashrc"; fi;;
	*) rc="$${HOME}/.profile";;
	esac
	case ":$${PATH}:" in
	*":$(BINARY_INSTALL_DIR):"*) exit 0;;
	esac
	if ! grep -qsx 'export PATH="$(BINARY_INSTALL_DIR):$$PATH"' "$$rc"; then
		echo 'export PATH="$(BINARY_INSTALL_DIR):$$PATH"' >> "$$rc"
		echo "Added $(BINARY_INSTALL_DIR) to PATH in $$rc: restart the shell or run 'source $$rc'."
	fi

make_usr_bin: bin_dir

# ------------------------------------------------------- the programs ----

install_diamond: bin_dir
	@$(SH_FUNCTIONS)
	wanted diamond diamond || exit 0
	echo "Installing DIAMOND $(DIAMOND_VERSION) in $(BINARY_INSTALL_DIR)"
	case "$$platform" in
	Linux-x86_64) download $(DIAMOND_RELEASE)/diamond-linux64.tar.gz "$$tmp/diamond.tar.gz";;
	Darwin-*) download $(DIAMOND_RELEASE)/diamond-macos.tar.gz "$$tmp/diamond.tar.gz";;
	esac
	if [ -f "$$tmp/diamond.tar.gz" ]; then
		tar -xzf "$$tmp/diamond.tar.gz" -C "$$tmp"
	else   # no binary released for it: from the source (needs cmake)
		clone $(DIAMOND_REPO) v$(DIAMOND_VERSION) "$$tmp/src"
		mkdir -p "$$tmp/src/build" && cd "$$tmp/src/build"
		run cmake .. || fail "could not configure DIAMOND (is cmake installed?)"
		run make -j$(JOBS) || fail "could not build DIAMOND"
		cp diamond "$$tmp/diamond"
	fi
	put "$$tmp/diamond" diamond

install_mcl: bin_dir
	@$(SH_FUNCTIONS)
	wanted mcl mcl || exit 0
	echo "Installing MCL $(MCL_VERSION) in $(BINARY_INSTALL_DIR)"
	download $(MCL_SOURCE)/cimfomfa-$(CIMFOMFA_VERSION).tar.gz "$$tmp/cimfomfa.tar.gz"
	download $(MCL_SOURCE)/mcl-$(MCL_VERSION).tar.gz "$$tmp/mcl.tar.gz"
	tar -xzf "$$tmp/cimfomfa.tar.gz" -C "$$tmp" && tar -xzf "$$tmp/mcl.tar.gz" -C "$$tmp"
	cd "$$tmp/cimfomfa-$(CIMFOMFA_VERSION)"
	run ./configure --prefix="$$tmp/prefix" --disable-shared || fail "could not configure cimfomfa (for MCL)"
	run make -j$(JOBS) && run make install || fail "could not build cimfomfa (for MCL)"
	cd "$$tmp/mcl-$(MCL_VERSION)"
	run ./configure CFLAGS="-O2 -I$$tmp/prefix/include" LDFLAGS="-L$$tmp/prefix/lib $(STATIC_FLAGS)" --prefix="$$tmp/prefix" || fail "could not configure MCL"
	run make -j$(JOBS) && run make install || fail "could not build MCL"
	put "$$tmp/prefix/bin/mcl" mcl

install_famsa: bin_dir
	@$(SH_FUNCTIONS)
	wanted famsa famsa || exit 0
	echo "Installing FAMSA $(FAMSA_VERSION) in $(BINARY_INSTALL_DIR)"
	case "$$platform" in
	Linux-x86_64) asset=x64_linux;;
	Linux-aarch64|Linux-arm64) asset=arm64_linux;;
	Darwin-x86_64) asset=x64_mac;;
	Darwin-arm64) asset=arm64_mac;;
	*) asset=;;
	esac
	if [ -n "$$asset" ]; then
		download $(FAMSA_RELEASE)/famsa-v$(FAMSA_VERSION)-$$asset.tar.gz "$$tmp/famsa.tar.gz"
		tar -xzf "$$tmp/famsa.tar.gz" -C "$$tmp"
	else   # no binary released for it: from the source
		clone $(FAMSA_REPO) v$(FAMSA_VERSION) "$$tmp/src"
		cd "$$tmp/src" && run make -j$(JOBS) || fail "could not build FAMSA"
		cp famsa "$$tmp/famsa"
	fi
	put "$$tmp/famsa" famsa

install_fasttree: bin_dir
	@$(SH_FUNCTIONS)
	wanted FastTree FastTree || exit 0
	echo "Installing FastTree $(FASTTREE_VERSION) in $(BINARY_INSTALL_DIR)"
	download $(FASTTREE_SOURCE) "$$tmp/FastTree.c"
	cd "$$tmp"
	run $(CC) $(FASTTREE_CFLAGS) $(NATIVE_FLAGS) FastTree.c -lm $(STATIC_FLAGS) -o FastTree || fail "could not build FastTree"
	put "$$tmp/FastTree" FastTree

install_aster: bin_dir
	@$(SH_FUNCTIONS)
	wanted astral-pro astral-pro || exit 0
	echo "Installing ASTRAL-Pro (ASTER $(ASTER_VERSION)) in $(BINARY_INSTALL_DIR)"
	clone $(ASTER_REPO) v$(ASTER_VERSION) "$$tmp/src"
	cd "$$tmp/src"
	run $(CXX) -std=gnu++11 -O3 $(NATIVE_FLAGS) -pthread src/astral-pro.cpp $(STATIC_FLAGS) -o astral-pro || fail "could not build ASTRAL-Pro"
	put "$$tmp/src/astral-pro" astral-pro

install_fastme: bin_dir
	@$(SH_FUNCTIONS)
	wanted fastme fastme || exit 0
	echo "Installing FastME $(FASTME_VERSION) in $(BINARY_INSTALL_DIR)"
	clone $(FASTME_REPO) v$(FASTME_VERSION) "$$tmp/src"
	# its release tarball: the build files of the repository are links to the developers' automake
	tar -xzf "$$tmp/src/tarball/fastme-$(FASTME_VERSION).tar.gz" -C "$$tmp"
	cd "$$tmp/fastme-$(FASTME_VERSION)"
	run ./configure LDFLAGS="$(STATIC_FLAGS)" || fail "could not configure FastME"
	run make -j$(JOBS) || fail "could not build FastME"
	put "$$tmp/fastme-$(FASTME_VERSION)/src/fastme" fastme

# MAFFT is a script and the programs it runs (libexec): in src/orthofinder/bin as
# mafft/bin/mafft and mafft/libexec/*, elsewhere as <dir>/mafft and <dir>/../libexec/mafft/*
install_mafft: bin_dir
	@$(SH_FUNCTIONS)
	if [ "$(PACKAGE_BIN)" = "true" ]; then wanted mafft/bin/mafft mafft || exit 0; else wanted mafft mafft || exit 0; fi
	echo "Installing MAFFT $(MAFFT_VERSION) in $(BINARY_INSTALL_DIR)"
	clone $(MAFFT_REPO) v$(MAFFT_VERSION) "$$tmp/src"
	cd "$$tmp/src/core"
	if [ "$(PACKAGE_BIN)" = "true" ]; then
		run make -j$(JOBS) PREFIX="$$tmp/prefix" LDFLAGS="$(STATIC_FLAGS)" || fail "could not build MAFFT"
		run make install PREFIX="$$tmp/prefix" || fail "could not install MAFFT"
		rm -rf "$(BINARY_INSTALL_DIR)/mafft"
		mkdir -p "$(BINARY_INSTALL_DIR)/mafft/bin" "$(BINARY_INSTALL_DIR)/mafft/libexec"
		cp "$$tmp/prefix/bin/mafft" "$(BINARY_INSTALL_DIR)/mafft/bin/mafft"
		cp "$$tmp/prefix/libexec/mafft/"* "$(BINARY_INSTALL_DIR)/mafft/libexec/"
		echo "Installed $(BINARY_INSTALL_DIR)/mafft"
	else
		prefix="$$(dirname "$(BINARY_INSTALL_DIR)")"
		run make -j$(JOBS) PREFIX="$$prefix" BINDIR="$(BINARY_INSTALL_DIR)" || fail "could not build MAFFT"
		run $(SUDO_PREFIX) make install PREFIX="$$prefix" BINDIR="$(BINARY_INSTALL_DIR)" || fail "could not install MAFFT"
		echo "Installed $(BINARY_INSTALL_DIR)/mafft (its programs in $$prefix/libexec/mafft)"
	fi

install_iqtree3: bin_dir
	@$(SH_FUNCTIONS)
	wanted iqtree3 iqtree3 || exit 0
	echo "Installing IQ-TREE $(IQTREE_VERSION) in $(BINARY_INSTALL_DIR)"
	case "$$platform" in
	Linux-x86_64) url=$(IQTREE_RELEASE)/iqtree-$(IQTREE_VERSION)-Linux-intel.tar.gz;;   # any x86_64 (its "Linux" release runs this one)
	Linux-aarch64|Linux-arm64) url=$(IQTREE_RELEASE)/iqtree-$(IQTREE_VERSION)-Linux-arm.tar.gz;;
	Darwin-arm64) url=$(IQTREE_RELEASE)/iqtree-$(IQTREE_VERSION)-macOS-arm.zip;;
	Darwin-x86_64) url=$(IQTREE_RELEASE)/iqtree-$(IQTREE_VERSION)-macOS-intel.zip;;
	*) fail "no IQ-TREE release for $$platform";;
	esac
	download "$$url" "$$tmp/iqtree.archive"
	case "$$url" in
	*.zip) run unzip -o "$$tmp/iqtree.archive" -d "$$tmp";;
	*) tar -xzf "$$tmp/iqtree.archive" -C "$$tmp";;
	esac
	binary=$$(find "$$tmp" -type f -path "*/bin/iqtree3" | head -1)   # the binary itself (each release for one architecture)
	[ -n "$$binary" ] || fail "iqtree3 not found in $$url"
	put "$$binary" iqtree3

install_mmseqs: bin_dir
	@$(SH_FUNCTIONS)
	wanted mmseqs mmseqs || exit 0
	echo "Installing MMseqs2 $(MMSEQS_VERSION) in $(BINARY_INSTALL_DIR)"
	case "$$platform" in
	Linux-x86_64)
		if grep -qsw avx2 /proc/cpuinfo && [ "$(PORTABLE)" != "true" ]; then asset=linux-avx2; else asset=linux-sse41; fi;;
	Linux-aarch64|Linux-arm64) asset=linux-arm64;;
	Darwin-*) asset=osx-universal;;
	*) fail "no MMseqs2 release for $$platform";;
	esac
	download $(MMSEQS_RELEASE)/mmseqs-$$asset.tar.gz "$$tmp/mmseqs.tar.gz"
	tar -xzf "$$tmp/mmseqs.tar.gz" -C "$$tmp"
	put "$$tmp/mmseqs/bin/mmseqs" mmseqs

install_raxml: bin_dir
	@$(SH_FUNCTIONS)
	wanted $(RAXML_BINARY_VERSION) $(RAXML_BINARY_VERSION) || exit 0
	echo "Installing RAxML $(RAXML_VERSION) ($(RAXML_BINARY_VERSION)) in $(BINARY_INSTALL_DIR)"
	clone $(RAXML_REPO) v$(RAXML_VERSION) "$$tmp/src"
	cd "$$tmp/src"
	run make -j$(JOBS) -f $(RAXML_MAKEFILE) || fail "could not build RAxML"
	put "$$tmp/src/$(RAXML_BUILT_NAME)" $(RAXML_BINARY_VERSION)

install_raxmlng: bin_dir
	@$(SH_FUNCTIONS)
	wanted raxml-ng raxml-ng || exit 0
	echo "Installing RAxML-NG $(RAXMLNG_VERSION) in $(BINARY_INSTALL_DIR)"
	case "$$platform" in
	Linux-x86_64) asset=linux_x86_64;;
	Darwin-*) asset=macos;;
	*) asset=;;
	esac
	if [ -n "$$asset" ]; then
		download $(RAXMLNG_RELEASE)/raxml-ng_v$(RAXMLNG_VERSION)_$$asset.zip "$$tmp/raxml-ng.zip"
		run unzip -o "$$tmp/raxml-ng.zip" -d "$$tmp/raxml-ng"
		binary=$$(find "$$tmp/raxml-ng" -type f -name raxml-ng | head -1)
	else   # no binary released for it: from the source (needs cmake)
		clone $(RAXMLNG_REPO) $(RAXMLNG_VERSION) "$$tmp/src"
		mkdir -p "$$tmp/src/build" && cd "$$tmp/src/build"
		run cmake .. || fail "could not configure RAxML-NG (is cmake installed?)"
		run make -j$(JOBS) || fail "could not build RAxML-NG"
		binary=$$(find "$$tmp/src" -type f -name raxml-ng | head -1)
	fi
	[ -n "$$binary" ] || fail "raxml-ng not found"
	put "$$binary" raxml-ng

install_muscle: bin_dir
	@$(SH_FUNCTIONS)
	wanted muscle muscle || exit 0
	echo "Installing MUSCLE $(MUSCLE_VERSION) in $(BINARY_INSTALL_DIR)"
	case "$$platform" in
	Linux-x86_64) asset=linux-x86;;
	Linux-aarch64|Linux-arm64) asset=aarch64;;
	Darwin-arm64) asset=osx-arm64;;
	Darwin-x86_64) asset=osx-x86;;
	*) fail "no MUSCLE release for $$platform";;
	esac
	download $(MUSCLE_RELEASE)/muscle-$$asset.v$(MUSCLE_VERSION) "$$tmp/muscle"
	put "$$tmp/muscle" muscle

install_blast: bin_dir
	@$(SH_FUNCTIONS)
	wanted blastp blastp || exit 0
	echo "Installing BLAST+ $(BLAST_VERSION) ($(BLAST_PROGRAMS)) in $(BINARY_INSTALL_DIR)"
	case "$$platform" in
	Linux-x86_64) asset=x64-linux;;
	Linux-aarch64|Linux-arm64) asset=aarch64-linux;;
	Darwin-arm64) asset=aarch64-macosx;;
	Darwin-x86_64) asset=x64-macosx;;
	*) fail "no BLAST+ release for $$platform";;
	esac
	download $(BLAST_RELEASE)/ncbi-blast-$(BLAST_VERSION)+-$$asset.tar.gz "$$tmp/blast.tar.gz"
	tar -xzf "$$tmp/blast.tar.gz" -C "$$tmp"
	for program in $(BLAST_PROGRAMS); do put "$$tmp/ncbi-blast-$(BLAST_VERSION)+/bin/$$program" $$program; done

# BLAST+ built from the source (takes a long time)
install_blast_src: bin_dir
	@$(SH_FUNCTIONS)
	wanted blastp blastp || exit 0
	echo "Building BLAST+ $(BLAST_VERSION) from the source, for $(BINARY_INSTALL_DIR)"
	download $(BLAST_RELEASE)/ncbi-blast-$(BLAST_VERSION)+-src.tar.gz "$$tmp/blast.tar.gz"
	tar -xzf "$$tmp/blast.tar.gz" -C "$$tmp"
	cd "$$tmp/ncbi-blast-$(BLAST_VERSION)+-src/c++"
	run ./configure --prefix="$$tmp/prefix" || fail "could not configure BLAST+"
	run make -j$(JOBS) || fail "could not build BLAST+"
	run make install || fail "could not install BLAST+"
	for program in $(BLAST_PROGRAMS); do put "$$tmp/prefix/bin/$$program" $$program; done

# The programs OrthoFinder uses by default (those in src/orthofinder/bin)
install_tools: install_diamond install_mcl install_famsa install_fasttree install_aster install_fastme install_mafft
	@echo "The programs OrthoFinder uses by default are installed in $(BINARY_INSTALL_DIR)."

# ... and the optional ones (-S, -A, -T)
install_all_tools: install_tools install_iqtree3 install_raxml install_raxmlng install_muscle install_mmseqs install_blast
	@echo "All the programs are installed in $(BINARY_INSTALL_DIR)."

# Replace the programs in src/orthofinder/bin, the ones OrthoFinder runs
update_package_bin:
	@$(MAKE) --no-print-directory install_tools PACKAGE_BIN=true
	$(MAKE) --no-print-directory show_versions PACKAGE_BIN=true

# The versions of the programs in the install directory
show_versions:
	@d="$(BINARY_INSTALL_DIR)"
	echo "The programs in $$d:"
	show() {   # name, its file in $$d, the pattern of its version line, the command that prints it
		name=$$1; file=$$2; pattern=$$3; shift 3
		[ -x "$$d/$$file" ] || return 0
		out=$$({ "$$@" 2>&1 < /dev/null || true; } | grep -m1 -iE "$$pattern" | cut -c1-70) || true
		printf "  %-12s %s\n" "$$name" "$${out:-?}"
	}
	show diamond diamond "version" "$$d/diamond" version
	show mcl mcl "^mcl " "$$d/mcl" --version
	show famsa famsa "version|ver\." "$$d/famsa"
	show FastTree FastTree "^FastTree" "$$d/FastTree"
	show astral-pro astral-pro "^version" "$$d/astral-pro"
	show fastme fastme "FastME" "$$d/fastme" --version
	if [ -x "$$d/mafft/bin/mafft" ]; then show mafft mafft/bin/mafft "^v[0-9]" env MAFFT_BINARIES="$$d/mafft/libexec" "$$d/mafft/bin/mafft" --version; else show mafft mafft "^v[0-9]" "$$d/mafft" --version; fi
	show iqtree3 iqtree3 "version" "$$d/iqtree3" --version
	show raxml-ng raxml-ng "RAxML-NG" "$$d/raxml-ng" --version
	show muscle muscle "muscle" "$$d/muscle" -version
	show mmseqs mmseqs "version" "$$d/mmseqs"
	show blastp blastp "blastp" "$$d/blastp" -version

# ------------------------------------------------------------- conda ----

check_conda:
	@if ! command -v conda > /dev/null; then
		echo "Error: Conda is not installed. Please install Conda first."
		exit 1
	fi

create_conda_env: check_conda
	@. "$$(conda info --base)/etc/profile.d/conda.sh"
	if conda env list | grep -q "^$(ENV_NAME)[[:space:]]"; then
		if [ "$(FORCE)" != "true" ]; then
			echo "Conda environment $(ENV_NAME) already exists."
			exit 0
		fi
		echo "Recreating the conda environment $(ENV_NAME)..."
		conda env remove -n $(ENV_NAME) -y || { echo "Error: Failed to remove the conda environment $(ENV_NAME)."; exit 1; }
	fi
	echo "Creating the conda environment $(ENV_NAME) with Python $(CONDA_PYTHON_VERSION)..."
	conda create -n $(ENV_NAME) python=$(CONDA_PYTHON_VERSION) -y || { echo "Error: Failed to create the conda environment $(ENV_NAME)."; exit 1; }

# $(call conda_install,<bioconda package>,<name>,<command>,<version>)
define conda_install
@if [ "$(FORCE)" != "true" ] && command -v $(3) > /dev/null 2>&1; then
echo "$(2) is already installed: $$(command -v $(3)). Skipping it (FORCE=true installs it again)."
exit 0
fi
echo "Installing $(2) $(or $(4),(the latest version)) in the conda environment $(ENV_NAME)..."
. "$$(conda info --base)/etc/profile.d/conda.sh"
conda activate $(ENV_NAME)
conda install -y bioconda::$(1)$(if $(4),=$(4)) || { echo "Error: Failed to install $(2)."; exit 1; }
endef

conda_install_diamond: create_conda_env
	$(call conda_install,diamond,DIAMOND,diamond,$(DIAMOND_VERSION))

conda_install_mcl: create_conda_env
	$(call conda_install,mcl,MCL,mcl,$(CONDA_MCL_VERSION))

conda_install_famsa: create_conda_env
	$(call conda_install,famsa,FAMSA,famsa,$(CONDA_FAMSA_VERSION))

conda_install_fasttree: create_conda_env
	$(call conda_install,fasttree,FastTree,FastTree,$(FASTTREE_VERSION))

conda_install_fastme: create_conda_env
	$(call conda_install,fastme,FastME,fastme,$(FASTME_VERSION))

conda_install_aster: create_conda_env
	$(call conda_install,aster,ASTER (ASTRAL-Pro),astral-pro,$(ASTER_VERSION))

conda_install_mafft: create_conda_env
	$(call conda_install,mafft,MAFFT,mafft,$(CONDA_MAFFT_VERSION))

conda_install_iqtree3: create_conda_env
	$(call conda_install,iqtree,IQ-TREE,iqtree3,$(IQTREE_VERSION))

conda_install_raxml: create_conda_env
	$(call conda_install,raxml,RAxML,raxmlHPC,$(RAXML_VERSION))

conda_install_raxmlng: create_conda_env
	$(call conda_install,raxml-ng,RAxML-NG,raxml-ng,$(RAXMLNG_VERSION))

conda_install_muscle: create_conda_env
	$(call conda_install,muscle,MUSCLE,muscle,$(MUSCLE_VERSION))

conda_install_mmseqs: create_conda_env
	$(call conda_install,mmseqs2,MMseqs2,mmseqs,$(CONDA_MMSEQS_VERSION))

conda_install_blast: create_conda_env
	$(call conda_install,blast,BLAST+,blastp,$(BLAST_VERSION))

conda_install_orthofinder: create_conda_env
	$(call conda_install,orthofinder,OrthoFinder,orthofinder,$(ORTHOFINDER_VERSION))

conda_install_tools: conda_install_diamond conda_install_mcl conda_install_famsa conda_install_aster conda_install_fastme conda_install_fasttree conda_install_mafft
	@echo "The programs OrthoFinder uses by default are installed in $(ENV_NAME)."

conda_install: conda_install_orthofinder conda_install_tools
	@echo "You have now installed OrthoFinder$(if $(ORTHOFINDER_VERSION), $(ORTHOFINDER_VERSION)) and its dependencies in $(ENV_NAME)!"

clean_conda_env: check_conda
	@. "$$(conda info --base)/etc/profile.d/conda.sh"
	if conda env list | grep -q "^$(ENV_NAME)[[:space:]]"; then
		echo "Removing the conda environment $(ENV_NAME)..."
		conda remove -n $(ENV_NAME) --all -y || { echo "Error: Failed to remove $(ENV_NAME)."; exit 1; }
	else
		echo "Conda environment $(ENV_NAME) does not exist."
	fi

# ------------------------------------------------- OrthoFinder itself ----

venv:
	@echo "Checking for existing virtual environment $(ENV_NAME)..."
	if [ -d "$(ENV_NAME)" ] && [ "$(FORCE)" != "true" ]; then
		echo "Virtual environment $(ENV_NAME) already exists."
	else
		echo "Creating virtual environment $(ENV_NAME) using $(PYTHON_VERSION)..."
		if ! command -v $(PYTHON_VERSION) > /dev/null; then
			echo "Error: $(PYTHON_VERSION) is not installed or not in PATH."
			exit 1
		fi
		rm -rf $(ENV_NAME)
		$(PYTHON_VERSION) -m venv $(ENV_NAME) || { echo "Error: Failed to create virtual environment."; exit 1; }
		echo "Virtual environment $(ENV_NAME) created."
	fi

install_dependencies: venv
	$(PIP) install -r requirements.txt

install: install_tools venv bin_dir
	@echo "Checking global paths for OrthoFinder..."
	if [ "$(FORCE)" != "true" ] && command -v orthofinder > /dev/null 2>&1; then
		echo "OrthoFinder already exists at: $$(command -v orthofinder). Skipping installation."
		exit 0
	fi
	echo "Installing OrthoFinder..."
	if [ "$(QUIET)" = "true" ]; then
		$(PIP) install -e . > /dev/null 2>&1 || { echo "Error: Failed to install OrthoFinder."; exit 1; }
	else
		$(PIP) install -e . || { echo "Error: Failed to install OrthoFinder."; exit 1; }
	fi
	mkdir -p $(USER_INSTALL_DIR)
	cp $(VENV_BIN)/orthofinder $(USER_INSTALL_DIR)/ || { echo "Error: Failed to copy OrthoFinder to $(USER_INSTALL_DIR)."; exit 1; }
	echo "OrthoFinder copied to $(USER_INSTALL_DIR)."

run: install
	@echo "Running OrthoFinder..."
	if [ -f "$(VENV_BIN)/orthofinder" ]; then
		$(VENV_BIN)/orthofinder -f ExampleData
	else
		$(ORTHOFINDER_DIR)/orthofinder -f ExampleData
	fi

clean:
	@echo "Cleaning build environment..."
	find . -name __pycache__ -type d -prune -not -path "./$(ENV_NAME)/*" -exec rm -rf {} +
	rm -rf ./build ./dist ./*.egg-info ./src/*.egg-info

purge_orthofinder:
	@echo "Removing OrthoFinder, its environment $(ENV_NAME) and build files..."
	rm -rf $(ENV_NAME) ./build ./dist ./*.egg-info ./src/*.egg-info
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
	if [ -f "$(ORTHOFINDER_DIR)/orthofinder" ]; then
		rm -f "$(ORTHOFINDER_DIR)/orthofinder"
		echo "$(ORTHOFINDER_DIR)/orthofinder has been removed."
	else
		echo "OrthoFinder not found in $(ORTHOFINDER_DIR)."
	fi

purge_binaries:
	@if [ -d "$(USER_INSTALL_DIR)" ]; then
		echo "WARNING: This will remove the directory $(USER_INSTALL_DIR) and all its contents."
		read -p "Are you sure you want to proceed? (y/n): " confirm
		if [ "$$confirm" = "y" ] || [ "$$confirm" = "Y" ]; then
			rm -rf $(USER_INSTALL_DIR) || { echo "Error: Failed to remove $(USER_INSTALL_DIR)."; exit 1; }
			echo "$(USER_INSTALL_DIR) has been removed."
		else
			echo "Skipping removal of $(USER_INSTALL_DIR)."
		fi
	else
		echo "$(USER_INSTALL_DIR) does not exist. Nothing to remove."
	fi

# tools/view_results as a zip file to download: it works on its own next to the
# results (Python and ete4, not OrthoFinder). The copies of the page code are
# made current first.
view_results_zip:
	python3 tools/view_results/make_copies.py --zip dist/view_results.zip

.PHONY: bin_dir make_usr_bin install_diamond install_mcl install_famsa install_fasttree install_aster \
	install_fastme install_mafft install_iqtree3 install_mmseqs install_raxml install_raxmlng install_muscle \
	install_blast install_blast_src install_tools install_all_tools update_package_bin show_versions \
	check_conda create_conda_env clean_conda_env conda_install conda_install_tools conda_install_orthofinder \
	conda_install_diamond conda_install_mcl conda_install_famsa conda_install_fasttree conda_install_fastme \
	conda_install_aster conda_install_mafft conda_install_iqtree3 conda_install_raxml conda_install_raxmlng \
	conda_install_muscle conda_install_mmseqs conda_install_blast venv install_dependencies install run \
	clean purge_orthofinder purge_binaries view_results_zip
