.ONESHELL:
# SHELL := /bin/bash

.DEFAULT_GOAL := install

ifndef MAKE_VERSION
$(error This Makefile requires GNU Make. Please use gmake instead of make.)
endif

QUIET ?= true
FORCE ?= false
CONDA_PYTHON_VERSION ?= 3.10
PYTHON_VERSION ?= python3
ORTHOFINDER_ENV_DEFAULT := OF3_venv
ENV_NAME ?= $(ORTHOFINDER_ENV_DEFAULT)

ORTHOFINDER_DEFAULT_VERSION := 3.0.1b1
ORTHOFINDER_VERSION ?= $(ORTHOFINDER_DEFAULT_VERSION)

DIAMOND_DEFAULT_VERSION := 2.1.10
DIAMOND_VERSION ?= $(DIAMOND_DEFAULT_VERSION)

MCL_DEFAULT_VERSION := 22.282
MCL_VERSION ?= $(MCL_DEFAULT_VERSION)

FAMSA_DEFAULT_VERSION := 2.2.3
FAMSA_VERSION ?= $(FAMSA_DEFAULT_VERSION)

FASTME_DEFAULT_VERSION := 2.1.6.3
FASTME_VERSION ?= $(FASTME_DEFAULT_VERSION)

FASTTREE_DEFAULT_VERSION := 2.1.11
FASTTREE_VERSION ?= $(FASTTREE_DEFAULT_VERSION)

ASTER_DEFAULT_VERSION := 1.19-0
ASTER_VERSION ?= $(ASTER_DEFAULT_VERSION)

MAFFT_DEFAULT_VERSION := 7.525
MAFFT_VERSION ?= $(MAFFT_DEFAULT_VERSION)

RAXML_DEFAULT_VERSION := 8.2.1
RAXML_VERSION ?= $(RAXML_DEFAULT_VERSION)

RAXMLNG_DEFAULT_VERSION := 1.2.2
RAXMLNG_VERSION ?= $(RAXMLNG_DEFAULT_VERSION)

MUSCLE_DEFAULT_VERSION := 5.3
MUSCLE_VERSION ?= $(MUSCLE_DEFAULT_VERSION)

IQTREE_DEFAULT_VERSION := 2.3.6
IQTREE_VERSION ?= $(IQTREE_DEFAULT_VERSION)


SYSTEM_WIDE ?= 0
HOME_DIR := $(if $(HOME),$(HOME),$(shell echo ~))

PROMPT_USER_INSTALL_DIR = \
	@if [ -d "$(HOME_DIR)/local/bin" ]; then \
		read -p "Directory $(HOME_DIR)/local/bin exists. Do you want to use it? (y/n) " choice; \
		case $$choice in \
			[yY]*) echo "$(HOME_DIR)/local/bin";; \
			[nN]*) \
				read -p "Enter a new directory name (relative to $(HOME_DIR)): " new_dir; \
				USER_INSTALL_DIR=$(HOME_DIR)/$$new_dir; \
				echo $$USER_INSTALL_DIR; \
				mkdir -p $$USER_INSTALL_DIR || { echo "Error creating directory $$USER_INSTALL_DIR. Exiting."; exit 1; }; \
				;; \
			*) echo "Invalid choice. Exiting."; exit 1; \
		esac; \
	else \
		echo "$(HOME_DIR)/local/bin"; \
	fi

USER_INSTALL_DIR := $(shell $(PROMPT_USER_INSTALL_DIR))
SYSTEM_INSTALL_DIR := /usr/local/bin
ORTHOFINDER_DIR := $(USER_INSTALL_DIR)

ifeq ($(SYSTEM_WIDE),1)
    BINARY_INSTALL_DIR := $(SYSTEM_INSTALL_DIR)
    SUDO_PREFIX := $(if $(shell [ "$(USER)" = "root" ] || echo 1),sudo)
else
    BINARY_INSTALL_DIR := $(USER_INSTALL_DIR)
    SUDO_PREFIX :=
endif

# URLs for OrthoFinder repo
ORTHOFINDER_REPO := https://github.com/OrthoFinder/OrthoFinder.git
ORTHOFINDER_BINARY := $(BINARY_INSTALL_DIR)/orthofinder

# URLs for DIAMOND repo
DIAMOND_REPO := https://github.com/bbuchfink/diamond.git
DIAMOND_BINARY := $(BINARY_INSTALL_DIR)/diamond

#URLs for MCL repo
MCL_INSTALL_SCRIPT := https://raw.githubusercontent.com/micans/mcl/main/install-this-mcl.sh
MCL_BINARY := $(BINARY_INSTALL_DIR)/mcl

# URLs for FAMSA repo
FAMSA_REPO := https://github.com/refresh-bio/FAMSA.git
FAMSA_BINARY := $(BINARY_INSTALL_DIR)/famsa

# URLs for ASTER repo
ASTER_REPO := https://github.com/chaoszhang/ASTER.git
ASTRALPRO_BINARY := $(BINARY_INSTALL_DIR)/astral-pro

# URLs for FastME repo
FASTME_REPO := https://gite.lirmm.fr/atgc/FastME.git
FASTME_BINARY := $(BINARY_INSTALL_DIR)/fastme

# URLs for FastTree repo
FASTTREE_REPO := https://github.com/morgannprice/fasttree.git
FASTTREE_BINARY := $(BINARY_INSTALL_DIR)/fasttree

# URLs for MAFFT repo
MAFFT_REPO := https://gitlab.com/sysimm/mafft.git
MAFFT_BINARY := $(BINARY_INSTALL_DIR)/mafft

# URLs for RAXML repo
RAXML_REPO := https://github.com/stamatak/standard-RAxML.git
RAXML_BINARY := $(BINARY_INSTALL_DIR)/raxml

# URLs for RAXML-NG repo
RAXMLNG_REPO := https://github.com/amkozlov/raxml-ng.git
RAXMLNG_BINARY := $(BINARY_INSTALL_DIR)/raxml-ng

# URLs for MUSCLE repo
MUSCLE_REPO := https://github.com/rcedgar/muscle.git
MUSCLE_BINARY := $(BINARY_INSTALL_DIR)/muscle

# URLs for IQ-TREE urlS
IQTREE_LINUX_INTEL := https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-Linux-intel.tar.gz
IQTREE_LINUX_ARM := https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-Linux-arm.tar.gz

IQTREE_MACOS_UNIVERSAL := https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-macOS.zip
IQTREE_MACOS_INTEL := https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-macOS-intel.zip
IQTREE_MACOS_ARM := https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-macOS-arm.zip

IQTREE_BINARY := $(BINARY_INSTALL_DIR)/iqtree2


# URLs for BLAST urlS
BLAST_LINUX=https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz
BLAST_MACOS=https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-macosx.tar.gz
BLAST_MACOS_UNIVERSAL=https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-universal-macosx.tar.gz

BLAST_SRC_URL = https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-src.tar.gz
BLAST_SRC_ARCHIVE = ncbi-blast-2.16.0+-src.tar.gz
BLASTN_BINARY := $(BINARY_INSTALL_DIR)/blastn 
BLASTP_BINARY := $(BINARY_INSTALL_DIR)/blastp 

PYTHON := ./$(ENV_NAME)/bin/python3
PIP := ./$(ENV_NAME)/bin/pip
VENV_BIN := ./$(ENV_NAME)/bin


check_conda:
	@if ! command -v conda > /dev/null; then \
		echo "Error: Conda is not installed. Please install Conda first."; \
		exit 1; \
	fi
	@echo "Conda is installed. Activating Conda..."; \
	. $(shell conda info --base)/etc/profile.d/conda.sh || { echo "Error: Failed to initialize Conda. Exiting."; exit 1; }; \
	echo "Conda activated successfully."

create_conda_env: check_conda
	@echo "Checking if Conda environment $(ENV_NAME) exists..."; \
	. $(shell conda info --base)/etc/profile.d/conda.sh && \
	if conda env list | grep -q "^$(ENV_NAME)[[:space:]]" && [ "$(FORCE)" = "0" ]; then \
		echo "Conda environment $(ENV_NAME) already exists. Skipping creation."; \
	else \
		if conda env list | grep -q "^$(ENV_NAME)[[:space:]]"; then \
			echo "Forcing recreation of Conda environment $(ENV_NAME)..."; \
			conda env remove -n $(ENV_NAME) -y || { echo "Error: Failed to remove existing Conda environment. Exiting."; exit 1; }; \
		fi; \
		echo "Creating Conda environment: $(ENV_NAME) with Python $(CONDA_PYTHON_VERSION)..."; \
		conda create -n $(ENV_NAME) python=$(CONDA_PYTHON_VERSION) -y || { echo "Error: Failed to create Conda environment. Exiting."; exit 1; }; \
		echo "Conda environment $(ENV_NAME) with Python $(CONDA_PYTHON_VERSION) created successfully."; \
	fi

conda_install_diamond: create_conda_env
	@echo "Checking global paths for DIAMOND..."; \
	diamond_exists=$$(command -v diamond > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$diamond_exists" = "0" ]; then \
		echo "Installing DIAMOND version $(DIAMOND_VERSION) in $(ENV_NAME)..."; \
		. $(shell conda info --base)/etc/profile.d/conda.sh && \
		conda activate $(ENV_NAME) && \
		conda install bioconda::diamond=$(DIAMOND_VERSION) -y || { echo "Error: Failed to install DIAMOND. Exiting."; exit 1; }; \
		echo "DIAMOND version $(DIAMOND_VERSION) installed successfully."
	elif [ "$$diamond_exists" = "1" ]; then \
		echo "DIAMOND already exist globally. Skipping installation."; \
	fi; \


conda_install_mcl: create_conda_env
	@echo "Checking global paths for MCL..."; \
	mcl_exists=$$(command -v mcl > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$mcl_exists" = "0" ]; then \
		echo "Installing MCL version $(MCL_VERSION) in $(ENV_NAME)..."; \
		. $(shell conda info --base)/etc/profile.d/conda.sh && \
		conda activate $(ENV_NAME) && \
		conda install bioconda::mcl=$(MCL_VERSION) -y || { echo "Error: Failed to install MCL. Exiting."; exit 1; }; \
		echo "MCL version $(MCL_VERSION) installed successfully."
	elif [ "$$mcl_exists" = "1" ]; then \
		echo "MCL already exist globally. Skipping installation."; \
	fi; \

conda_install_famsa: create_conda_env
	@echo "Checking global paths for FAMSA..."; \
	famsa_exists=$$(command -v famsa > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$famsa_exists" = "0" ]; then \
		echo "Installing FAMSA version $(FAMSA_VERSION) in $(ENV_NAME)..."; \
		. $(shell conda info --base)/etc/profile.d/conda.sh && \
		conda activate $(ENV_NAME) && \
		conda install bioconda::famsa=$(FAMSA_VERSION) -y || { echo "Error: Failed to install FAMSA. Exiting."; exit 1; }; \
		echo "FAMSA version $(FAMSA_VERSION) installed successfully."
	elif [ "$$famsa_exists" = "1" ]; then \
		echo "FAMSA already exist globally. Skipping installation."; \
	fi; \

conda_install_fasttree: create_conda_env
	@echo "Checking global paths for FastTree..."; \
	fasttree_exists=$$(command -v fasttree > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$fasttree_exists" = "0" ]; then \
		echo "Installing FastTree version $(FASTTREE_VERSION) in $(ENV_NAME)..."; \
		. $(shell conda info --base)/etc/profile.d/conda.sh && \
		conda activate $(ENV_NAME) && \
		conda install bioconda::fasttree=$(FASTTREE_VERSION) -y || { echo "Error: Failed to install FastTree. Exiting."; exit 1; }; \
		echo "FastTree version $(FASTTREE_VERSION) installed successfully."
	elif [ "$$fasttree_exists" = "1" ]; then \
		echo "FastTree already exist globally. Skipping installation."; \
	fi; \

conda_install_fastme: create_conda_env
	@echo "Checking global paths for FastME..."; \
	fastme_exists=$$(command -v fastme > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$fastme_exists" = "0" ]; then \
		echo "Installing FastME version $(FASTME_VERSION) in $(ENV_NAME)..."; \
		. $(shell conda info --base)/etc/profile.d/conda.sh && \
		conda activate $(ENV_NAME) && \
		conda install bioconda::fastme=$(FASTME_VERSION) -y || { echo "Error: Failed to install FastME. Exiting."; exit 1; }; \
		echo "FastME version $(FASTME_VERSION) installed successfully."
	elif [ "$$fastme_exists" = "1" ]; then \
		echo "FastME already exist globally. Skipping installation."; \
	fi; \

conda_install_aster: create_conda_env
	@echo "Checking global paths for ASTRAL-Pro..."; \
	astralpro_exists=$$(command -v astral-pro > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$astralpro_exists" = "0" ]; then \
		echo "Installing ASTER version $(ASTER_VERSION) in $(ENV_NAME)..."; \
		. $(shell conda info --base)/etc/profile.d/conda.sh && \
		conda activate $(ENV_NAME) && \
		conda install aster=$(ASTER_VERSION) -y || { echo "Error: Failed to install ASTER. Exiting."; exit 1; }; \
		echo "ASTER version $(ASTER_VERSION) installed successfully."
	elif [ "$$astralpro_exists" = "1" ]; then \
		echo "ASTRAL-Pro already exist globally. Skipping installation."; \
	fi; \

conda_install_mafft: create_conda_env
	@echo "Checking global paths for MAFFT..."; \
	mafft_exists=$$(command -v mafft > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$mafft_exists" = "0" ]; then \
		echo "Installing MAFFT version $(mafft_VERSION) in $(ENV_NAME)..."; \
		. $(shell conda info --base)/etc/profile.d/conda.sh && \
		conda activate $(ENV_NAME) && \
		conda install bioconda::mafft=$(MAFFT_VERSION) -y || { echo "Error: Failed to install MAFFT. Exiting."; exit 1; }; \
		echo "MAFFT version $(MAFFT_VERSION) installed successfully."
	elif [ "$$mafft_exists" = "1" ]; then \
		echo "MAFFT already exist globally. Skipping installation."; \
	fi; \


conda_install_iqtree2: create_conda_env
	@echo "Checking global paths for IQ-TREE2..."; \
	iqtree2_exists=$$(command -v iqtree2 > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$iqtree2_exists" = "0" ]; then \
		echo "Installing IQ-TREE2 version $(IQTREE_VERSION) in $(ENV_NAME)..."; \
		. $(shell conda info --base)/etc/profile.d/conda.sh && \
		conda activate $(ENV_NAME) && \
		conda install bioconda::iqtree=$(IQTREE_VERSION) -y || { echo "Error: Failed to install IQ-TREE2. Exiting."; exit 1; }; \
		echo "IQ-TREE2 version $(IQTREE_VERSION) installed successfully."
	elif [ "$$iqtree2_exists" = "1" ]; then \
		echo "IQ-TREE2 already exist globally. Skipping installation."; \
	fi

conda_install_raxml: create_conda_env
	@echo "Checking global paths for RAxML..."; \
	raxml_exists=$$(command -v raxml > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$raxml_exists" = "0" ]; then \
		echo "Installing RAxML version $(RAXML_VERSION) in $(ENV_NAME)..."; \
		. $(shell conda info --base)/etc/profile.d/conda.sh && \
		conda activate $(ENV_NAME) && \
		conda install bioconda::raxml=$(RAXML_VERSION) -y || { echo "Error: Failed to install RAxML. Exiting."; exit 1; }; \
		echo "RAxML version $(RAXML_VERSION) installed successfully."
	elif [ "$$raxml_exists" = "1" ]; then \
		echo "RAxML already exist globally. Skipping installation."; \
	fi; \

conda_install_raxmlng: create_conda_env
	@echo "Checking global paths for RAXML-NG..."; \
	raxmlng_exists=$$(command -v raxml-ng > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$raxmlng_exists" = "0" ]; then \
		echo "Installing RAXML-NG version $(RAXMLNG_VERSION) in $(ENV_NAME)..."; \
		. $(shell conda info --base)/etc/profile.d/conda.sh && \
		conda activate $(ENV_NAME) && \
		conda install bioconda::raxml-ng=$(RAXMLNG_VERSION) -y || { echo "Error: Failed to install RAXML-NG. Exiting."; exit 1; }; \
		echo "RAXML-NG version $(RAXMLNG_VERSION) installed successfully."
	elif [ "$$raxmlng_exists" = "1" ]; then \
		echo "RAXML-NG already exist globally. Skipping installation."; \
	fi; \

conda_install_muscle: create_conda_env
	@echo "Checking global paths for MUSCLE..."; \
	muscle_exists=$$(command -v muscle > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$muscle_exists" = "0" ]; then \
		echo "Installing MUSCLE version $(MUSCLE_VERSION) in $(ENV_NAME)..."; \
		. $(shell conda info --base)/etc/profile.d/conda.sh && \
		conda activate $(ENV_NAME) && \
		conda install bioconda::muscle=$(MUSCLE_VERSION) -y || { echo "Error: Failed to install MUSCLE. Exiting."; exit 1; }; \
		echo "MUSCLE version $(MUSCLE_VERSION) installed successfully."
	elif [ "$$muscle_exists" = "1" ]; then \
		echo "MUSCLE already exist globally. Skipping installation."; \
	fi; \

conda_install_orthofinder: create_conda_env
	@echo "Checking global paths for OrthoFinder..."; \
	orthofinder_exists=$$(command -v orthofinder > /dev/null && echo 1 || echo 0); \

	if [ "$$orthofinder_exists" = "1" ]; then \
		echo " OrthoFinder already exist globally. Skipping installation."; \
	elif [ "$(FORCE)" = "true" ] || [ "$$orthofinder_exists" = "0" ]; then \
		echo "Installing  OrthoFinder version $(ORTHOFINDER_VERSION) in $(ENV_NAME)..."; \
		. $(shell conda info --base)/etc/profile.d/conda.sh && \
		conda activate $(ENV_NAME) && \
		conda install bioconda::orthofinder=$(ORTHOFINDER_VERSION) -y || { echo "Error: Failed to install OrthoFinder. Exiting."; exit 1; }; \
		echo " OrthoFinder version $(ORTHOFINDER_VERSION) installed successfully."
	fi; \

conda_install_tools: conda_install_diamond conda_install_mcl conda_install_famsa conda_install_aster conda_install_fastme conda_install_fasttree
	@echo "External Dependencies Installation Complete!"

conda_install: conda_install_orthofinder conda_install_tools
	@echo "You have now installed OrthoFinder $(ORTHOFINDER_VERSION) and its dependencies in $(ENV_NAME)!"

clean_conda_venv:
	@echo "Checking if Conda environment $(ENV_NAME) exists..."; \
	. $(shell conda info --base)/etc/profile.d/conda.sh && \
	if conda env list | grep -q "^$(ENV_NAME)[[:space:]]"; then \
		echo "Removing $(ENV_NAME) from Conda..."; \
		conda remove -n $(ENV_NAME) --all -y || { echo "Error: Failed to remove $(ENV_NAME) from Conda. Exiting."; exit 1; }; \
		echo "You have successfully removed $(ENV_NAME) from Conda!"; \
	else \
		echo "Conda environment $(ENV_NAME) does not exist. Skipping removal."; \
	fi

make_usr_bin:
	@if [ ! -d "$(USER_INSTALL_DIR)" ]; then \
		echo "Directory $(USER_INSTALL_DIR) does not exist. Creating it..."; \
		mkdir -p $(USER_INSTALL_DIR); \
	fi; \
	echo "Checking if $(USER_INSTALL_DIR) is already in the PATH..."; \
	if ! grep -qx 'export PATH="$(USER_INSTALL_DIR):$$PATH"' ~/.bashrc; then \
		echo "Adding $(USER_INSTALL_DIR) to the PATH permanently."; \
		echo 'export PATH="$(USER_INSTALL_DIR):$$PATH"' >> ~/.bashrc || { echo "Error: Failed to update PATH in ~/.bashrc. Exiting."; exit 1; }; \
		echo "PATH update added to ~/.bashrc. Please restart your shell or run 'source ~/.bashrc' to apply changes."; \
	else \
		echo "$(USER_INSTALL_DIR) is already in the PATH. Skipping addition to ~/.bashrc."; \
	fi

install_diamond: make_usr_bin
	@echo "Checking global paths for DIAMOND..."; \
	diamond_exists=$$(command -v diamond > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$diamond_exists" = "0" ]; then \
		echo "Creating a temporary directory for DIAMOND installation..."; \
		temp_dir=$$(mktemp -d); \
		echo "Cloning the latest version of DIAMOND from GitHub..."; \
		if [ "$(QUIET)" = "false" ]; then \
			git clone --depth 1 $(DIAMOND_REPO) $$temp_dir/diamond-src || { echo "Error: Failed to clone DIAMOND repository."; rm -rf $$temp_dir; exit 1; }; \
		else \
			git clone --depth 1 $(DIAMOND_REPO) $$temp_dir/diamond-src > /dev/null 2>&1 || { echo "Error: Failed to clone DIAMOND repository."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		cd $$temp_dir/diamond-src && \
		echo "Building DIAMOND..."; \
		mkdir -p build && cd build && \
		if [ "$(QUIET)" = "false" ]; then \
			cmake .. && make || { echo "Error: Failed to build DIAMOND."; rm -rf $$temp_dir; exit 1; }; \
		else \
			cmake .. > /dev/null 2>&1 && make > /dev/null 2>&1 || { echo "Error: Failed to build DIAMOND."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		echo "Moving DIAMOND binary to $(BINARY_INSTALL_DIR)..."; \
		$(SUDO_PREFIX) mv diamond $(BINARY_INSTALL_DIR) || { echo "Error: Failed to move DIAMOND binary."; rm -rf $$temp_dir; exit 1; }; \
		rm -rf $$temp_dir; \
		echo "DIAMOND installation completed successfully."; \
	else \
		diamond_path=$$(command -v diamond); \
		echo "DIAMOND already exists at: $$diamond_path. Skipping installation."; \
	fi


install_mcl: make_usr_bin
	@echo "Checking global paths for MCL..."; \
	mcl_exists=$$(command -v mcl > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$mcl_exists" = "0" ]; then \
		echo "Creating a temporary directory for MCL installation..."; \
		temp_dir=$$(mktemp -d); \
		cd $$temp_dir && \
		echo "Downloading the official MCL installation script..."; \
		if [ "$(QUIET)" = "true" ]; then \
			wget $(MCL_INSTALL_SCRIPT) -O install-this-mcl.sh > /dev/null 2>&1 || { echo "Error: Failed to download MCL installation script."; rm -rf $$temp_dir; exit 1; }; \
		else \
			wget $(MCL_INSTALL_SCRIPT) -O install-this-mcl.sh || { echo "Error: Failed to download MCL installation script."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		chmod u+x install-this-mcl.sh && \
		echo "Running the MCL installation script..."; \
		if [ "$(QUIET)" = "true" ]; then \
			./install-this-mcl.sh > /dev/null 2>&1 || { echo "Error: Failed to install MCL."; rm -rf $$temp_dir; exit 1; }; \
		else \
			./install-this-mcl.sh || { echo "Error: Failed to install MCL."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		rm -rf $$temp_dir $(HOME)/local/lib $(HOME)/local/share $(HOME)/local/include; \
		echo "MCL installation completed successfully."; \
	else \
		mcl_path=$$(command -v mcl); \
		echo "MCL already exists at: $$mcl_path. Skipping installation."; \
	fi


install_famsa: make_usr_bin
	@echo "Checking global paths for FAMSA..."; \
	famsa_exists=$$(command -v famsa > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$famsa_exists" = "0" ]; then \
		echo "Creating a temporary directory for FAMSA installation..."; \
		temp_dir=$$(mktemp -d); \
		echo "Cloning the latest version of FAMSA repository from GitHub..."; \
		if [ "$(QUIET)" = "true" ]; then \
			git clone --depth 1 $(FAMSA_REPO) --recursive $$temp_dir/famsa-src > /dev/null 2>&1 || { echo "Error: Failed to clone FAMSA repository."; rm -rf $$temp_dir; exit 1; }; \
		else \
			git clone --depth 1 $(FAMSA_REPO) --recursive $$temp_dir/famsa-src || { echo "Error: Failed to clone FAMSA repository."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		cd $$temp_dir/famsa-src && \
		echo "Building FAMSA from source..."; \
		if [ "$(QUIET)" = "true" ]; then \
			make > /dev/null 2>&1 || { echo "Error: Failed to build FAMSA."; rm -rf $$temp_dir; exit 1; }; \
		else \
			make || { echo "Error: Failed to build FAMSA."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		echo "Moving FAMSA binary to $(BINARY_INSTALL_DIR)..."; \
		$(SUDO_PREFIX) mv famsa $(BINARY_INSTALL_DIR) || { echo "Error: Failed to move FAMSA binary."; rm -rf $$temp_dir; exit 1; }; \
		rm -rf $$temp_dir; \
		echo "FAMSA installation completed successfully."; \
	else \
		famsa_path=$$(command -v famsa); \
		echo "FAMSA already exists at: $$famsa_path. Skipping installation."; \
	fi


install_fasttree: make_usr_bin
	@echo "Checking global paths for FastTree..."; \
	fasttree_exists=$$(command -v FastTree > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$fasttree_exists" = "0" ]; then \
		echo "Creating a temporary directory for FastTree installation..."; \
		temp_dir=$$(mktemp -d); \
		echo "Cloning the latest version of FastTree repository from GitHub..."; \
		if [ "$(QUIET)" = "true" ]; then \
			git clone --depth 1 $(FASTTREE_REPO) $$temp_dir/fasttree-src > /dev/null 2>&1 || { echo "Error: Failed to clone FastTree repository."; rm -rf $$temp_dir; exit 1; }; \
		else \
			git clone --depth 1 $(FASTTREE_REPO) $$temp_dir/fasttree-src || { echo "Error: Failed to clone FastTree repository."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		cd $$temp_dir/fasttree-src && \
		echo "Building FastTree from source..."; \
		if [ "$(QUIET)" = "true" ]; then \
			gcc -O3 -funroll-loops -fopenmp -DOPENMP FastTree.c -lm -o FastTree > /dev/null 2>&1 || { echo "Error: Failed to build FastTree."; rm -rf $$temp_dir; exit 1; }; \
		else \
			gcc -O3 -funroll-loops -fopenmp -DOPENMP FastTree.c -lm -o FastTree || { echo "Error: Failed to build FastTree."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		echo "Moving FastTree binary to $(BINARY_INSTALL_DIR)..."; \
		$(SUDO_PREFIX) mv FastTree $(BINARY_INSTALL_DIR) || { echo "Error: Failed to move FastTree binary."; rm -rf $$temp_dir; exit 1; }; \
		rm -rf $$temp_dir; \
		echo "FastTree installation completed successfully."; \
	else \
		fasttree_path=$$(command -v FastTree); \
		echo "FastTree already exists at: $$fasttree_path. Skipping installation."; \
	fi


install_aster: make_usr_bin
	@echo "Checking global paths for Astral-Pro..."; \
	astralpro_exists=$$(command -v astral-pro > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$astralpro_exists" = "0" ]; then \
		echo "Creating a temporary directory for Astral-Pro installation..."; \
		temp_dir=$$(mktemp -d); \
		echo "Cloning the latest version of ASTER repository from GitHub..."; \
		if [ "$(QUIET)" = "true" ]; then \
			git clone --depth 1 $(ASTER_REPO) $$temp_dir/aster-src > /dev/null 2>&1 || { echo "Error: Failed to clone ASTER repository."; rm -rf $$temp_dir; exit 1; }; \
		else \
			git clone --depth 1 $(ASTER_REPO) $$temp_dir/aster-src || { echo "Error: Failed to clone ASTER repository."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		cd $$temp_dir/aster-src && \
		echo "Building ASTER from source..."; \
		if [ "$(QUIET)" = "true" ]; then \
			make > /dev/null 2>&1 || { echo "Error: Failed to build ASTER."; rm -rf $$temp_dir; exit 1; }; \
		else \
			make || { echo "Error: Failed to build ASTER."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		echo "Moving Astral-Pro binaries to $(BINARY_INSTALL_DIR)..."; \
		$(SUDO_PREFIX) cp -r bin/* $(BINARY_INSTALL_DIR) || { echo "Error: Failed to move ASTER binaries."; rm -rf $$temp_dir; exit 1; }; \
		rm -rf $$temp_dir; \
		echo "Astral-Pro installation completed successfully."; \
	else \
		astralpro_path=$$(command -v astral-pro); \
		echo "Astral-Pro already exists at: $$astralpro_path. Skipping installation."; \
	fi


install_fastme: make_usr_bin
	@echo "Checking global paths for FastME..."; \
	fastme_exists=$$(command -v fastme > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$fastme_exists" = "0" ]; then \
		echo "Creating a temporary directory for FastME installation..."; \
		temp_dir=$$(mktemp -d); \
		echo "Cloning the latest version of FastME repository from GitLab..."; \
		if [ "$(QUIET)" = "true" ]; then \
			git clone --depth 1 $(FASTME_REPO) $$temp_dir/fastme-src > /dev/null 2>&1 || { echo "Error: Failed to clone FastME repository."; rm -rf $$temp_dir; exit 1; }; \
		else \
			git clone --depth 1 $(FASTME_REPO) $$temp_dir/fastme-src || { echo "Error: Failed to clone FastME repository."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		cd $$temp_dir/fastme-src && \
		echo "Preparing FastME build system..."; \
		if [ "$(QUIET)" = "true" ]; then \
			./configure > /dev/null 2>&1 || { echo "Error: Failed to configure FastME."; rm -rf $$temp_dir; exit 1; }; \
		else \
			./configure || { echo "Error: Failed to configure FastME."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		echo "Building FastME..."; \
		if [ "$(QUIET)" = "true" ]; then \
			make > /dev/null 2>&1 || { echo "Error: Failed to build FastME."; rm -rf $$temp_dir; exit 1; }; \
		else \
			make || { echo "Error: Failed to build FastME."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		echo "Moving FastME binary to $(BINARY_INSTALL_DIR)..."; \
		$(SUDO_PREFIX) mv src/fastme $(BINARY_INSTALL_DIR) || { echo "Error: Failed to move FastME binary."; rm -rf $$temp_dir; exit 1; }; \
		rm -rf $$temp_dir; \
		echo "FastME installation completed successfully."; \
	else \
		fastme_path=$$(command -v fastme); \
		echo "FastME already exists at: $$fastme_path. Skipping installation."; \
	fi


install_mafft: make_usr_bin
	@echo "Checking global paths for MAFFT..."; \
	mafft_exists=$$(command -v mafft > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$mafft_exists" = "0" ]; then \
		echo "Creating a temporary directory for MAFFT installation..."; \
		temp_dir=$$(mktemp -d); \
		echo "Cloning the latest version of MAFFT repository from GitLab..."; \
		if [ "$(QUIET)" = "true" ]; then \
			git clone --depth 1 $(MAFFT_REPO) $$temp_dir/mafft-src > /dev/null 2>&1 || { echo "Error: Failed to clone MAFFT repository."; rm -rf $$temp_dir; exit 1; }; \
		else \
			git clone --depth 1 $(MAFFT_REPO) $$temp_dir/mafft-src || { echo "Error: Failed to clone MAFFT repository."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		cd $$temp_dir/mafft-src/core && \
		echo "Updating PREFIX and BINDIR in the Makefile..."; \
		if [ "$(QUIET)" = "true" ]; then \
			sed -i 's|^PREFIX = /usr/local|PREFIX = $(USER_INSTALL_DIR)|' Makefile > /dev/null 2>&1 && \
			sed -i 's|^BINDIR = .*|BINDIR = $(BINARY_INSTALL_DIR)|' Makefile > /dev/null 2>&1 || { echo "Error: Failed to update Makefile."; rm -rf $$temp_dir; exit 1; }; \
		else \
			sed -i 's|^PREFIX = /usr/local|PREFIX = $(USER_INSTALL_DIR)|' Makefile && \
			sed -i 's|^BINDIR = .*|BINDIR = $(BINARY_INSTALL_DIR)|' Makefile || { echo "Error: Failed to update Makefile."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		echo "Building MAFFT from source..."; \
		if [ "$(QUIET)" = "true" ]; then \
			make clean > /dev/null 2>&1 && make > /dev/null 2>&1 || { echo "Error: Failed to build MAFFT."; rm -rf $$temp_dir; exit 1; }; \
		else \
			make clean && make || { echo "Error: Failed to build MAFFT."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		echo "Installing MAFFT to $(BINARY_INSTALL_DIR)..."; \
		if [ "$(QUIET)" = "true" ]; then \
			make install > /dev/null 2>&1 || { echo "Error: Failed to install MAFFT."; rm -rf $$temp_dir; exit 1; }; \
		else \
			make install || { echo "Error: Failed to install MAFFT."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		rm -rf $$temp_dir; \
		echo "MAFFT installation completed successfully."; \
	else \
		mafft_path=$$(command -v mafft); \
		echo "MAFFT already exists at: $$mafft_path. Skipping installation."; \
	fi


install_iqtree2: make_usr_bin
	@echo "Checking global paths for IQ-TREE2..."; \
	iqtree2_exists=$$(command -v iqtree2 > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$iqtree2_exists" = "0" ]; then \
		echo "Detecting system architecture..."; \
		OS=$$(uname -s); ARCH=$$(uname -m); \
		if [ "$$OS" = "Linux" ]; then \
			if [ "$$ARCH" = "x86_64" ]; then \
				IQTREE_URL=$(IQTREE_LINUX_INTEL); \
			elif [ "$$ARCH" = "aarch64" ]; then \
				IQTREE_URL=$(IQTREE_LINUX_ARM); \
			else \
				echo "Error: Unsupported Linux architecture: $$ARCH"; exit 1; \
			fi; \
		elif [ "$$OS" = "Darwin" ]; then \
			if [ "$$ARCH" = "arm64" ]; then \
				IQTREE_URL=$(IQTREE_MACOS_ARM); \
			elif [ "$$ARCH" = "x86_64" ]; then \
				IQTREE_URL=$(IQTREE_MACOS_INTEL); \
			else \
				IQTREE_URL=$(IQTREE_MACOS_UNIVERSAL); \
			fi; \
		else \
			echo "Error: Unsupported operating system: $$OS"; exit 1; \
		fi; \
		echo "Downloading IQ-TREE2 from $$IQTREE_URL..."; \
		temp_dir=$$(mktemp -d); \
		download_path=$$temp_dir/iqtree2-src; \
		if [ "$(QUIET)" = "true" ]; then \
			wget -O $$download_path $$IQTREE_URL > /dev/null 2>&1 || { echo "Error: Failed to download IQ-TREE2."; rm -rf $$temp_dir; exit 1; }; \
		else \
			wget -O $$download_path $$IQTREE_URL || { echo "Error: Failed to download IQ-TREE2."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		echo "Extracting IQ-TREE..."; \
		if echo "$$IQTREE_URL" | grep -q '.tar.gz'; then \
			if [ "$(QUIET)" = "true" ]; then \
				tar -xzf $$download_path -C $$temp_dir > /dev/null 2>&1 || { echo "Error: Failed to extract IQ-TREE2 tar.gz file."; rm -rf $$temp_dir; exit 1; }; \
			else \
				tar -xzf $$download_path -C $$temp_dir || { echo "Error: Failed to extract IQ-TREE2 tar.gz file."; rm -rf $$temp_dir; exit 1; }; \
			fi; \
		elif echo "$$IQTREE_URL" | grep -q '.zip'; then \
			if [ "$(QUIET)" = "true" ]; then \
				unzip -o $$download_path -d $$temp_dir > /dev/null 2>&1 || { echo "Error: Failed to extract IQ-TREE2 zip file."; rm -rf $$temp_dir; exit 1; }; \
			else \
				unzip -o $$download_path -d $$temp_dir || { echo "Error: Failed to extract IQ-TREE2 zip file."; rm -rf $$temp_dir; exit 1; }; \
			fi; \
		else \
			echo "Error: Unknown file format for IQ-TREE2."; rm -rf $$temp_dir; exit 1; \
		fi; \
		echo "Locating extracted IQ-TREE2 binary..."; \
		iqtree2_binary=$$(find $$temp_dir -type f -name "iqtree*" -executable | head -1); \
		if [ -z "$$iqtree2_binary" ]; then \
			echo "Error: IQ-TREE2 binary not found after extraction."; rm -rf $$temp_dir; exit 1; \
		fi; \
		echo "Moving IQ-TREE2 binary to $(BINARY_INSTALL_DIR)..."; \
		$(SUDO_PREFIX) mv $$iqtree2_binary $(BINARY_INSTALL_DIR) || { echo "Error: Failed to move IQ-TREE2 binary."; rm -rf $$temp_dir; exit 1; }; \
		rm -rf $$temp_dir; \
		echo "IQ-TREE2 installation completed successfully."; \
	else \
		iqtree2_path=$$(command -v iqtree2); \
		echo "IQ-TREE2 already exists at: $$iqtree2_path. Skipping installation."; \
	fi


install_raxml: make_usr_bin
	@echo "Checking global paths for RAxML..."; \
	raxml_exists=$$(command -v raxmlHPC > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$raxml_exists" = "0" ]; then \
		echo "Creating a temporary directory for RAxML installation..."; \
		temp_dir=$$(mktemp -d); \
		echo "Cloning the latest version of RAxML repository from GitHub..."; \
		if [ "$(QUIET)" = "true" ]; then \
			git clone --depth 1 $(RAXML_REPO) $$temp_dir/raxml-src > /dev/null 2>&1 || { echo "Error: Failed to clone RAxML repository."; rm -rf $$temp_dir; exit 1; }; \
		else \
			git clone --depth 1 $(RAXML_REPO) $$temp_dir/raxml-src || { echo "Error: Failed to clone RAxML repository."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		cd $$temp_dir/raxml-src && \
		echo "Building RAxML from source..."; \
		if [ "$(QUIET)" = "true" ]; then \
			make -f Makefile.gcc > /dev/null 2>&1 || { echo "Error: Failed to build RAxML."; rm -rf $$temp_dir; exit 1; }; \
		else \
			make -f Makefile.gcc || { echo "Error: Failed to build RAxML."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		echo "Moving RAxML binaries to $(BINARY_INSTALL_DIR)..."; \
		$(SUDO_PREFIX) mv raxmlHPC* $(BINARY_INSTALL_DIR)/ || { echo "Error: Failed to move RAxML binaries."; rm -rf $$temp_dir; exit 1; }; \
		rm -rf $$temp_dir; \
		echo "RAxML installation completed successfully."; \
	else \
		raxml_path=$$(command -v raxmlHPC); \
		echo "RAxML already exists at: $$raxml_path. Skipping installation."; \
	fi


install_raxmlng: make_usr_bin
	@echo "Checking global paths for RAxML-NG..."; \
	raxmlng_exists=$$(command -v raxml-ng > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$raxmlng_exists" = "0" ]; then \
		echo "Creating a temporary directory for RAxML-NG installation..."; \
		temp_dir=$$(mktemp -d); \
		echo "Detecting the best RAxML-NG version for your system..."; \
		AVX=$$(cat /proc/cpuinfo | grep -m 1 avx > /dev/null && echo 1 || echo 0); \
		MPI=$$(command -v mpirun > /dev/null && echo 1 || echo 0); \
		if [ "$$MPI" = "1" ]; then \
			echo "MPI detected. Building MPI version of RAxML-NG..."; \
			cmake_flags="-DUSE_MPI=ON"; \
		elif [ "$$AVX" = "0" ]; then \
			echo "No AVX support detected. Building Portable PTHREADS version of RAxML-NG..."; \
			cmake_flags="-DSTATIC_BUILD=ON -DENABLE_RAXML_SIMD=OFF -DENABLE_PLLMOD_SIMD=OFF"; \
		else \
			echo "Building PTHREADS version of RAxML-NG..."; \
			cmake_flags=""; \
		fi; \
		echo "Cloning the RAxML-NG repository from GitHub..."; \
		if [ "$(QUIET)" = "true" ]; then \
			git clone --recursive --depth 1 $(RAXMLNG_REPO) $$temp_dir/raxml-ng-src > /dev/null 2>&1 || { echo "Error: Failed to clone RAxML-NG repository."; rm -rf $$temp_dir; exit 1; }; \
		else \
			git clone --recursive --depth 1 $(RAXMLNG_REPO) $$temp_dir/raxml-ng-src || { echo "Error: Failed to clone RAxML-NG repository."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		cd $$temp_dir/raxml-ng-src && \
		mkdir build && cd build && \
		echo "Configuring RAxML-NG build system..."; \
		if [ "$(QUIET)" = "true" ]; then \
			cmake -DCMAKE_INSTALL_PREFIX=$(BINARY_INSTALL_DIR) $$cmake_flags .. > /dev/null 2>&1 || { echo "Error: Failed to configure RAxML-NG."; rm -rf $$temp_dir; exit 1; }; \
		else \
			cmake -DCMAKE_INSTALL_PREFIX=$(BINARY_INSTALL_DIR) $$cmake_flags .. || { echo "Error: Failed to configure RAxML-NG."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		echo "Building RAxML-NG..."; \
		if [ "$(QUIET)" = "true" ]; then \
			make > /dev/null 2>&1 || { echo "Error: Failed to build RAxML-NG."; rm -rf $$temp_dir; exit 1; }; \
		else \
			make || { echo "Error: Failed to build RAxML-NG."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		echo "Installing RAxML-NG..."; \
		if [ "$(QUIET)" = "true" ]; then \
			make install > /dev/null 2>&1 || { echo "Error: Failed to install RAxML-NG."; rm -rf $$temp_dir; exit 1; }; \
		else \
			make install || { echo "Error: Failed to install RAxML-NG."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		echo "Moving RAxML-NG binary to $(BINARY_INSTALL_DIR)..."; \
		$(SUDO_PREFIX) mv $(BINARY_INSTALL_DIR)/bin/raxml-ng $(BINARY_INSTALL_DIR)/raxml-ng || { echo "Error: Failed to move RAxML-NG binary."; rm -rf $$temp_dir; exit 1; }; \
		rm -rf $(BINARY_INSTALL_DIR)/bin; \
		rm -rf $$temp_dir; \
		echo "RAxML-NG installation completed successfully."; \
	else \
		raxmlng_path=$$(command -v raxml-ng); \
		echo "RAxML-NG already exists at: $$raxmlng_path. Skipping installation."; \
	fi


install_muscle: make_usr_bin
	@echo "Checking global paths for MUSCLE..."; \
	muscle_exists=$$(command -v muscle > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$muscle_exists" = "0" ]; then \
		echo "Creating a temporary directory for MUSCLE installation..."; \
		temp_dir=$$(mktemp -d); \
		echo "Cloning the latest version of MUSCLE repository from GitHub..."; \
		git clone --depth 1 $(MUSCLE_REPO) $$temp_dir/muscle-src $(if $(QUIET),> /dev/null 2>&1) || { echo "Error: Failed to clone MUSCLE repository."; rm -rf $$temp_dir; exit 1; }; \
		cd $$temp_dir/muscle-src/src && \
		if [ "$$(uname)" = "Linux" ]; then \
			echo "Building MUSCLE for Linux..."; \
			chmod u+x build_linux.bash; \
			if [ "$(QUIET)" = "true" ]; then \
				./build_linux.bash > /dev/null 2>&1 || { echo "Error: Failed to build MUSCLE for Linux."; rm -rf $$temp_dir; exit 1; }; \
			else \
				./build_linux.bash || { echo "Error: Failed to build MUSCLE for Linux."; rm -rf $$temp_dir; exit 1; }; \
			fi; \
		elif [ "$$(uname)" = "Darwin" ]; then \
			echo "Building MUSCLE for macOS..."; \
			chmod u+x build_osx.bash; \
			if [ "$(QUIET)" = "true" ]; then \
				./build_osx.bash > /dev/null 2>&1 || { echo "Error: Failed to build MUSCLE for macOS."; rm -rf $$temp_dir; exit 1; }; \
			else \
				./build_osx.bash || { echo "Error: Failed to build MUSCLE for macOS."; rm -rf $$temp_dir; exit 1; }; \
			fi; \
		else \
			echo "Error: Unsupported operating system. MUSCLE can only be built on Linux or macOS."; \
			rm -rf $$temp_dir; \
			exit 1; \
		fi; \
		echo "Moving MUSCLE binary to $(BINARY_INSTALL_DIR)..."; \
		$(SUDO_PREFIX) mv $$temp_dir/muscle-src/bin/muscle $(BINARY_INSTALL_DIR) || { echo "Error: Failed to move MUSCLE binary."; rm -rf $$temp_dir; exit 1; }; \
		rm -rf $$temp_dir; \
		echo "MUSCLE installation completed successfully."; \
	else \
		muscle_path=$$(command -v muscle); \
		echo "MUSCLE already exists at: $$muscle_path. Skipping installation."; \
	fi


install_blast: make_usr_bin
	@echo "Checking global paths for BLAST..."; \
	blastn_exists=$$(command -v blastn > /dev/null && echo 1 || echo 0); \
	blastp_exists=$$(command -v blastp > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$blastn_exists" = "0" ] || [ "$$blastp_exists" = "0" ]; then \
		echo "Detecting system architecture..."; \
		OS=$$(uname -s); ARCH=$$(uname -m); \
		if [ "$$OS" = "Linux" ]; then \
			if [ "$$ARCH" = "x86_64" ]; then \
				BLAST_URL=$(BLAST_LINUX); \
			else \
				echo "Error: Unsupported Linux architecture: $$ARCH"; exit 1; \
			fi; \
		elif [ "$$OS" = "Darwin" ]; then \
			if [ "$$ARCH" = "arm64" ]; then \
				BLAST_URL=$(BLAST_MACOS_UNIVERSAL); \
			elif [ "$$ARCH" = "x86_64" ]; then \
				BLAST_URL=$(BLAST_MACOS); \
			else \
				echo "Error: Unsupported macOS architecture: $$ARCH"; exit 1; \
			fi; \
		else \
			echo "Error: Unsupported operating system: $$OS"; exit 1; \
		fi; \
		echo "Downloading BLAST from $$BLAST_URL..."; \
		temp_dir=$$(mktemp -d); \
		blast_archive=$$temp_dir/blast.tar.gz; \
		if [ "$(QUIET)" = "true" ]; then \
			wget -O $$blast_archive $$BLAST_URL > /dev/null 2>&1 || { echo "Error: Failed to download BLAST."; rm -rf $$temp_dir; exit 1; }; \
		else \
			wget -O $$blast_archive $$BLAST_URL || { echo "Error: Failed to download BLAST."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		echo "Extracting BLAST archive..."; \
		if [ "$(QUIET)" = "true" ]; then \
			tar -xzf $$blast_archive -C $$temp_dir > /dev/null 2>&1 || { echo "Error: Failed to extract BLAST."; rm -rf $$temp_dir; exit 1; }; \
		else \
			tar -xzf $$blast_archive -C $$temp_dir || { echo "Error: Failed to extract BLAST."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		blast_binary_dir=$$(find $$temp_dir -mindepth 1 -maxdepth 1 -type d)/bin; \
		if [ -z "$$blast_binary_dir" ]; then \
			echo "Error: BLAST binaries not found after extraction."; rm -rf $$temp_dir; exit 1; \
		fi; \
		echo "Moving BLAST binaries to $(BINARY_INSTALL_DIR)..."; \
		$(SUDO_PREFIX) cp -r $$blast_binary_dir/* $(BINARY_INSTALL_DIR) || { echo "Error: Failed to move BLAST binaries."; rm -rf $$temp_dir; exit 1; }; \
		rm -rf $$temp_dir; \
		echo "BLAST installation completed successfully."; \
	else \
		blastn_path=$$(command -v blastn); \
		blastp_path=$$(command -v blastp); \
		echo "BLASTN already exists at: $$blastn_path."; \
		echo "BLASTP already exists at: $$blastp_path."; \
		echo "Skipping installation."; \
	fi

install_blast_src: make_usr_bin
	@echo "Checking if BLAST is already installed..."; \
	blastn_exists=$$(command -v blastn > /dev/null && echo 1 || echo 0); \
	blastp_exists=$$(command -v blastp > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$blastn_exists" = "0" ] || [ "$$blastp_exists" = "0" ]; then \
		echo "Creating a temporary directory for BLAST installation..."; \
		temp_dir=$$(mktemp -d); \
		echo "Downloading BLAST source..."; \
		blast_archive_path=$$temp_dir/$(BLAST_SRC_ARCHIVE); \
		if [ "$(QUIET)" = "true" ]; then \
			wget $(BLAST_SRC_URL) -O $$blast_archive_path > /dev/null 2>&1 || { echo "Error: Failed to download BLAST source."; rm -rf $$temp_dir; exit 1; }; \
		else \
			wget $(BLAST_SRC_URL) -O $$blast_archive_path || { echo "Error: Failed to download BLAST source."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		echo "Extracting BLAST source..."; \
		if [ "$(QUIET)" = "true" ]; then \
			tar -xvzf $$blast_archive_path -C $$temp_dir > /dev/null 2>&1 || { echo "Error: Failed to extract BLAST source."; rm -rf $$temp_dir; exit 1; }; \
		else \
			tar -xvzf $$blast_archive_path -C $$temp_dir || { echo "Error: Failed to extract BLAST source."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		blast_src_dir=$$(find $$temp_dir -mindepth 1 -maxdepth 1 -type d); \
		cd $$blast_src_dir/c++ && \
		echo "Preparing BLAST build system..."; \
		if [ "$(QUIET)" = "true" ]; then \
			./configure > /dev/null 2>&1 || { echo "Error: Failed to configure BLAST build."; rm -rf $$temp_dir; exit 1; }; \
		else \
			./configure || { echo "Error: Failed to configure BLAST build."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		echo "Building BLAST..."; \
		if [ "$(QUIET)" = "true" ]; then \
			make > /dev/null 2>&1 || { echo "Error: Failed to build BLAST."; rm -rf $$temp_dir; exit 1; }; \
		else \
			make || { echo "Error: Failed to build BLAST."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		echo "Installing BLAST binaries..."; \
		if [ "$(QUIET)" = "true" ]; then \
			$(SUDO_PREFIX) make install > /dev/null 2>&1 || { echo "Error: Failed to install BLAST binaries."; rm -rf $$temp_dir; exit 1; }; \
		else \
			$(SUDO_PREFIX) make install || { echo "Error: Failed to install BLAST binaries."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		rm -rf $$temp_dir; \
		echo "BLAST installation completed successfully."; \
	else \
		blastn_path=$$(command -v blastn); \
		blastp_path=$$(command -v blastp); \
		echo "BLASTN already exists at: $$blastn_path. Skipping installation."; \
		echo "BLASTP already exists at: $$blastp_path. Skipping installation."; \
	fi


install_tools: install_diamond install_mcl install_aster install_famsa install_fasttree install_famsa
	@echo "All the essential external dependencies required by OrthoFinder have all been installed!"

venv:
	@echo "Checking for existing virtual environment $(ENV_NAME)..."
	@if [ -d "$(ENV_NAME)" ] && [ "$(FORCE)" = "0" ]; then \
		echo "Virtual environment $(ENV_NAME) already exists."; \
		echo "Activating $(ENV_NAME)..."; \
		. $(ENV_NAME)/bin/activate; \
		echo "Virtual environment $(ENV_NAME) activated successfully."; \
	else \
		echo "Creating virtual environment $(ENV_NAME) using $(PYTHON_VERSION)..."; \
		if ! command -v $(PYTHON_VERSION) > /dev/null; then \
			echo "Error: $(PYTHON_VERSION) is not installed or not in PATH. Exiting."; \
			exit 1; \
		fi; \
		rm -rf $(ENV_NAME); \
		$(PYTHON_VERSION) -m venv $(ENV_NAME) || { echo "Error: Failed to create virtual environment. Exiting."; exit 1; }; \
		chmod +x $(ENV_NAME)/bin/activate; \
		. $(ENV_NAME)/bin/activate; \
		echo "Virtual environment $(ENV_NAME) created and activated successfully."; \
	fi

install_dependencies: venv 
	$(PIP) install -r requirements.txt

install: install_tools venv make_usr_bin
	@echo "Checking global paths for OrthoFinder..."; \
	orthofinder_exists=$$(command -v orthofinder > /dev/null 2>&1 && echo 1 || echo 0); \
	if [ "$(FORCE)" = "true" ] || [ "$$orthofinder_exists" = "0" ]; then \
		echo "Installing OrthoFinder..."; \
		if [ "$(QUIET)" = "true" ]; then \
			$(PIP) install -e . > /dev/null 2>&1 || { echo "Error: Failed to install OrthoFinder. Exiting."; exit 1; }; \
		else \
			$(PIP) install -e . || { echo "Error: Failed to install OrthoFinder. Exiting."; exit 1; }; \
		fi; \
		mkdir -p $(USER_INSTALL_DIR) || { echo "Error: Failed to create directory $(USER_INSTALL_DIR). Exiting."; exit 1; }; \
		echo "Copying OrthoFinder to $(USER_INSTALL_DIR)..."; \
		cp $(VENV_BIN)/orthofinder $(USER_INSTALL_DIR)/ && \
		echo "OrthoFinder successfully copied to $(USER_INSTALL_DIR)." || \
		{ echo "Error: Failed to copy OrthoFinder to $(USER_INSTALL_DIR). Exiting."; exit 1; }; \
	elif [ "$$orthofinder_exists" = "1" ]; then \
		orthofinder_path=$$(command -v orthofinder); \
		echo "OrthoFinder already exists at: $$orthofinder_path. Skipping installation."; \
	fi


run: install_tools install
	@echo "Running OrthoFinder..."
	@if [ -f "$(VENV_BIN)/orthofinder" ]; then \
		$(VENV_BIN)/orthofinder -f ExampleData 
	elif [ -f "$(ORTHOFINDER_DIR)/orthofinder" ]; then \
		$(ORTHOFINDER_DIR)/orthofinder -f ExampleData 
	else \
		echo "OrthoFinder not found. Installing..."; \
		$(MAKE) install; \
	fi

clean:
	@echo "Cleaning build environment..."
	rm -rf **/__pycache__
	rm -rf ./build ./dist ./orthofinder.egg-info

purge_orthofinder:
	@echo "Remove OrthoFinder and it's environment and dependencies..."
	rm -rf $(ENV_NAME) 
	rm -rf **/__pycache__
	rm -rf ./build ./dist ./orthofinder.egg-info

	if [ -f "$(ORTHOFINDER_DIR)/orthofinder" ]; then \
		rm -f "$(ORTHOFINDER_DIR)/orthofinder"
		echo "$(ORTHOFINDER_DIR)/orthofinder including it's environment have been removed..."; \
	else \
		echo "OrthoFinder not found..."; \
	fi

purge_binaries:
	@if [ -d "$(USER_INSTALL_DIR)" ]; then \
		echo "WARNING: This will remove the directory $(USER_INSTALL_DIR) and all its contents."; \
		read -p "Are you sure you want to proceed? (y/n): " confirm; \
		if [ "$$confirm" = "y" ] || [ "$$confirm" = "Y" ]; then \
			echo "Removing $(USER_INSTALL_DIR)..."; \
			rm -rf $(USER_INSTALL_DIR) || { echo "Error: Failed to remove $(USER_INSTALL_DIR)."; exit 1; }; \
			echo "$(USER_INSTALL_DIR) has been removed."; \
		else \
			echo "Skipping removal of $(USER_INSTALL_DIR)."; \
		fi; \
	else \
		echo "$(USER_INSTALL_DIR) does not exist. Nothing to remove."; \
	fi

.PHONY: make_usr_bin clean clean_conda_venv purge_orthofinder purge_binaries
