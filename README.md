# OrthoFinder [![CI](https://github.com/OrthoFinder/OrthoFinder/actions/workflows/orthofinder-python-main.yml/badge.svg?branch=main)](https://github.com/OrthoFinder/OrthoFinder/actions/workflows/orthofinder-python-main.yml)[![GitHub release](https://img.shields.io/github/v/release/OrthoFinder/OrthoFinder)](https://github.com/OrthoFinder/OrthoFinder/releases)![License](https://img.shields.io/github/license/OrthoFinder/OrthoFinder)![Stars](https://img.shields.io/github/stars/OrthoFinder/OrthoFinder?style=social)


<p align="center">
  <img src="assets/images/coreassignfigure-images-2.png" alt="concept_figure" width="800"/>
</p>

OrthoFinder identifies orthogroups, infers gene trees for all orthogroups, and analyzes the gene trees to identify the rooted species tree. The method subsequently identifies all gene duplication events in the complete set of gene trees, and analyses them at both gene tree and species tree level. OrthoFinder further analyzes all of this phylogenetic information to identify the complete set of orthologs between all species, and provides extensive comparative genomics statistics.

## Table of contents

- [What's new in v3.2](#whats-new-in-v32)
- [Installation](#installation)
- [Simple Usage](#simple-usage)
- [Advanced Usage - Scaling to Thousands of Species](#advanced-usage---scaling-to-thousands-of-species)
- [Command line Options](#command-line-options)
- [Output files](#output-files)
- [Latest additions](#latest-additions)
- [Citation](#citation)
- [System Requirements](#system-requirements)

For more information please visit [our website](https://orthofinder.github.io/OrthoFinder/).

## What's new in v3.2

These changes are in OrthoFinder v3.2 and later. Some of them change the names and files of the results: please check them if you have scripts that read OrthoFinder results.

1. **"Orthologue" is now "ortholog"**, in the names of the results and in the messages:

   | Before v3.2 | From v3.2 |
   |-------------|-----------|
   | `Orthologues/` | `Orthologs/` |
   | `Orthologues/Orthologues_<species>/` | `Orthologs/Orthologs_<species>/` (with `-pof`) |
   | `Single_Copy_Orthologue_Sequences/` | `Single_Copy_Ortholog_Sequences/` |
   | `Comparative_Genomics_Statistics/OrthologuesStats_*.tsv` | `Comparative_Genomics_Statistics/OrthologsStats_*.tsv` |
   | `Orthogroups/Orthogroups_SingleCopyOrthologues.txt` | `Orthogroups/Orthogroups_SingleCopyOrthologs.txt` |

   The contents of the files are unchanged. A restart (`-b`, `-fg`, `-fgt`, `-fst`) in results from an earlier version first renames the old names to the new ones.

2. **With `-ST` (e.g. `-ST astral-pro`), the concatenated alignment is no longer written.** The species tree is then inferred from the gene trees, and the concatenated alignment of the orthogroups is not needed for it. To write it as well, add `--concat-alignment` (it is then in `WorkingDirectory/Alignments_ids/SpeciesTreeAlignment.fa`, with the orthogroups used listed in `Species_Tree/Orthogroups_for_concatenated_alignment.txt`), or for the results of a finished run use `species_tree_alignment <Results_dir>`. Without `-ST` (the default) the concatenated alignment is written as before.

3. **Scoring matrices and gap penalties for DIAMOND searches**: `--matrix` (e.g. `BLOSUM45`, default `BLOSUM62`), `--custom-matrix <file>` for a matrix file of your own, and `-go`/`--gapopen` and `-ge`/`--gapextend` for the gap penalties. A restart or `--assign` from the results uses the same scoring. See [Sequence search scoring options (DIAMOND)](#command-line-options).

4. **Pairwise ortholog files only with `-pof`.** By default, the orthologs are written as one compressed file per species, `Orthologs/<species>.tsv.gz`, which takes far less space and far fewer files with many species. With `-pof`, OrthoFinder writes a file for each pair of species as before, `Orthologs/Orthologs_<species>/<species>__v__<other species>.tsv`. For the results of a run without `-pof`, `split_ortholog_files <Results_dir>` writes the pairwise files afterwards.

5. **The results in a web browser**: `--tree-view` (the species tree, the resolved gene trees with their alignments, and a species tree explorer), `--stats-view` (the statistics, with interactive plots) and `--groups-view` (the genes of each orthogroup and HOG, searchable). `--view` writes all three, linked to each other, and `OrthoFinder_Results.html` to open them, a tab for each. For the results of a finished run: `view_results <Results_dir>`. The pages work offline in any recent browser, in light or dark colours. See [Viewing the results in a web browser](#command-line-options).

6. **Stopping and restarting an analysis at each stage**, also for `--assign`:
   - Stop after a stage: `-op` (the search commands, written to `WorkingDirectory/blast_commands.txt`, e.g. to run on a cluster), `-og` (orthogroups), `-ogt` (gene trees) or `-ost` (species tree). The end of the run prints the command to continue it.
   - Restart from a stage: `-b` (search results), `-fg` (orthogroups), `-fgt` (gene trees) or `-fst` (species tree), given the results directory.
   - A restart continues with the settings of the run (search program, scoring, `-M`, `-A`, `-T`, `-ST`, `-I` etc., kept in `WorkingDirectory/run_state.json`), unless they are given again.
   - With `-b`, OrthoFinder runs only the searches whose results are missing or incomplete, e.g. after `-op`, or after a run that was interrupted (Ctrl-C, a job time limit) or failed. The finished searches are recorded in `WorkingDirectory/checkpoint.txt`; the results of an interrupted search are moved aside (`.incomplete`) and the search is run again.
   - With `-op`, the search commands are no longer printed to the terminal: they are written to `WorkingDirectory/blast_commands.txt`. A restart with `-b` reuses the search results that are there, and runs only the missing ones, with the search program of the run. It does not create a new results directory: the results are written to the one created by `-op`.

Earlier change, from v3.1: `N0.tsv` is no longer in `Phylogenetic_Hierarchical_Orthogroups/`; `Orthogroups/Orthogroups.tsv` contains its orthogroups.

## Installation

### Install in conda (recommended) 

The simplest way to install OrthoFinder is through [conda](https://anaconda.org/bioconda/orthofinder). If you're unfamiliar with conda, [this tutorial](https://www.machinelearningplus.com/deployment/conda-create-environment-and-everything-you-need-to-know-to-manage-conda-virtual-environment/) offers a beginner-friendly introduction.

```bash
conda create -n of3_env python=3.12
conda activate of3_env
conda install orthofinder
```

Alternatively, you could install via github, or download the source code and install locally.

### Install via github

```bash
python3 -m venv of3_env 
. of3_env/bin/activate
pip install git+https://github.com/OrthoFinder/OrthoFinder.git
```

### Install locally from source code

The following commands provide three ways to download the source code of OrthoFinder locally into a directory named `OrthoFinder`.
```bash
# Download via git 
git clone https://github.com/OrthoFinder/OrthoFinder.git

# or download the orthofinder-linux-intel-3.2.tar.gz and unzip it into OrthoFinder if you are on a Linux Intel machine
mkdir OrthoFinder && \
  wget -qO- https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.2/orthofinder-linux-intel-3.2.tar.gz | \
  tar -xz --strip-components=1 -C OrthoFinder
```

Next, you can run the following commands to install OrthoFinder inside the of3_env virtural environment.
```bash
cd OrthoFinder
python3 -m venv of3_env # Create an virtural environment named of3_env
. of3_env/bin/activate # Activate of3_env
pip install .
```

Whether you've installed OrthoFinder directly from GitHub or downloaded and set it up locally, the OrthoFinder package will only be available within the `of3_env` virtual environment. This avoids potential conflicts with Python dependencies.

To deactivate the virtual environment when you are finished, run:
```bash
deactivate
```
To activate the virtual environment you have created, run:

```bash
. of3_env/bin/activate
```
### Test your installation

Once you have installed OrthoFinder, you can print the help information and version, and test it on the [example data](https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.2/ExampleData.tar.gz).

```bash
orthofinder --help # Print out help informatioin
orthofinder --version # Check the version
orthofinder -f ExampleData # Test OrthoFinder on an example dataset - this should take a few minutes to run. 
```

### Uninstalling

To uninstall on conda:
```bash
conda deactivate
conda remove -n of3_env --all
```
To remove the virtual environment where OrthoFinder is installed:
```bash
deactivate
cd ..
rm -rf OrthoFinder
```

## Simple Usage

Run OrthoFinder on FASTA format proteomes in `<dir>`

```python
orthofinder [options] -f <dir>
```

OrthoFinder requires one FASTA file for each species. Each file should contain the complete set of protein sequences from that species' genome,  with a single representative sequence for each gene.


If your proteomes contain multiple transcript isoforms per gene, we provide the script `primary_transcript.py` to retain only the longest isoform for each gene. This script should be run before running OrthoFinder.

If you are using OrthoFinder v3.2 or later, and installed OrthoFinder into your environment, the script is installed automatically and can be run directly as `primary_transcript`. It accepts:

- a single FASTA file
- a directory containing FASTA files
- a single .tar, .tar.gz, .tgz, or .zip archive
- a directory containing archive files (.tar, .tar.gz, .tgz, .zip)

For example:
```bash
primary_transcript proteomes/
```
or 
```bash
primary_transcript proteins.fa
```
For older OrthoFinder versions, run the script manually for each FASTA file:
```bash
for f in *.fa; do
    python primary_transcript.py "$f"
done
```

## Advanced Usage - Scaling to Thousands of Species

If you are analysing >100 species, we recommend that you use the scalable implementation. 

Add the files for 64 species into one directory `<core>`
Add the remaining files into another directory `<additional>`

First, run OrthoFinder on the subset of 64 species

```python
orthofinder [options] -f <core>
```

Then, add the additional species to the results of the core run

```python
orthofinder [options] --assign <additional> --core <Results_Dir>
```

To choose which 64 species to include in the core, aim to capture a broad range of the evolutionary diversity of your species.

<p align="center">
  <img src="assets/images/core_assign_figure.jpg" alt="core_assign_figure" width="800"/>
</p>

**Note that this alternative way of running OrthoFinder requires that the core species are run using the multiple sequence alignment option. You cannot add additional species to OrthoFinder results that were run with the `-M dendroblast` option, which was the default for OrthoFinder2**


## Command-line options

Command-line options for OrthoFinder

- **Adding additional species**<br>
  | Parameter | Description                               |
  |-----------|-------------------------------------------|
  | `--assign <dir1> --core <dir2>`      | Assign species from `<dir1>` to existing orthogroups in `<dir2>`.                   |

- **Method choices**<br>
  | Parameter | Description                               | Default   | Options                                                                                     |
  |-----------|-------------------------------------------|-----------|---------------------------------------------------------------------------------------------|
  | `-M`      | Method for gene tree inference.           | `msa`     | `dendroblast`, `msa`                                                                        |
  | `-S`      | Sequence search program                   | `diamond` | `diamond`, `diamond_ultra_sens`, `blastp`, `mmseqs`, `blastn` |
  | `-A`      | MSA program, requires `-M msa`            | `famsa`   | `famsa`, `mafft`, `muscle`,                                                          |
  | `-T`      | Tree inference method, requires `-M msa`  | `fasttree`| `fasttree`, `fasttree_fastest`, `raxml`, `iqtree`                               |
  | `-ST`     | Species tree program, from the gene trees, requires `-M msa` | a tree of the concatenated alignment | `astral-pro`; others can be added to `config.json` (`"program_type": "species_tree"`) |
  | `-I`      | MCL inflation parameter                   | `1.2`     | `1-10`                                                                                         |

- **Sequence search scoring options (DIAMOND)**<br>
  | Parameter        | Description                                                               |
  |------------------|---------------------------------------------------------------------------|
  | `--matrix`       | Scoring matrix of the DIAMOND searches (default `BLOSUM62`): `BLOSUM45`, `BLOSUM50`, `BLOSUM62`, `BLOSUM80`, `BLOSUM90`, `PAM250`, `PAM70` or `PAM30`. |
  | `--custom-matrix` | A scoring matrix file for the DIAMOND searches (see below); requires `--gapopen` and `--gapextend`. |
  | `-go`, `--gapopen` | Gap open penalty of the DIAMOND searches (default: the matrix's, `11` for `BLOSUM62`), as DIAMOND allows for the matrix. |
  | `-ge`, `--gapextend` | Gap extend penalty of the DIAMOND searches (default: the matrix's, `1` for `BLOSUM62`), as DIAMOND allows for the matrix. |

  **Custom scoring matrices (`--custom-matrix`)**: a text file in the format of NCBI's matrices, e.g. any file from https://ftp.ncbi.nlm.nih.gov/blast/matrices/ as it is:
  ```
  # comment lines start with #
     A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
  A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
  R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
  ...
  ```
  - The first line is the one-letter amino acid codes; then, for each code in that order, a line of the code and its integer score with each code of the first line.
  - The 20 standard amino acids are needed (DIAMOND would otherwise find no hits); others, such as `B`, `Z`, `X` and `*`, are optional.
  - DIAMOND has no default gap penalties for a custom matrix, so `--gapopen` and `--gapextend` are needed. Not every pair works with every matrix: if DIAMOND cannot compute its statistics for them ("The ALP library failed to compute the statistical parameters for this matrix"), try others.
  - The file can be anywhere: give its path with `--custom-matrix`. OrthoFinder checks it before the run starts, and keeps a copy with the results (`WorkingDirectory/custom_scoring_matrix.txt`) that the searches use. A restart from those results (e.g. `-b` or `--assign`) uses the same matrix without `--custom-matrix` being given again, even if the original file has since been moved or changed.

- **Input options**<br>
  | Parameter | Description                               |
  |-----------|-------------------------------------------|
  | `-d`      | Input is DNA sequences.                   |
  | `-s`      | User-specified rooted species tree.        |

- **Output options**<br>
  | Parameter   | Description                                                                 |
  |-----------  |-----------------------------------------------------------------------------|
  | `-X`      | Don’t add species names to sequence IDs.                                    |
  | `-n <txt>`      | Name to append to the results directory.                                    |
  | `-o <txt>`      | Specify a non-default results directory.                                    |
  | `-pof`    | Save the pairwise ortholog files.                                           |
  | `--orthoxml` | Also write the orthogroups and HOGs in OrthoXML format, for other tools. |
  | `--concat-alignment` | With `-ST`, also write the concatenated alignment of the orthogroups chosen for a species tree, as a run without `-ST` does (`WorkingDirectory/Alignments_ids/SpeciesTreeAlignment.fa`, `Species_Tree/Orthogroups_for_concatenated_alignment.txt`). For the results of a finished run: `species_tree_alignment <Results_dir>` (`-o <dir>` to write them elsewhere, `--species-names` to name the rows by species). |

- **Viewing the results in a web browser**<br>
  | Parameter | Description |
  |-----------|-------------|
  | `--tree-view` | Also write a page to view the species tree and the resolved gene trees in a web browser, with each gene tree's alignment, and a species tree explorer: the HOGs at each node, an orthogroup's duplications and losses on the species tree, each species' genes; and all genes on the species tree, each orthogroup's gene tree where it was gained (`Resolved_Gene_Trees/Resolved_Gene_Trees.html`). |
  | `--stats-view` | Also write a page of the run's statistics with interactive plots (`Comparative_Genomics_Statistics/Statistics.html`). |
  | `--groups-view` | Also write a page of the genes of each orthogroup and HOG, searchable (`Orthogroups/Groups.html`). |
  | `--view` | All three pages, linked to each other, and `OrthoFinder_Results.html` to open them, a tab for each (also with two of them). |

  - For the results of a finished run: `view_results <Results_dir>` (same options). With many species the pages can be large: each can also be written on its own.
  - The pages work in any recent web browser (Safari, Edge, Chrome, Firefox) on macOS, Windows or Linux, opened from disk. To look at results from a server on your own computer, download the whole results directory, as each page reads its data from the folders next to it.
  - On Linux, open them with `xdg-open <Results_dir>/OrthoFinder_Results.html` (the default browser) or `google-chrome <Results_dir>/OrthoFinder_Results.html 2>/dev/null` (with Chrome already open, it prints errors such as "Broken pipe" that are not about the page).
  - Light or dark colours: **Colours**, at the top right of each page (Auto follows the system).

- **Parallel processing options**<br>
  | Parameter | Description                                 | Default |
  |-----------|---------------------------------------------|---------|
  | `-t`      | Number of parallel sequence search threads. | `All available`    |
  | `-a`      | Number of parallel analysis threads.        | `16 or t/8 (whichever lower)`     |

- **Workflow stopping options**<br>
  | Parameter | Description                                                                 |
  |-----------|-----------------------------------------------------------------------------|
  | `-op`     | Stop after preparing input files for BLAST.                                 |
  | `-og`     | Stop after MCL Orthogroups.                                                 |
  | `-ogt`    | Stop after gene trees.                                                      |
  | `-ost`    | Stop after species trees.                                                   |

- **Workflow restart options**<br>
  | Parameter    | Description                                                   |
  |--------------|---------------------------------------------------------------|
  | `-b <dir>`   | Restart OrthoFinder from pre-computed sequence search results in `<dir>`. |
  | `-fg <dir>`  | Restart OrthoFinder from MCL Orthogroups results in `<dir>`.            |
  | `-fgt <dir>` | Restart OrthoFinder from gene trees results in `<dir>`.                 |
  | `-fst <dir>` | Restart OrthoFinder from species trees results in `<dir>`.              |

  > *`<dir>` refers to the OrthoFinder results directory*

- **Other options**<br>
  | Parameter        | Description                                                               |
  |------------------|---------------------------------------------------------------------------|
  | `-1`             | Only perform one-way sequence search.                                     |
  | `-z`             | Don’t trim MSAs (columns >= 90% gap, min. alignment length 500).          |
  | `-y`             | Split paralogous clades below the root of a HOG into separate HOGs.        |
  | `-h`             | Print this help text.                                                     |
  | `-v`             | Print version.                                                     |

## Output files

A standard OrthoFinder run produces a set of files describing the orthogroups, orthologs, gene trees, resolve gene trees, the rooted species tree, gene duplication events, and comparative genomic statistics for the set of species being analysed. These files are located in an intuitive directory structure.

Full details on the output files and directories can be found [here](http://orthofinder.github.io/OrthoFinder/tutorials/guide-to-results/). The directories that are useful for most users are

**```/Orthogroups```**
- `Orthogroups.tsv` is the main orthogroup file. Each row contains the genes belonging to a single orthogroup. The genes from each orthogroup are organized into columns, one per species. 
- `Orthogroups.txt` is a text file with each line showing the genes in a single orthogroup. It differs from Orthogroups.tsv in that it doesn’t show the species which each gene belongs to.
- `Orthogroups.GeneCount.tsv` is a tab separated text file that contains counts of the number of genes for each species in each orthogroup.
- `Orthogroups_SingleCopyOrthologs.txt` is a list of orthogroups that contain exactly one gene per species
- `Orthogrouops_UnassignedGenes.tsv` is a tab separated text file that contains all of the genes that were not assigned to any orthogroup.

**```/Phylogenetic_Hierarchical_Orthogroups```**
- Each file is a phylogenetic hierarchical orthogroup (HOG) for a different node of the species tree.
- Each row of a file contain the genes belonging to a single orthogroup.
- Each species is represented by a single column.
- `N0.tsv` from the old version is now `Orthogroups/Orthogroups.tsv`

**```/Orthologs```**
- `<species>.tsv.gz` for each species lists its orthologs in all the other species.
- With `-pof`, each species has a sub-directory, `Orthologs_<species>`, that in turn contains a file for each pairwise species comparison, listing the orthologs between that species pair (`split_ortholog_files` writes them for the results of a run without `-pof`).

**```/Comparative_Genomics_Statistics```**
- Files containing summary statistics across all orthogroups, as well as comparisons between each pair of species.

**```/Resolved_Gene_Trees```**
- A rooted phylogenetic tree inferred for each orthogroup with 4 or more sequences and resolved using the OrthoFinder hybrid species-overlap/duplication-loss coalescent model.

**```/Species_Tree```**
- `SpeciesTree_rooted.txt` is a species tree inferred using STAG or ASTRAL-Pro.
- `SpeciesTree_rooted_node_labels.txt` is the same tree, but with nodes labels instead of support values. This labelled version is useful for interpreting and analysing the results of the gene duplication analyses.

```/Gene_Duplication_Events```
- `Duplications.tsv` has a row for each gene duplication event, with information on orthogroup in which it occured, the species that contain the duplicated gene, the node in the species tree on which the gene duplication event occured, and the support score for the gene duplication event.
- `SpeciesTree_Gene_Duplications_0.5_Support.txt` provides a summation of the above duplications over the branches of the species tree.

**```/Orthogroup_Sequences```**
- A FASTA file for each orthogroup giving the amino acid sequences for each gene in the orthogroup.


## Latest additions
The current version of OrthoFinder has several major changes compared to OrthoFinder version 2 ([Emms & Kelly 2019]((https://img.shields.io/badge/DOI-10.1186%2Fs13059--019--1832--y-blue))).

**New workflow for scalability**

The ``--core --assign`` workflow uses the [SHOOT](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02652-8) algorithm to create profiles for previously computed orthogroups, and adds new genes to these orthogroups without requiring a costly all-versus-all sequence search. Genes that cannot be assigned using the SHOOT approach are analysed using a standard OrthoFinder workflow. 

**Phylogenetic Hierarchical Orthogroups**

OrthoFinder has now extended its phylogenetic approach to orthogroups, allowing orthogroups to be defined for each node within the species tree. This significantly increases the accuracy of orthogroups, and enables users to perform orthogroup analyses for any clade of species in the species tree. 

<p align="center">
  <img src="assets/images/HOG_figure.png" alt="HOGs" width="400"/>
</p>

## Citation

- Latest<br>
  [1] *David M. Emms, Yi Liu, Laurence Belcher, Jonathan Holmes, Steven Kelly, 2026.* **OrthoFinder: improved phylogenetic orthology inference with enhanced accuracy and scalability.** *Nature Methods* **23**, 1327–1333. [![DOI:10.1038/s41592-026-03126-6](https://img.shields.io/badge/DOI-10.1038%2Fs41592--026--03126--6-blue)](https://doi.org/10.1038/s41592-026-03126-6)

- Introduced the SHOOT method to perform phylogenetic gene search<br>
  [2] *Emms, D.M., Kelly, S*. **SHOOT: phylogenetic gene search and ortholog inference**. Genome Biol 23, 85 (2022). [![DOI:10.1186/s13059-022-02652-8](https://img.shields.io/badge/DOI-10.1186%2Fs13059--022--02652--8-blue)](https://doi.org/10.1186/s13059-022-02652-8)

- Introduced the phylogenetic inference of orthologs, including rooted gene and species trees, and gene duplication events<br>
  [3] *Emms, D.M., Kelly, S*. **OrthoFinder: phylogenetic orthology inference for comparative genomics**. Genome Biol 20, 238 (2019). [![DOI:10.1186/s13059-019-1832-y](https://img.shields.io/badge/DOI-10.1186%2Fs13059--019--1832--y-blue)](https://doi.org/10.1186/s13059-019-1832-y)

- Introduced the STRIDE method to root an unrooted species tree.<br>
  [4] *Emms DM, Kelly S*. **STRIDE: Species Tree Root Inference from Gene Duplication Events**. Mol Biol Evol. 2017 Dec 1;34(12):3267-3278. [![DOI:10.1093/molbev/msx259](https://img.shields.io/badge/DOI-10.1093%2Fmolbev%2Fmsx259-blue)](https://doi.org/10.1093/molbev/msx259)

- Introduced the STAG method of species tree inference<br>
  [5] *D.M. Emms, S. Kelly, 2017*. **STAG: Species Tree Inference from All Genes** bioRxiv. [![DOI:10.1101/267914](https://img.shields.io/badge/DOI-10.1101%2F267914-blue)](https://doi.org/10.1101/267914)

- Introduced the orthogroup inference method<br>
  [6] *Emms, D.M., Kelly, S*. **OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy**. Genome Biol 16, 157 (2015). [![DOT:10.1186/s13059-015-0721-2](https://img.shields.io/badge/DOI-10.1186%2Fs13059--015--0723--2-blue)](https://doi.org/10.1186/s13059-015-0721-2)

## System Requirements 

**Operating system** 

OrthoFinder was designed to run on Linux (including WSL2). 

We have tested OrthoFinder v3.1 on debian 12.9, centOS v8, macOS 14.4.1, macOS 13.2.1. 

**Dependencies**
- Python `>=3.11`
- Diamond `>=2.1.7,<2.2`
- Famsa `>=2.2.3`
- Fasttree `>=2.1.11`
- Numpy `>=2.3.2`
- Scipy `>=1.16`
- Biopython `>=1.85`
- Rich `>=14.1.0`
- Scikit-learn `>=1.7.1`

## Meet the team

OrthoFinder was developed by David Emms & Steve Kelly

Current members of the OrthoFinder team:

Yi Liu, Jonathan Holmes, Laurie Belcher

