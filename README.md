# OrthoFinder
<p align="center">
  <img src="assets/images/coreassignfigure-images-2.JPG" alt="concept_figure" width="800"/>
</p>
OrthoFinder identifies orthogroups, infers gene trees for all orthogroups, and analyzes the gene trees to identify the rooted species tree. The method subsequently identifies all gene duplication events in the complete set of gene trees, and analyses them at both gene tree and species tree level. OrthoFinder further analyzes all of this phylogenetic information to identify the complete set of orthologs between all species, and provide extensive comparative genomics statistics.

## Table of contents
- [Installation](#installation)
- [Simple Usage](#simple-usage)
- [Advanced Usage - Scaling to Thousands of Species](#advanced-usage---scaling-to-thousands-of-species)
- [Command line Options](#command-line-options)
- [Output files](#output-files)
- [Latest additions](#latest-additions)
- [Citation](#citation)

A single PDF with all documentation and tutorial is available [here](https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/OrthoFinder3_tutorial.pdf). For more information please visit [our website](https://orthofinder.github.io/OrthoFinder/).

## Installation

### Install in conda (recommended) 

The simplest way to install OrthoFinder is through [conda](https://anaconda.org/bioconda/orthofinder). If you're unfamiliar with conda, [this tutorial](https://www.machinelearningplus.com/deployment/conda-create-environment-and-everything-you-need-to-know-to-manage-conda-virtual-environment/) offers a beginner-friendly introduction.

```bash
conda create -n of3_env python=3.10
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

# or download the orthofinder-linux-intel-3.1.0.tar.gz and unzip it into OrthoFinder if you are on a Linux Intel machine
mkdir OrthoFinder && wget -qO- https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/orthofinder-linux-intel-3.1.0.tar.gz | tar -xz --strip-components=1 -C OrthoFinder

# or download the orthofinder-macos-intel-3.1.0.tar.gz and unzip it into OrthoFinder a MacOS Intel machine
mkdir OrthoFinder && wget -qO- https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/orthofinder-macos-intel-3.1.0.tar.gz | tar -xz --strip-components=1 -C OrthoFinder
```

Next, you can run the following commands to install OrthoFinder inside the of3_env virtural environment.
```bash
cd OrthoFinder
python3 -m venv of3_env && . of3_env/bin/activate
pip install .
```

Whether you've installed OrthoFinder directly from GitHub or downloaded and set it up locally, the OrthoFinder package will only be available within the `of3_env` virtual environment. This avoids potential conflicts with Python dependencies.

To deactivate the virtual environment when you are finished, run:
```bash
deactivate
```

### Test your installation

Once you have installed OrthoFinder, you can print the help information and version, and test it on the [example data](https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/ExampleData.tar.gz).

```bash
orthofinder --help # Print out help informatioin
orthofinder --version # Check the version
orthofinder -f ExampleData # Test OrthoFinder on an example dataset
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

If your files have multiple transcript variants for each gene, then we provide a script `primary_transcripts.py` to extract the longest variant per gene. This script should be run on your files prior to running OrthoFinder;

```bash
for f in *fa ; do python primary_transcript.py $f ; done
```

## Advanced Usage - Scaling to Thousands of Species

If you are analysing >100 species, we reccommend that you use the scalable implementation. 

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

**Adding additional species**
| Parameter | Description                               |
|-----------|-------------------------------------------|
| `--assign <dir1> --core <dir2>`      | Assign species from `<dir1>` to existing orthogroups in `<dir2>`.                   |

**Method choices**
| Parameter | Description                               | Default   | Options                                                                                     |
|-----------|-------------------------------------------|-----------|---------------------------------------------------------------------------------------------|
| `-M`      | Method for gene tree inference.           | `msa`     | `dendroblast`, `msa`                                                                        |
| `-S`      | Sequence search program                   | `diamond` | `blast`, `diamond`, `diamond_ultra_sens`, `blastp`, `mmseqs`, `blastn` |
| `-A`      | MSA program, requires `-M msa`            | `famsa`   | `famsa`, `mafft`, `muscle`,                                                          |
| `-T`      | Tree inference method, requires `-M msa`  | `fasttree`| `fasttree`, `fasttree_fastest`, `raxml`, `iqtree`                               |
| `-I`      | MCL inflation parameter                   | `1.2`     | `1-10`                                                                                         |

**Input options**
| Parameter | Description                               |
|-----------|-------------------------------------------|
| `-d`      | Input is DNA sequences.                   |
| `-s`      | User-specified rooted species tree.        |

**Output options**
| Parameter   | Description                                                                 |
|-----------  |-----------------------------------------------------------------------------|
| `-X`      | Don’t add species names to sequence IDs.                                    |
| `-n <txt>`      | Name to append to the results directory.                                    |
| `-o <txt>`      | Specify a non-default results directory.                                    |

**Parallel processing options**
| Parameter | Description                                 | Default |
|-----------|---------------------------------------------|---------|
| `-t`      | Number of parallel sequence search threads. | `All available`    |
| `-a`      | Number of parallel analysis threads.        | `16 or t/8 (whichever lower)`     |

**Workflow stopping options**
| Parameter | Description                                                                 |
|-----------|-----------------------------------------------------------------------------|
| `-op`     | Stop after preparing input files for BLAST.                                 |

**Workflow restart options**
| Parameter  | Description                                                  |
|------------|--------------------------------------------------------------|
| `-b <dir>` | Start OrthoFinder from pre-computed BLAST results in `<dir>`. |

**Other options**
| Parameter        | Description                                                               |
|------------------|---------------------------------------------------------------------------|
| `-1`             | Only perform one-way sequence search.                                     |
| `-z`             | Don’t trim MSAs (columns >= 90% gap, min. alignment length 500).          |
| `-y`             | Split paralogous clades below the root of a HOG into separate HOGs.        |
| `-h`             | Print this help text.                                                     |
| `-v`             | Print version.                                                     |


## Output files

A standard OrthoFinder run produces a set of files describing the orthogroups, orthologs, gene trees, resolve gene trees, the rooted species tree, gene duplication events, and comparative genomic statistics for the set of species being analysed. These files are located in an intuitive directory structure.

Full details on the output files and directories can be found in the PDF [here](link). The directories that are useful for most users are;

```/Orthogroups```
- Orthogroups.tsv is the main orthogroup file. Each row contains the genes belonging to a single orthogroup. The genes from each orthogroup are organized into columns, one per species. 
- Orthogroups.txt is a text file with each line showing the genes in a single orthogroup. It differs from Orthogroups.tsv in that it doesn’t show the species which each gene belongs to.
- Orthogroups.GeneCount.tsv is a tab separated text file that contains counts of the number of genes for each species in each orthogroup.
- Orthogroups_SingleCopyOrthologues.txt is a list of orthogroups that contain exactly one gene per species
- Orthogrouops_UnassignedGenes.tsv is a tab separated text file that contains all of the genes that were not assigned to any orthogroup.

```/Phylogenetic_Hierarchical_Orthogroups```
- Each file is a phylogenetic hierarchical orthogroup (HOG) for a different node of the species tree
- Each row of a file contain the genes belonging to a single orthogroup
- Each species is represented by a single column

```/Orthologues```
- Each species has a sub-directory that in turn contains a file for each pairwise species comparison, listing the orthologs between that species pair.

```/Comparative_Genomics_Statistics```
- Files containing summary statistics across all orthogroups, as well as comparisons between each pair of species

```/Resolved_Gene_Trees```
- A rooted phylogenetic tree inferred for each orthogroup with 4 or more sequences and resolved using the OrthoFinder hybrid species-overlap/duplication-loss coalescent model.

```/Species_Tree```
- SpeciesTree_rooted.txt = A species tree inferred using ASTRAL-Pro.
- SpeciesTree_rooted_node_labels.txt = The same tree, but with nodes labels instead of support values. This labelled version is useful for interpreting and analysing the results of the gene duplication analyses.

```/Gene_Duplication_Events```
- `Duplications.tsv` has a row for each gene duplication event, with information on orthogroup in which it occured, the species that contain the duplicated gene, the node in the species tree on which the gene duplication event occured, and the support score for the gene duplication event.
- `SpeciesTree_Gene_Duplications_0.5_Support.txt` provides a summation of the above duplications over the branches of the species tree.

```/Orthogroup_Sequences```
- A FASTA file for each orthogroup giving the amino acid sequences for each gene in the orthogroup.


## Latest additions
The current version of OrthoFinder has several major changes comapred to OrthoFinder version 2 (Emms & Kelly 2019)

**New workflow for scalability**

The ``--core --assign`` workflow uses [SHOOT](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02652-8) to create profiles for previously computed orthogroups, and adds new genes to these orthogroups without requiring a costly all-versus-all sequence search. Genes that cannot be assigned using the SHOOT approach are analysed using a standard OrthoFinder workflow. 

**Phylogenetic Hierarchical Orthogroups**

OrthoFinder has now extended its phylogenetic approach to orthogroups, allowing orthogroups to be defined for each node within the species tree. This significantly increases the accuracy of orthogroups, and enables users to perform orthogroup analyses for any clade of species in the species tree. 

<p align="center">
  <img src="assets/images/hog.png" alt="HOGs" width="300"/>
</p>

## Citation

- Latest<br>
  [1] *David M Emms, Yi Liu, Laurence Belcher, Jonathan Holmes, Steven Kelly, 2025.* **OrthoFinder: scalable phylogenetic orthology inference for comparative genomics** bioRxiv. 

- Introduced the orthogroup inference method<br>
  [2] *Emms, D.M., Kelly, S*. **OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy**. Genome Biol 16, 157 (2015). [![DOT:10.1186/s13059-015-0721-2](https://img.shields.io/badge/DOI-10.1186%2Fs13059--015--0723--2-blue)](https://doi.org/10.1186/s13059-015-0721-2)

- Introduced the phylogenetic inference of orthologs, including rooted gene and species trees, and gene duplication events<br>
  [3] *Emms, D.M., Kelly, S*. **OrthoFinder: phylogenetic orthology inference for comparative genomics**. Genome Biol 20, 238 (2019).  
  [![DOI:10.1186/s13059-019-1832-y](https://img.shields.io/badge/DOI-10.1186%2Fs13059--019--1832--y-blue)](https://doi.org/10.1186/s13059-019-1832-y)

- Introduced the STRIDE method to root an unrooted species tree.<br>
  [4] *Emms DM, Kelly S*. **STRIDE: Species Tree Root Inference from Gene Duplication Events**. Mol Biol Evol. 2017 Dec 1;34(12):3267-3278.  
  [![DOI:10.1093/molbev/msx259](https://img.shields.io/badge/DOI-10.1093%2Fmolbev%2Fmsx259-blue)](https://doi.org/10.1093/molbev/msx259)

- Introduced the STAG method of species tree inference<br>
  [5] *D.M. Emms, S. Kelly, 2017*. **STAG: Species Tree Inference from All Genes** bioRxiv.  
  [![DOI:10.1101/267914](https://img.shields.io/badge/DOI-10.1101%2F267914-blue)](https://doi.org/10.1101/267914)



## Meet the team

OrthoFinder was developed by David Emms & Steve Kelly

Current members of the OrthoFinder team:

Yi Liu, Jonathan Holmes, Laurie Belcher

