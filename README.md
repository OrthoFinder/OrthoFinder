# OrthoFinder3

<img src="assets/OF2_Workflow.png" alt="OF2_Workflow" width="750"/>

OrthoFinder identifies orthogroups, infers gene trees for all orthogroups, and analyzes these gene trees to identify the rooted species tree. The method subsequently identifies all gene duplication events in the complete set of gene trees, and analyzes this information in the context of the species tree to provide both gene tree and species tree-level analysis of gene duplication events. OrthoFinder further analyzes all of this phylogenetic information to identify the complete set of orthologs between all species and provide extensive comparative genomics statistics.

## Table of contents
- [Installation](#Installation)
- [Simple Usage](#Simple-Usage)
- [Advanced Usage](#Advanced-Usage)
- [Command line Options](#Command-line-Options)
- [Output files](#Output-files)
- [Latest additions](#Latest-additions)
- [Citation](#Citation)

Tutorials and further documentation can be found on our [website](https://sites.google.com/view/orthofinder3)

A single PDF with all documentation, tutorials, and this README is available [here](https://github.com/lauriebelch/OrthoDocs/blob/main/of3_v2.pdf)

## Installation

OrthoFinder3 provides a user-freindly way to use. It can be run with or without installation. Either way, the first step requires a user to download a copy of the package to their local machine. Whichever way a user choose to run OrthoFinder3, it is recommended that the software is ran inside a virtual environment or a name space in conda. 

> Run OrthoFinder3 with installation

When running OrthoFinder3 with installation, we recommend the user to run it inside a conda environment and miniconda is suffice. The installation instructions for miniconda can be found [here](https://www.anaconda.com/docs/getting-started/miniconda/install#macos-linux-installation).

- Install OrthoFinder3 in conda

  OrthoFinder3 is a listed software in bioconda. More information can be found [here](https://www.machinelearningplus.com/deployment/conda-create-environment-and-everything-you-need-to-know-to-manage-conda-virtual-environment/).

  ```bash
    conda create -n of3_env python=3.10
    conda activate of3_env
    conda install orthofinder
  ```
- [Install OrthoFinder3 via src code](#install-via-src-code)

  Generally speaking, more steps is evolved in installing OrthoFinder3 via the src code. In the meantime, more OS related issues can occur when running orthofinder via src code. The main issue is the binaries shipped with OrthoFinder3 might not be comptable with your underlying OS. Therefore, this is not a recommended way to run OrthoFinder3. Nevertheless, for those Linux users who prefer running OrthoFinder3 in such way, you can running the following commands to install the software and run it on your local machine.

  - Download via `git`
      
    ```bash
      git clone https://github.com/OrthoFinder/OrthoFinder.git
      cd OrthoFinder
    ```

  - Download zipfile or tar.gz file
    ```bash
      wget https://github.com/OrthoFinder/OrthoFinder/archive/refs/heads/master.zip
      ## link need to be updated later
      cd OrthoFinder
    ```

  - Install OrthoFinder 
    ```bash
      python3 -m venv of3_env && . of3_env/bin/activate
      pip install .
    ```
  Once you have installed OrthoFinder3, you can run the following commands to check the version, print out the help informatioin or test if it has been installed sucessfully.

    ```bash
      python3 -m venv of3_env && . of3_env/bin/activate
      pip install .
    ```

> Run OrthoFinder3 without installation

Running OrthoFinder3 without installation also requires the user to download a copy of the src code from GitHub. Flease follow the instruction in the <a name="install-via-src-code">Install OrthoFinder3 via src code</a> section. 

Once you have package src code available you can running the following commands to installed the required Python depencies inside a virtual environment. 

  ```bash
    python3 -m venv of3_env && . of3_env/bin/activate
    pip install -r requirements.txt
  ```

Having done so, you then test and run OrthoFinder3 without actaully install the software. To test it, please run 

  ```bash
    python3 orthofinder.py --version # Check the version
    python3 orthofinder.py --help # Print out help informatioin
    python3 orthofinder.py -f ExampleData # Test the OrthoFinder on an example dataset shipped with OrthoFinder3
  ```




    





The easiest way to install OrthoFinder3 is using [conda](https://www.machinelearningplus.com/deployment/conda-create-environment-and-everything-you-need-to-know-to-manage-conda-virtual-environment/).
```bash
conda install orthofinder
orthofinder -h
```

If you are on a mac that has an M1/M2/M3 chip, you might have to adjust your conda architecture. Instructions can be found [here](https://towardsdatascience.com/how-to-manage-conda-environments-on-an-apple-silicon-m1-mac-1e29cb3bad12).

Alternatively, you can also download OrthoFinder3 directly from GitHub. 

Via GitHub, you have the following options to download OrthoFinder to your local machine.
- via `git`
    
  ```bash
  git clone https://github.com/OrthoFinder/OrthoFinder.git
  python OrthoFinder/orthofinder -h
  ```

A docker image will also available [here](###linktodocker###)

**Installing dependencies**

Add some info on how to manually install dependencies

## Simple Usage

OrthoFinder requires one FASTA format file for each species being analysed. Each file should contain the complete set of protein sequences encoded by genes present in that species genome,  with a single representative protein sequence for each gene.

OrthoFinder does not have the ability to distinguish proteins derived from transcript variants from proteins derived from separate genes. If your files have protein sequences derived from multiple transcript variants for each gene then we provide a script `primary_transcripts.py` to extract the longest variant per gene that should be run on your input files prior to running OrthoFinder;

```bash
for f in *fa ; do python primary_transcript.py $f ; done
```

Run OrthoFinder3 on FASTA format proteomes in `<dir>`

```python
orthofinder [options] -f <dir> -M dendroblast
```


## Advanced Usage

If you have large numbers of species to analyse OrthoFinder has a scalable implementation that can be used to add new species to an already inferred set of orthogroups for a smaller, core group of species.

<!-- ![OrthoFinder3 workflow](of3.png) -->
<img src="assets/OF3_Workflow.png" alt="OrthoFinder3 workflow" width="750"/>

We reccomend that if you are analysing >32 species, you use the this core/assign method. You will need to partition your input dataset into 16 “core” and the remainder “additional” species. 

The core set of species should be the subset that maximally captures the evolutionary uniqueness of the species in the input dataset.

You need to make a folder `core` with the core species, and a folder `additional`  with the additional species

You can then run OrthoFinder3 on the core species

```python
orthofinder [options] -f <dir_core>
```

Next, you can add the additional species

```python
orthofinder [options] --assign <dir_additional> --core <dir_core>
```

**Note that this alternative way of running OrthoFinder requires that the core species set is run using the multiple sequence alignment option. You cannot add additional species to OrthoFinder results that were run with the `-M dendroblast` option, which was the default for OrthoFinder2**

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
| `-S`      | Sequence search program                   | `diamond` | `blast`, `diamond`, `diamond_ultra_sens`, `diamond_custom`, `diamond_ultra_sens_custom`, `blast_gz`, `mmseqs`, `blast_nucl` |
| `-A`      | MSA program, requires `-M msa`            | `famsa`   | `famsa`, `mafft`, `muscle`, `mafft_memsave`                                                          |
| `-T`      | Tree inference method, requires `-M msa`  | `fasttree`| `fasttree`, `fasttree_fastest`, `raxml`, `raxml-ng`, `iqtree`                               |
| `-I`      | MCL inflation parameter                   | `1.2`     | `1-10`                                                                                         |

**Input options**
| Parameter | Description                               |
|-----------|-------------------------------------------|
| `-d`      | Input is DNA sequences.                   |
| `-s`      | User-specified rooted species tree.        |

**Output options**
| Parameter   | Description                                                                 |
|-----------  |-----------------------------------------------------------------------------|
| `-x <file>`      | Info for outputting results in OrthoXML format.                             |
| `-p <dir>`      | Write the temporary pickle files to `<dir>`.                                |
| `-X`      | Don’t add species names to sequence IDs.                                    |
| `-n <txt>`      | Name to append to the results directory.                                    |
| `-o <txt>`      | Specify a non-default results directory.                                    |
| `-efn`    | Extend the output directory name with the name of the scoring matrix, gap penalties, search program, MSA program, and tree program. |

**Parallel processing options**
| Parameter | Description                                 | Default |
|-----------|---------------------------------------------|---------|
| `-t`      | Number of parallel sequence search threads. | `11`    |
| `-a`      | Number of parallel analysis threads.        | `1`     |

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
| `--matrix`       | Scoring matrix allowed by DIAMOND.                                        |
| `--custom-matrix`| Custom scoring matrix.                                                    |
| `-z`             | Don’t trim MSAs (columns >= 90% gap, min. alignment length 500).          |
| `--save-space`   | Only create one compressed orthologs file per species.                    |
| `-y`             | Split paralogous clades below the root of a HOG into separate HOGs.        |
| `-h`             | Print this help text.                                                     |
| `-v`             | Print version.                                                     |


## Output files

A standard OrthoFinder run produces a set of files describing the orthogroups, orthologs, gene trees, resolve gene trees, the rooted species tree, gene duplication events, and comparative genomic statistics for the set of species being analysed. These files are located in an intuitive directory structure.

Full details on the output files and directories can be found [here](https://sites.google.com/view/orthofinder3/guide-to-results-files). The directories that are useful for most users are;

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

OrthoFinder now provides a ``--core --assign`` workflow to assign additional species to a previously computed OrthoFinder run. This workflow makes use of SHOOT to create profiles for the previously computed orthogroups, and new genes are assigned to these orthogroups without requiring a costly all-versus-all sequence search. Genes that cannot be assigned using the [SHOOT](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02652-8) approach are analysed on a clade-by-clade basis using a standard OrthoFinder workflow. 

**Phylogenetic Hierarchical Orthogroups**

OrthoFinder3 uses phylogenetic trees to determine orthologs. This is in contrast to most other methods that use only sequence similarity to infer orthologs. We have now extended our phylogenetic approach to orthogroups. This means that the orthogroups are more accurate than before, and also that orthogroups are now provided for each node within the species tree.

This change significantly increases the accuracy of the orthogroups identified by OrthoFinder. It also provides new functionality by enabling users to perform orthogroup-level analyses for any clade of species in the species tree. Orthogroups for every internal node in the species tree are provided in `/Comparative_Genomics_Statistics` (e.g. `N3.tsv` is the orthogroups for node N3 in the species tree).

<img src="assets/hog.png" alt="HOGs" width="500"/>

**Performance improvements**

(4x quicker runtime, 2.5x lower RAM usage, 15% more accurate orthogroups)


## Citation

The manuscript "OrthoFinder is the best" is now published in *Nature Methods*
[link here](https://uk.linkedin.com/in/jonathan-holmes-ab2794294).

[Emms & Kelly (2015)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0721-2]) introduced the orthogroup inference method.

[Emms & Kelly (2019)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y) introduced the phylogenetic inference of orthologs, including rooted gene and species trees, and gene duplication events.

[Emms & Kelly (2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5850722/) introduced the STRIDE method to root an unrooted species tree.

[Emms & Kelly (2017)](https://www.biorxiv.org/content/10.1101/267914v1) introduced the STAG method of species tree inference.

## Meet the team

OrthoFinder was developed by David Emms & Steve Kelly

Current members of the OrthoFinder team:

Yi Liu, Jonathan Holmes, Laurie Belcher


