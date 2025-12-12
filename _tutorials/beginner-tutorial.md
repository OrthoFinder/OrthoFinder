---
layout: post
title:  "Beginner Tutorial"
date:   2025-05-01 09:01:16 +0000
updated: 2025-05-08
docs: /assets/docs/beginner-tutorial.pdf
pdfs: https://GitHub.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/beginner-tutorial.pdf
---


This tutorial will cover:

- [Downloading OrthoFinder](#downloading-orthofinder)
- [Running OrthoFinder](#running-orthofinder)
- [Exploring the results of OrthoFinder](#exploring-the-results-of-orthofinder)

OrthoFinder requires as input the amino acid sequences for all the protein coding genes in
your species of interest. We provide a separate tutorial for [getting input files]({ site.baseurl}/tutorials/gettting-input-data/) for OrthoFinder.

All these steps will be done on the command line so that you can just copy and
paste the commands yourself. If you are not familiar with the command line there
are many online tutorials and reference pages, [here](https://www.techspot.com/guides/835-linux-command-line-basics/) is a nice short one that covers
the basics.

### Downloading OrthoFinder

There are two main ways of getting OrthoFinder. You can either use conda, or you can
install it directly from GitHub. Installing directly from GitHub will always give you the latest
version, but you might have to manually install other software that OrthoFinder is
dependent on, and it can be trickier to troubleshoot if you aren’t familiar with the
command line. Conda automates the installation process and handles all dependencies,
making it very beginner-friendly.
To install via conda, we first need to install miniconda. Follow the instructions [here](https://docs.anaconda.com/miniconda/).

We then need to run these commands

```
conda config --add channels defaults conda config --add channels bioconda conda config
--add channels conda-forge
conda create -n orthofinder
conda activate orthofinder
conda install orthofinder
```

If you are on one of the newer Macs with the new chips (M1/M2/M3), you might need to
follow a few extra steps to use conda. Check out the guide [here](https://towardsdatascience.com/how-to-manage-conda-environments-on-an-applesilicon-
m 1-mac-1e29cb3bad12).

To install directly from GitHub, we need to run these commands
```
python3 -m venv of3_env
. of3_env/bin/activate
pip install git+https://GitHub.com/OrthoFinder/OrthoFinder.git
```

You can test that OrthoFinder has been installed by printing its help file
```orthofinder -h```, which will print all of the command line options.

You can test that OrthoFinder is working correctly by running it on the example
dataset, which you can download from our [GitHub](https://GitHub.com/OrthoFinder/OrthoFinder/)

```orthofinder -f ExampleData/ ```

OrthoFinder will print lots of information to the command line as it runs. If you get an
error message, the best way to troubleshoot is to just google the error message. You
can also ask a question on our [GitHub](https://GitHub.com/OrthoFinder/OrthoFinder/issues).

When OrthoFinder has finished running, it will generate a folder containing the output,
with the folder named according to today’s date, for example, `ExampleData/OrthoFinder/Results_Dec08`. We’ll discuss how to interpret and analyse these files and folders later on, in the
[Exploring the results](#exploring-the-results-of-orthofinder) section of the tutorial.

### Running OrthoFinder

You can now run OrthoFinder
First, you have to open a terminal and navigate to the directory where your files
are. You can now run OrthoFinder on your proteomes.

```bash
orthofinder -f primary_transcripts 
```

That’s it! OrthoFinder will print updates on its progress to the terminal, and tell you when
it’s finished. To see what options you might want to adjust for your own data, check out the [GitHub](https://GitHub.com/OrthoFinder/OrthoFinder/), or the [Advanced Tutorial]({ site.baseurl }/tutorials/advanced-tutorial/) page

### Exploring the results of OrthoFinder

OrthoFinder creates a results directory named OrthoFinder inside the proteome
directory, and puts the results here. 

Results_Dec12/
   ├── Citation.txt
   ├── Comparative_Genomics_Statistics/
   │   ├── Duplications_per_Orthogroup.tsv
   │   ├── OrthologuesStats_many-to-many.tsv
   │   ├── OrthologuesStats_many-to-one.tsv
   │   ├── OrthologuesStats_one-to-many.tsv
   │   ├── OrthologuesStats_one-to-one.tsv
   │   ├── OrthologuesStats_Totals.tsv
   │   ├── Statistics_Overall.tsv
   ├── Gene_Duplication_Events/
   │   ├── Duplications.tsv
   ├── Log.txt
   ├── MultipleSequenceAlignments/
   │   ├── OG0000000.fa
   │   ├── OG000####.fa
   │   └── OG0000598.fa
   ├── Orthogroup_Sequences/
   │   ├── OG0000000.fa
   │   ├── OG000####.fa
   │   └── OG0001116.fa
   ├── Orthogroups/
   │   ├── Orthogroups.GeneCount.tsv
   │   ├── Orthogroups.tsv
   │   ├── Orthogroups.txt
   │   ├── Orthogroups_SingleCopyOrthologues.txt
   │   └── Orthogroups_UnassignedGenes.tsv
   ├── Orthologues/
   │   ├── Species_1.tsv
   │   ├── Species_#.tsv
   ├── Phylogenetic_Hierarchical_Orthogroups/
   │   └── N0.tsv
   │   └── N#.tsv
   ├── Phylogenetically_Misplaced_Genes/
   │   ├── Species_1.tsv
   │   ├── Species_#.tsv
   ├── Putative_Xenologs/
   │   ├── Species_1.tsv
   │   ├── Species_#.tsv
   ├── Resolved_Gene_Trees/
   │   └── Resolved_Gene_Trees.txt
   ├── Single_Copy_Orthologue_Sequences/
   │   ├── OG0000000.fa
   │   ├── OG000####.fa
   │   └── OG0000325.fa
   │   ├── Orthogroups_for_concatenated_alignment.txt
   └── WorkingDirectory/
       ├── Alignments_ids/
       ├── Blast0_0.txt.gz
       ├── Blast#_#.txt.gz
       ├── clusters_OrthoFinder_I1.2.txt
       ├── clusters_OrthoFinder_I1.2.txt_id_pairs.txt
       ├── dependencies/
       │   ├── SimpleTest.fa
       │   ├── Species0.fa
       │   ├── Species0.fa.output.txt
       │   └── test_search_results.txt.gz
       ├── N0.ids.tsv
       ├── OrthoFinder_graph.txt
       ├── pickle/
       ├── SequenceIDs.txt
       ├── Sequences_ids/
       │   ├── OG0000000.fa
       │   ├── OG000####.fa
       │   └── OG0001087.fa
       ├── Species0.fa
       ├── Species#.fa
       └── Trees_ids/
           ├── OG0000000.txt
           ├── OG000####.txt
           └── OG0000407.txt

#### Step 1: Quality Control

Before we start diving into the orthogroups, it would behoove us to check the quality
of the OrthoFinder run. We want to make sure that most genes across all species
have been assigned to orthogroups, and that the species tree looks realistic.

Open the file Statistics_Overall.tsv from the folder `Comparative_Genomics_Statistics`. This file can be opened in spreadsheet software
like Microsoft Excel, or in a text editor like Notepad.

On the 5th line, we can see the `Percentage of genes in orthogroups`, which in my case
is `81.0
`.

|------------------------------------|------|
| Number of species                  |    4 |
| Number of genes                    | 2733 |
| Number of genes in orthogroups     | 2215 |
| Number of unassigned genes         |  518 |
| Percentage of genes in orthogroups |   81 |
| Percentage of unassigned genes     |   19 |
| Number of orthogroups              |  599 |

A good rule of thumb is that this number should be `>80%`. If not, you are likely missing
some orthology relationships that actually exist. The best way to fix this would be better
species sampling.

Now open the file `Statistics_PerSpecies`, from the same folder. This file gives us the
`%` of genes in each species that are assigned to orthogroups, rather than the
percentage for all genes across species.

You can see here that we capture most genes across all species.

|                                    |   Mycoplasma_agalactiae |   Mycoplasma_gallisepticum |   Mycoplasma_genitalium |   Mycoplasma_hyopneumoniae |
|------------------------------------|-------------------------|----------------------------|-------------------------|----------------------------|
| Number of genes                    |                   820   |                      763   |                   476   |                      674   |
| Number of genes in orthogroups     |                   650   |                      596   |                   417   |                      552   |
| Number of unassigned genes         |                   170   |                      167   |                    59   |                      122   |
| Percentage of genes in orthogroups |                    79.3 |                       78.1 |                    87.6 |                       81.9 |
| Percentage of unassigned genes     |                    20.7 |                       21.9 |                    12.4 |                       18.1 |

The lowest percentage is the *`lowest_percent_species`*, but we still managed to assign `78.1` of its
genes to orthogroups. The key message here is that it’s always a good idea to look at
this information before you start interpreting your results. If the numbers were too low for
one species, we might want to consider sampling more species to fill in the long
evolutionary divergence between species (e.g. something in between a Kiwi and a
Kakapo, such as a Hoatzin).

One more useful thing to do before we really start to dive in is to look at the species
tree. Go to the [iTOL website](https://itol.embl.de/), and click `Upload a tree`.
You can then drag and drop the tree file, which is in `Species_Tree/SpeciesTree_rooted.txt`
You will now see the phylogenetic tree that OrthoFinder has produced. I have annotated
my version with icons PhyloPic, so that we can see what is going on

 ╭─┬╴Mycoplasma_hyopneumoniae
─┤ ╰╴Mycoplasma_agalactiae
 ╰─┬╴Mycoplasma_gallisepticum
   ╰╴Mycoplasma_genitalium

We now want to do some common-sense checking that everything appears to be in
order, and we aren’t rewriting the history of life on earth. With our six species, this tree
looks exactly as we would expect.
If the tree doesn’t look correct, then this won’t impact orthogroup inference, but will affect
our measures of gene duplication, and might affect our assignment of orthologs and
paralogs within an orthogroup. If you need to, you can run OrthoFinder with your own
species tree (use the `-s` option).

#### Step 2: Interpreting results

Now that we are happy with our OrthoFinder run, we can start diving into the results.

- ***Orthologues***<br>
    We will start by finding orthologues of a gene that we are interested in. We will focus on the gene `ENSVURG00010002700.1` in wombats, which is an
    olfactory receptor. Let’s find out what its orthologues are in the Tammar wallaby.
    In the Orthologues directory there is a sub-directory for each species. 

    Open `Orthologues/Orthologues_Mycoplasma_hyopneumoniae/Mycoplasma_hyopneumoniae__v__Mycoplasma_agalactiae.tsv`, in a spreadsheet program (specifying that it’s tab-delimited
    if necessary). The file has three columns, `Orthogroup`, `Mycoplasma_hyopneumoniae`, and
    `Mycoplasma_agalactiae`. Find `gi|71851854|gb|AAZ44462.1|` in the table, I can see that
    the gene is in orthogroup `OG0000014` and that it has three orthologues in wallabies:
    `['gi|290752976|emb|CBH40952.1|', ' gi|290752482|emb|CBH40454.1|', ' gi|290752494|emb|CBH40466.1|']`

- ***Gene trees***<br>
    Next, we are going to look at the gene tree to see how these orthologues arose.
    OrthoFinder infers orthlologues from `resolved` gene trees using a Duplication-Loss-
    Coalescence analysis to identify the more parsimonious interpretation of the tree (see
    the OrthoFinder2 [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y) for more details).

    All of the gene trees are in one file (`Resolved_Gene_Trees/Resolved_Gene_Trees.txt`).
    Each line of the file contains the ID of an orthogroup (e.g. `OG0000008`:), followed by the
    gene tree for that orthogroup. To find the tree for certain orthogroup, just search for the
    orthogroup ID.

    We are going to view the tree for `OG0000008` on [iTOL](https://itol.embl.de/).

     ╭─┬╴Mycoplasma_hyopneumoniae_gi|144227629|gb|AAZ44523.2|
─┤ ╰╴Mycoplasma_hyopneumoniae_gi|144227727|gb|AAZ44698.2|
 │ ╭─┬╴Mycoplasma_gallisepticum_gi|284812289|gb|AAP57058.2|
 ╰─┤ ╰─┬╴Mycoplasma_gallisepticum_gi|284812274|gb|AAP57044.2|
   │   ╰╴Mycoplasma_gallisepticum_gi|284812288|gb|AAP57057.2|
   ╰─┬╴Mycoplasma_gallisepticum_gi|31541755|gb|AAP57054.1|
     ╰─┬╴Mycoplasma_gallisepticum_gi|31541760|gb|AAP57059.1|
       ╰─┬╴Mycoplasma_gallisepticum_gi|284812286|gb|AAP57055.2|
         ╰─┬╴Mycoplasma_gallisepticum_gi|284812285|gb|ADB96899.1|
           ╰─┬╴Mycoplasma_gallisepticum_gi|284812284|gb|ADB96898.1|
             ╰╴Mycoplasma_gallisepticum_gi|284812290|gb|AAP57060.2|


    Looking at the gene tree, we can see that there have been several gene duplications in
    the lineage leading to wallabies (*`Notamacropus`*). This has resulted in a one-to-three
    orthology relationship, i.e. all three of the wallaby genes are equally related to the
    wombat gene `ENSVURG00010002700.1`. It’s often the case that orthology relationships
    aren’t one-to-one, and it’s important to know this—you don’t want to spend months doing
    experiments on ‘the orthologue’ only to find out later there are actually three!

- ***Gene duplications***<br>
    Having the gene trees means that OrthoFinder can identify all gene duplication events
    that occurred. There is a folder called `Gene_Duplication_Events` that has two files that
    allow us to explore duplications. Let’s first open
    `Gene_Duplication_Events/SpeciesTree_Gene_Duplications_0.5_Support.txt` in iTOL.
    Go into the `Advanced` tab on the Control Panel and select `Display` next to `Node IDs`
    to see the node labels

          ╭╴N1_6╶┬╴Mycoplasma_hyopneumoniae_69
╴N0_0╶┤      ╰╴Mycoplasma_agalactiae_123
      ╰╴N2_11╶┬╴Mycoplasma_gallisepticum_119
              ╰╴Mycoplasma_genitalium_10

    This gives a summary of gene duplication events. Each node shows the node name
    followed by an underscore and then the number of well-supported gene duplication
    events mapped to each node in the species tree. Gene-duplication events are
    considered `well-supported` if at least `50%` of the descendant species have retained
    both copies of the duplicated gene. For the common ancestor of the mammals, `N2`, 
    there were `11` of these well-supported gene duplication events. The numbers after the
    species names are the number of `terminal` duplications that map to that species, rather
    than an internal node of the species tree. 

    We can see the full list of gene duplication
    events in the file `Gene_Duplication_Events/Duplications.tsv`. Here are just a few lines
    from the file:

    | Orthogroup   | Species Tree Node    | Gene Tree Node   |   Support | Type     | Genes 1              | Genes 2              |
|--------------|----------------------|------------------|-----------|----------|----------------------|----------------------|
| OG0000000    | Mycoplasma_gallisept | n0               |         1 | Terminal | Mycoplasma_gallisept | Mycoplasma_gallisept |
| OG0000000    | Mycoplasma_gallisept | n1               |         1 | Terminal | Mycoplasma_gallisept | Mycoplasma_gallisept |
| OG0000000    | Mycoplasma_gallisept | n2               |         1 | Terminal | Mycoplasma_gallisept | Mycoplasma_gallisept |
| OG0000000    | Mycoplasma_gallisept | n3               |         1 | Terminal | Mycoplasma_gallisept | Mycoplasma_gallisept |
| OG0000000    | Mycoplasma_gallisept | n4               |         1 | Terminal | Mycoplasma_gallisept | Mycoplasma_gallisept |

    Each gene duplication event is cross-referenced to the species tree node, and the node
    in the gene tree. It also lists the genes descended from each of the two copies arising
    from the gene duplication event. We can check this out for our wombat olfactory
    receptor orthologues.

    # <p align="center" class="figure-wrapper">
    # <img src="{ site.baseurl }/assets/images/Beginner_image9.png" alt="Beginner_img9" width="1200"/>
    # </p>

    These events are also summarised by orthogroup and by species tree node in the files
    Duplications_per_Orthogroup.tsv and Duplications_per_Species_Tree_Node.tsv which are
    both in the directory `Comparative_Genomics_Statistics/`.

- ***Orthogroups***<br>
    Often we’re interested in group-wise species comparisons, that is comparisons across a
    clade of species rather than between a pair of species. The generalisation of orthology
    to multiple species is the orthogroup. Just like orthologues are the genes descended
    from a single gene in the last common ancestor of a pair of species an orthogroup is the
    set of genes descended from a single gene in a group of species. Each gene tree from
    OrthoFinder, for example the one above, is for one orthogroup. The orthogroup gene
    tree is the tree we need to look at if we want it to include all pairwise orthologues. And
    even though some of the genes within an orthogroup can be paralogs of one another, if
    we tried to take any genes out then we would also be removing orthologs too.

    So if we want to do a comparison of the `equivalent` genes in a set of species, we need
    to do the comparison across the genes in an orthogroup. The orthogroups are in the file
    `Orthogroups/Orthogroups.tsv`. This table has one orthogroup per line and one species
    per column and is ordered from largest orthogroup to smallest.

- ***Hierarchical Orthogroups***<br>
    OrthoFinder3 also infers hierarchical orthogroups for each node in the species tree. A
    file equivalent to `Orthogroups/Orthogroups.tsv` is available for each node in
    `/Phylogenetic_Hierarchical_Orthogroups`. You can compare the node number (e.g. `N3`)
    to the species tree, to see which species will be included.

- ***Orthogroup sequences***<br>
    For each orthogroup there is a FASTA file in Orthogroup_Sequences/ which contains the
    sequences for the genes in that orthogroup.

- ***Other results files***<br>
    We have now covered all of the main output files that will be useful to most users, but
    OrthoFinder also outputs much more useful information! A full description of the output
    files is available below.

    There are also some useful community tools that allow interactive viewing of results,
    such as [OrthoBrowser](https://orthobrowserexamples.netlify.app/)

