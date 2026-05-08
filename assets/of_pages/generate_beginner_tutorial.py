import os
import argparse
from datetime import datetime, timezone
from ete4 import Tree
from tabulate import tabulate


TUTORIAL = """---
layout: post
title:  "Beginner Tutorial"
date:   {tutorial_datetime}
updated: {tutorial_update_date}
docs: /assets/docs/beginner-tutorial.pdf
pdfs: https://github.com/OrthoFinder/OrthoFinder/blob/gh-pages/assets/docs/beginner-tutorial.pdf
---


This tutorial will cover:

- [Downloading OrthoFinder](#downloading-orthofinder)
- [Running OrthoFinder](#running-orthofinder)
- [Exploring the results of OrthoFinder](#exploring-the-results-of-orthofinder)

OrthoFinder requires as input the amino acid sequences for all the protein coding genes in
your species of interest. We provide a separate tutorial for
[getting input files]({{{{ site.baseurl}}}}/tutorials/gettting-input-data/) for OrthoFinder.

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

We then need to run these commands.

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
dataset, which you can download from our [GitHub](https://GitHub.com/OrthoFinder/OrthoFinder/).

```orthofinder -f ExampleData/ ```

OrthoFinder will print lots of information to the command line as it runs. If you get an
error message, the best way to troubleshoot is to just google the error message. You
can also ask a question on our [GitHub](https://GitHub.com/OrthoFinder/OrthoFinder/issues).

When OrthoFinder has finished running, it will generate a folder containing the output,
with the folder named according to today’s date, for example, `ExampleData/OrthoFinder/Results_Dec18`. We’ll discuss how to interpret and analyse these files and folders later on, in the
[Exploring the results](#exploring-the-results-of-orthofinder) section of the tutorial.

### Running OrthoFinder

You can now run OrthoFinder
First, you have to open a terminal and navigate to the directory where your files
are. You can now run OrthoFinder on your proteomes.

```bash
orthofinder -f primary_transcripts
```

That’s it! OrthoFinder will print updates on its progress to the terminal, and tell you when
it’s finished. To see what options you might want to adjust for your own data, check out the [GitHub](https://GitHub.com/OrthoFinder/OrthoFinder/),
or the [Advanced Tutorial]({{{{ site.baseurl }}}}/tutorials/advanced-tutorial/) page.

### Exploring the results of OrthoFinder

OrthoFinder creates a results directory named OrthoFinder inside the proteome
directory, and puts the results here.
My results directory `{target_dir}` looks like this:

{{% include_relative of-output-tree.md %}}


#### Step 1: Quality Control

Before we start diving into the orthogroups, it would behoove us to check the quality
of the OrthoFinder run. We want to make sure that most genes across all species
have been assigned to orthogroups, and that the species tree looks realistic.

Open the file `Statistics_Overall.tsv` from the folder `Comparative_Genomics_Statistics`. This file can be opened in spreadsheet software
like Microsoft Excel, or in a text editor like Notepad.

On the 5th line, we can see the `Percentage of genes in orthogroups`, which in my case
is `{percentage_of_genes_in_orthogroup}`.

{stats_overall}

A good rule of thumb is that this number should be `>80%`. If not, you are likely missing
some orthology relationships that actually exist. The best way to fix this would be better
species sampling.

Now open the file `Statistics_PerSpecies.tsv`, from the same folder. This file gives us the
`%` of genes in each species that are assigned to orthogroups, rather than the
percentage for all genes across species.

You can see here that we capture most genes across all species.

{statisitcs_per_species}

The lowest percentage is the *`{lowest_percent_species}`*, but we still managed to assign `{lowest_percent_value}` of its
genes to orthogroups. The key message here is that it’s always a good idea to look at
this information before you start interpreting your results. If the numbers were too low for
one species, we might want to consider sampling more species to fill in the long
evolutionary divergence between species.

One more useful thing to do before we really start to dive in is to look at the species
tree. You can do this by opening the tree in [iTOL](https://itol.embl.de/upload.cgi) by either copy and pasting the file content
or uploading the file directly.

<p align="center" class="figure-wrapper">
  <img src="{{{{ site.baseurl }}}}/assets/images/species_tree.png" alt="species_tree" width="500"/>
</p>

We now want to do some common-sense checking that everything appears to be in
order, and we aren’t rewriting the history of life on earth. With our species, this tree
looks exactly as we would expect.
If the tree doesn’t look correct, then this won’t impact orthogroup inference, but will affect
our measures of gene duplication, and might affect our assignment of orthologs and
paralogs within an orthogroup. If you need to, you can run OrthoFinder with your own
species tree (use the `-s` option).

#### Step 2: Interpreting results

Now that we are happy with our OrthoFinder run, we can start diving into the results.

- ***Orthologues***<br>
    We will start by finding orthologues of a gene that we are interested in.
    In the Orthologues directory there is a sub-directory for each species.

    Open `{ortholog_path}`, in a spreadsheet program (specifying that it’s tab-delimited
    if necessary). The file has three columns, `Orthogroup`, `{species1}`, and
    `{species2}`. Find `{gene_of_interest}` in the table, I can see that
    the gene is in orthogroup `{orthogroup}` and that its orthologs are:
    `{orthologs}`.

- ***Gene trees***<br>
    Next, we are going to look at the gene tree to see how these orthologues arose.
    OrthoFinder infers orthlologues from `resolved` gene trees using a Duplication-Loss-
    Coalescence analysis to identify the more parsimonious interpretation of the tree (see
    the OrthoFinder2 [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y) for more details).

    All of the gene trees are in one file (`Resolved_Gene_Trees/Resolved_Gene_Trees.txt`).
    Each line of the file contains the ID of an orthogroup (e.g. `{gene_tree_orthogroup}`:), followed by the
    gene tree for that orthogroup. To find the tree for certain orthogroup, just search for the
    orthogroup ID.

    We are going to view the tree for `{gene_tree_orthogroup}`.

    <p align="center" class="figure-wrapper">
        <img src="{{{{ site.baseurl }}}}/assets/images/gene_tree.png" alt="gene_tree" width="1000"/>
    </p>


    Looking at the gene tree, we can see if there are any gene duplications.

- ***Gene duplications***<br>
    Having the gene trees means that OrthoFinder can identify all gene duplication events
    that occurred. There is a folder called `Gene_Duplication_Events` that has two files that
    allow us to explore duplications. Let’s first open
    `Gene_Duplication_Events/SpeciesTree_Gene_Duplications_0.5_Support.txt` in iTOL.
    Go into the `Advanced` tab on the Control Panel and select `Display` next to `Node IDs`
    to see the node labels.

    <p align="center" class="figure-wrapper">
        <img src="{{{{ site.baseurl }}}}/assets/images/duplication_tree.png" alt="duplication_tree" width="600"/>
    </p>


    This gives a summary of gene duplication events. Each node shows the node name
    followed by an underscore and then the number of well-supported gene duplication
    events mapped to each node in the species tree. Gene-duplication events are
    considered `well-supported` if at least `50%` of the descendant species have retained
    both copies of the duplicated gene. For the node `{node_example}`,
    there were `{number_of_duplications_for_node}` of these well-supported gene duplication events. The numbers after the
    species names are the number of `terminal` duplications that map to that species, rather
    than an internal node of the species tree.

    We can see the full list of gene duplication
    events in the file `Gene_Duplication_Events/Duplications.tsv`. Here are just a few lines
    from the file:

    {duplication_table}

    Each gene duplication event is cross-referenced to the species tree node, and the node
    in the gene tree. It also lists the genes descended from each of the two copies arising
    from the gene duplication event.

    These events are also summarised by orthogroup and by species tree node in the files
    `Duplications_per_Orthogroup.tsv` and `Duplications_per_Species_Tree_Node.tsv` which are
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

    So if we want to do a comparison of the equivalent genes in a set of species, we need
    to do the comparison across the genes in an orthogroup. The orthogroups are in the file
    `Orthogroups/Orthogroups.tsv`. This table has one orthogroup per line and one species
    per column and is ordered from largest orthogroup to smallest.

- ***Hierarchical Orthogroups***<br>
    OrthoFinder3 also infers hierarchical orthogroups for each node in the species tree. A
    file equivalent to `Orthogroups/Orthogroups.tsv` is available for each node in
    `/Phylogenetic_Hierarchical_Orthogroups`. You can compare the node number (e.g. `N2`)
    to the species tree, to see which species will be included.

- ***Orthogroup sequences***<br>
    For each orthogroup there is a FASTA file in Orthogroup_Sequences/ which contains the
    sequences for the genes in that orthogroup.

- ***Other results files***<br>
    We have now covered all of the main output files that will be useful to most users, but
    OrthoFinder also outputs much more useful information! A full description of the output
    files is available below.

    There are also some useful community tools that allow interactive viewing of results,
    such as [OrthoBrowser](https://orthobrowserexamples.netlify.app/).

"""


def get_Statistics_Overall(target_dir):
    #
    # Returns string containing the first six lines of Statistics_Overall.tsv" and the Percentage of genes in orthogroups
    percentage_of_genes_in_orthogroup = 0
    Statistics_Overall_Lines = []
    max_len = 0
    comparaitve_genome_statisitcs_overall = "/".join([target_dir,"Comparative_Genomics_Statistics","Statistics_Overall.tsv"])
    with open(comparaitve_genome_statisitcs_overall) as statisitcs_overall:
        for lines in range(0,7):
            line = next(statisitcs_overall)
            #sys.exit()
            Statistics_Overall_Lines.append(line.split("\t"))
            if max_len < len(line.split("\t")[0]):
                max_len = len(line.split("\t")[0])
            if line.startswith("Percentage of genes in orthogroups"):
                percentage_of_genes_in_orthogroup = line.split("\t")[1]

    table = []
    for line in Statistics_Overall_Lines:
        new_space = (" ")*(max_len - len(line[0]) + 3)
        table.append([line[0],line[1]])
        #sys.exit()

    return tabulate(table,tablefmt="github"), percentage_of_genes_in_orthogroup


def get_Statistics_PerSpecies(target_dir):
    comparaitve_genome_statisitcs_perSpecies = "/".join([target_dir,"Comparative_Genomics_Statistics","Statistics_PerSpecies.tsv"])
    first_6_lines = []
    lowest_percent_species = ""
    lowest_percent_value = 100

    with open(comparaitve_genome_statisitcs_perSpecies) as statisitcs_per_species:
        for lines in range(0,6):
            line = next(statisitcs_per_species)
            first_6_lines.append(line.rstrip())
    lengths = []
    for species in first_6_lines[0].split("\t"):
        lengths.append(len(species))

    first_col_lengths = []
    for first_col in first_6_lines:
        first_col_lengths.append(len(first_col.split("\t")[0]))

    lengths[0] = max(first_col_lengths)
    reformatted_statisitcs_per_species = []
    for line in first_6_lines:
        new_line = []
        for pos,element in enumerate(line.split("\t")):
            if len(element) < lengths[pos]:
                new_line.append(element + " "*(lengths[pos] - len(element)))
            else:
                new_line.append(element)
        reformatted_statisitcs_per_species.append("   ".join(new_line))

    table = []

    for i, line in enumerate(reformatted_statisitcs_per_species):
        line = line.split()
        if i == 0:
            header = line
        else:
            name = " ".join(line[:-4])
            table.append([name] + line[-4:])

    for species,percent_of_genes_in_og in zip(first_6_lines[0].split("\t")[1:], first_6_lines[4].split("\t")[1:]):
        if float(percent_of_genes_in_og) < lowest_percent_value:
            lowest_percent_value = float(percent_of_genes_in_og)
            lowest_percent_species = species

    return tabulate(table, header, tablefmt="github"), lowest_percent_species, lowest_percent_value


def get_ortholog(target_dir, species1, species2):
    # pick some thing??Mycoplasma_hyopneumoniae v  Mycoplasma_agalactiae
    orthologs_target_path = "/".join([
        target_dir,
        "Orthologues",
        "_".join(("Orthologues", species1)),
        "__v__".join((species1, species2)) + ".tsv"
    ])
    ## first three columns will also be the same here...Orthogroup	Mycoplasma_hyopneumoniae	Mycoplasma_agalactiae
    with open(orthologs_target_path) as ortholog_file:
        for line in ortholog_file:
            no_genes = 0
            for i in line.split("\t"):
                no_genes += 1
                no_genes += i.count(",")
            no_genes -= 1
            if no_genes == 4 and line.split("\t")[0].count(",") == 0:
                gene_of_interest = line.split("\t")[1]
                orthologs = line.rstrip().split("\t")[2].split(",")
                orthogroup = line.split("\t")[0]
                break
    return gene_of_interest,orthologs,orthogroup



def get_random_gene_tree(target_dir):
    Resolved_genes_tree_path = "/".join([target_dir,"Resolved_Gene_Trees","Resolved_Gene_Trees.txt"])
    with open(Resolved_genes_tree_path) as gene_tree_file:
        for line in gene_tree_file:
            OG_Tree_Pair = line.rstrip().split(": ")
            ete3_tree = Tree(OG_Tree_Pair[1], parser=1)

            if len(list(ete3_tree.leaves())) < 12 and 4 < len(list(ete3_tree.leaves())):
                orthogroup = OG_Tree_Pair[0]
                break
    return orthogroup



def gene_duplications(target_dir):
    gene_duplications = "/".join([target_dir,"Gene_Duplication_Events"])
    gene_duplication_tree = Tree(open("/".join([gene_duplications,"SpeciesTree_Gene_Duplications_0.5_Support.txt"])).read(),parser=1)

    ## Example NODE for now N2
    for node in gene_duplication_tree.traverse():
        if node.name.startswith("N2"):
            node_example = node.name.split("_")[0]
            number_of_duplications_for_node = node.name.split("_")[1]


    #print(gene_duplication_tree)
    first_duplication_lines = []
    gene_duplication_table = "/".join([gene_duplications,"Duplications.tsv"])
    with open(gene_duplication_table) as duplication_table:
        for lines in range(0,6):
            line = next(duplication_table)
            first_duplication_lines.append(line.rstrip())

    gap = 0
    for col_name in first_duplication_lines[0].split("\t"):
        if gap < len(col_name):
            gap = len(col_name)
    reformatted_table = []

    for line in first_duplication_lines:
        new_line = []
        for element in line.split("\t"):
            if len(element) < gap + 3:
                new_line.append(element)
            else:
                new_line.append(element[:gap + 3])
        reformatted_table.append(new_line)
    return tabulate(reformatted_table[1:], reformatted_table[0], tablefmt="github"),\
        node_example,\
        number_of_duplications_for_node



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("path")
    args = parser.parse_args()
    target_dir = args.path
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S %z")
    update_date = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    if not os.path.exists(target_dir):
        print("The target directory doesn't exist")
        raise SystemExit(1)

    if not os.path.isdir(target_dir):
        print("The input path is not a directory")
        raise SystemExit(1)

    stats_overall, percentage_of_genes_in_orthogroup = get_Statistics_Overall(target_dir)
    statisitcs_per_species, lowest_percent_species, lowest_percent_value = get_Statistics_PerSpecies(target_dir)

    species1 = "Mycoplasma_hyopneumoniae"
    species2 = "Mycoplasma_agalactiae"

    gene_of_interest,orthologs,ortholog_orthogroup = get_ortholog(target_dir, species1, species2)
    gene_tree_orthogroup = get_random_gene_tree(target_dir)
    duplication_table, node_example, number_of_duplications_for_node = gene_duplications(target_dir)

    tutorial = TUTORIAL.format(
        tutorial_datetime=ts,
        tutorial_update_date=update_date,
        target_dir=target_dir,
        stats_overall=stats_overall,
        percentage_of_genes_in_orthogroup=percentage_of_genes_in_orthogroup,
        statisitcs_per_species=statisitcs_per_species,
        lowest_percent_species=lowest_percent_species,
        lowest_percent_value=lowest_percent_value,
        ortholog_path=os.path.join("Orthologues",
                                    "_".join(("Orthologues", species1)),
                                    "__v__".join((species1, species2)) + ".tsv"),
        species1=species1,
        species2=species2,
        gene_of_interest=gene_of_interest,
        orthologs=", ".join(orthologs),
        orthogroup=ortholog_orthogroup,
        gene_tree_orthogroup=gene_tree_orthogroup,
        duplication_table=duplication_table,
        node_example=node_example,
        number_of_duplications_for_node=number_of_duplications_for_node
    )
    print(tutorial)
