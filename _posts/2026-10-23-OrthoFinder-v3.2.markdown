---
layout: post
title:  "OrthoFinder-v3.2"
date:   2026-10-23 18:00:00 +0000
categories: release
permalink: /:categories/:day/:month/:year/:title.html
---

<div class="download-item">
  📦
    <a href="https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.2/orthofinder-linux-intel-3.2.tar.gz" class="btn">
        orthofinder-linux-intel-3.2.tar.gz
    </a>
    <span class="download-count" data-asset="orthofinder-linux-intel-3.2.tar.gz">
    <i class="fa fa-download" aria-hidden="true"></i>
    </span>
</div>
<!--more-->

### Changes to the results files

Please check these if you have scripts that read OrthoFinder results.

- **"Orthologue" is now "ortholog"** in the names of the results files and directories, e.g. `Orthologs/`, `Single_Copy_Ortholog_Sequences/`, `Comparative_Genomics_Statistics/OrthologsStats_*.tsv` and `Orthogroups/Orthogroups_SingleCopyOrthologs.txt`. The contents of the files are unchanged. Restarting from the results of an earlier OrthoFinder version renames these files and directories to the new names.
- **Pairwise ortholog files only with `-pof`.** By default, the orthologs of each species are in one compressed file, `Orthologs/<species>.tsv.gz`. This reduces the number of files on large runs. The new `split_ortholog_files` tool writes the pairwise files afterwards, for the species you need.
- **Gene and species names are the same in every results file**, on every operating system. Each name is cleaned once from the FASTA header (or file name).

### New

- **Stop and restart an analysis at each stage**, also with `--assign`. Stop with `-op`, `-og`, `-ogt` or `-ost`, and restart with `-b`, `-fg`, `-fgt` or `-fst`. A restart continues in the same results directory, with the same run settings as the original, and `-b` runs only the sequence searches that are missing or incomplete (e.g. after a job timeout or error). With `-op`, the search commands are written to `WorkingDirectory/blast_commands.txt`.
- **New tools** `split_ortholog_files`, `species_tree_alignment` and `mcl_orthogroup_sequences`. `species_tree_alignment` writes the concatenated alignment that OrthoFinder infers the species tree from, for a finished run, e.g. one whose species tree was inferred from the gene trees with `-ST astral-pro`.
- **The species tree.** A core run infers it from the concatenated alignment of the single-copy orthogroups. An `--assign` run, or a core run with `-ST astral-pro`, infers it from the gene trees with ASTRAL-Pro.
- **More analysis threads by default.** `-a` is a quarter of `-t` (it was an eighth, at most 16). `orthofinder -h` shows the defaults for your computer.
- **Clearer messages.** When OrthoFinder stops because of an error, it says where to get help and which log file to send.
- **The tutorial** (attached, `OrthoFinder3_tutorial.pdf`) is rewritten for v3.2, with a step-by-step walk-through of the example data.

### Fixes

- `-b` on a results directory whose `Log.txt` names its working directory without a final `/` (e.g. results of an earlier version that were moved) failed with a "No such file or directory" error for `…WorkingDirectorySpeciesIDs.txt`.
- Installs from the release archives (e.g. Bioconda) report the right version (the v3.1.5 archive reported 3.1.4).
- The RAxML entry of the example user config (`user_config.json`) named a program that conda does not install, and RAxML's trees were not found with it.
### Requirements

Python 3.11 or later (tested with 3.11 to 3.14) and DIAMOND 2.1.7 or later. See [System Requirements](https://github.com/OrthoFinder/OrthoFinder#system-requirements). The programs bundled with OrthoFinder (`orthofinder-linux-intel-3.2.tar.gz`) are for Linux on Intel/AMD (x86-64) processors. On other systems, install OrthoFinder with conda.

**Full details** in [What's new in v3.2](https://github.com/OrthoFinder/OrthoFinder#whats-new-in-v32).
