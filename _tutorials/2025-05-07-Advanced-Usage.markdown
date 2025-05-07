---
layout: post
title:  "Advanced Usage"
date:   2025-05-07 09:16:16 +0000
categories: tutorials
permalink: /:categories/:title.html
---



## Advanced usage
### Manually installing dependencies

To run OrthoFinder3 on default setting, you will need diamond, famsa, fasttree, MCL
If you install OrthoFinder using a recommended method, you shouldn’t need to install
any of the dependencies manually, but in case you do;
```bash
Diamond
https://github.com/bbuchfink/diamond
FAMSA
https://github.com/refresh-bio/FAMSA
FastTree
http://www.microbesonline.org/fasttree/#Install
MCL
https://github.com/micans/mcl

```
### Running BLAST searches separately
The `-op` option will prepare the files in the format required by
OrthoFinder and print the set of BLAST commands to run.
```bash
orthofinder -f fasta_files_directory -op
```
This is useful if you want to manage the BLAST searches yourself. For example, you
may want to distribute them across multiple machines. Once the BLAST searches have 
been completed the orthogroups can be calculated using the `-b` command (see below)

### Use pre-computed BLAST
It is possible to run OrthoFinder with pre-computed BLAST results provided they are in 
the correct format. They can be prepared in the correct format using the `-op` command 
and, equally, the files from a previous OrthoFinder run are also in the correct format to rerun using the `-b` option. The command is simply:
```bash
orthofinder -b directory_with_processed_fasta_and_blast_results
```
If you are running the BLAST searches yourself it is strongly recommended that you
use the `-op` option to prepare the files first (see Section "Running BLAST Searches 
Separately")