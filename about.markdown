---
layout: page
title: About
permalink: /about/
---


This is a software developed by the [Steven Kelly Lab](http://www.stevekellylab.com/)




You can find the source code for OrthoFinder at GitHub:
[OrthoFinder](https://github.com/OrthoFinder/OrthoFinder)

The Kelly Lab develops software for analysis of biological data (including, genes genomes and quantitative data). This software provides the foundation of our experimental research. Details on individual methods can be found at the links below.
- [RECUR](https://github.com/OrthoFinder/RECUR): Identifying recurrent amino acid substitutions from multiple sequence alignments

- [OGtoberfest](https://github.com/OrthoFinder/OGtoberfest): Benchmarking Orthogroups and Orthologs

- [Clust](https://github.com/BaselAbujamous/clust): Automatic extraction of optimal co-expressed gene clusters from gene expression data

- [TransRate](http://hibberdlab.com/transrate): Reference free quality assesment of de novo assembled transcriptomes.

MergeAlign: Higher accuracy multiple sequence alignments by alignment averaging.

CodonMuSe: A method for analysing the effect of translation selection and selection acting on biosynthetic cost on nucleotide sequences.

OrthoFiller: An automated method for detecting and correcting false negative genome annotation errors.

OMGene: a method for detecting and correcting erroneous gene models by optimizing evolutionary consistency

STRIDE: Species Tree Root Inference from gene Duplicaiton Events

STAG: A method for inferring species trees from paralogous gene families

DendroBLAST:  Fast and accurate phylogenetic trees from BLAST scores.

SLaP mapper:  Identification and and analysis of Spliced-Leader Addition and Polyadenylation sites in kinetoplastid genomes.

Conditional Orthology Assignment: a webserver for identifying homologous transcripts from a whole de novo transcriptome assembly.


## Meet&nbsp;the&nbsp;Team {#team}

### Current&nbsp;Team&nbsp;Member {#team-current}

<div class="team-grid">
{% assign current_team = site.data.team | where: "status", "Current" %}
{% for member in current_team %}
  {% include team-card.html member=member %}
{% endfor %}
</div>

---

### Former&nbsp;Team&nbsp;Member {#team-former}

<div class="team-grid">
{% assign former_team = site.data.team | where: "status", "Former" %}
{% for member in former_team %}
  {% include team-card.html member=member %}
{% endfor %}
</div>
