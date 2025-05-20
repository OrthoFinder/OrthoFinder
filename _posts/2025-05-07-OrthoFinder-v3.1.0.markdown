---
layout: post
title:  "OrthoFinder-v3.1.0"
date:   2025-05-07 09:16:16 +0000
categories: release
permalink: /:categories/:day/:month/:year/:title.html
---

<div class="download-item">
  📦 
    <a href="https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/orthofinder-linux-intel-3.1.0.tar.gz" class="btn">
        orthofinder-linux-intel-3.1.0.tar.gz
    </a>
    <span class="download-count" data-asset="orthofinder-linux-intel-3.1.0.tar.gz">
    <i class="fa fa-download" aria-hidden="true"></i>
    </span>
</div>
<div class="download-item">
  📦 
    <a href="https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/orthofinder-linux-arm-3.1.0.tar.gz" class="btn">
        orthofinder-linux-arm-3.1.0.tar.gz
    </a>
    <span class="download-count" data-asset="orthofinder-linux-arm-3.1.0.tar.gz">
    <i class="fa fa-download" aria-hidden="true"></i>
    </span>
</div>
<div class="download-item">
  📦 
    <a href="https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/orthofinder-macos-arm-3.1.0.tar.gz" class="btn">
        orthofinder-macos-arm-3.1.0.tar.gz
    </a>
    <span class="download-count" data-asset="orthofinder-macos-arm-3.1.0.tar.gz">
    <i class="fa fa-download" aria-hidden="true"></i>
    </span>
</div>
<div class="download-item">
  📦 
    <a href="https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/orthofinder-macos-intel-3.1.0.tar.gz" class="btn">
        orthofinder-macos-intel-3.1.0.tar.gz
    </a>
    <span class="download-count" data-asset="orthofinder-macos-intel-3.1.0.tar.gz">
    <i class="fa fa-download" aria-hidden="true"></i>
    </span>
</div>
<div class="download-item">
  📦 
    <a href="https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/orthofinder-macos-universal-3.1.0.tar.gz" class="btn">
        orthofinder-macos-universal-3.1.0.tar.gz
    </a>
    <span class="download-count" data-asset="orthofinder-macos-universal-3.1.0.tar.gz">
    <i class="fa fa-download" aria-hidden="true"></i>
    </span>
</div>
<!--more-->


### 🧬 Input
- Protein or codon MSA (FASTA format)
- Defined outgroup species or clade
- *(Optional)* Species tree in Newick format — used to constrain phylogenetic inference if provided

---

### 📤 Output
- `*.OrthoFinder.tsv` — list of OrthoFinderrent amino acid substitutions inferred from the tree
- `.OrthoFinder/` — directory containing intermediate files:
  - Model selection output
  - Inferred phylogenetic tree
  - Ancestral state reconstructions
  - Site substitution matrices

---

### ⚙️ Dependencies
- OrthoFinder v3.1.0 includes the **IQ-TREE2** binary (bundled in `OrthoFinder/bin/`) and uses it for:
  - Phylogenetic tree inference
  - Ancestral sequence reconstruction
  - Simulation of alignments for null models

No external installation of **IQ-TREE2** is required.

---
### ⬇️  Dowload 

- orthofinder-linux-intel-3.1.0.tar.gz
```bash
mkdir OrthoFinder && \
    wget -qO- https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/orthofinder--linux-intel-3.1.0.tar.gz | \
    tar -xz --strip-components=1 -C OrthoFinder
```

- orthofinder-linux-arm-3.1.0.tar.gz
```bash
mkdir OrthoFinder && \
    wget -qO- https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/orthofinder-linux-arm-3.1.0.tar.gz | \
    tar -xz --strip-components=1 -C OrthoFinder
```

- orthofinder-macos-arm-3.1.0.tar.gz
```bash
mkdir OrthoFinder && \
    wget -qO- https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/orthofinder-macos-arm-3.1.0.tar.gz | \
    tar -xz --strip-components=1 -C OrthoFinder
```

- orthofinder-macos-intel-3.1.0.tar.gz
```bash
mkdir OrthoFinder && \
    wget -qO- https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/orthofinder-macos-intel-3.1.0.tar.gz | \
    tar -xz --strip-components=1 -C OrthoFinder
```

- orthofinder-macos-universal-3.1.0.tar.gz
```bash
mkdir OrthoFinder && \
    wget -qO- https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/OrthoFinder-iqtree2-macos-universal-3.1.0.tar.gz | \
    tar -xz --strip-components=1 -C OrthoFinder
```

### 📦 Install 

- Linux and MacOS
    ```bash
    cd OrthoFinder
    python3 -m venv of3_env 
    . of3_env/bin/activate
    pip install .
    OrthoFinder --version
    OrthoFinder -f ExampleData/example_alignments.aln -st AA --outgroups ExampleData/example_alignments.outgroups.txt 
    ```

    To deactivate the virtual environment, run
    ```bash
    deactivate
    ```
    Having deactivated the virtual environment, to remove OrthoFinder, please run
    ```bash
    cd ..
    rm -rf OrthoFinder
    ```


