---
layout: post
title:  "OrthoFinder-v3.1.0"
date:   2025-05-07 09:16:16 +0000
categories: release
permalink: /:categories/:day/:month/:year/:title.html
---

<div class="download-item">
    <a href="https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/orthofinder-linux-intel-3.1.0.tar.gz" class="btn">
        <svg class="icon" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" width="16" height="16">
            <path d="M12 16l4-4H8l4 4zm1-13h2v6h-2V3zm-4 6h8v2H9V9zM4 22h16v-2H4v2z" />
        </svg>
        orthofinder-linux-intel-3.1.0.tar.gz
    </a>
</div>
<div class="download-item">
    <a href="https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/orthofinder-linux-arm-3.1.0.tar.gz" class="btn">
        <svg class="icon" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" width="16" height="16">
            <path d="M12 16l4-4H8l4 4zm1-13h2v6h-2V3zm-4 6h8v2H9V9zM4 22h16v-2H4v2z" />
        </svg>
        orthofinder-linux-arm-3.1.0.tar.gz
    </a>
</div>
<div class="download-item">
    <a href="https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/orthofinder-macos-arm-3.1.0.tar.gz" class="btn">
        <svg class="icon" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" width="16" height="16">
            <path d="M12 16l4-4H8l4 4zm1-13h2v6h-2V3zm-4 6h8v2H9V9zM4 22h16v-2H4v2z" />
        </svg>
        orthofinder-macos-arm-3.1.0.tar.gz
    </a>
</div>
<div class="download-item">
    <a href="https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/orthofinder-macos-intel-3.1.0.tar.gz" class="btn">
        <svg class="icon" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" width="16" height="16">
            <path d="M12 16l4-4H8l4 4zm1-13h2v6h-2V3zm-4 6h8v2H9V9zM4 22h16v-2H4v2z" />
        </svg>
        orthofinder-macos-intel-3.1.0.tar.gz
    </a>
</div>
<div class="download-item">
    <a href="https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/orthofinder-macos-universal-3.1.0.tar.gz" class="btn">
        <svg class="icon" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" width="16" height="16">
            <path d="M12 16l4-4H8l4 4zm1-13h2v6h-2V3zm-4 6h8v2H9V9zM4 22h16v-2H4v2z" />
        </svg>
        orthofinder-macos-universal-3.1.0.tar.gz
    </a>
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
