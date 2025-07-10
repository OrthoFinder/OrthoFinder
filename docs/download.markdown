---
layout: post
title: Download & Install
permalink: /download_and_install/
---

<!-- <h2 style="margin-bottom: 1.2rem;"><strong>Download & Install</strong></h2> -->

### Download

Please pick the right `.tar.gz` file to download from the list.



<div class="download-item">
  📦 
    <a href="https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/orthofinder-linux-intel-3.1.0.tar.gz" class="btn">
        orthofinder-linux-intel-3.1.0.tar.gz
    </a>
    <span class="download-count" data-asset="orthofinder-linux-intel-3.1.0.tar.gz">
    <i class="fa fa-download" aria-hidden="true"></i>
    </span>
</div>
<!--more-->

Older versions of OrthoFinder can be downloaded from [here]({{ site.baseurl }}/older-versions/).  
For a full list of OrthoFinder before v3.0.1b1, see the [GitHub Releases](https://github.com/davidemms/OrthoFinder/releases).


### Installation

#### Install in conda (recommended) 

The simplest way to install OrthoFinder is through [conda](https://anaconda.org/bioconda/orthofinder). If you're unfamiliar with conda, [this tutorial](https://www.machinelearningplus.com/deployment/conda-create-environment-and-everything-you-need-to-know-to-manage-conda-virtual-environment/) offers a beginner-friendly introduction.

```bash
conda create -n of3_env python=3.10
conda activate of3_env
conda install orthofinder
```

Alternatively, you could install via github, or download the source code and install locally.

#### Install via github
```bash
python3 -m venv of3_env 
. of3_env/bin/activate
pip install git+https://github.com/OrthoFinder/OrthoFinder.git
```

#### Install locally from source code

The following commands provide three ways to download the source code of OrthoFinder locally into a directory named `OrthoFinder`.
```bash
# Download via git 
git clone https://github.com/OrthoFinder/OrthoFinder.git
# or download the orthofinder-linux-intel-3.1.0.tar.gz and unzip it into OrthoFinder if you are on a Linux Intel machine
mkdir OrthoFinder && \
  wget -qO- https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.0/orthofinder-linux-intel-3.1.0.tar.gz | \
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


#### Test your installation

Once you have installed OrthoFinder, you can print the help information and version, and test it on the [example data](https://github.com/OrthoFinder/OrthoFinder/tree/master/ExampleData).

```bash
orthofinder --help # Print out help informatioin
orthofinder --version # Check the version
orthofinder -f ExampleData # Test OrthoFinder on an example dataset
```

#### Uninstalling
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