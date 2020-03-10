Installation
============

List of required software/packages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- snakemake
- bowtie2
- samtools
- vcflib
- freebayes
- picard
- numpy

Some tips for installing software
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The easiest and recommended way to install this pipeline is through conda in an isolated environment.
Bellow an example of how to install conda and the pipeline in Linux is presented:

1. Installing conda (Linux)

Download miniconda::

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

Then, execute the script::

    bash Miniconda3-latest-Linux-x86_64.sh

Answer "yes" to question. Then, close and then re-open your terminal window. Now the command conda will work.

It may be necessary to source, with the command::

    source ~/.bashrc

2. Download the pipiline and Install the pipeline software

Clone the repositiory from github::

    git clone <pipeline repository>

Now, go to the directory Input_POGENOM::

    cd Input_POGENOM

This directory contains:

1. Empty RAW_DATA/Reads and RAW_DATA/MAGs folders
2. Config_file, snakefiles and src directories
3. The Input_POGENOM.sh script.

Create the virtual environment using the command::

    conda env create -f config_files/Input_POGENOM_conda_env_setup.yaml
