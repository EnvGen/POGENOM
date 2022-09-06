Installation
============

List of required software/packages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- `snakemake <https://snakemake.readthedocs.io/en/stable/>`_
- `samtools <http://www.htslib.org/>`_
- `vcflib <https://github.com/vcflib/vcflib>`_
- `freebayes <https://github.com/ekg/freebayes>`_
- `numpy <https://numpy.org/>`_
- `python <https://www.python.org/>`_

Aligners:
- `bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/>`_
- `bwa-mem2 <https://github.com/bwa-mem2/bwa-mem2>`_

bam filter:
- `coverm filter <https://wwood.github.io/CoverM/coverm-filter.html>`_


Installing the pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The easiest and recommended way to install this pipeline is through conda in an isolated environment.
Below an example of how to install Miniconda3 (on Linux) and how to set up the pipeline (including its required software) in an environment:

1. Installing miniconda (Linux)

Download miniconda::

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

Then, execute the script::

    bash Miniconda3-latest-Linux-x86_64.sh

Answer "yes" to question. Then, close and then re-open your terminal window. Now the command conda will work.

It may be necessary to source, with the command::

    source ~/.bashrc

2. Download the pipeline and Install the pipeline software

Clone the repository from GitHub::

    git clone https://github.com/EnvGen/POGENOM

Now, go to the directory Input_POGENOM::

    cd POGENOM/Input_POGENOM

This directory contains:

1. Empty RAW_DATA/Reads and RAW_DATA/Genomes folders
2. Config_file, snakefiles and src directories
3. The Input_POGENOM.sh script.

Create the virtual environment (with name ip_env) using the command::

    conda env create -f config_files/Input_POGENOM_conda_env_setup.yaml
