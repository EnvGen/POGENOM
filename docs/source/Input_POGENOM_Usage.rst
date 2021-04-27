Usage
=====

Organise your data
^^^^^^^^^^^^^^^^^^

The pipeline is capable of analysing several datasets with samples (paired fastq read files) and several genomes (fasta files). The minimum data required is one dataset with one pair of read files (forward and reverse) and one genome.

Input data (reads and genomes) must be stored in the directory ``RAW_DATA/``, as explained below:

1. READS

Go to the directory Input_POGENOM/RAW_DATA/Reads/::

    cd Input_POGENOM/RAW_DATA/Reads/

Create a directory for each dataset (with an optional name)::

    mkdir <dataset_name>

There must be a least one dataset.

Copy or link (recommended) reads to this directory::

    cp path/to/reads/* <dataset_name>/. or ln -s path/to/reads/* .

Make sure that read file names follow the syntax:

forward reads::

    <sample_name><fwd_index><reads_ext> e.g., P6071_505_R1.fq.gz, where sample_name = P6071_505, fwd_index = _R1 , and reads_ext = .fq.gz

reverse reads::

    <sample_name><rev_index><reads_ext> e.g., P6071_505_R1.fq.gz, where sample_name = P6071_505, rev_index = _R2 , and reads_ext = .fq.gz

2. GENOMES

Go to the directory cd Input_POGENOM/RAW_DATA/Genomes/::

    cd Input_POGENOM/RAW_DATA/Genomes/

Create a directory for each dataset::

    mkdir <dataset_name>

There must be a least one dataset.

Copy or link (recommended) genomes in FASTA format to this file::

    cp path/to/Genomes/* <dataset_name>/. or ln -s path/to/Genomes/* <dataset_name>/.

Make sure that genome names follow the syntax::

    <genome_name><genome_ext>, e.g., P6071_521_bin125.fa, where genome_name = P6071_521_bin125, and genome_ext = .fa


Run the pipeline
^^^^^^^^^^^^^^^^
If you are using conda, activate the pipeline environment by typing::

    conda activate ip_env

If you are not in the working directory, go there using the command::

    cd path/to/Input_POGENOM

1. Set parameters in the config file

In the "Input_POGENOM_config.json" file, set the parameters to be used. It contains the pipeline parameters. Below an example:

"workdir":"absolute/path/to/Input_POGENOM",
  It must be an absolute path.

"dataset": "name of the dataset to be analysed",
  It cannot be empty.

"mode_prefilt": "TRUE",
  To activate mode_prefilt set option to ``"TRUE"``.
  The pipeline will do a quick prescreening by mapping a subset of the reads from each sample, to estimate the     coverage of the samples and determine which should be included.
  If no prescreening (prefilt) is required, set option to ``"FALSE"``.

"fraction": "0.15",
  Fraction of reads to be subsampled when running the pipeline using "mode_prefilt".
  Lowering the fraction increases the uncertaintly in the coverage estimates.
  Increasing the fraction increases the size of the directory ``<temp_sub_Reads_dir>/Reads/`` and the runtime.
  Required when "mode_prefilt" used.

"temp_sub_Reads_dir": "PREFILT",
  Directory storing the subsampled reads when running the pipeline using "mode_prefilt". The size of this directory will be    "fraction" * the size of "dataset".
  Required when mode "mode_prefilt" used.

"remove_subreads": "FALSE",
  If ``"TRUE"`` the directory of subsampled reads (i.e., ``<temp_sub_Reads_dir>/Reads/``) will be removed after usage. This directory is created during sample prescreening, when "mode_prefilt" is used.

"min_coverage": 20,
  Minimum median coverage depth per sample per genome. Integer. Sample/genome combinations below this threshold will not be included in the subsequent analysis.
  It cannot be empty.

"min_breadth": 40,
  Minimum coverage breadth (percentage of genome covered) per sample per genome. Integer.
  It cannot be empty.

"subsample": "TRUE",
  Sample/genome combinations fulfilling the user-defined thresholds (i.e., min_coverage, min_breadth) will, per default (``"TRUE"``), be downsampled to the user-defined min_coverage median coverage. If the user do not want to downsample, set option to ``"FALSE"``. It cannot be empty.

"min_bsq_for_cov_median_calculation": 15,
  Minimum base quality when counting the coverage depth per genome position during coverage calculation. Integer. It cannot be empty.

"threads": 15,
  Number of threads. Integer. It cannot be empty.

"genomes_ext": ".fa",
  Extention used on your genome files.

"reads_ext": ".fq.gz",
  Extention used on your read files. For instance, ".fq.gz" if files are named "sample_R1.fq.gz & sample_R2.fq.gz".

"fwd_index": "_R1",
  Index used to define forward reads.

"rev_index": "_R2",
  Index used to define reverse reads.

"bowtie2_params": "--ignore-quals --mp 1,1 --np 1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.05",
  Bowtie2 mapping parameters. The –score-min then gives the minimum score that is allowed to report an alignment.
  Here, it represents a 95% identity threshold.
  For more information visit http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

"mapqual": 20,
  Read mapping quality threshold in BAM files. Integer. Parameter used in samtools view -q {}. It cannot be empty.

"samtools_view_alignment_extra_filters": "-f 2 -F 1024",
  Filters used for selecting mapped reads to be included in the BAM file.
  Here it selects only paired reads (-f 2) and avoids optical duplicates (-F 1024).
  If no filters are required, then set an empty string ("samtools_view_alignment_extra_filters": "",)

"freebayes_parameters": "-C 4 -p 1 --pooled-continuous --read-max-mismatch-fraction 0.05 --min-alternate-fraction 0.01 -q 15",
  Parameters used during variant calling.
  By default, freebayes exclude duplicates marked as such in alignments.
  If you want to include duplicates, use the tag ``--use-duplicate-reads`` and remove "-F 1024" in "samtools_view_alignment_extra_filters".
  The flag ``-q --min-base-quality Q``, exclude alleles from analysis if their supporting base quality is less than Q.

"vcffilter_qual": "'QUAL > 20'"
  Filtering variant calling.
  Here it removes any sites with an estimated probability of not being polymorphic less than Phred 20 (corresponding to 99% probability of being a real SNP).
  It cannot be empty.

"snakemake_extra_params": "<command line 1>, <command line 2>"
  Snakemake extra command line options (comma-separated) to be used. If you don't want to use any extra command line, set an empty string, "snakemake_extra_params": "".


To access and modify this file, you can use the following command::

    nano config_files/Input_POGENOM_config.json

Modify the required items and save the file. Use Ctrl +x and answer y, to save the modifications and exit the file.

2. Run

The workflow is run with the following command::

    bash Input_POGENOM.sh

If you need to set a different path to the config file ( flag -d=<absolute path to configFile> ), please do not use relative paths (~/ nor ./)

If you are using conda, before exiting the workflow, the environment needs to be deactivated using the following command::

    conda deactivate
