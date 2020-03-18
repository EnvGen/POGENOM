Usage
=====

Organise your data
^^^^^^^^^^^^^^^^

The pipeline is capable of analyzing several datasets with samples (paired fastq read files) and several genomes (fasta files). The minimum data required is one dataset with one pair of read files (forward and reverse) and one genome.

Input data (reads and genomes) must be stored in the directory RAW_DATA/, as explained below:

1. READS

Go to the directory Input_POGENOM/RAW_DATA/Reads/::

    cd Input_POGENOM/RAW_DATA/Reads/

Create a directory for each dataset::

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

Go to the directory cd Input_POGENOM/RAW_DATA/MAGs/::

    cd Input_POGENOM/RAW_DATA/MAGs/

Create a directory for each dataset::

    mkdir <dataset_name>

There must be a least one dataset.

Copy or link (recommended) genomes in FASTA format to this file::

    cp path/to/MAGs/* <dataset_name>/. or ln -s path/to/MAGs/* <dataset_name>/.

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

"threads": 15,
  Number of threads.

"fraction": "0.15",
  Fraction of reads to be subsampled when running the pipeline using dataset: "prefilt".
  Lowering the fraction increases the uncertaintly in the coverage estimates.
  Increasing the fraction increases the size of the directory <temp_sub_Reads_dir>/Reads/ and the runtime.

"temp_sub_Reads_dir": "PREFILT",
  Directory storing the subsampled reads when running the pipeline using dataset: "prefilt". The size of this directory will be "fraction" * the size of "dataset".

"mags_ext": ".fa",
  Extention used on your genome files.

"reads_ext": ".fq.gz",
  Extention used on your read files. For instance, ".fq.gz" if files are named "sample_R1.fq.gz & sample_R2.fq.gz".

"fwd_index": "_R1",
  Index used to define forward reads.

"rev_index": "_R2",
  Index used to define reverse reads.

"bowtie2_params": "--ignore-quals --mp 1,1 --np 1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.05",
  Bowtie2 mapping parameters. The –score-min then gives the minimum score that is allowed to report an alignment.
  Here, It represents a 95% identity threshold.
  For more information visit http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

"mapqual": 20,
  Read mapping quality threshold in BAM files. Integer. Parameter used in samtools view -q {}.

"samtools_view_alignment_extra_filters": "-f 2 -F 1024",
  Filters used for selecting mapped reads to be included in the BAM file.
  Here it selects only paired reads (-f 2) and avoids optical duplicates (-F 1024).
  If no filters are required, then set an empty string ("samtools_view_alignment_extra_filters": "",)

"min_coverage": 10,
  Minimum genome coverage per sample per genome (higher than, not included), integer.
  When dataset: "prefilt", a min_coverage value lower than 10 will select all samples, and the prefilter will be obsolete.

"min_breadth": 40,
  Minimum genome breadth per sample per genome. Integer.

"min_bsq_for_cov_median_calculation": 15,
  Minimum base quality when counting the number of bases per genome position during coverage calculation. Integer number.

"freebayes_parameters": "-C 4 -p 1 --pooled-continuous --read-max-mismatch-fraction 0.05 --min-alternate-fraction 0.01 -q 15",
  Parameters used during variant calling.
  By default, freebayes exclude duplicates marked as such in alignments.
  If you want to include use tag --use-duplicate-reads and remove "-F 1024" in "samtools_view_alignment_extra_filters".
  The flag '-q --min-base-quality Q', exclude alleles from analysis if their supporting base quality is less than Q.

"vcffilter_qual": "'QUAL > 20'"
  Filtering variant calling.
  Here it removes any sites with an estimated probability of not being polymorphic less than Phred 20 (corresponding to 99% probability of being a real SNP).


To access and modify this file, you can use the following command::

    nano config_files/Input_POGENOM_config.json

Modify the required items and save the file. Use Ctrl +x and answer y, to save the modifications and exit the file.

2. Run

The workflow is run with the following command::

    bash Input_POGENOM.sh

If you need to set a different path to the config file ( flag -d=<absolute path to configFile> ), please do not use relative paths (~/ nor ./)

2.1) A dataset

If you want to run the pipeline on one dataset, please set the corresponding name in the config_file, "dataset": <dataset_name>

2.2) Several datasets

2.2.1) "prefilt"

If you want to run the pipeline on the entire sampling dataset, and only on those genomes and their corresponding samples with median coverage higher than a certain threshold (i.e., min_coverage), please set "prefilt" in the config_file, "dataset": "prefilt" and "fraction": "<your float value, default=0.15>."

When running the pipeline with dataset "prefilt", the created RAW_DATA/Reads/prefilt and RAW_DATA/Mags/prefilt folders contains symbolic links files.

If you are using conda, before exiting the workflow, the environment needs to be deactivated using the following command::

    conda deactivate

