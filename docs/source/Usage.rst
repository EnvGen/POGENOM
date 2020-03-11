Usage
=====

Import your data
^^^^^^^^^^^^^^^^

The pipeline is capable of analyzing several datasets, mags, and samples. The minimum data required is one dataset, one sample Reads (forward and reverse), and one MAG.

Input data (Reads and MAGs) must be stored in the directory RAW_DATA/, as explained below:

1. READS

Go to the directory Input_POGENOM/RAW_DATA/Reads/::

    cd Input_POGENOM/RAW_DATA/Reads/

Create a directory for each dataset::

    mkdir <dataset_name>

There must be a least one dataset.

Copy or link (recommended) reads to this directory::
                                                                                                                                                           
    cp path/to/reads/* <dataset_name>/. or ln -s path/to/reads/* .

MAke sure that reads name follows the syntax:

forward reads::

    <sample_name><fwd_index><reads_ext> e.g., P6071_505_R1.fq.gz, where sample_name = P6071_505, fwd_index = _R1 , and reads_ext = .fq.gz

reverse reads::

    <sample_name><rev_index><reads_ext> e.g., P6071_505_R1.fq.gz, where sample_name = P6071_505, rev_index = _R2 , and reads_ext = .fq.gz

2. MAGS

Go to the directory cd Input_POGENOM/RAW_DATA/MAGs/::

    cd Input_POGENOM/RAW_DATA/MAGs/

Create a directory for each dataset::

    mkdir <dataset_name>

There must be a least one dataset.

Copy or link (recommended) MAGs to this file::

    cp path/to/MAGs/* <dataset_name>/. or ln -s path/to/MAGs/* <dataset_name>/.

MAke sure that MAG name follows the  syntax::

    <MAG_name><mags_ext>, e.g., P6071_521_bin125.fa, where MAG_name = P6071_521_bin125, and mags_ext = .fa


Run the pipeline
^^^^^^^^^^^^^^^^

If you are using conda, activate the pipeline environment by typing::

    conda activate ip_env

If you are not in the working directory, go there using the command::

    cd ~/Input_POGENOM

1. Set parameters in the config file

In the "Input_POGENOM_config.json" file, set the parameters to be used. It contains the pipeline parameters. Below an example:

{

"workdir":"absolute/path/to/working_directory",
  #it must be an absolute path

"dataset": "name of the dataset to be analyzed",
  #it cannot be empty

"threads": 15,
  #Number of threads

"fraction": "0.15",
  #Fraction of reads to be subsampled when running the pipeline using dataset: "prefilt". 
  Values lower then 15% will difficult the selection of samples with coverage close to the min_coverage 10. 
  Values higher then 25% will considerable increase the size of the directory PREFILT/Reads/ and the running time.  

"mags_ext": ".fa",
  #extention used on your MAGs

"reads_ext": ".fq.gz",
  #extention used on your Reads, after R1. For instance, ".fq.gz" if reads are named "sample_R1.fq.gz, sample_R2.fq.gz"

"fwd_index": "_R1",
  #index used to define foward reads in your sample dataset
"rev_index": "_R2",
  #index used to define reverse reads in your dataset

"bowtie2_params": "--ignore-quals --mp 1,1 --np 1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.05",
  #Mapping parameters. The –score-min then gives the minimum score that is allowed to report an alignment.
  Here, It represents a 95% identity threshold.
  For more information visit http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

"mapqual": 20,
  #Reads mapping quality threshold in BAM files. integer number. Parameter used in samtools view -q {}.

"samtools_view_alignment_extra_filters": "-f 2 -F 1024",
  #filters used after mapping reads.
  Here it selects only paired reads (-f 2) and avoids optical duplicates (-F 1024).                                                                                                                                  If no filters are required, then set an empty string ("samtools_view_alignment_extra_filters": "",)

"min_coverage": 10,
  #minimum genome coverage per sample per MAG (higher than, not included), integer number.
  When dataset: "prefilt", a min_coverage value lower than 10 will select all samples, and the prefilter will be obsolete.

"min_breadth": 40,
  #minimum genome breadth per sample per MAG, integer number.

"min_bsq_for_cov_median_calculation": 15,
  #minimum base quality when counting the number of bases per genome position during coverage calculation. Integer number.

"freebayes_parameters": "-C 4 -p 1 --pooled-continuous --read-max-mismatch-fraction 0.05 --min-alternate-fraction 0.01 -q 15",
  #parameters used during variant calling.
  By default, freebayes exclude duplicates marked as such in alignments.
  If you want to include use tag --use-duplicate-reads and remove "-F 1024" in "samtools_view_alignment_extra_filters".
  # -q --min-base-quality Q. Exclude alleles from analysis if their supporting base quality is less than Q

"vcffilter_qual": "'QUAL > 20'"
  #filtering variant calling.
  Here it removes any sites with an estimated probability of not being polymorphic less than Phred 20 (corresponding to 99% probability of being a real SNP)

}

To access and modify this file, you can use the following command::

    nano config_files/Input_POGENOM_config.json

Modify the required items and save the file. Use Ctrl +x and answer y, to save the modifications and exit the file.

2. Run

The workflow is run with the following command::

    bash Input_POGENOM.sh -d=<absolute path to configFile. Default, d=/from/root/to/config_files/Input_POGENOM_config.json>

If you use the default parameters, the pipeline can be run with the command::

    bash Input_POGENOM.sh

If you need to set a different path to the config file, please do not use relative paths (~/ nor ./)

2.1) A dataset

If you want to run the pipeline on one dataset, please set the corresponding name in the config_file, "dataset": <dataset_name>

2.2) Several datasets

2.2.1) "prefilt"

If you want to run the pipeline on the entire sampling dataset, and only on those MAGs and their corresponding samples with Median coverage higher than a certain threshold (i.e., min_coverage),
please set "prefilt" in the config_file, "dataset": "prefilt" and "fraction": "<your float value, default=0.10>."

When running the pipeline with dataset "prefilt", the created RAW_DATA/Reads/prefilt and RAW_DATA/Mags/prefilt folders contains symbolic links files.

If you are using conda, before exiting the workflow, the environment needs to be deactivated using the following command::

    conda deactivate
