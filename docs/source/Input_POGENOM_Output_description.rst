Output description
==================

Intermediate files
^^^^^^^^^^^^^^^^^^

A) 01_INDEXING.
 The indexed genome file(s) are stored in this directory. A subdirectory per dataset and genome is created.  Example of filename::

    01_INDEXING/<dataset>/<genome_name>/

B) 02_MAPPING.
 A subdirectory per dataset and genome is created. The BAM files of mapped reads corresponding to each genome and sample are stored here.
 The read group (@RG) information for each BAM file corresponds to the sample name.
 Example of filename::

    02_MAPPING/<dataset>/<genome_name>/<sample_name>_<genome_name>_mpq_<min mapping_quality>_RG_sorted_position.bam

 The Bowtie2 log files are stored in this directory. Example of filename::

    02_MAPPING/log_bowtie2/<dataset>/<genome_name>/<sample_name>_<genome_name>_mpq_<min mapping quality>.log

C) 03_MPILEUP.
 A subdirectory per dataset and genome is created. It contains the samtools mpileup file, in which the coverage per genome position is  stored.
 Example of filename::

    03_MPILEUP/<dataset>/<genome_name>/<sample_name>_<genome_name>_mpq_<min mapping_quality>_bq_<min base_quality>_mpileup_<subsampling_fraction>

D) 04_mergeable.
 Files passing the filter (minimum coverage, minimum breadth) are subsampled (using samtools view) up to user-defined coverage (defined by "subsampling_fraction"), sorted by position and stored in this directory.
 Example of filename::

    04_mergeable/<dataset>/params_<parameters>/<genome_name>/<sample_name>_<genome_name>_RG_sorted_position_subsampled.bam
    where <parameters> is a string of the key parameters used, for instance 'cov_20_bdth_40_subsamp_TRUE_mpq_20_bq_15'.

The corresponding log file for these steps is(are)::

    log_files/<dataset>_genomes_coverage_breadth.log or log_files/<dataset>.<genome_name>.coverage_breadth.log (when "mode_prefilt" used)

E) 05_BAM_merged.
 When the number of BAM files in 04_mergeable/ directory is more than 1, the files are merged into one BAM file per genome. The read group (@RG) information from each BAM file, corresponding to the sample name, is kept.
 Example of filename::

    05_BAM_merged/<dataset>/params_<parameters>/<genome_name>_merged_sorted_position.bam


Intermediate files when "mode_prefilt" is used
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When "mode_prefilt" is used, the suffix "_prefilt" will be added to <dataset> in intermediate files B-E, e.g., 05_BAM_merged/<dataset>_prefilt/<genome_name>_merged_sorted_position.bam

Additionally, the directory ``PREFILT/<dataset>`` is created and contains the subdirectories:

``02_MAPPING``, ``03_MPILEUP``, and ``params_<parameters>``, where <parameters> is a string of the key parameters used, for instance 'cov_20_mpq_20_bq_15_fr_0.15'.

In the subdirectory ``params_<parameters>``, the following files are stored:

Estimated_median_cov_per_sample.tsv, and Selected_samples_genomes.txt (The file names describe their contents).

The directories, 02_MAPPING, and 03_MPILEUP have the same format as described above (Intermediate files). Example of filenames::

    PREFILT/<dataset>/02_MAPPING/params_<parameters>/<genome_name>/<sample_name>_<genome_name>_sorted_position.bam
    where <parameters> is a string of the key parameters used, for instance 'mpq_20_fr_0.15'.

    PREFILT/<dataset>/03_MPILEUP/params_<parameters>/<genome_name>/<sample_name>_<genome_name>_mpileup
    where <parameters> is a string of the key parameters used, for instance 'mpq_20_bq_15_fr_0.15'.

The reads used to generated those files are the Reads subsets, which are stored in the folder ``<temp_sub_Reads_dir>/Reads/<fraction>/``.

The corresponding log file for these steps is ``log_files/samples_filter.log``
The Bowtie2 log files generated when mapping Reads subset, are stored in ``PREFILT/<dataset>/02_MAPPING``. Example of filename::

    PREFILT/<dataset>/02_MAPPING/params_<parameters>/log_bowtie2/<genome_name>/<sample_name>_<genome_name>.log


VCF files
^^^^^^^^^

Variant calling files per genome (input for POGENOM) are stored in the directory 06_VCF.
Example of filename::

    06_VCF/<dataset>/params_<parameters>/<genome_name>.vcf
    where <parameters> is a string of the key parameters used, for instance 'cov_20_bdth_40_subsamp_TRUE_mpq_20_bq_15'.

The list of samples used for the generation of the vcf files can be found in the files ``06_VCF/<dataset>/params_<parameters>/<genome_name>_samples.txt``

When no BAM file passes the filter (coverage and breadth), a vcf file cannot be created.
In this case, the corresponding <genome_name>_samples.txt file will contain the following statement: "The genome <genome_name> has not BAM file that passes the filter breadth and coverage. A vcf file cannot be created."

When "mode_prefilt" is used, the suffix "_prefilt" will be added to <dataset> in VCF files, e.g.,
06_VCF/<dataset>_prefilt/params_<parameters>/<genome_name>.vcf

The corresponding log file for these steps is (are)::

    log_files/<dataset>_genomes_vcf_files.log or log_files/<dataset>.<genome_name>_vcf_files.log (when "mode_prefilt" is used)

Genome size files
^^^^^^^^^^^^^^^^^

The size of the genome (number of bases) is stored in file ``<genome_name>.size``. This value may be used later as input for POGENOM.
