Output description
==================

Intermediate files
^^^^^^^^^^^^^^^^^^

A) 01_INDEXING.
The indexed genome file(s) are stored in this directory. A subdirectory per dataset and genome is created, e.g., "01_INDEXING/<dataset>/<genome_name>/"

B) 02_MAPPING.
A subdirectory per dataset and genome is created. The BAM files of mapped reads corresponding to each genome and sample are stored here.
The read group (@RG) information for each BAM file corresponds to the sample name.
Example of filename::

    02_MAPPING/<dataset>/<genome_name>/<sample_name>_<genome_name>_mpq_<mapping_quality>_RG_sorted_position.bam

C) 03_MPILEUP.
A subdirectory per dataset and genome is created. It contains the samtools mpileup file, in which the coverage per genome position is stored.
Example of filename::

    03_MPILEUP/<dataset>/<genome_name>/<sample_name>_<genome_name>_mpq_<mapping_quality>_mpileup

D) 04_mergeable.
Files passing the filter (minimum coverage, minimum breadth) are subsampled (using samtools view) up to minimum coverage, sorted by position and stored in this directory.
Example of filename::

    04_mergeable/<dataset>/<genome_name>/<sample_name>_<genome_name>_mpq_<mapping_quality>_RG_sorted_position_subsampled.bam

The corresponding log file for these steps is(are)::

    log_files/<dataset>_genomes_coverage_breadth.log or log_files/coverage_breadth_<genome_name>.log (when dataset: "prefilt")

E) 05_BAM_merged.
When the number of BAM files in 04_mergeable/ directory is more than 1, the files are merged into one BAM file per genome. The read group (@RG) information from each BAM file, corresponding to the sample name, is kept.
Example of filename::

    05_BAM_merged/<dataset>/<genome_name>_merged_sorted_position.bam


Intermediate files when dataset: "prefilt"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The directory PREFILT/ is created and contains:

02_MAPPING, 03_MPILEUP

and the files:

Estimated_median_cov_per_sample.tsv, and Selected_samples_genomes.tmp (The file names describe their contents)

The directories, 02_MAPPING, and 03_MPILEUP have the same format as described above (Intermediate files).
The reads used to generated those files are the Reads subsets, which are stored in the folder <temp_sub_Reads_dir>/Reads/.

The corresponding log file for these steps is log_files/samples_filter.log


VCF files
^^^^^^^^^

Variant calling files per MAG are stored in the directory 06_VCF.
Example of filename::

    06_VCF/<dataset>/<genome_name>_mpq_<mapping_quality>.vcf

The list of samples used for the generation of the vcf files can be found in the files <genome_name>_samples.txt

When no BAM file passes the filter (coverage and breadth), a vcf file cannot be created. In this case, the corresponding <genome_name>_samples.txt file will contain the following statement:
"The mag <MAG_name>.fa has not BAM file that passes the filter breadth and coverage. A .vcf file cannot be created."

The corresponding log file for these steps is (are)::

    log_files/<dataset>_genomes_vcf_files.log or log_files/vcf_<genome_name>.log (when dataset: "prefilt")

Genome size files
^^^^^^^^^^^^^^
The size of the genome (number of bases) is stored in file <genome_name>.size. This value may be used later as input for POGENOM.

