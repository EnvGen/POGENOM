Output description
==================

Intermediate files
^^^^^^^^^^^^^^^^^^

A) 01_INDEXING.
The indexed genomes (MAGs) files are stored in this directory. A subdirectory per dataset and MAG is created, e.g., "01_INDEXING/<dataset>/<MAG_name>/"

B) 02_MAPPING.
A subdirectory per dataset and MAG is created. The Bam files corresponding to each MAG and sample reads mapping are stored here.
The reads group (@RG) information for each BAM file corresponds to the sample name.
Example of filename: "02_MAPPING/<dataset>/<MAG_name>/<sample_name>_<MAG_name>_mpq_<mapping_quality>_RG_sorted_position.bam"

C) 03_MPILEUP.
A subdirectory per dataset and MAG is created. It contains the samtools mpileup file, in which the coverage per genome position is stored.
Example of filename: "03_MPILEUP/<dataset>/<MAG_name>/<sample_name>_<MAG_name>_mpq_<mapping_quality>_mpileup"

D) 04_mergeable.
Files passing the filter (minimum coverage, minimum breadth) are subsampled (using samtools view) up to minimum coverage, sorted by position and stored in this directory.
Example of filename: 04_mergeable/<dataset>/<MAG_name>/<sample_name>_<MAG_name>_mpq_<mapping_quality>_RG_sorted_position_subsampled.bam

The corresponding log file for these steps is(are) log_files/<dataset>_MAGs_coverage_breadth.log or log_files/coverage_breadth_<mag_name>.log (when dataset: "prefilt")

E) 05_BAM_merged.
When the number of BAM files in 04_mergeable/ directory is more than 1, the files are merged into one BAM file per MAG. The read group (@RG) information from each BAM file, corresponding to the sample name, is $Example of filename: 05_BAM_merged/<dataset>/<MAG_name>_merged_sorted_position.bam                                                                                                                                 


Intermediate files when dataset: "prefilt"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The directory PREFILT/ is created and contains:

02_MAPPING, 03_MPILEUP, Reads

and the files:

Estimated_median_cov_per_sample.tsv, and Selected_samples_MAGs.tmp (The names of the files describe their containing)

The directories, 02_MAPPING, and 03_MPILEUP have the same format as described above (Intermediate files).
The reads used to generated those files are the Reads subsets, which are stored in the folder PREFILT/Reads/.

The corresponding log file for these steps is log_files/samples_filter.log


VCF files
^^^^^^^^^

Variant calling files per MAG are stored in the directory 06_VCF.
Example of filename: 06_VCF/<dataset>/<MAG_name>_mpq_<mapping_quality>.vcf

The list of samples used for the generation of the vcf files can be found in the files <mag_name>_samples.txt

When no BAM file passes the filter (coverage and breadth), a vcf file cannot be created. In this case, the corresponding <mag_name>_samples.txt file will contain the following statement:
"The mag <MAG_name>.fa has not BAM file that passes the filter breadth and coverage. A .vcf file cannot be created."


MAG_size files
^^^^^^^^^^^^^^
The size of the MAG is stored in file <mag_name>.size (genome size: #bases). This value may be used later as input for POGENOM.

The corresponding log file for these steps is(are) log_files/<dataset>_MAGs_vcf_files.log or log_files/vcf_<mag_name>.log (when dataset: "prefilt")
