Usage
=====

Minimum input:

Either::

perl pogenom.pl --vcf_file VCF_FILE --out OUTPUT_FILES_PREFIX --genome_size GENOME_SIZE [--help]

Or::

`perl pogenom.pl --vcf_file VCF_FILE --out OUTPUT_FILES_PREFIX --gff_file GFF_FILE [--help]`

Or::

`perl pogenom.pl --vcf_file VCF_FILE --out OUTPUT_FILES_PREFIX --fasta_file FASTA_FILE [--help]`


Required arguments
^^^^^^^^^^^^^^^^^^


--vcf_file VCF_FILE                   
 Specify vcf file with data from a single or multiple samples.

--out OUTPUT_FILES_PREFIX             
 Specify the prefix of the output file name(s) (overwrites existing files with same names).

--genome_size GENOME_SIZE             
 Specify genome size (in bp; integer). Not required if `--gff_file` or `--fasta_file` with genome sequence is given.


Optional arguments
^^^^^^^^^^^^^^^^^^


--gff_file GFF_FILE                   
 Specify gff file. Either this, `--genome_size` or `--fasta_file` must be given.

--fasta_file FASTA_FILE
 Specify fasta file. Either this, `--genome_size` or `--gff_file` must be given

--genetic_code_file GENETIC_CODE_FILE
 Specify genetic code file. E.g. `standard_genetic_code.txt` in the POGENOM distribution.

--loci_file LOCI_FILE
 Specify file with ids of loci to include.

--min_count MIN_COUNT
 Specify minimum coverage for a locus to be included for the sample.

--min_found MIN_FOUND_IN
 Specify minimum number samples that a locus needs to be present in to be included.

--subsample SUBSAMPLE
 Specify coverage level at which to subsample.

--keep_haplotypes
 If this is used, POGENOM will not split haplotypes into single-nucleotide variants, which is otherwise the default behaviour.

--vcf_format
 Specify VCF file format version. Can be set to freebayes (default) or GATK.

--help
 To print help message on screen.

