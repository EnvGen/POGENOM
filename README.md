# POGENOM
Population genomics from metagenomes



USAGE (minimum input)

perl pogenom.pl -vcf_file VCF_FILE -o OUTPUT_FILES_PREFIX -genome_size GENOME_SIZE [-h]

OR:

perl pogenom.pl -vcf_file VCF_FILE -o OUTPUT_FILES_PREFIX -gff_file GFF_FILE [-h]



POSITIONAL ARGUMENTS

-vcf_file VCF_FILE          Specify vcf file with data from a single or multiple samples

-o OUTPUT_FILES_PREFIX		Specify the prefix of the output file name(s) (overwrites existing files with same names)

-genome_size GENOME_SIZE    Specify genome size (in bp; integer). Not required if -gff_file is given

-gff_file GFF_FILE          Specify gff file. Either this or genome_size must be given.



OPIONAL ARGUMENTS

-min_count MIN_COUNT        Specify minimum coverage for a locus to be included for the sample

-min_found MIN_FOUND_IN     Specify minimum number samples that a locus need to be present in to be included

-loci_file LOCI_FILE        Specify file with ids for loci to include

-subsample SUBSAMPLE        Specify coverage level at which to subsample

-h							To print help message on screen


