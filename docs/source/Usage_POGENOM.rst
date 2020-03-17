Usage
=====

Minimum input (only genome-wide 𝜋 and FST will be calculated)::

    perl pogenom.pl --vcf_file VCF_FILE --out OUTPUT_FILES_PREFIX --genome_size GENOME_SIZE

or::

    perl pogenom.pl --vcf_file VCF_FILE --out OUTPUT_FILES_PREFIX --gff_file GFF_FILE

or::

    perl pogenom.pl --vcf_file VCF_FILE --out OUTPUT_FILES_PREFIX --fasta_file FASTA_FILE

If a GFF file is provided, gene-wise 𝜋 and gene-wise FST will also be calculated::

    perl pogenom.pl --vcf_file VCF_FILE --out OUTPUT_FILES_PREFIX --gff_file GFF_FILE
    
And if a genetic code file is provided (such as ``standard_genetic_code.txt`` in the POGENOM distribution), amino acid frequencies will be calculated for each codon position in each gene and sample, and gene-wise 𝜋 and FST will be calculated also at the amino acid level. Now also non-synonymous to synonymous polymorphism rates (pN/pS) will be calculated for each gene in each sample::

    perl pogenom.pl --vcf_file VCF_FILE --out OUTPUT_FILES_PREFIX --gff_file GFF_FILE --genetic_code_file GENETIC_CODE_FILE
    

Required arguments
^^^^^^^^^^^^^^^^^^

--vcf_file VCF_FILE                   
Specify VCF file with data from a single or multiple samples.

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

--sample_file SAMPLE_FILE
 Specify file with ids of samples to be included.

--min_count MIN_COUNT
 Specify minimum coverage (integer) for a locus to be included for the sample.

--min_found MIN_FOUND_IN
 Specify minimum number samples (integer) that a locus needs to be present in to be included. If set to 0, it will be set to the number of samples of the VCF file.

--subsample SUBSAMPLE
 Specify coverage level (integer) at which to subsample.

--keep_haplotypes
 If this is used, POGENOM will not split haplotypes into single-nucleotide variants, which is otherwise the default behaviour.

--vcf_format
 Specify VCF file format version. Can be set to freebayes (default) or GATK.
 
--fst_perm FST_PERM         
 Specify number of permutations (integer) for making randomised gene-wise Fst. Without setting this randomised Fst are not generated. Warning: use with care, output files can become huge.

--pi_only                   
 Set this to make POGENOM only calculate and output genome-wide pi (fast).

--help
 To print help message on screen.

