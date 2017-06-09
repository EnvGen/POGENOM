#!/usr/bin/perl -w

# pi per gene: jfr med randomiserat dataset: 1) shuffle pi mellan locus (bara loci som har snps), eller 2) shuflle pi över hela genomet (alla loci).

=head1 NAME

pogenom.pl - Calculates nucleotide diveristy (pi) and fixation indices (FST) from VCF file(s)

=head1 USAGE (minimum input)

perl pogenom.pl -vcf_file VCF_FILE -o OUTPUT_FILES_PREFIX -genome_size GENOME_SIZE [-h]

 OR:
 
perl pogenom.pl -vcf_file VCF_FILE -o OUTPUT_FILES_PREFIX -gff_file GFF_FILE [-h]


=head1 POSITIONAL ARGUMENTS

-vcf_file VCF_FILE          Specify vcf file with data from a single or multiple samples

-o OUTPUT_FILES_PREFIX		Specify the prefix of the output file name(s) (overwrites existing files with same names)


=head1 OPIONAL ARGUMENTS

-genome_size GENOME_SIZE    Specify genome size (in bp)
 
-gff_file GFF_FILE          Specify gff file. Either this or genome_size must be given.

-min_count MIN_COUNT        Specify minimum coverage for a locus to be included for the sample
 
-min_found MIN_FOUND_IN     Specify minimum number samples that a locus need to be present in to be included
 
-loci_file LOCI_FILE        Specify file with ids for loci to include

-subsample SUBSAMPLE        Specify level at which to subsmaple

-h							Print this help message

[Press q to close this help message]

=cut

use Getopt::Long;
use List::MoreUtils qw(uniq);

$genome_size = undef;
$min_found_in = 1;
$min_count = 2;
$vcf_file = undef;
$outprefix = undef;
$subsample_level = undef;
$reference = undef;
$keep_haplotypes = undef;
$use_pseudocounts = undef;
$subsample = undef;
$vcf_version = "4.1";
$na_if_missing_loci = 1;

&GetOptions('vcf_file=s' => \$vcf_file, 'vcf_version=s' => \$vcf_version, 'gff_file=s' => \$gff_file, 'genetic_code_file=s' => \$genetic_code_file, 'o=s' => \$outprefix, 'min_count=i' => \$min_count, 'min_found=i' => \$min_found_in, 'sub=i' => \$subsample_level, 'ref=s' => \$reference, 'genome_size=i' => \$genome_size, 'keep_haplotypes!' => \$keep_haplotypes, 'loci_file=s' => \$loci_file, 'subsample=s' => \$subsample, 'use_pseudocounts' => \$use_pseudocounts, 'h!' => \$help);

if (!$outprefix) {
    system ('perldoc', $0);
    exit;
}
if ($help) {
    system ('perldoc', $0);
    exit;
}
if (!$vcf_file) {
    system ('perldoc', $0);
    exit;
}
if (!$genome_size and !$gff_file) {
    system ('perldoc', $0);
    exit;
}
if ($min_count < 2) {
    print"Error: min_count cannot be set to a number < 2\n";
    exit;
}
if ($vcf_version eq "4.1") {
    #fileformat=VCFv4.1
    $tot_count_ix = 1;
    $ref_count_ix = 2;
    $alt_count_ix = 4;
} elsif ($vcf_version eq "4.2") {
    #fileformat=VCFv4.2
    $tot_count_ix = 1;
    $ref_count_ix = 3;
    $alt_count_ix = 5;
} else {
    print"Error: Unrecognized vcf_version. Should be either 4.1 (default) or 4.2\n";
    exit;
}


####################

print"\n### Running pogenom ###\n";
if ($genome_size) {
    print"genome_size set to $genome_size\n";
} else {
    print"genome_size will be calculated from GFF\n";
}
print"min_count set to $min_count\n";
print"min_found set to $min_found_in\n";
if ($subsample) {
    print"Subsampling will be done with $subsample reads per sample\n";
}
if ($loci_file) {
    print"Analysis restricted to loci in file: $loci_file\n";
    &read_loci_to_include;
}
print"\n### Read variant data ###\n";
if ($keep_haplotypes) {
    &get_snp_data_combined_vcf;
} else {
    &get_snp_data_combined_vcf_split_haplotypes;
}
if ($subsample) {
    &subsample_allele_counts;
    #&subsample_allele_counts_dupl;
}
if ($gff_file) {
    print"\n### Reading GFF file ###\n";
    &read_gff;
}
if ($genetic_code_file) {
    print"\n### Reading Genetic Code file ###\n";
    &read_genetic_code;
}
print"\n### Calculating Nucleotide Diversity (pi) ###\n";
&calc_pi;
if ($gff_file) {
    print"\n### Calculating Gene-wise Nucleotide Diversity (pi) ###\n";
    &calc_per_gene_pi;
    if ($genetic_code_file) {
        print"\n### Calculating Gene-wise Aminoacid Diversity (aa-pi) ###\n";
        &calc_per_gene_aminoacid_pi;
        print"\n### Calculating Gene-wise pN/pS ###\n";
        &calc_pN_pS;
    }
}
if (@samples > 1) {
    print"\n### Calculating Fixation Index (FST) ###\n";
    &calc_fst;
    if ($gff_file) {
        print"\n### Calculating Gene-wise Fixation Index (FST) ###\n";
        &calc_per_gene_fst;
        if ($genetic_code_file) {
            print"\n### Calculating Gene-wise Aminoacid Fixation Index (aa-FST) ###\n";
            &calc_per_gene_aminoacid_fst;
        }
    }
}
print"\n### Printing results to files ###\n";
&print_output_to_file;
print"\n### Finished pogenom succesfully ###\n\n";

####################

sub read_loci_to_include {
    open (INFILE, "$loci_file") || die ("Can't open $loci_file");
    while (<INFILE>) {
        chomp;
        $row = $_;
        @fields = split(/\t/, $row);
        $locus = $fields[0]."|".$fields[1];
        $include_locus{$locus} = 1;
        #print"$locus\n";
    }
    $temp = keys %include_locus;
    print"Number of loci specified in file: $temp\n";
}

sub read_gff {
    $fasta_started = 0;
    $seq = "";
    open (INFILE, "$gff_file") || die ("Can't open $gff_file");
    print"Reading $gff_file\n";
    while (<INFILE>) {
        chomp;
        $row = $_;
        #print"$row\n";
        if ($fasta_started == 1) {
            if (substr($row, 0, 1) eq ">") {
                if ($seq ne "") {
                    $contig_seq{$contig} = $seq;
                }
                $contig = $row;
                substr($contig, 0, 1) = "";
                $seq = "";
            } else {
                $seq = $seq.$_;
                $contig_seq{$contig} = $seq;
            }
        }
        elsif (substr($row, 0, 7) eq "##FASTA") {
            $fasta_started = 1;
        } else {
            next if (substr($row, 0, 1) eq "#");
            @fields = split(/\t/, $row);
            next if (@fields == 1);
            $contig = $fields[0];
            $start = $fields[3];
            $end = $fields[4];
            $strand = $fields[6];
            $end = $end - 3 if ($strand eq "+"); # to exclude the stop codon
            $start = $start + 3 if ($strand eq "-"); # to exclude the stop codon
            $annotation = $fields[8];
            @subfields = split(/;/, $annotation);
            $gene = $subfields[0];
            if (defined $gene_start{$gene}) {
                print"Error: Non-unique gene identifiers in GFF file. The program will exit without finishing.\n";
            }
            $gene_start{$gene} = $start;
            $gene_end{$gene} = $end;
            $gene_strand{$gene} = $strand;
            $gene_length{$gene} = $gene_end{$gene} - $gene_start{$gene} + 1;
            $gene_contig{$gene} = $contig;
            if (!defined $contig_genes{$contig}) {
                push(@contigs, $contig);
            }
            push(@{ $contig_genes{$contig} }, $gene);
            for ($pos = $start; $pos < $end; $pos++) {
                $locus = $contig."|".$pos;
                #print"$locus#\n";
                if (defined $locus_found{$locus}) {
                    if ($locus_found{$locus} >= $min_found_in) {
                        $gene_locus{$gene}{$locus} = 1;
                    }
                }
            }
        }
    }
    close (INFILE);
    if ($fasta_started == 1) {
        $temp_genome_size = 0;
        foreach $contig (@contigs) {
            $temp_genome_size = $temp_genome_size + length($contig_seq{$contig});
            @genes = @{ $contig_genes{$contig} };
            foreach $gene (@genes) {
                if ($gene_strand{$gene} eq "+") {
                    $gene_seq{$gene} = substr($contig_seq{$contig}, ($gene_start{$gene} - 1), $gene_length{$gene});
                    #$peptide = "";
                    #for ($j = 0; $j < $gene_length{$gene}; $j=$j+3) {
                    #    $codon = substr($gene_seq{$gene}, $j, 3);
                    #    $peptide = $peptide.$codon_aminoacid{$codon};
                    #}
                    #print">$gene\n$gene_seq{$gene}\n";
                    #print">$gene\n$peptide\n";
                }
                if ($gene_strand{$gene} eq "-") {
                    $gene_seq{$gene} = substr($contig_seq{$contig}, ($gene_start{$gene} - 1), $gene_length{$gene});
                    $gene_seq{$gene} = &make_revcomp($gene_seq{$gene});
                    #$peptide = "";
                    #for ($j = 0; $j < $gene_length{$gene}; $j=$j+3) {
                    #    $codon = substr($gene_seq{$gene}, $j, 3);
                    #    $peptide = $peptide.$codon_aminoacid{$codon};
                    #}
                    #print">$gene\n$gene_seq{$gene}\n";
                    #print">$gene\n$peptide\n";
                }
            }
        }
        if (!defined $genome_size) {
            $genome_size = $temp_genome_size;
            print"Genome size calculated from GFF to $genome_size bp\n";
        }
    }
}

sub read_genetic_code {
    open(INFILE, "$genetic_code_file") || die ("could not open $genetic_code_file");
    print"Reading $genetic_code_file\n";
    while (<INFILE>) {
        chomp;
        @fields = split(/\t/);
        $codon_aminoacid{$fields[0]} = $fields[1];
        #print"#$fields[0]#$fields[1]#\n";
    }
    close (INFILE);
}

sub get_snp_data_combined_vcf {
    %nt = ('A', 1, 'T', 1, 'C', 1, 'G', 1);
    @samples = ();
    @samples_plus = ();
    open (INFILE, "$vcf_file") || die ("Can't open $vcf_file");
    print"Reading $vcf_file\n";
    while (<INFILE>) {
        #chomp;
        $_ =~ s/\R//g;
        $row = $_;
        @fields = split(/\t/, $row);
        if (substr($row, 0, 6) eq "#CHROM") {
            for ($i = 9; $i < @fields; $i++) {
                $samples[$i - 9] = $fields[$i];
            }
            @samples_plus = (@samples, 'All_samples_combined');
        }
        next if (substr($row, 0, 1) eq "#");
        #next if (@fields == 1);
        $contig = $fields[0];
        $pos = $fields[1];
        $locus = $contig."|".$pos;
        #print"$locus $locus_found{$locus} $min_found_in\n";
        if ($loci_file) {
            #print"$locus\n";
            next if (!defined $include_locus{$locus});
        }
        $locus_found{$locus} = 0;
        for ($i = 9; $i < @fields; $i++) {
            @subfields = split(/:/, $fields[$i]);
            if (@subfields == 7) {
                $tot_count = $subfields[1];
                if ($tot_count ne ".") {
                    if ($tot_count >= $min_count) {
                        $locus_found{$locus}++
                        #print"$locus_found{$locus}\n";
                    }
                }
            }
        }
        next if ($locus_found{$locus} < $min_found_in);
        $contig_pos{$contig}{$pos} = 1;
        $ref = $fields[3];
        @alt = split(/,/, $fields[4]);
        $sample_locus_totcount{'All_samples_combined'}{$locus} = 0;
        $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$ref} = 0;
        for ($i = 9; $i < @fields; $i++) {
            $sample = $samples[$i - 9];
            @subfields = split(/:/, $fields[$i]);
            next if (@subfields == 1);
            $ref_count = $subfields[$ref_count_ix];
            $tot_count = $subfields[$tot_count_ix];
            next if ($tot_count eq ".");
            next if ($tot_count < $min_count);
            @alt_count = split(/,/, $subfields[$alt_count_ix]);
            $sample_locus_ref{$sample}{$locus} = $ref;
            $sample_locus_totcount{$sample}{$locus} = $tot_count;
            $sample_locus_allel_counts{$sample}{$locus}{$ref} = $ref_count;
            $sample_locus_totcount{'All_samples_combined'}{$locus} = $sample_locus_totcount{'All_samples_combined'}{$locus} + $tot_count;
            $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$ref} = $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$ref} + $ref_count;
            #print"Ref: #$ref# $ref_count\n";
            for ($j = 0; $j < @alt; $j++) {
                $sample_locus_allel_counts{$sample}{$locus}{$alt[$j]} = $alt_count[$j];
                $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$alt[$j]} = 0 if (!defined $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$alt[$j]});
                $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$alt[$j]} = $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$alt[$j]} + $alt_count[$j];
                #print"Alt: #$alt[$j]# $alt_count[$j]\n";
                #print"$sample_locus_allel_counts{$sample}{$locus}{$alt[$j]}\n";
            }
        }
        
        $contig_pos{$contig}{$pos} = 0 if (!defined $nt{$ref});
        foreach $string (@alt) {
            $contig_pos{$contig}{$pos} = 0 if (!defined $nt{$string});
        }
    }
    foreach $locus (keys %locus_found) {
        $found[$locus_found{$locus}]++;
    }
    for ($i = 1; $i < @found; $i++) {
        if (defined $found[$i]) {
            print"Number of loci found $i times: $found[$i]\n"
        } else {
            print"Number of loci found $i times: 0\n"
        }
    }
}

sub get_snp_data_combined_vcf_split_haplotypes {
    %nt = ('A', 1, 'T', 1, 'C', 1, 'G', 1);
    @samples = ();
    @samples_plus = ();
    open (INFILE, "$vcf_file") || die ("Can't open $vcf_file");
    print"Reading $vcf_file\n";
    print"Haplotypes will be split into individual bases\n";
    while (<INFILE>) {
        #chomp;
        $_ =~ s/\R//g;
        $row = $_;
        @fields = split(/\t/, $row);
        if (substr($row, 0, 6) eq "#CHROM") {
            for ($i = 9; $i < @fields; $i++) {
                $samples[$i - 9] = $fields[$i];
                #print"$fields[$i]\n";
            }
            @samples_plus = (@samples, 'All_samples_combined');
        }
        next if (substr($row, 0, 1) eq "#");
        #next if (@fields == 1);
        $contig = $fields[0];
        $pos = $fields[1];
        $ref = $fields[3];
        @alt = split(/,/, $fields[4]);
        $indel = undef;
        for ($i = 0; $i < @alt; $i++) {
            if (length($alt[$i]) != length($ref)) {
                $indel = 1;
            }
        }
        next if $indel;
        #print"Ref: $ref Alt: @alt\n";
        for ($i = 0; $i < length($ref); $i++) {
            $modpos = $pos + $i;
            $locus = $contig."|".$modpos;
            #print"#$locus#\n";
            if ($loci_file) {
                next if (!defined $include_locus{$locus});
                #print"$locus\n";
            }
            $locus_found{$locus} = 0;
            for ($j = 9; $j < @fields; $j++) {
                @subfields = split(/:/, $fields[$j]);
                if (@subfields > 1) {
                    $tot_count = $subfields[1];
                    if ($tot_count ne ".") {
                        if ($tot_count >= $min_count) {
                            $locus_found{$locus}++;
                            #print"$locus_found{$locus}\n";
                        }
                    }
                }
            }
            next if ($locus_found{$locus} < $min_found_in);
            $contig_pos{$contig}{$modpos} = 1;
            $subref = substr($ref, $i, 1);
            $contig_pos{$contig}{$modpos} = 0 if (!defined $nt{$subref});
            $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$subref} = 0;
            $sample_locus_totcount{'All_samples_combined'}{$locus} = 0;
            for ($j = 9; $j < @fields; $j++) {
                $sample = $samples[$j - 9];
                @subfields = split(/:/, $fields[$j]);
                next if (@subfields == 1);
                $ref_count = $subfields[$ref_count_ix];
                $tot_count = $subfields[$tot_count_ix];
                next if ($tot_count eq ".");
                next if ($tot_count < $min_count);
                @alt_count = split(/,/, $subfields[$alt_count_ix]);
                $sample_locus_ref{$sample}{$locus} = $subref;
                $sample_locus_totcount{$sample}{$locus} = $tot_count;
                $sample_locus_allel_counts{$sample}{$locus}{$subref} = $ref_count;
                $sample_locus_totcount{'All_samples_combined'}{$locus} = $sample_locus_totcount{'All_samples_combined'}{$locus} + $tot_count;
                $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$subref} = $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$subref} + $ref_count;
                #print"    Sample: $sample Subpos: $i Modpos: $modpos\n";
                #print"        Subref: $subref $ref_count\n";
                for ($k = 0; $k < @alt; $k++) {
                    $subalt = substr($alt[$k], $i, 1);
                    if (defined $sample_locus_allel_counts{$sample}{$locus}{$subalt}) {
                        $sample_locus_allel_counts{$sample}{$locus}{$subalt} = $sample_locus_allel_counts{$sample}{$locus}{$subalt} + $alt_count[$k];
                    } else {
                        $sample_locus_allel_counts{$sample}{$locus}{$subalt} = $alt_count[$k];
                    }
                    $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$subalt} = 0 if (!defined $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$subalt});
                    $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$subalt} = $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$subalt} + $alt_count[$k];
                    $contig_pos{$contig}{$modpos} = 0 if (!defined $nt{$subalt});
                    #print"        Subalt: $subalt $alt_count[$k] $sample_locus_allel_counts{$sample}{$locus}{$subalt}\n";
                }
                
            }
        }
    }
    foreach $locus (keys %locus_found) {
        $found[$locus_found{$locus}]++;
    }
    $temp = 0;
    for ($i = 1; $i < @found; $i++) {
        if (defined $found[$i]) {
            $temp = 1;
            print"Number of loci found $i times: $found[$i]\n"
        } else {
            print"Number of loci found $i times: 0\n"
        }
    }
    #if ($temp == 0) {
    if ((keys %sample_locus_allel_counts) == 0) {
        print"Zero loci found fulfilling criteria\n";
        print"No files will be generated\n";
        exit;
    }
}

sub subsample_allele_counts {
    foreach $sample (@samples_plus) {
        foreach $locus (keys %{$sample_locus_allel_counts{$sample}}) {
            if ($sample_locus_totcount{$sample}{$locus} > $subsample) {
                @alleles = (keys %{$sample_locus_allel_counts{$sample}{$locus}});
                @temp = ();
                foreach $allele (@alleles) {
                    $count = $sample_locus_allel_counts{$sample}{$locus}{$allele};
                    for ($i = 0; $i < $count; $i++) {
                        push(@temp, $allele);
                    }
                    $sample_locus_allel_counts{$sample}{$locus}{$allele} = 0;
                }
                for ($i = 0; $i < $subsample; $i++) {
                    $ix = int(rand(@temp));
                    $sample_locus_allel_counts{$sample}{$locus}{$temp[$ix]}++;
                }
                $sample_locus_totcount{$sample}{$locus} = $subsample;
            }
        }
    }
}

sub subsample_allele_counts_dupl {
    foreach $sample (@samples_plus) {
        $sample_dupl = $sample."_duplicate";
        push(@samples_plus_dupl, $sample);
        push(@samples_plus_dupl, $sample_dupl);
        foreach $locus (keys %{$sample_locus_allel_counts{$sample}}) {
            if ($sample_locus_totcount{$sample}{$locus} > $subsample) {
                @alleles = (keys %{$sample_locus_allel_counts{$sample}{$locus}});
                @temp = ();
                foreach $allele (@alleles) {
                    $count = $sample_locus_allel_counts{$sample}{$locus}{$allele};
                    for ($i = 0; $i < $count; $i++) {
                        push(@temp, $allele);
                    }
                    $sample_locus_allel_counts{$sample}{$locus}{$allele} = 0;
                    $sample_locus_allel_counts{$sample_dupl}{$locus}{$allele} = 0;
                }
                
                for ($i = 0; $i < $subsample; $i++) {
                    $ix = int(rand(@temp));
                    $sample_locus_allel_counts{$sample}{$locus}{$temp[$ix]}++;
                }
                for ($i = 0; $i < $subsample; $i++) {
                    $ix = int(rand(@temp));
                    $sample_locus_allel_counts{$sample_dupl}{$locus}{$temp[$ix]}++;
                }
                $sample_locus_totcount{$sample}{$locus} = $subsample;
                $sample_locus_totcount{$sample_dupl}{$locus} = $subsample;
            }
        }
    }
    @samples = @samples_plus = @samples_plus_dupl;
}

sub calc_pi {
    print "Sample\tpi\tTotal_num_loci\tTotal_num_alleles\tAverage_depth\n";
    foreach $sample (@samples_plus) {
        $intra_pi = 0;
        $tot_alleles = 0;
        $av_count = 0;
        @loci = (keys %{$sample_locus_allel_counts{$sample}});
        foreach $locus (@loci) {
            $tot_count = $sample_locus_totcount{$sample}{$locus};
            #print"Locus totcount $locus $tot_count\n";
            @alleles = (keys %{$sample_locus_allel_counts{$sample}{$locus}});
            #print"Alleles: @alleles\n";
            $av_count = $av_count + $tot_count;
            $locus_intra_pi = 0;
            for ($i = 0; $i < @alleles; $i++) {
                #print"Allele 1: $alleles[$i]\n";
                $counts_1 = $sample_locus_allel_counts{$sample}{$locus}{$alleles[$i]};
                #print"Counts1: $counts_1\n";
                $tot_alleles++ if ($counts_1 > 0);
                for ($j = 0; $j < @alleles; $j++) {
                    next if ($alleles[$i] eq $alleles[$j]);
                    #print"Allele 2: $alleles[$j]\n";
                    $counts_2 = $sample_locus_allel_counts{$sample}{$locus}{$alleles[$j]};
                    #print"Counts2: $counts_2\n";
                    $locus_intra_pi = $locus_intra_pi + ($counts_1/$tot_count)*($counts_2/($tot_count - 1)); # According to Schloising
                    #$locus_intra_pi = $locus_intra_pi + ($counts_1/$tot_count)*($counts_2/($tot_count)); # According to our initial version
                }
            }
            #print"$locus $locus_intra_pi\n\n";
            $sample_locus_pi{$sample}{$locus} = $locus_intra_pi;
            $intra_pi = $intra_pi + $locus_intra_pi;
        }
        $intra_pi = $intra_pi/$genome_size;
        $sample_pi{$sample} = $intra_pi;
        $av_count = $av_count/@loci;
        $sample_avcount{$sample} = $av_count;
        $sample_totalleles{$sample} = $tot_alleles;
        $num_loci{$sample} = @loci;
        print "$sample\t$intra_pi\t$num_loci{$sample}\t$tot_alleles\t$av_count\n";
    }
}

sub calc_per_gene_pi {
    foreach $contig (@contigs) {
        #print "$contig\n";
        @genes = @{ $contig_genes{$contig} };
        #print "@genes\n"; die;
        foreach $gene (@genes) {
            foreach $sample (@samples_plus) {
                $intra_pi = 0;
                $missing_loci = 0;
                foreach $locus (keys %{$gene_locus{$gene}}) {
                    #print "$locus/n";
                    $missing_loci = 1 if (!defined $sample_locus_pi{$sample}{$locus});
                    next if (!defined $sample_locus_pi{$sample}{$locus});
                    $intra_pi = $intra_pi + $sample_locus_pi{$sample}{$locus};
                }
                $intra_pi = $intra_pi/$gene_length{$gene};
                $sample_gene_pi{$sample}{$gene} = $intra_pi;
                if ($na_if_missing_loci == 1 && $missing_loci > 0) {
                    $sample_gene_pi{$sample}{$gene} = "NA";
                }
            }
        }
    }
}

sub calc_per_gene_aminoacid_pi {
    foreach $sample (@samples_plus) {
        $tot_peptides = 0;
        foreach $contig (@contigs) {
            @genes = @{ $contig_genes{$contig} };
            foreach $gene (@genes) {
                $intra_pi = 0;
                $missing_loci = 0;
                #print"\n$gene\n";
                foreach $locus (keys %{$gene_locus{$gene}}) {
                    $missing_loci = 1 if (!defined $sample_locus_pi{$sample}{$locus});
                    next if (!defined $sample_locus_pi{$sample}{$locus});
                    @fields = split(/\|/, $locus);
                    $contig = $fields[0];
                    $pos = $fields[1];
                    $tot_count = $sample_locus_totcount{$sample}{$locus};
                    @alleles = (keys %{$sample_locus_allel_counts{$sample}{$locus}});
                    #print"Alleles: @alleles\n";
                    for ($i = 0; $i < @alleles; $i++) {
                        $mod_contig_seq = $contig_seq{$contig};
                        substr($mod_contig_seq, ($pos-1), length($alleles[$i])) = $alleles[$i];
                        $mod_gene_seq = substr($mod_contig_seq, ($gene_start{$gene} - 1), $gene_length{$gene}); # OBS!!!
                        $mod_gene_seq = &make_revcomp($mod_gene_seq) if ($gene_strand{$gene} eq "-");
                        $peptide = "";
                        for ($j = 0; $j < $gene_length{$gene}; $j=$j+3) {
                            $codon = substr($mod_gene_seq, $j, 3);
                            $peptide = $peptide.$codon_aminoacid{$codon};
                        }
                        if (defined $sample_locus_gene_peptide_counts{$sample}{$locus}{$gene}{$peptide}) {
                            $sample_locus_gene_peptide_counts{$sample}{$locus}{$gene}{$peptide} = $sample_locus_gene_peptide_counts{$sample}{$locus}{$gene}{$peptide} + $sample_locus_allel_counts{$sample}{$locus}{$alleles[$i]};
                        } else {
                            $sample_locus_gene_peptide_counts{$sample}{$locus}{$gene}{$peptide} = $sample_locus_allel_counts{$sample}{$locus}{$alleles[$i]};
                        }
                        #print"\nLocus: $locus\nGene: $gene\nSample: $sample\nAllele: $alleles[$i]\nStrand: $gene_strand{$gene}\n>Gene_seq:\n$gene_seq{$gene}\n\n>Mod_seq:\n$mod_gene_seq\n\n>Mod_peptide:\n$peptide\n";
                    }
                    $locus_intra_pi = 0;
                    @peptides = (keys %{$sample_locus_gene_peptide_counts{$sample}{$locus}{$gene}});
                    for ($i = 0; $i < @peptides; $i++) {
                        #print"\nPeptide 1: $peptides[$i]\n";
                        $counts_1 = $sample_locus_gene_peptide_counts{$sample}{$locus}{$gene}{$peptides[$i]};
                        #print"Counts1: $counts_1\n";
                        $tot_peptides++ if ($counts_1 > 0);
                        for ($j = 0; $j < @peptides; $j++) {
                            next if ($peptides[$i] eq $peptides[$j]);
                            #print"Peptide 2: $peptides[$j]\n";
                            $counts_2 = $sample_locus_gene_peptide_counts{$sample}{$locus}{$gene}{$peptides[$j]};
                            #print"Counts2: $counts_2\n";
                            $locus_intra_pi = $locus_intra_pi + ($counts_1/$tot_count)*($counts_2/($tot_count - 1)); # According to Schloising
                            #$locus_intra_pi = $locus_intra_pi + ($counts_1/$tot_count)*($counts_2/($tot_count)); # According to our initial version
                        }
                    }
                    #print"Sample: $sample $locus locus_intra_pi: $locus_intra_pi\n";
                    $sample_locus_aminoacid_pi{$sample}{$locus} = $locus_intra_pi;
                    $intra_pi = $intra_pi + $locus_intra_pi;
                    #print"Sample: $sample $locus intra_pi: $intra_pi\n\n";
                }
                $intra_pi = $intra_pi/$gene_length{$gene};
                $sample_gene_aminoacid_pi{$sample}{$gene} = $intra_pi;
                #$sample_totpeptides{$sample} = $tot_peptides;
                #print"Sample: $sample $locus intra_pi: $intra_pi\n\n";
                if ($na_if_missing_loci == 1 && $missing_loci == 1) {
                    $sample_gene_aminoacid_pi{$sample}{$gene} = "NA";
                }
            }
        }
    }
}

sub calc_fst { # is it more logical to only include the loci where any inter-pi exist for the sample pair, and calculate mean fst for these loci?
    print "Sample1\tSample2\tpi_1\tpi_2\tpi_1-2\tfst\n";
    for ($ix1 = 0; $ix1 < @samples; $ix1++) {
        $sample1 = $samples[$ix1];
        #print"Sample1 $sample1\n";
        for ($ix2 = $ix1 + 1; $ix2 < @samples; $ix2++) {
        #for ($ix2 = $ix1 + 0; $ix2 < @samples; $ix2++) {
            $sample2 = $samples[$ix2];
            #print"Sample2 $sample2\n";
            $inter_pi = 0;
            $sample1_pi = $sample2_pi = 0; # intra-pi based on only loci shared by both samples
            @loci = (keys %{$sample_locus_allel_counts{$sample1}});
            #$temp = @loci;
            #print"sample1: $sample1 num_loci: $temp\n"; die;
            #print"sample1: $sample1 loci: @loci\n"; die;
            foreach $locus (@loci) {
                next if (!defined $sample_locus_allel_counts{$sample2}{$locus});
                #print"Locus $locus\n";
                @alleles1 = (keys %{$sample_locus_allel_counts{$sample1}{$locus}});
                @alleles2 = (keys %{$sample_locus_allel_counts{$sample2}{$locus}});
                $tot_count1 = $sample_locus_totcount{$sample1}{$locus};
                $tot_count2 = $sample_locus_totcount{$sample2}{$locus};
                #print"Alleles: @alleles1 @alleles2\n";
                for ($i = 0; $i < @alleles1; $i++) {
                    #print"Allele 1: $alleles1[$i]\n";
                    $counts_1 = $sample_locus_allel_counts{$sample1}{$locus}{$alleles1[$i]};
                    #print"Counts1: $counts_1\n";
                    for ($j = 0; $j < @alleles2; $j++) {
                        next if ($alleles1[$i] eq $alleles2[$j]);
                        #print"Allele 2: $alleles2[$j]\n";
                        $counts_2 = $sample_locus_allel_counts{$sample2}{$locus}{$alleles2[$j]};
                        #print"Counts2: $counts_2\n";
                        $inter_pi = $inter_pi + ($counts_1/$tot_count1)*($counts_2/$tot_count2);
                    }
                }
                $sample1_pi = $sample1_pi + $sample_locus_pi{$sample1}{$locus};
                $sample2_pi = $sample2_pi + $sample_locus_pi{$sample2}{$locus};
                
                #$sample1_pi = $sample1_pi;
                #$sample2_pi = $sample2_pi;
            }
            $inter_pi = $inter_pi/$genome_size;
            $sample_sample_pi{$sample1}{$sample2} = $inter_pi;
            $sample_sample_pi{$sample2}{$sample1} = $inter_pi;
            $sample1_pi = $sample1_pi/$genome_size;
            $sample2_pi = $sample2_pi/$genome_size;
            #$fst = 1 - 0.5*($sample_pi{$sample1} + $sample_pi{$sample2})/$inter_pi; # intra pi values only based on all loci
            $fst = 1 - 0.5*($sample1_pi + $sample2_pi)/$inter_pi; # intra pi values only based on shared loci
            $fst = sprintf("%.4f", $fst);
            $sample_sample_fst{$sample1}{$sample2} = $fst;
            $sample_sample_fst{$sample2}{$sample1} = $fst;
            print"$sample1\t$sample2\t$sample_pi{$sample1}\t$sample_pi{$sample2}\t$inter_pi\t$fst\n";
        }
    }
}

sub calc_per_gene_fst {
    foreach $contig (@contigs) {
        #print "Contig $contig\n";
        @genes = @{ $contig_genes{$contig} };
        foreach $gene (@genes) {
            #print "Gene $gene\n";
            @loci = (keys %{$gene_locus{$gene}});
            next if (@loci == 0);
            for ($ix1 = 0; $ix1 < @samples; $ix1++) {
                $sample1 = $samples[$ix1];
                #print"Sample1 $sample1\n";
                for ($ix2 = $ix1 + 1; $ix2 < @samples; $ix2++) {
                    $sample2 = $samples[$ix2];
                    #print"Sample2 $sample2\n";
                    $inter_pi = 0;
                    $sample1_pi = $sample2_pi = 0;
                    foreach $locus (@loci) {
                        #print"Locus $locus\n";
                        next if (!defined $sample_locus_allel_counts{$sample1}{$locus});
                        next if (!defined $sample_locus_allel_counts{$sample2}{$locus});
                        @alleles1 = (keys %{$sample_locus_allel_counts{$sample1}{$locus}});
                        @alleles2 = (keys %{$sample_locus_allel_counts{$sample2}{$locus}});
                        $tot_count1 = $sample_locus_totcount{$sample1}{$locus};
                        $tot_count2 = $sample_locus_totcount{$sample2}{$locus};
                        #print"Alleles1: @alleles1 Alleles2: @alleles2\n";
                        for ($i = 0; $i < @alleles1; $i++) {
                            #print"Allele 1: $alleles1[$i]\n";
                            $counts_1 = $sample_locus_allel_counts{$sample1}{$locus}{$alleles1[$i]};
                            #print"Counts1: $counts_1\n";
                            for ($j = 0; $j < @alleles2; $j++) {
                                next if ($alleles1[$i] eq $alleles2[$j]);
                                #print"Allele 2: $alleles2[$j]\n";
                                $counts_2 = $sample_locus_allel_counts{$sample2}{$locus}{$alleles2[$j]};
                                #print"Counts2: $counts_2\n";
                                $inter_pi = $inter_pi + ($counts_1/$tot_count1)*($counts_2/$tot_count2);
                                #print"inter_pi $inter_pi\n";
                            }
                        }
                        $sample1_pi = $sample1_pi + $sample_locus_pi{$sample1}{$locus};
                        $sample2_pi = $sample2_pi + $sample_locus_pi{$sample2}{$locus};
                    }
                    $inter_pi = $inter_pi/$gene_length{$gene};
                    $sample_sample_gene_pi{$sample1}{$sample2}{$gene} = $inter_pi;
                    $sample_sample_gene_pi{$sample2}{$sample1}{$gene} = $inter_pi;
                    $sample1_pi = $sample1_pi/$gene_length{$gene};
                    $sample2_pi = $sample2_pi/$gene_length{$gene};
                    if ($inter_pi > 0) {
                        #$fst = 1 - 0.5*($sample_gene_pi{$sample1}{$gene} + $sample_gene_pi{$sample2}{$gene})/$inter_pi; # intra pi based on all loci
                        $fst = 1 - 0.5*($sample1_pi + $sample2_pi)/$inter_pi; # intra pi based on only shared loci
                        $fst = sprintf("%.4f", $fst);
                    } else {
                        $fst = NA; # i.e. no intra-pi in any of the two samples and consequently no inter-pi.
                    }
                    $sample_sample_gene_fst{$sample1}{$sample2}{$gene} = $fst;
                    $sample_sample_gene_fst{$sample2}{$sample1}{$gene} = $fst;
                    #print"$gene $sample1 $sample2 $sample_gene_pi{$sample1}{$gene} $sample_gene_pi{$sample2}{$gene} $inter_pi $fst\n";
                }
            }
        }
    }
}

# double check!
sub calc_per_gene_aminoacid_fst {
    foreach $contig (@contigs) {
        #print "Contig $contig\n";
        @genes = @{ $contig_genes{$contig} };
        foreach $gene (@genes) {
            #print "Gene $gene\n";
            @loci = (keys %{$gene_locus{$gene}});
            next if (@loci == 0);
            for ($ix1 = 0; $ix1 < @samples; $ix1++) {
                $sample1 = $samples[$ix1];
                #print"Sample1 $sample1\n";
                for ($ix2 = $ix1 + 1; $ix2 < @samples; $ix2++) {
                    $sample2 = $samples[$ix2];
                    #print"Sample2 $sample2\n";
                    $inter_pi = 0;
                    $sample1_pi = $sample2_pi = 0;
                    foreach $locus (@loci) {
                        #print"Locus $locus\n";
                        next if (!defined $sample_locus_allel_counts{$sample1}{$locus});
                        next if (!defined $sample_locus_allel_counts{$sample2}{$locus});
                        @peptides1 = (keys %{$sample_locus_gene_peptide_counts{$sample1}{$locus}{$gene}});
                        @peptides2 = (keys %{$sample_locus_gene_peptide_counts{$sample2}{$locus}{$gene}});
                        $tot_count1 = $sample_locus_totcount{$sample1}{$locus};
                        $tot_count2 = $sample_locus_totcount{$sample2}{$locus};
                        #print"Peptides1: @peptides1 Peptides2: @peptides2\n";
                        for ($i = 0; $i < @peptides1; $i++) {
                            #print"Peptide 1: $peptides1[$i]\n";
                            $counts_1 = $sample_locus_gene_peptide_counts{$sample1}{$locus}{$gene}{$peptides1[$i]};
                            #print"Counts1: $counts_1\n";
                            for ($j = 0; $j < @peptides2; $j++) {
                                next if ($peptides1[$i] eq $peptides2[$j]);
                                #print"Allele 2: $peptides2[$j]\n";
                                $counts_2 = $sample_locus_gene_peptide_counts{$sample2}{$locus}{$gene}{$peptides2[$j]};
                                #print"Counts2: $counts_2\n";
                                $inter_pi = $inter_pi + ($counts_1/$tot_count1)*($counts_2/$tot_count2);
                            }
                        }
                        $sample1_pi = $sample1_pi + $sample_locus_aminoacid_pi{$sample1}{$locus};
                        $sample2_pi = $sample2_pi + $sample_locus_aminoacid_pi{$sample2}{$locus};
                    }
                    $inter_pi = $inter_pi/$gene_length{$gene};
                    #print"$gene $sample1 $sample2 $inter_pi\n";
                    $sample_sample_gene_aminoacid_pi{$sample1}{$sample2}{$gene} = $inter_pi;
                    $sample_sample_gene_aminoacid_pi{$sample2}{$sample1}{$gene} = $inter_pi;
                    if ($inter_pi > 0) {
                        #$fst = 1 - 0.5*($sample_gene_aminoacid_pi{$sample1}{$gene} + $sample_gene_aminoacid_pi{$sample2}{$gene})/$inter_pi;
                        $fst = 1 - 0.5*($sample1_pi + $sample2_pi)/$inter_pi; # intra pi based on only shared loci
                        $fst = sprintf("%.4f", $fst);
                    } else {
                        $fst = NA; # i.e. no intra-pi in any of the two samples and consequently no inter-pi.
                    }
                    $sample_sample_gene_aminoacid_fst{$sample1}{$sample2}{$gene} = $fst;
                    $sample_sample_gene_aminoacid_fst{$sample2}{$sample1}{$gene} = $fst;
                    #print"$gene $sample1 $sample2 $sample_gene_aminoacid_pi{$sample1}{$gene} $sample_gene_aminoacid_pi{$sample2}{$gene} $inter_pi $fst\n";
                }
            }
        }
    }
}

sub calc_pN_pS {
    foreach $sample (@samples_plus) {
        foreach $contig (@contigs) {
            @genes = @{ $contig_genes{$contig} };
            foreach $gene (@genes) {
                $PN = 0; # non-synonomous polymorphisms
                $PS = 0; # synonomous polymorphisms
                $TN = 0; # possible non-synonomous polymorphisms
                $TS = 0; # possible synonomous polymorphisms
                # $pNpS = ($PN/$TN) / ($PS/$TS) # NA if $PS = 0;
                # pN equals the fraction of polymorphic nonsynonymous sites, pS equals the fraction of polymorphic synonymous sites (Simmons et al, 2008)
                $contig = $gene_contig{$gene};
                $contig_seq = $contig_seq{$contig};
                $gene_seq = substr($contig_seq, ($gene_start{$gene} - 1), $gene_length{$gene});
                # Gene on "+" strand
                
                #@loci = (keys %{$sample_locus_allel_counts{$sample}});
                #$temp = @loci;
                #print"sample: $sample gene: $gene num_loci: $temp\n";
                
                if ($gene_strand{$gene} eq "+") {
                    #print"\n$contig $gene $sample\n$gene_seq\n";
                    for ($i = 3; $i < $gene_length{$gene}; $i = $i + 3) { # skipping the start codon
                        $codon = substr($gene_seq, $i, 3);
                        # modify codon to represent majority sequence
                        for ($j = 0; $j < 3; $j++) {
                            $pos = $gene_start{$gene} + $i + $j;
                            $locus = $contig."|".$pos;
                            $base = substr($gene_seq, $i + $j, 1);
                            if (defined $sample_locus_totcount{$sample}{$locus}) {
                                #print"$locus $gene_start{$gene} $i  $j\n";
                                #print"codon before: $codon aa before: $codon_aminoacid{$codon}\n";
                                @alleles = (keys %{$sample_locus_allel_counts{$sample}{$locus}});
                                #print"Alleles: @alleles\n";
                                #print"Base: $base\n";
                                if (defined $sample_locus_allel_counts{$sample}{$locus}{$base}) {
                                    $base_counts = $sample_locus_allel_counts{$sample}{$locus}{$base};
                                } else {
                                    $base_counts = 0;
                                }
                                for ($k = 0; $k < @alleles; $k++) {
                                    if (length($alleles[$k]) != 1) {
                                        print "\nError: Variant > 1 bp. Violates pN/pS calculation.\n\n"; die;
                                    }
                                    if ($sample_locus_allel_counts{$sample}{$locus}{$alleles[$k]} > $base_counts) {
                                        $base_counts = $sample_locus_allel_counts{$sample}{$locus}{$alleles[$k]};
                                        $base = $alleles[$k];
                                    }
                                }
                                substr($codon, $j, 1) = $base;
                                #print"codon after: $codon aa after: $codon_aminoacid{$codon}\n";
                            }
                        }
                        next if ($codon_aminoacid{$codon} eq "*"); # CHECK!!!
                        # update $PN $TN $PS $TS
                        for ($j = 0; $j < 3; $j++) {
                            $pos = $gene_start{$gene} + $i + $j;
                            $locus = $contig."|".$pos;
                            $base = substr($codon, $j, 1);
                            foreach $nt (keys %nt) {
                                next if ($base eq $nt);
                                $mod_codon = $codon;
                                substr($mod_codon, $j, 1) = $nt;
                                next if ($codon_aminoacid{$mod_codon} eq "*"); # CHECK!!!
                                if ($codon_aminoacid{$mod_codon} eq $codon_aminoacid{$codon}) {
                                    $TS++;
                                    if (defined $sample_locus_allel_counts{$sample}{$locus}) {
                                        if (defined $sample_locus_allel_counts{$sample}{$locus}{$nt}) {
                                            #print"PS $j $nt codon $codon aa $codon_aminoacid{$codon} mod_codon $mod_codon mod_aa $codon_aminoacid{$mod_codon}\n";
                                            $PS++ if ($sample_locus_allel_counts{$sample}{$locus}{$nt} > 0);
                                        }
                                    }
                                } else {
                                    $TN++;
                                    if (defined $sample_locus_allel_counts{$sample}{$locus}) {
                                        if (defined $sample_locus_allel_counts{$sample}{$locus}{$nt}) {
                                            #print"PN $j $nt codon $codon aa $codon_aminoacid{$codon} mod_codon $mod_codon mod_aa $codon_aminoacid{$mod_codon}\n";
                                            $PN++ if ($sample_locus_allel_counts{$sample}{$locus}{$nt} > 0);
                                        }
                                    }
                                }
                            }

                        }
                    }
                }
                
                # Gene on "-" strand
                if ($gene_strand{$gene} eq "-") {
                    $gene_seq = &make_revcomp($gene_seq);
                    #print"$gene_strand{$gene}\n$gene_seq\n\n";
                    #print"\n$contig $gene $sample\n";
                    for ($i = 3; $i < $gene_length{$gene}; $i = $i + 3) { # skipping the start codon
                        $codon = substr($gene_seq, $i, 3);
                        # modify codon to represent majority sequence
                        for ($j = 0; $j < 3; $j++) {
                            $pos = $gene_end{$gene} - ($i + $j); # ok?
                            $locus = $contig."|".$pos;
                            $base = substr($gene_seq, $i + $j, 1);
                            if (defined $sample_locus_totcount{$sample}{$locus}) {
                                #print"$codon\n";
                                $rc_base = &make_revcomp($base);
                                @alleles = (keys %{$sample_locus_allel_counts{$sample}{$locus}});
                                #print"Alleles: @alleles\n";
                                #print"RC_base: $rc_base\n\n";
                                if (defined $sample_locus_allel_counts{$sample}{$locus}{$rc_base}) {
                                    $base_counts = $sample_locus_allel_counts{$sample}{$locus}{$rc_base};
                                } else {
                                    $base_counts = 0;
                                }
                                for ($k = 0; $k < @alleles; $k++) {
                                    if (length($alleles[$k]) != 1) {
                                        print "\nError: Variant > 1 bp. Violates pN/pS calculation.\n\n"; die;
                                    }
                                    if ($sample_locus_allel_counts{$sample}{$locus}{$alleles[$k]} > $base_counts) {
                                        $base_counts = $sample_locus_allel_counts{$sample}{$locus}{$alleles[$k]};
                                        $rc_base = $alleles[$k];
                                    }
                                }
                                $base = &make_revcomp($rc_base);
                                substr($codon, $j, 1) = $base;
                            }
                        }
                        #print"$codon\n";
                        next if ($codon_aminoacid{$codon} eq "*"); # CHECK!!!
                        # upgrade $PN $TN $PS $TS
                        for ($j = 0; $j < 3; $j++) {
                            $pos = $gene_end{$gene} - ($i + $j); # ok?
                            $locus = $contig."|".$pos;
                            $base = substr($codon, $j, 1);
                            foreach $nt (keys (%nt)) {
                                next if ($base eq $nt);
                                $mod_codon = $codon;
                                substr($mod_codon, $j, 1) = $nt;
                                next if ($codon_aminoacid{$mod_codon} eq "*"); # CHECK!!!
                                $rc_nt = &make_revcomp($nt);
                                if ($codon_aminoacid{$mod_codon} eq $codon_aminoacid{$codon}) {
                                    $TS++;
                                    if (defined $sample_locus_allel_counts{$sample}{$locus}) {
                                        if (defined $sample_locus_allel_counts{$sample}{$locus}{$rc_nt}) {
                                            $PS++ if ($sample_locus_allel_counts{$sample}{$locus}{$rc_nt} > 0);
                                        }
                                    }
                                } else {
                                    $TN++;
                                    if (defined $sample_locus_allel_counts{$sample}{$locus}) {
                                        if (defined $sample_locus_allel_counts{$sample}{$locus}{$rc_nt}) {
                                            $PN++ if ($sample_locus_allel_counts{$sample}{$locus}{$rc_nt} > 0);
                                        }
                                    }
                                }
                                
                            }
                        }
                    }
                }

                if ($use_pseudocounts) {
                    if ($TS > 0) {
                        $pNpS = (($PN + 1)/$TN) / (($PS + 1)/$TS);
                    } else {
                        $pNpS = "NA";
                    }
                } elsif ($PS > 0) {
                    $pNpS = ($PN/$TN) / ($PS/$TS);
                } else {
                    $pNpS = "NA";
                }
                $sample_gene_pNpS{$sample}{$gene} = $pNpS;
                #print"$gene_strand{$gene} PN: $PN PS: $PS TN: $TN TS: $TS pNpS: $pNpS\n";
            }
        }
    }
}

sub make_revcomp {
    local($seq) = $_[0];
    local($modseq) = "";
    local($i);
    local($nt);
    local($i);
    for ($i = 0; $i < length($seq); $i++) {
        $nt = substr($seq, $i, 1);
        if ($nt eq "A") {
            $modseq = T.$modseq;
        } elsif ($nt eq "T") {
            $modseq = A.$modseq;
        } elsif ($nt eq "C") {
            $modseq = G.$modseq;
        } elsif ($nt eq "G") {
            $modseq = C.$modseq;
        } else {
            $modseq = $nt.$modseq;
        }
    }
    return ($modseq);
}

sub print_output_to_file {
    print"$outprefix.intradiv.txt\n";
    open (OUT, ">$outprefix.intradiv.txt");
    print OUT "Sample\tIntra_pi\tNum_loci\tTotal_num_alleles\tAverage_read_depth\n";
    foreach $sample (@samples_plus) {
        print OUT "$sample\t$sample_pi{$sample}\t$num_loci{$sample}\t$sample_totalleles{$sample}\t$sample_avcount{$sample}\n";
    }
    close(OUT);
    ###
    print"$outprefix.intradiv-per-locus.txt\n";
    open (OUT, ">$outprefix.intradiv-per-locus.txt");
    print OUT "Contig\tPosition";
    foreach $sample (@samples_plus) {
        print OUT "\t$sample Intra_locus_pi\t$sample Locus_Depth";
    }
    print OUT "\n";
    foreach $contig (keys %contig_pos) {
        @positions = (keys %{$contig_pos{$contig}});
        @positions = sort {$a <=> $b} @positions;
        foreach $pos (@positions) {
            print OUT "$contig\t$pos";
            $locus = $contig."|".$pos;
            foreach $sample (@samples_plus) {
                if (defined $sample_locus_pi{$sample}{$locus}) {
                    print OUT "\t$sample_locus_pi{$sample}{$locus}\t$sample_locus_totcount{$sample}{$locus}";
                } else {
                    print OUT "\tNA\tNA";
                }
            }
            print OUT "\n";
        }
    }
    close(OUT);
    ###
    print"$outprefix.allele-freqs.txt\n";
    @nts = ('A', 'T', 'C', 'G');
    open (OUT, ">$outprefix.allele-freqs.txt");
    print OUT "Contig\tPosition";
    foreach $sample (@samples_plus) {
        print OUT "\t$sample A\t$sample T\t$sample C\t$sample G";
    }
    print OUT "\n";
    foreach $contig (keys %contig_pos) {
        @positions = (keys %{$contig_pos{$contig}});
        @positions = sort {$a <=> $b} @positions;
        foreach $pos (@positions) {
            next if ($contig_pos{$contig}{$pos} == 0);
            #print "$contig\t$pos\t$contig_pos{$contig}{$pos}\n";
            print OUT "$contig\t$pos";
            $locus = $contig."|".$pos;
            foreach $sample (@samples_plus) {
                if (defined $sample_locus_pi{$sample}{$locus}) {
                    foreach $nt (@nts) {
                        $counts = 0;
                        if (defined $sample_locus_allel_counts{$sample}{$locus}{$nt}) {
                            $counts = $sample_locus_allel_counts{$sample}{$locus}{$nt};
                        }
                        print OUT "\t$counts";
                    }
                } else {
                    print OUT "\tNA\tNA\tNA\tNA";
                }
            }
            print OUT "\n";
        }
    }
    close(OUT);
    ###
    if ($gff_file) {
        print"$outprefix.intradiv-per-gene.txt\n";
        open (OUT, ">$outprefix.intradiv-per-gene.txt");
        print OUT "Contig\tGene\tLength\tStart\tEnd\tStrand\tNum_loci";
        foreach $sample (@samples_plus) {
            #print OUT "\t$sample Intra_gene_pi\t$sample Intra_gene_pi_p-value";
            print OUT "\t$sample Intra_gene_pi";
        }
        print OUT "\n";
        foreach $contig (@contigs) {
            #print "$contig\n";
            @genes = @{ $contig_genes{$contig} };
            #print "@genes\n"; die;
            foreach $gene (@genes) {
                $num_loci = (keys %{$gene_locus{$gene}});
                print OUT "$contig\t$gene\t$gene_length{$gene}\t$gene_start{$gene}\t$gene_end{$gene}\t$gene_strand{$gene}\t$num_loci";
                foreach $sample (@samples_plus) {
                    #print OUT "\t$sample_gene_pi{$sample}{$gene}\tNA";
                    print OUT "\t$sample_gene_pi{$sample}{$gene}";
                }
                print OUT "\n";
            }
        }
        close(OUT);
    }
    ###
    if ($gff_file and $genetic_code_file) {
        print"$outprefix.intradiv-aminoacid-per-gene.txt\n";
        open (OUT, ">$outprefix.intradiv-aminoacid-per-gene.txt");
        print OUT "Contig\tGene\tLength\tStart\tEnd\tStrand\tNum_loci";
        foreach $sample (@samples_plus) {
            #print OUT "\t$sample Intra_gene_aminoacid_pi\t$sample Intra_gene_aminoacid_pi_p-value";
            print OUT "\t$sample Intra_gene_aminoacid_pi";
        }
        print OUT "\n";
        foreach $contig (@contigs) {
            #print "$contig\n";
            @genes = @{ $contig_genes{$contig} };
            #print "@genes\n"; die;
            foreach $gene (@genes) {
                $num_loci = (keys %{$gene_locus{$gene}});
                print OUT "$contig\t$gene\t$gene_length{$gene}\t$gene_start{$gene}\t$gene_end{$gene}\t$gene_strand{$gene}\t$num_loci";
                foreach $sample (@samples_plus) {
                    #print OUT "\t$sample_gene_aminoacid_pi{$sample}{$gene}\tNA";
                    print OUT "\t$sample_gene_aminoacid_pi{$sample}{$gene}";
                }
                print OUT "\n";
            }
        }
        close(OUT);
        
        print"$outprefix.pNpS-per-gene.txt\n";
        open (OUT, ">$outprefix.pNpS-per-gene.txt");
        print OUT "Contig\tGene\tLength\tStart\tEnd\tStrand\tNum_loci";
        foreach $sample (@samples_plus) {
            #print OUT "\t$sample Intra_gene_aminoacid_pi\t$sample Intra_gene_aminoacid_pi_p-value";
            print OUT "\t$sample pNpS";
        }
        print OUT "\n";
        foreach $contig (@contigs) {
            #print "$contig\n";
            @genes = @{ $contig_genes{$contig} };
            #print "@genes\n"; die;
            foreach $gene (@genes) {
                $num_loci = (keys %{$gene_locus{$gene}});
                print OUT "$contig\t$gene\t$gene_length{$gene}\t$gene_start{$gene}\t$gene_end{$gene}\t$gene_strand{$gene}\t$num_loci";
                foreach $sample (@samples_plus) {
                    if (defined $sample_gene_pNpS{$sample}{$gene}) {
                        print OUT "\t$sample_gene_pNpS{$sample}{$gene}";
                    } else {
                        print OUT "\tNA";
                    }
                }
                print OUT "\n";
            }
        }
        close(OUT);

    }
    ###
    if (@samples > 1) {
        print"$outprefix.fst.txt\n";
        open (OUT, ">$outprefix.fst.txt");
        foreach $sample (@samples) {
            print OUT "\t$sample";
        }
        print OUT "\n";
        for ($ix1 = 0; $ix1 < @samples; $ix1++) {
            print OUT "$samples[$ix1]";
            for ($ix2 = 0; $ix2 < @samples; $ix2++) {
                if ($ix1 == $ix2) {
                    print OUT "\tNA";
                } else {
                    print OUT "\t$sample_sample_fst{$samples[$ix1]}{$samples[$ix2]}";
                }
            }
            print OUT "\n";
        }
        close(OUT);
    }
    ###
    if ((@samples > 1) and $gff_file) {
        print"$outprefix.fst-per-gene.txt\n";
        open (OUT, ">$outprefix.fst-per-gene.txt");
        print OUT "Contig\tGene\tLength\tStart\tEnd\tStrand\tNum_loci";
        for ($ix1 = 0; $ix1 < @samples; $ix1++) {
            for ($ix2 = ($ix1 + 1); $ix2 < @samples; $ix2++) {
                print OUT "\t$samples[$ix1] $samples[$ix2]";
            }
        }
        print OUT "\n";
        foreach $contig (@contigs) {
            @genes = @{ $contig_genes{$contig} };
            foreach $gene (@genes) {
                $num_loci = (keys %{$gene_locus{$gene}});
                print OUT "$contig\t$gene\t$gene_length{$gene}\t$gene_start{$gene}\t$gene_end{$gene}\t$gene_strand{$gene}\t$num_loci";
                for ($ix1 = 0; $ix1 < @samples; $ix1++) {
                    for ($ix2 = ($ix1 + 1); $ix2 < @samples; $ix2++) {
                        #print"$ix1 $ix2\n";
                        if (defined $sample_sample_gene_fst{$samples[$ix1]}{$samples[$ix2]}{$gene}) {
                            print OUT "\t$sample_sample_gene_fst{$samples[$ix1]}{$samples[$ix2]}{$gene}";
                        } else {
                            print OUT "\tNA";
                        }
                    }
                }
                print OUT "\n";
            }
        }
    }
    ###
    if ((@samples > 1) and $gff_file and $genetic_code_file) {
        print"$outprefix.fst-aminoacid-per-gene.txt\n";
        open (OUT, ">$outprefix.fst-aminoacid-per-gene.txt");
        print OUT "Contig\tGene\tLength\tStart\tEnd\tStrand\tNum_loci";
        for ($ix1 = 0; $ix1 < @samples; $ix1++) {
            for ($ix2 = ($ix1 + 1); $ix2 < @samples; $ix2++) {
                print OUT "\t$samples[$ix1] $samples[$ix2]";
            }
        }
        print OUT "\n";
        foreach $contig (@contigs) {
            @genes = @{ $contig_genes{$contig} };
            foreach $gene (@genes) {
                $num_loci = (keys %{$gene_locus{$gene}});
                print OUT "$contig\t$gene\t$gene_length{$gene}\t$gene_start{$gene}\t$gene_end{$gene}\t$gene_strand{$gene}\t$num_loci";
                for ($ix1 = 0; $ix1 < @samples; $ix1++) {
                    for ($ix2 = ($ix1 + 1); $ix2 < @samples; $ix2++) {
                        #print"$ix1 $ix2\n";
                        if (defined $sample_sample_gene_aminoacid_fst{$samples[$ix1]}{$samples[$ix2]}{$gene}) {
                            print OUT "\t$sample_sample_gene_aminoacid_fst{$samples[$ix1]}{$samples[$ix2]}{$gene}";
                        } else {
                            print OUT "\tNA";
                        }
                    }
                }
                print OUT "\n";
            }
        }
    }
}


