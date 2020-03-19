#!/usr/bin/env python3

import os
import argparse
import gzip
import io
import numpy as np

usage = 'parse_subsample_median_cov.py -i -o -f -t'
description = 'This program creates a table reporting the estimated median coverage per sample per Genome in the input dataset. Outputfile Estimated_median_cov_per_sample.tsv'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument ('-i', dest= 'inf', help='input Genome directory', required = True)
parser.add_argument ('-o', dest= 'out', help='output directory', required = True)
parser.add_argument ('-f', dest= 'f', help='fraction used when subsampling, e.g., "0.10"', required = True)
parser.add_argument ('-t', dest= 't', help='Median coverage threshold, e.g., "10"', required = True)

args = parser.parse_args()

MAG_list = [ os.path.basename(d.path) for d in os.scandir(args.inf) if d.is_dir() ]

with open(os.path.join(args.out,"Estimated_median_cov_per_sample.tsv"), "w") as fout1, open(os.path.join(args.out,"Selected_samples_Genomes.tmp"), "w") as fout2:
    print("#Genome", "Sample", "Estimated_median_cov",sep="\t", file=fout1)
    print("Selected Genomes, and samples")
    for m in MAG_list:
        selected_sample_list = []
        path1=os.path.join(args.inf,m)
        file_list = [ f for f in os.listdir(path1) if f.endswith("mpileup")]
        for file in file_list:
            sample = file.split("_in_")[0]
            with open(os.path.join(args.inf,m,file), "r") as fin:
                value = []
                for line in fin:
                    line = line.rstrip()
                    col = int(line.split()[3])
                    if col != 0:
                        value.append(col)
                print(m, sample, float(np.median(value)/float(args.f)),sep="\t", file=fout1)
                if float(np.median(value)/float(args.f)) > float(args.t):
                    selected_sample_list.append(sample)

        if len(selected_sample_list) != 0:
           print("Genome: {} SAMPLE(S): {}".format(m," ".join(map(str,selected_sample_list))))
           print(m,",".join(map(str,selected_sample_list)),sep="\t", file=fout2)
        else:
            print("#{} has not sample with Estimated Median coverage higher then threshold".format(m), file=fout2)
