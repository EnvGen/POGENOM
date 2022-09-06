#!/bin/bash -l

err_report() {
  echo "Error on line $1 - script freebayes_parallel.sh"
  exit 2
}

trap 'err_report $LINENO' ERR
samtools faidx $1 -o $1.fai
freebayes-parallel <(fasta_generate_regions.py $1.fai 100000) $2 -f $1 $4 $3 > $5
