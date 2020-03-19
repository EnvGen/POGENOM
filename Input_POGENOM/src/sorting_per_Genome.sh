
#!/bin/bash -l

err_report() {
    echo "Error on line $1 - script sorting_per_Genome.sh"
    exit 1
}

trap 'err_report $LINENO' ERR

#---arguments

i=$1
wdir=$2
MAGs_path=$3
threads=$4
dts=$5

if [ "$(ls -A $wdir/$MAGs_path/$i)" ]
  then
    files=($(ls $wdir/$MAGs_path/$i/*.bam))
    if (( $(echo "${#files[@]} > 1" | bc -l) ))
       then
          snakemake -s snakefiles/step2_prefilt_pogenom_input --config mag_name="$i" -j $threads -F --quiet

    elif (( $(echo "${#files[@]} == 1" | bc -l) ))
        then
          snakemake -s snakefiles/step2_B_prefilt_pogenom_input --config mag_name="$i" -j $threads -F --quiet
    fi
else
  mkdir -p $wdir/06_VCF
  mkdir -p $wdir/06_VCF/$dts
  echo "The genome $i has not BAM file that passes the filter breadth and coverage. A vcf file cannot be created" > $wdir/06_VCF/$dts/$i"_samples.txt"
  echo "The genome $i has not BAM file that passes the filter breadth and coverage. Directory $MAGs_path/$i is empty, and a vcf file cannot be created"
fi
