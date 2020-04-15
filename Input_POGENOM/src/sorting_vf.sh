#!/bin/bash -l

err_report() {
    echo "Error on line $1 - script sorting_vf.sh"
    exit 1
}

trap 'err_report $LINENO' ERR

#---arguments

wdir=$1
MAGs_path=$2
threads=$3
dts=$4
pdir=$5

#---- Main
declare -a mags
mags=()
cd $wdir/$MAGs_path/
mags=($(ls -d *))
cd $wdir

for i in "${mags[@]}"
do
if [ "$(ls -A $wdir/$MAGs_path/$i)" ]
  then
    files=($(ls $wdir/$MAGs_path/$i/*.bam))
    if (( $(echo "${#files[@]} > 1" | bc -l) ))
       then
            snakemake -s snakefiles/step2 --config mag_name="$i" -j $threads -F

    elif (( $(echo "${#files[@]} == 1" | bc -l) ))
        then
          snakemake -s snakefiles/step2_B --config mag_name="$i" -j $threads -F
    fi
else
  mkdir -p $wdir/06_VCF/$dts/$pdir
  echo "The genome $i has not BAM file that passes the filter breadth and coverage. A vcf file cannot be created" > $wdir/06_VCF/$dts/$pdir/$i"_samples.txt"
  echo "The genome $i has not BAM file that passes the filter breadth and coverage. Directory $MAGs_path/$i is empty, and a vcf file cannot be created"
fi
done
