#!/bin/bash -l
#usage: bash src/create_prefilt_Reads_MAGs_subdir.sh
err_report() {
  echo "Error on line $1 - script create_prefilt_Reads_MAGs_subdir.sh"
  exit 1
}

trap 'err_report $LINENO' ERR

#---arguments---

wd=$(pwd)
fract=$1
mags_ext=$2
reads_ext=$3
dataset=$4

#--- Main ----

mkdir -p $dataset
mkdir -p $dataset/Reads

cd RAW_DATA/Reads
Lgs=($(ls -d *))
cd $wd
for n in "${Lgs[@]}"
  do
    if [[ "$n" != prefilt ]]; then
     Rds=($(ls $wd/RAW_DATA/Reads/$n/*$reads_ext))
     for r in "${Rds[@]}"
        do
           read_file=$(basename $r)
           if [ ! -f $wd/$dataset/Reads/$read_file ] #If file doesn't exit
            then seqtk sample -s100 $r $fract > $wd/$dataset/Reads/$read_file
                 gzip -q $wd/$dataset/Reads/$read_file
           fi
        done
    fi
  done

cd RAW_DATA/MAGs
Mgs=($(ls -d *))
cd $wd
mkdir -p RAW_DATA/MAGs/prefilt
for f in "${Mgs[@]}"
  do
   if [[ "$f" != prefilt ]]; then
     Mds=($(ls $wd/RAW_DATA/MAGs/$f/*$mags_ext))
     for m in "${Mds[@]}"
        do
           mag_file=$(basename $m)
           if [ ! -f $wd/RAW_DATA/MAGs/prefilt/$mag_file ] #If file doesn't exit
            then ln -s $m $wd/RAW_DATA/MAGs/prefilt/.
           fi
        done
   fi
  done
