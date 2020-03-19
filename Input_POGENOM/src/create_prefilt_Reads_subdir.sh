#!/bin/bash -l

err_report() {
  echo "Error on line $1 - script create_prefilt_Reads_subdir.sh"
  exit 1
}

trap 'err_report $LINENO' ERR

#---arguments---

wd=$(pwd)
fract=$1
mags_ext=$2
reads_ext=$3
dataset=$4
Dts=$5

#--- Main ----
mkdir -p $dataset
mkdir -p $dataset/Reads

Rds=($(ls $wd/RAW_DATA/Reads/$Dts/*$reads_ext))

for r in "${Rds[@]}"
   do
      read_file=$(basename $r)
       if [ ! -f $wd/$dataset/Reads/$read_file ] #If file doesn't exit
          then seqtk sample -s100 $r $fract > $wd/$dataset/Reads/$read_file
               gzip -q $wd/$dataset/Reads/$read_file
       fi
   done
