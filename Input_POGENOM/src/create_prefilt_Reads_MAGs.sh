#!/bin/bash -l
#usage: bash src/create_prefilt_Reads_MAGs_subdir.sh
err_report() {
  echo "Error on line $1 - script create_prefilt_Reads_MAGs.sh"
  exit 1
}

trap 'err_report $LINENO' ERR

#-----arguments
wd=$(pwd)
dataset="prefilt"
reads_ext=$1

#---- Main
cd RAW_DATA/Reads
Lgs=($(ls -d *))
cd $wd
mkdir -p RAW_DATA/Reads/$dataset
for n in "${Lgs[@]}"
  do
   if [[ "$n" != prefilt ]]; then
     Rds=($(ls $wd/RAW_DATA/Reads/$n/*$reads_ext))
     for r in "${Rds[@]}"
        do
           read_file=$(basename $r)
           if [ ! -f $wd/RAW_DATA/Reads/$dataset/$read_file ] #If file doesn't exit
            then ln -s $r $wd/RAW_DATA/Reads/$dataset/.
           fi
        done
   fi
  done
