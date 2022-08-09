#!/bin/bash -l

err_report() {
  echo "Error on line $1 - script create_prefilt_Reads_subdir.sh"
  exit 2
}

trap 'err_report $LINENO' ERR

start=`date +%s`
#---arguments---
wd=$(pwd)
fract=$1
reads_ext=$2
dataset=$3
Dts=$4
maxjobs=$5
#--- Main ----

mkdir -p $dataset/Reads/fraction_$fract

subsample_reads() {
  r=$1
  wd=$2
  dataset=$3
  fract=$4
  read_file=$(basename $r)
  echo "INFO: Working on $read_file"

   if ! test -s $wd/$dataset/Reads/fraction_$fract/$read_file  #If file doesn't exit or if it exist but it is empty
      then
           all_lines_in=$(gzip -cd $r | wc -l)
           all_number_reads=$( echo $all_lines_in/4 | bc )
           subsample=$( echo $all_number_reads*$fract | bc )
           seqtk sample -s100 $r $subsample > $wd/$dataset/Reads/fraction_$fract/temponame_$read_file
           lines_in=$(cat $wd/$dataset/Reads/fraction_$fract/temponame_$read_file | wc -l)
           number_reads=$( echo $lines_in/4 | bc )
           gzip -c $wd/$dataset/Reads/fraction_$fract/temponame_$read_file > $wd/$dataset/Reads/fraction_$fract/$read_file
           rm $wd/$dataset/Reads/fraction_$fract/temponame_$read_file
           echo "      Subset $read_file created - Number of reads in subset: $number_reads"
   else
        echo "      Subset $read_file already created"
   fi

}

Rds=($(ls $wd/RAW_DATA/Reads/$Dts/*$reads_ext))
system_="linux"
cat /proc/meminfo >/dev/null 2>/dev/null || system_="Mac"

for r in "${Rds[@]}"
   do
      if [ "$system_" == "Mac" ]; then
        memavble=$(memory_pressure | tail -1 | cut -d ":" -f2 | sed s/%// | tr -d " ")
      else
        memavble=$(mf=$(cat /proc/meminfo | grep "MemFree" | cut -d":" -f2 | sed s/"kB"// | tr -d " ");\
                  mt=$(cat /proc/meminfo | grep "MemTotal" | cut -d":" -f2 | sed s/"kB"// | tr -d " ");\
                  echo $(echo $mf*100/$mt | bc))
      fi
    subsample_reads "$r" "$wd" "$dataset" "$fract" &
    if [ "$memavble" -le 2 ] || [ "$(jobs -r -p | wc -l)" -ge "$maxjobs" ]; then wait; fi
   done
wait

end=`date +%s`
runtimes=$( echo "$end - $start" | bc -l )
runtimem=$( printf "%.2f \n" $(echo $runtimes/60 | bc -l ) )
echo "INFO: Subsampling reads done in (s) $runtimes | (min) $runtimem"
