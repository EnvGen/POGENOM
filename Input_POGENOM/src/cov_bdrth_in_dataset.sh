
#!/bin/bash -l

err_report() {
    echo "Error on line $1 - script cov_bdrth_in_dataset.sh"
    exit 1
}

trap 'err_report $LINENO' ERR

# ---- arguments
mpileupfile=$1
bamfile=$2
outbamfile=$3
mag=$4
mincov=$5
minbreadth=$6
threads=$7
dataset=$8
samplename=$9
pdir="${10}"
subsamp="${11}"
m_ext="${12}"
wkd=$(pwd)

#--- Median coverage
cov=$(cut -f4 $mpileupfile | grep -vw "0" | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')

#---size
direct=$(echo $dataset | sed s/"_prefilt"//)

positions=$(awk 'BEGIN{i=0}; !/^>/ {i=i+length($0)} END {print i}' $wkd/RAW_DATA/Genomes/$direct/$mag$m_ext )


mkdir -p Genome_sizes
echo "genome size:" $positions > Genome_sizes/$mag.size

#---breadth
non_zero=$(cut -f4 $mpileupfile | grep -cvw "0")
breadth=$(echo $non_zero*100/$positions | bc -l )

echo "Genome:" $mag "- Sample:" $samplename "Median_coverage:" $cov " breadth %:" $breadth

mkdir -p 04_mergeable/$dataset/$pdir/$mag

#---selection of BAM files and subsample
if (( $(echo "$breadth >= $minbreadth" | bc -l) )) && (( $(echo "$cov >= $mincov" | bc -l) )); then
  if [[ $subsamp == FALSE ]]; then
     echo "        You have selected not to downsample selected BAM file for - Genome: $mag - Sample: $samplename "
     samtools sort -o $outbamfile --threads $threads $bamfile

  else
    echo "        Downsampling coverage to $mincov - Genome: $mag - Sample: $samplename "
    limite=$(echo "scale=3; $mincov/$cov" | bc )
    samp=$(echo "scale=3; ($limite)+10" | bc)
    samtools view -Sbh --threads $threads -s $samp $bamfile | samtools sort -o $outbamfile --threads $threads
  fi
fi
