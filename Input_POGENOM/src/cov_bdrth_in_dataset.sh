
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

wkd=$(pwd)

#--- Median coverage
cov=$(cut -f4 $mpileupfile | grep -v "0" | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')

#---size
check=$(echo "$dataset" | cut -d "_" -f2)

if [[ "$check" == "prefilt" ]]; then
   direct=$(echo $dataset | sed s/"_prefilt"//)
else
    direct="$dataset"
fi

positions=$(bowtie2-inspect -s $wkd/01_INDEXING/$direct/$mag/$mag | awk '{ FS = "\t" } ; BEGIN{L=0}; {L=L+$3}; END{print L}')

mkdir -p Genome_sizes
echo "genome size:" $positions > Genome_sizes/$mag.size

#---breadth
non_zero=$(cut -f4 $mpileupfile | grep -cv "0")
breadth=$(echo $non_zero*100/$positions | bc -l )

echo "Genome:" $mag "- Sample:" $samplename "Median_coverage:" $cov " breadth %:" $breadth

mkdir -p 04_mergeable/$dataset/$pdir/$mag

#---selection of BAM files and subsample
if (( $(echo "$breadth > $minbreadth" | bc -l) )) && (( $(echo "$cov > $mincov" | bc -l) )); then
  limite=$(echo "scale=3; $mincov/$cov" | bc )
  samp=$(echo "scale=3; ($limite)+10" | bc)
  samtools view -Sbh --threads $threads -s $samp $bamfile | samtools sort -o $outbamfile --threads $threads
fi
