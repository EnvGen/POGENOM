#!/bin/bash -l

wd=$(pwd)

err_report() {
    echo "Error on line $1 - script Input_POGENOM.sh"
    cd $wd
    if test -f "temporal"; then rm temporal; fi
    exit 1
}

trap 'err_report $LINENO' ERR

#Default options
configFile="here"

#----Argument parse----------------------

for a in "$@"
do
case $a in
  -d=*|--path_to_config_file=*)
  if [ -z "${a#*=}" ];  then
  echo "value to argument -d No supplied"
  exit 0
  else configFile="${a#*=}"
  fi
  shift # past argument
  ;;

  *)
  echo -e '\nUsage: bash Input_POGENOM.sh [options]\n -d=<absolute path to configFile. Default=/<current>/<directory>/<to>/config_files/Input_POGENOM_config.json>\n'
  echo -e 'Description:\nThis program executes a pipeline that generates the required input files for POGENOM.\nThe aim of this pipeline is to increase the reproducibility of the data analysis, and to simplify the use of POGENOM.\nPOGENOM is a computer program that calculates several population genetic parameters for a genome in relation to a set of samples (https://github.com/EnvGen/POGENOM).'
  exit 0
  ;;

esac
done

if [[ "$configFile" == here ]]; then configFile=$wd/config_files/Input_POGENOM_config.json; fi


if [[ "$configFile" != /* ]] || [ -z "$configFile" ]; then
    echo "Please provide an absoltute path to configfile e.g., bash Input_POGENOM.sh '/absolute/path/to/configfile' "
    exit 0
fi

cat $configFile | sed s/"[{|}]"//g | sed s/":"/"="/g | sed s/",$"//g | sed s/" ="/"="/g | sed s/"= "/"="/g | sed s/'"'//g | sed s/" "//g > temporal
. temporal

if [[ "$workdir" != /* ]] || [ -z "$workdir" ]; then
    echo "Please provide an absoltute path to the working directory in configfile e.g., 'workdir': '/absolute/path/to/working_directory/' "
    rm temporal
    exit 0
fi


#----Using prefilt mode - full workflow
mkdir -p $workdir/log_files
if  [[ "$mode" == prefilt ]]; then
           cd $workdir
           echo "INFO: Generating Reads subsets - Fraction used $fraction"
           bash src/create_prefilt_Reads_subdir.sh $fraction $genomes_ext $reads_ext $temp_sub_Reads_dir $dataset
           echo "INFO: Calculating Genome Median coverage - sub-samples - Median coverage threshold $min_coverage"

   snakemake -s snakefiles/step_filter -j $threads --quiet 2>log_files/samples_filter_$dataset.log

           if [[ "$remove_subreads" == yes ]] && test -d "$temp_sub_Reads_dir/Reads"; then
                echo "WARNING: You have chosen to remove $temp_sub_Reads_dir/Reads/"
                rm -rf $temp_sub_Reads_dir/Reads/
           fi

           file_empty=$(grep -v "#" PREFILT/$dataset/Selected_samples_Genomes.tmp | wc -l)
             if [ "$file_empty" -eq 0 ]; then
                echo -e "INFO: With the current parameter setting: Dataset $dataset - Fraction $fraction - Median coverage threshold $min_coverage - Min-base quality $min_bsq_for_cov_median_calculation - Mapping quality $mapqual\n      There is no Genome - sample with Estimated Median Coverage higher than threshold.\n      A vcf file cannot be created\n"
             else
                cat PREFILT/$dataset/Selected_samples_Genomes.tmp | grep -v "#" | while read line
                do
                  echo "INFO: Calculating Genome Median coverage and breadth - Dataset $dataset - Median coverage threshold $min_coverage - Breadth threshold $min_breadth %"
                  mag=$(echo $line | cut -d " " -f1)
                  samples=$(echo $line | cut -d " " -f2)

   snakemake -s snakefiles/step_pogenom_input step1_all --config my_mag="$mag" my_samples="$samples" -j $threads --quiet 2> log_files/$dataset.$mag.coverage_breadth.log

                  echo "INFO: Generating VCF files - Genome $mag"

   snakemake -s snakefiles/step_pogenom_input vcf --config my_mag="$mag" my_samples="$samples" -j $threads --quiet 2> log_files/$dataset.$mag"_vcf_files.log"

                done
             fi

            no_genome=$(grep "#" PREFILT/$dataset/Selected_samples_Genomes.tmp | wc -l)
             if [ "$no_genome" -ne 0 ]; then
                 echo "**********************************************"
                 echo "The following Genome(s) has(have) not been analysed"
                 cat PREFILT/$dataset/Selected_samples_Genomes.tmp | grep "#" | while read line; do echo $line; done
                 echo -e "**********************************************\n"
             fi

           echo "INFO: Input_POGENOM pipeline is done !!!"

rm temporal

exit 0
fi
#---End of prefilt mode

#---Option when analysing a dataset without prefilt

cd $workdir
    echo "INFO: Calculating Genome Median coverage and breadth - Dataset $dataset - Median coverage threshold $min_coverage - Breadth threshold $min_breadth %"

snakemake -s snakefiles/step1_pogenom_input step1_all -j $threads --quiet 2> log_files/$dataset"_Genomes_coverage_breadth.log"

    echo "INFO: Generating VCF files"
snakemake -s snakefiles/step1_pogenom_input vcf -j $threads --quiet 2> log_files/$dataset"_Genomes_vcf_files.log"

    rm temporal

echo 'INFO: Input_POGENOM pipeline is ready !!!'
