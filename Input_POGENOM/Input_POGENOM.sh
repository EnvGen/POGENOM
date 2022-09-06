#!/bin/bash -l
wd=$(pwd)
err_report() {
    echo "Error on line $1 - script Input_POGENOM.sh"
    cd $wd
    if [ -z "$mag" ]; then mag=$(fullname=$(basename $(ls RAW_DATA/Genomes/$dataset/*$genomes_ext | head -1)) ; echo ${fullname%.*}); fi
    if [ -z "$samples" ]; then samples=$(fulln=$(basename $(ls RAW_DATA/Reads/$dataset/*$fwd_index$reads_ext | head -1)) ; echo ${fulln%$fwd_index$reads_ext}); fi
    mess1="Check if the 'ip_env' has been activated - command:\n    'conda activate ip_env'"
    mess2="Use the command:\n    'snakemake -s snakefiles/step1_pogenom_input --unlock'"
    mess3="Use the command:\n    'snakemake -s snakefiles/step_pogenom_input --config my_mag='"$mag"' my_samples='"$samples"' --unlock'"
    mess="or\n    'snakemake -s snakefiles/step2 --config mag_name='"$mag"' --unlock' or\n    'snakemake -s snakefiles/step2_B --config mag_name='"$mag"' --unlock'"
    if [ "$1" == 79 ]; then echo -e "if you are using conda, check if the 'ip_env' has been activated - command:\n    'conda activate ip_env' "; fi
    if [ "$1" == 81 ]; then echo -e "TIP 1 - Look at $wd/log_files/samples_filter_$dataset.log\nTIP 2 - $mess1\nTIP 3 - Use the command:\n    'snakemake -s snakefiles/step_filter --unlock'\n     and run the pipeline again"; fi
    if [ "$1" == 96 ]; then echo -e "TIP 1 - Look at $wd/log_files/$dataset.$mag.coverage_breadth.log\nTIP 2 - $mess1\nTIP 3 - $mess3\n     and run the pipeline again"; fi
    if [ "$1" == 98 ]; then echo -e "TIP 1 - Look at $wd/log_files/$dataset.$mag"_vcf_files.log"\nTIP 2 - $mess1\nTIP 3 - $mess3 $mess\n     and run the pipeline again"; fi
    if [ "$1" == 110 ]; then echo -e "TIP 1 - Look at\n    $wd/log_files/$dataset.$mag.coverage_breadth.log or\n    $wd/log_files/$dataset.$mag"_vcf_files.log"\nTIP 2 - $mess1\nTIP 3 - $mess3\n     and run the pipeline again"; fi
    if [ "$1" == 114 ]; then echo -e "TIP 1 - Look at $wd/log_files/$dataset"_Genomes_coverage_breadth.log"\nTIP 2 - $mess1\nTIP 3 - $mess2\n     and run the pipeline again"; fi
    if [ "$1" == 116 ]; then echo -e "TIP 1 - Look at $wd/log_files/$dataset"_Genomes_vcf_files.log"\nTIP 2 - $mess1\nTIP 3 - $mess2 $mess\n     and run the pipeline again"; fi
    if test -f "temporal"; then rm temporal; fi
    exit 1
}
trap 'err_report $LINENO' ERR
start=`date +%s`
#Default options
configFile=$wd/config_files/Input_POGENOM_config.json
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
  echo -e "\nUsage: bash Input_POGENOM.sh [options]\n -d=<absolute path to configFile. Default=$configFile>\n"
  echo -e 'Description:\nThis program executes a pipeline that generates the required input files for POGENOM.\nThe aim of this pipeline is to increase the reproducibility of the data analysis, and to simplify the use of POGENOM.\nPOGENOM is a computer program that calculates several population genetic parameters for a genome in relation to a set of samples (https://github.com/EnvGen/POGENOM).'
  exit 0
  ;;
esac
done
if [[ "$configFile" != /* ]] || [ -z "$configFile" ]; then
    echo "Please provide an absoltute path to configfile e.g., bash Input_POGENOM.sh '/absolute/path/to/configfile' "
    exit 0
fi
cat $configFile | sed s/"[{|}]"//g | sed s/":"/"="/g | sed s/",$"//g | sed s/" ="/"="/g | sed s/"= "/"="/g | sed s/'"'//g | sed s/" "//g > temporal
. temporal

if [[ "$workdir" != /* ]] || [ -z "$workdir" ]; then
    echo "Please provide an absoltute path to the working directory in configfile e.g., 'workdir': '/absolute/path/to/working_directory/' "
    echo "Please double-check the absolute path to working directory $workdir and to configfile $configFile"
    rm temporal
    exit 0
fi
#Checking key parameters setting
options=("$dataset" "$min_coverage" "$min_breadth" "$min_bsq_for_cov_median_calculation" "$threads" "$genomes_ext" "$reads_ext" "$fwd_index" "$rev_index" "$mapper" "$mapqual" "$freebayes_parameters" "$vcffilter_qual" "$subsample")
for o in "${options[@]}"; do if [ -z "$o" ]; then echo "A key parameter is undefined, please check in the config_files/Input_POGENOM_config.json file the parameters used"; exit 1; fi; done
if [[ "$mapper" == "bwa" ]] && [ -z "$coverm_filter_params" ]; then echo "CoverM filter parameters need to be defined, please check in the config_files/Input_POGENOM_config.json file the parameters used"; exit 1; fi
if [[ "$mapper" != "bwa" ]] && [[ "$mapper" != "bowtie2" ]]; then echo "Aligner/mapper needs to be defined, options: bwa or bowtie2. Please check in the config_files/Input_POGENOM_config.json file the parameters used"; exit 1; fi
if [[ $snakemake_extra_params == *","* ]]; then extra_params=$( echo $snakemake_extra_params | sed s/","/" "/g); else extra_params=$snakemake_extra_params; fi
echo "INFO: Starting Input_POGENOM pipeline - Working directory: $workdir"
if [[ "$mapper" == "bowtie2" ]] && [[ "$coverm" == FALSE ]]; then echo -e "WARNING: bowtie2 params to be used: $bowtie2_params \n         CoverM Filtre won't be used\n         Minimum bowtie2 parameters suggested: --ignore-quals --mp 1,1 --np 1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.05"; fi
#----Using prefilt mode - full workflow
mkdir -p $workdir/log_files
if  [[ "$mode_prefilt" == TRUE ]]; then
#Checking key parameters setting
  options2=("$fraction" "$temp_sub_Reads_dir")
  for p in "${options2[@]}"; do if [ -z "$p" ]; then echo 'A key parameter in "mode_prefilt" is undefined, please check in the config_files/Input_POGENOM_config.json file the parameters used'; exit 1; fi; done
  # main - mode prefilt
      cd $workdir
      result_dir="PREFILT/"$dataset"/params_cov_"$min_coverage"_mpq_"$mapqual"_bq_"$min_bsq_for_cov_median_calculation"_fr_"$fraction
      file_ready=$result_dir/Selected_samples_Genomes.txt
      if [ ! -s "$file_ready" ]; then
           echo "INFO: Generating Reads subsets - Fraction used $fraction"
           bash src/create_prefilt_Reads_subdir.sh $fraction $reads_ext $temp_sub_Reads_dir $dataset $threads
           echo "INFO: Calculating Genome Median coverage - sub-samples - Median coverage threshold $min_coverage"
           snakemake -s snakefiles/step_filter -j $threads --use-conda $extra_params 2>log_files/samples_filter_$dataset.log
           if [[ "$remove_subreads" == TRUE ]] && test -d "$temp_sub_Reads_dir/Reads"; then
                echo "WARNING: You have chosen to remove $temp_sub_Reads_dir/Reads/"
                rm -rf $temp_sub_Reads_dir/Reads/
           fi
       else echo -e "INFO: The file $result_dir/Selected_samples_Genomes.txt already exists and it is not empty"; fi
               file_empty=$(grep -v "#" $result_dir/Selected_samples_Genomes.txt | wc -l)
               if [ "$file_empty" -eq 0 ]; then
                  echo -e "INFO: With the current parameter setting: Dataset $dataset - Fraction $fraction - Median coverage threshold $min_coverage - Min-base quality $min_bsq_for_cov_median_calculation - Mapping quality $mapqual\n      There is no Genome - sample with Estimated Median Coverage higher than threshold.\n      A vcf file cannot be created\n"
               else
                  echo "INFO: Calculating Genome Median coverage and breadth - Dataset: $dataset - Median coverage threshold: $min_coverage - Breadth threshold: $min_breadth %"
                  grep -v "#" $result_dir/Selected_samples_Genomes.txt | while read line
                  do
                    mag=$(echo $line | cut -d " " -f1)
                    samples=$(echo $line | cut -d " " -f2)
                    snakemake -s snakefiles/step_pogenom_input step1_all --config my_mag="$mag" my_samples="$samples" -j $threads --use-conda $extra_params 2> log_files/$dataset.$mag.coverage_breadth.log
                    echo "INFO: Generating VCF files - Genome $mag"
                    snakemake -s snakefiles/step_pogenom_input vcf --config my_mag="$mag" my_samples="$samples" -j $threads $extra_params 2> log_files/$dataset.$mag"_vcf_files.log"
                  done
               fi
               no_genome=$(grep "#" $result_dir/Selected_samples_Genomes.txt | wc -l )
               if [ "$no_genome" -ne 0 ]; then
                   echo "**********************************************"
                   echo "The following Genome(s) has(have) not been analysed"
                   grep "#" $result_dir/Selected_samples_Genomes.txt
                   echo -e "**********************************************\n"
               fi
 rm temporal; echo "INFO: Input_POGENOM pipeline is done !!!"; end=`date +%s`; runtimes=$( echo "$end - $start" | bc -l ); runtimem=$( printf "%.2f \n" $(echo $runtimes/60 | bc -l ) ); echo "***** Total time (s) $runtimes | (min) $runtimem"
  #---End of prefilt mode
else
#---Option when analysing a dataset without prefilt
  cd $workdir
      echo "INFO: Calculating Genome Median coverage and breadth - Dataset: $dataset - Median coverage threshold: $min_coverage - Breadth threshold: $min_breadth %"
      snakemake -s snakefiles/step1_pogenom_input step1_all -j $threads --use-conda $extra_params 2> log_files/$dataset"_Genomes_coverage_breadth.log"
      echo "INFO: Generating VCF files"
      snakemake -s snakefiles/step1_pogenom_input vcf -j $threads $extra_params 2> log_files/$dataset"_Genomes_vcf_files.log"
      rm temporal
echo 'INFO: Input_POGENOM pipeline is done !!!'; end=`date +%s`; runtimes=$( echo "$end - $start" | bc -l ); runtimem=$( printf "%.2f \n" $(echo $runtimes/60 | bc -l ) ); echo "***** Total time (s) $runtimes | (min) $runtimem"
fi
