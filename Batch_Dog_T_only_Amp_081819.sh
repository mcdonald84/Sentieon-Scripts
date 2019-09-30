#!/usr/bin/env bash


if [[ -z $1 ]]
then
    echo "the path to the input fastq files can not be empty"
    exit 1 
fi

#if [[ -z $2 ]]
#then
#    echo "the name of the output dir can not be empty "
#    exit 1 
#fi

sample_dir="$1"

#make the output directory

#mkdir "$output_dir"
#output_dir="/home/ec2-user/output/${2}"

#mkdir -p "$output_dir"

first_tumor_sample=""
second_tumor_sample=""

echo -e "\e[2mINPUT Directory:\e[0m \e[95m$sample_dir\e[0m"

#The fastq files for each sample are next to each other when listed alphabetically 
find "${sample_dir}" -name '*.fastq.gz' | grep -v "Undetermined" | sort | while read line
do
    
    if [ "$first_tumor_sample" == '' ]
    then
    first_tumor_sample="$line"
    continue
    fi
    
    if [ "$second_tumor_sample" == '' ]
    then
    second_tumor_sample="$line"
    fi

    #Parse the sample name out of the file name 
    sample=$(echo "${line}" | rev | cut -d '/' -f 1 | cut -d '.' -f 3 | cut -d '_' -f 4- | rev | sed -e 's/ /_/g' )

    echo -e "\e[2mSampleName:\e[0m \e[91m${sample}\e[0m"
    echo -e "\e[2mFirst  Tumor  Sample:\e[0m \e[95m$first_tumor_sample\e[0m"
    echo -e "\e[2mSecond Tumor  Sample:\e[0m \e[95m$second_tumor_sample\e[0m\n"

    
    
#    ./pipeline-example-tumor_normal.sh "$first_tumor_sample" "$second_tumor_sample" "$first_normal_sample" "$second_normal_sample" "$sample"
    ./Dog_T_only_Amp_081819.sh "$first_tumor_sample" "$second_tumor_sample" "$sample"
    
    
    
first_tumor_sample=""
second_tumor_sample=""
    
done;


