#!/usr/bin/env bash

number_of_threads=4
gold_indels="D:\User_Data\Sentieon\Windows_Sentieon\Windows_Sentieon\Secondary-Analysis\resources\Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
reference="D:\User_Data\Sentieon\Windows_Sentieon\Windows_Sentieon\Secondary-Analysis\resources\reference_uncompressed\human_g1k_v37_decoy.fasta"
dbsnp="D:\User_Data\Sentieon\Windows_Sentieon\Windows_Sentieon\Secondary-Analysis\resources\dbsnp_149.vcf.gz"
platform="ILLUMINA"
sentieon_dir="D:\User_Data\Sentieon\Windows_Sentieon\Windows_Sentieon\Secondary-Analysis\tools\sentieon"
vspipeline_dir="Set Me"
emit="variant"

export SENTIEON_LICENSE="D:\User_Data\Sentieon\Windows_Sentieon\Windows_Sentieon\Secondary-Analysis\GoldenHelix-GenesisDx_cluster.lic"

if ! "$sentieon_dir/bin/sentieon" licsrvr --ping ;
then
   echo "# ******************************************"
   echo "# 0. Starting the license server "
   echo "# ******************************************"
    "$sentieon_dir/bin/sentieon" licsrvr --start --log "D:\User_Data\Sentieon\Windows_Sentieon\Windows_Sentieon\Secondary-Analysis\sention_lic.log" "$SENTIEON_LICENSE"
else
   echo "# ******************************************"
   echo "# 0. The license server is already running "
   echo "# ******************************************"
fi

if [[ -z $1 ]]
then
    echo "the path to the input fastq files can not be empty"
    exit 1 
fi

if [[ -z $2 ]]
then
    echo "the name of the output dir can not be empty "
    exit 1 
fi

sample_dir="$1"

#make the output directory

#mkdir "$output_dir"
output_dir="./${2}"

mkdir -p "$output_dir"

CALL_VARIANTS () #sample #fq1 #fq2
{
    sample="$1"
    group="$sample"
    first_sample_file="$2"
    second_sample_file="$3"



    cd "$output_dir"

    mkdir -p "$sample"
    cd "$sample"

    local current_dir=$(pwd)

    echo "$sample"
    echo "$first_sample_file"
    echo "$second_sample_file"

    echo "# ******************************************"
    echo "# 1. Mapping reads with BWA-MEM, sorting"
    echo "# ******************************************"

    "$sentieon_dir/bin/bwa" mem \
    -M \
    -R "@RG\tID:$group\tSM:$sample\tPL:$platform" \
    -t $number_of_threads \
    -K 10000000 \
    "$reference" \
    "$first_sample_file" \
    "$second_sample_file" \
    | \
    "$sentieon_dir/bin/sentieon" util sort \
    -r "$reference" \
    -o "${current_dir}/${sample}_sorted.bam" \
    -t $number_of_threads \
    --sam2bam \
    -i -

    output_name="${current_dir}/${sample}_sorted.bam"
    
    cat << 'EOF' >> "$sample"

echo "# ******************************************"
echo "# 2. Metrics"
echo "# ******************************************"

"$sentieon_dir/bin/sentieon" driver \
    -r "$reference" \
    -t $number_of_threads \
    -i "${sample}_sorted.bam" \
    --algo MeanQualityByCycle mq_metrics.txt \
    --algo QualDistribution qd_metrics.txt \
    --algo GCBias --summary gc_summary.txt gc_metrics.txt \
    --algo AlignmentStat --adapter_seq '' aln_metrics.txt \
    --algo InsertSizeMetricAlgo is_metrics.txt

"$sentieon_dir/bin/sentieon" plot metrics \
    -o metrics-report.pdf \
    gc=gc_metrics.txt \
    qd=qd_metrics.txt \
    mq=mq_metrics.txt \
    isize=is_metrics.txt

EOF



    echo "# ******************************************"
    echo "# 3. Remove Duplicate Reads"
    echo "# ******************************************"

    "$sentieon_dir/bin/sentieon" driver \
    -t $number_of_threads \
    -i "$output_name" \
    --algo LocusCollector --fun score_info score.txt

    "$sentieon_dir/bin/sentieon" driver \
    -t $number_of_threads \
    -i "$output_name" \
    --algo Dedup --rmdup --score_info score.txt \
    --metrics dedup_metrics.txt "${current_dir}/${sample}_deduped.bam"

    output_name="${current_dir}/${sample}_deduped.bam"


    echo "# ******************************************"
echo "# 4. Indel realigner"
echo "# ******************************************"
"$sentieon_dir/bin/sentieon" driver \
    -r "$reference" \
    -t $number_of_threads \
    -i "$output_name" \
    --algo Realigner \
    -k "$gold_indels" \
    "${sample}_realigned.bam"


output_name="${sample}_realigned.bam"

echo "# ******************************************"
echo "# 5. Base recalibration"
echo "# ******************************************"

"$sentieon_dir/bin/sentieon" driver \
    -r "$reference" \
    -t $number_of_threads \
    -i "$output_name" \
    --algo QualCal \
    -k "$dbsnp" \
    -k "$gold_indels" \
     recal_data.table

"$sentieon_dir/bin/sentieon" driver \
    -r "$reference" \
    -t $number_of_threads \
    -i "$output_name" \
    -q recal_data.table \
    --algo QualCal \
    -k "$dbsnp" \
    -k "$gold_indels" \
    recal_data.table.post

"$sentieon_dir/bin/sentieon" driver \
    -t $number_of_threads \
    --algo QualCal \
    --plot  \
    --before recal_data.table \
    --after recal_data.table.post \
    recal.csv

"$sentieon_dir/bin/sentieon" plot QualCal -o recal_plots.pdf recal.csv

"$sentieon_dir/bin/sentieon" driver \
    -r "$reference" \
    -t $number_of_threads \
    -i "$output_name" \
    -q recal_data.table \
    --algo ReadWriter "${sample}_recal.bam"

output_name="${sample}_recal.bam"


echo "# ******************************************"
echo "# 6. HC Variant caller"
echo "# ******************************************"

"$sentieon_dir/bin/sentieon" driver \
    -r "$reference" \
    -t $number_of_threads \
    -i "$output_name" \
    --algo Haplotyper  \
    --emit_conf=10 \
    --emit_mode=$emit \
    -d "$dbsnp" \
    --call_conf=30 \
    "${sample}_output-hc.g.vcf.gz"
	

}


f_sample_file=""
r_sample_file=""

echo "$sample_dir"
script="$(pwd)"
#The fastq files for each sample are next to each other when listed alphabetically 
find "${sample_dir}" -name '*.fastq.gz' -exec readlink -f {} \; | grep -v "Undetermined" | sort | while read line

do
    if [ "$f_sample_file" == '' ]
    then
    f_sample_file="$line"
    continue
    fi
    r_sample_file="$line"

    #Parse the sample name out of the file name 
    sample=$(echo "${line}" | rev | cut -d '/' -f 1 | cut -d '.' -f 3 | cut -d '_' -f 3- | cut -d ' ' -f 1 | rev)

    echo "sampleName: ${sample}"
    echo "$f_sample_file"
    echo "$r_sample_file"
    
   CALL_VARIANTS "$sample" "$f_sample_file" "$r_sample_file"
   cd "$script"
    f_sample_file=""
	
done;


