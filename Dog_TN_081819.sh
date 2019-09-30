set -x

#!/bin/sh
# *******************************************
# Script to perform TN seq variant calling
# using a matched paired Tumor+normal sample with fastq
# files named normal_1.fastq.gz, normal_2.fastq.gz
# tumor_1.fastq.gz, tumor_2.fastq.gz
# *******************************************

# Update with the fullpath location of your sample fastq
fastq_folder="/home/goldenadmin/SECONDARY-ANALYSIS/NCSU_FASTQ/" #Change to local folder
tumor_fastq_1="$3"
tumor_fastq_2="$4"
tumor_sample="${5}_TS"
tumor_group="${5}_TG"
normal_fastq_1="$1"
normal_fastq_2="$2"
normal_sample="${5}_NS"
normal_group="${5}_NG"
platform="ILLUMINA"

# Update with the location of the reference data files
fasta="/home/goldenadmin/SECONDARY-ANALYSIS/resources/Reference_Assemblies/canFam3.fa"
#dbsnp="/home/pipeline/ref_hg19/dbsnp_138.hg19.vcf.gz"


# Update with the location of the Sentieon software package and license file
export SENTIEON_INSTALL_DIR="/home/goldenadmin/SECONDARY-ANALYSIS/tools/sentieon/" #Change to local folder
#export SENTIEON_LICENSE=/home/Licenses/Sentieon.lic

# Other settings
nt=32 #number of threads to use in computation
workdir="/home/goldenadmin/SECONDARY-ANALYSIS/NCSU_FASTQ/$5" #Determine where the output files will be stored


export SENTIEON_LICENSE="/home/goldenadmin/SECONDARY-ANALYSIS/GoldenHelix_cluster.lic" #Change to local folder



# ******************************************
# 0. Setup
# ******************************************
mkdir -p $workdir
logfile=$workdir/run.log
exec >$logfile 2>&1
cd $workdir

if ! "$SENTIEON_INSTALL_DIR/bin/sentieon" licsrvr --ping ;
then
   echo "# ******************************************"
   echo "# 0. Starting the license server "
   echo "# ******************************************"
    "$SENTIEON_INSTALL_DIR/bin/sentieon" licsrvr --start --log "/home/goldenadmin/SECONDARY-ANALYSIS/sentieon_lic.log" "$SENTIEON_LICENSE" #Change to local folder
else
   echo "# ******************************************"
   echo "# 0. The license server is already running "
   echo "# ******************************************"
fi



# ******************************************
# 1a. Mapping reads with BWA-MEM, sorting for tumor sample
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 
( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -M -R "@RG\tID:$tumor_group\tSM:$tumor_sample\tPL:$platform" -t $nt -K 10000000 "$fasta" "$fastq_folder/$tumor_fastq_1" "$fastq_folder/$tumor_fastq_2" || echo -n "error" ) | $SENTIEON_INSTALL_DIR/bin/sentieon util sort -o ${tumor_sample}_sorted.bam -t $nt --sam2bam -i -
# ******************************************



# 1b. Mapping reads with BWA-MEM, sorting for normal sample
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 
( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -M -R "@RG\tID:$normal_group\tSM:$normal_sample\tPL:$platform" -t $nt -K 10000000 "$fasta" "$fastq_folder/$normal_fastq_1" "$fastq_folder/$normal_fastq_2" || echo -n "error" ) | $SENTIEON_INSTALL_DIR/bin/sentieon util sort -o ${normal_sample}_sorted.bam -t $nt --sam2bam -i -



# ******************************************
# 2a. Metrics for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i ${tumor_sample}_sorted.bam --algo MeanQualityByCycle ${tumor_sample}_mq_metrics.txt --algo QualDistribution ${tumor_sample}_qd_metrics.txt --algo GCBias --summary ${tumor_sample}_gc_summary.txt ${tumor_sample}_gc_metrics.txt --algo AlignmentStat --adapter_seq '' ${tumor_sample}_aln_metrics.txt --algo InsertSizeMetricAlgo ${tumor_sample}_is_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o ${tumor_sample}_gc-report.pdf ${tumor_sample}_gc_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o ${tumor_sample}_qd-report.pdf ${tumor_sample}_qd_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o ${tumor_sample}_mq-report.pdf ${tumor_sample}_mq_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o ${tumor_sample}_is-report.pdf ${tumor_sample}_is_metrics.txt

# ******************************************
# 2b. Metrics for normal sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i ${normal_sample}_sorted.bam --algo MeanQualityByCycle ${normal_sample}_mq_metrics.txt --algo QualDistribution ${normal_sample}_qd_metrics.txt --algo GCBias --summary ${normal_sample}_gc_summary.txt ${normal_sample}_gc_metrics.txt --algo AlignmentStat --adapter_seq '' ${normal_sample}_aln_metrics.txt --algo InsertSizeMetricAlgo ${normal_sample}_is_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o ${normal_sample}_gc-report.pdf ${normal_sample}_gc_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o ${normal_sample}_qd-report.pdf ${normal_sample}_qd_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o ${normal_sample}_mq-report.pdf ${normal_sample}_mq_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o ${normal_sample}_is-report.pdf ${normal_sample}_is_metrics.txt




# ******************************************
# 3a. Remove Duplicate Reads for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i ${tumor_sample}_sorted.bam --algo LocusCollector --fun score_info ${tumor_sample}_score.txt
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i ${tumor_sample}_sorted.bam --algo Dedup --rmdup --score_info ${tumor_sample}_score.txt --metrics ${tumor_sample}_dedup_metrics.txt ${tumor_sample}_deduped.bam 
# ******************************************
# 3b. Remove Duplicate Reads for normal sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i ${normal_sample}_sorted.bam --algo LocusCollector --fun score_info ${normal_sample}_score.txt
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i ${normal_sample}_sorted.bam --algo Dedup --rmdup --score_info ${normal_sample}_score.txt --metrics ${normal_sample}_dedup_metrics.txt ${normal_sample}_deduped.bam 



# ******************************************
# 3. Corealignment of tumor and normal
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i ${tumor_sample}_sorted.bam -i ${normal_sample}_sorted.bam --algo Realigner tn_corealigned.bam

# ******************************************
# 4. Somatic Variant Calling
# ******************************************

$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i tn_corealigned.bam --algo TNhaplotyper --tumor_sample ${tumor_sample} --normal_sample $normal_sample ${5}-tnhaplotyper.vcf.gz