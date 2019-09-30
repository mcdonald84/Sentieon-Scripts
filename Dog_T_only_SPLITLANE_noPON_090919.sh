set -x
#!/bin/sh
# *******************************************
# Script to perform T_only seq variant calling
# using a Tumor sample with fastq
# files named # tumor_1.fastq.gz, tumor_2.fastq.gz, tumor_3.fastq.gz, tumor_4.fastq.gz, 
# and data from a DUMMY Panel of Normals
# *******************************************

# Update with the fullpath location of your sample fastq
#fastq_folder="/home/ec2-user/input"
tumor_fastq_1="Sample_5Q_S15_L001_R1_001.fastq.gz"
tumor_fastq_2="Sample_5Q_S15_L001_R2_001.fastq.gz"
tumor_fastq_3="Sample_5Q_S15_L002_R1_001.fastq.gz"
tumor_fastq_4="Sample_5Q_S15_L002_R2_001.fastq.gz"
tumor_sample="tumor_sample_name"
tumor_group="tumor_read_group_name"
platform="ILLUMINA"

# Update with the location of the reference data files
fasta="/home/ec2-user/Secondary-Analysis/resources/Reference_Assemblies/canFam3.fa"

# We recommend that you create the panel of normal file with the corresponding algorithm that you plan to use for the somatic mutation calling. 
panel_of_normal_TNhaplotyper="/home/ec2-user/Secondary-Analysis/resources/noPON_haplotyper.vcf.gz"

# Update with the location of the Sentieon software package and license file
export SENTIEON_INSTALL_DIR=/home/ec2-user/Secondary-Analysis/tools/sentieon/
#export SENTIEON_LICENSE=/home/Licenses/Sentieon.lic

# Other settings
nt=32 #number of threads to use in computation
workdir="/home/ec2-user/output/5Q" #Determine where the output files will be stored

# ******************************************
# 0. Setup
# ******************************************
mkdir -p $workdir
logfile=$workdir/run.log
exec >$logfile 2>&1
cd $workdir

# ******************************************
# 1. Mapping reads with BWA-MEM, sorting for tumor sample
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 
( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -M -R "@RG\tID:$tumor_group\tSM:$tumor_sample\tPL:$platform" -t $nt -K 10000000 "$fasta" "$fastq_folder/$tumor_fastq_1" "$fastq_folder/$tumor_fastq_2" || echo -n 'error' ) | $SENTIEON_INSTALL_DIR/bin/sentieon util sort -o ${tumor_sample}_sorted.bam -t $nt --sam2bam -i -
# ******************************************


# ******************************************
# 2. Metrics for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i ${tumor_sample}_sorted.bam --algo MeanQualityByCycle ${tumor_sample}_mq_metrics.txt --algo QualDistribution ${tumor_sample}_qd_metrics.txt --algo GCBias --summary ${tumor_sample}_gc_summary.txt ${tumor_sample}_gc_metrics.txt --algo AlignmentStat --adapter_seq '' ${tumor_sample}_aln_metrics.txt --algo InsertSizeMetricAlgo ${tumor_sample}_is_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o ${tumor_sample}_gc-report.pdf ${tumor_sample}_gc_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o ${tumor_sample}_qd-report.pdf ${tumor_sample}_qd_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o ${tumor_sample}_mq-report.pdf ${tumor_sample}_mq_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o ${tumor_sample}_is-report.pdf ${tumor_sample}_is_metrics.txt

# ******************************************
# 3. Remove Duplicate Reads for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i ${tumor_sample}_sorted.bam --algo LocusCollector --fun score_info ${tumor_sample}_score.txt
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i ${tumor_sample}_sorted.bam --algo Dedup --rmdup --score_info ${tumor_sample}_score.txt --metrics ${tumor_sample}_dedup_metrics.txt ${tumor_sample}_deduped.bam 
# ******************************************


# ******************************************
# 4. Somatic Variant Calling
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i ${tumor_sample}_deduped.bam --algo TNhaplotyper --tumor_sample ${tumor_sample} --pon $panel_of_normal_TNhaplotyper ${tumor_sample}_tnhaplotyper.vcf.gz