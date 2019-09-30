#!/bin/sh
# *******************************************
# Script to perform T_only seq variant calling
# using a Tumor sample with fastq
# files named # tumor_1.fastq.gz, tumor_2.fastq.gz
# and data from a Panel of Normals
# *******************************************

# Update with the fullpath location of your sample fastq
#fastq_folder="/home/ec2-user/input"
tumor_fastq_1="$1"
tumor_fastq_2="$2"
tumor_sample="3"
tumor_group="4"
platform="ILLUMINA"

# Update with the location of the reference data files
fasta="/home/ec2-user/Secondary-Analysis/resources/Reference_Assemblies/canFam3.fa"

# Update with the location of the panel of normal vcf file
panel_of_normal_TNsnv="/home/ec2-user/input/panel_of_normals_snv.vcf"
# We recommend that you create the panel of normal file with the corresponding algorithm that you plan to use for the somatic mutation calling. 
panel_of_normal_TNhaplotyper="/home/ec2-user/input/panel_of_normals_haplotyper.vcf"

# Update with the location of the Sentieon software package and license file
export SENTIEON_INSTALL_DIR=/home/ec2-user/Secondary-Analysis/tools/sentieon/
#export SENTIEON_LICENSE=/home/Licenses/Sentieon.lic

# Other settings
nt=32 #number of threads to use in computation
workdir="/home/ec2-user/output/$5" #Determine where the output files will be stored

# ******************************************
# 0. Setup
# ******************************************
mkdir -p $workdir
logfile=$workdir/run.log
exec >$logfile 2>&1
cd $workdir

# ******************************************
# 1a. Mapping reads with BWA-MEM, sorting for tumor sample
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 
( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -M -R "@RG\tID:$tumor_group\tSM:$tumor_sample\tPL:$platform" -t $nt -K 10000000 $fasta $fastq_folder/$tumor_fastq_1 $fastq_folder/$tumor_fastq_2 || echo -n 'error' ) | $SENTIEON_INSTALL_DIR/bin/sentieon util sort -o tumor_sorted.bam -t $nt --sam2bam -i -

# ******************************************
# 2a. Metrics for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i tumor_sorted.bam --algo MeanQualityByCycle tumor_mq_metrics.txt --algo QualDistribution tumor_qd_metrics.txt --algo GCBias --summary tumor_gc_summary.txt tumor_gc_metrics.txt --algo AlignmentStat --adapter_seq '' tumor_aln_metrics.txt --algo InsertSizeMetricAlgo tumor_is_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o tumor_gc-report.pdf tumor_gc_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o tumor_qd-report.pdf tumor_qd_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o tumor_mq-report.pdf tumor_mq_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o tumor_is-report.pdf tumor_is_metrics.txt

# ******************************************
# 3a. Remove Duplicate Reads for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i tumor_sorted.bam --algo LocusCollector --fun score_info tumor_score.txt
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i tumor_sorted.bam --algo Dedup --rmdup --score_info tumor_score.txt --metrics tumor_dedup_metrics.txt tumor_deduped.bam 

# ******************************************
# 4a. Indel realigner for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i tumor_deduped.bam --algo Realigner tumor_realigned.bam

# ******************************************
# 5a. Base recalibration for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i tumor_realigned.bam --algo QualCal tumor_recal_data.table
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i tumor_realigned.bam -q tumor_recal_data.table --algo QualCal  tumor_recal_data.table.post
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt --algo QualCal --plot --before tumor_recal_data.table --after tumor_recal_data.table.post tumor_recal.csv
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualCal -o tumor_recal_plots.pdf tumor_recal.csv

# ******************************************
# 6. Somatic Variant Calling
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i tumor_realigned.bam -q tumor_recal_data.table --algo TNsnv --tumor_sample $tumor_sample --pon $panel_of_normal__TNsnv --call_stats_out output-call.stats output-tnsnv.vcf.gz
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i tumor_realigned.bam -q tumor_recal_data.table --algo TNhaplotyper --tumor_sample $tumor_sample --pon $panel_of_normal__TNhaplotyper output-call.stats output-tnhaplotyper.vcf.gz
