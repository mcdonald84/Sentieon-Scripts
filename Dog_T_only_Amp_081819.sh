set -x
#!/bin/sh
# *******************************************
# Script to perform T_only seq variant calling
# using a Tumor sample with fastq
# files named # tumor_1.fastq.gz, tumor_2.fastq.gz
# and data from a Panel of Normals
# *******************************************

# Update with the fullpath location of your sample fastq
fastq_folder="/home/goldenadmin/SECONDARY-ANALYSIS/NCSU_FASTQ/" #Change to local folder
tumor_fastq_1="$1"
tumor_fastq_2="$2"
tumor_sample="${3}_TG"
tumor_group="${3}_TG"
platform="ILLUMINA"

# Update with the location of the reference data files
fasta="/home/goldenadmin/SECONDARY-ANALYSIS/resources/Reference_Assemblies/canFam3.fa"


# We recommend that you create the panel of normal file with the corresponding algorithm that you plan to use for the somatic mutation calling. 
panel_of_normal_TNhaplotyper="/home/goldenadmin/SECONDARY-ANALYSIS/resources/panel_of_normals_haplotyper.vcf.gz" #Change to local folder

# Update with the location of the Sentieon software package and license file
export SENTIEON_INSTALL_DIR="/home/goldenadmin/SECONDARY-ANALYSIS/tools/sentieon/" #Change to local folder
#export SENTIEON_LICENSE=/home/Licenses/Sentieon.lic

# Other settings

nt=48 #number of threads to use in computation
workdir="/home/goldenadmin/SECONDARY-ANALYSIS/NCSU_FASTQ/$3" #Determine where the output files will be stored

# ******************************************
# 0. Setup
# ******************************************
mkdir -p "$workdir"
logfile="$workdir"/run.log
exec >"$logfile" 2>&1
cd "$workdir"

export SENTIEON_LICENSE="/home/goldenadmin/SECONDARY-ANALYSIS/GoldenHelix_cluster.lic" #Change to local folder

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
( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -M -R "@RG\tID:$tumor_group\tSM:$tumor_sample\tPL:$platform" -t $nt -K 10000000 "$fasta" "$fastq_folder/$tumor_fastq_1" "$fastq_folder/$tumor_fastq_2" || echo -n 'error' ) | $SENTIEON_INSTALL_DIR/bin/sentieon util sort -o ${tumor_sample}_sorted.bam -t $nt --sam2bam -i -
# ******************************************


# ******************************************
# 2a. Metrics for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i ${tumor_sample}_sorted.bam --algo MeanQualityByCycle ${tumor_sample}_mq_metrics.txt --algo QualDistribution ${tumor_sample}_qd_metrics.txt --algo GCBias --summary ${tumor_sample}_gc_summary.txt ${tumor_sample}_gc_metrics.txt --algo AlignmentStat --adapter_seq '' ${tumor_sample}_aln_metrics.txt --algo InsertSizeMetricAlgo ${tumor_sample}_is_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o ${tumor_sample}_gc-report.pdf ${tumor_sample}_gc_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o ${tumor_sample}_qd-report.pdf ${tumor_sample}_qd_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o ${tumor_sample}_mq-report.pdf ${tumor_sample}_mq_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o ${tumor_sample}_is-report.pdf ${tumor_sample}_is_metrics.txt


# ******************************************
# 4. Somatic Variant Calling
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i ${tumor_sample}_sorted.bam --algo TNhaplotyper --tumor_sample ${tumor_sample} --pon $panel_of_normal_TNhaplotyper ${tumor_sample}_tnhaplotyper.vcf.gz
