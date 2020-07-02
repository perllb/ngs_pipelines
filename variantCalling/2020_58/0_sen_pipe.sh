#
# Create dirs and set paths/license
#

SENTIEON_INSTALL_DIR=/data/bnf/sw/sentieon/sentieon-genomics-201808.01
export SENTIEON_LICENSE=mtlucmds2.lund.skane.se:8990
REFERENCE="/data/bnf/ref/b37/human_g1k_v37_decoy.fasta"
MODELFILE="/data/bnf/ref/sentieon/Sentieon_GiAB_HighAF_LowFP_201711.05.model"

mkdir -p qual bam vcf

#
# Define samples IDs
#
Tumor=CTG-L18-T
Normal=CTG-L18-N

#
# Path to fastq
#
T_path=/trannel/proj/ctg/hedenfalk_twist-exome/cmd-pipe/bam/twist-exome/$Tumor.bwa.sort.bam
N_path=/trannel/proj/ctg/hedenfalk_twist-exome/cmd-pipe/bam/twist-exome/$Normal.bwa.sort.bam


#
# Align using standard BWA

#$Tumor/trannel/proj/ctg/hedenfalk_twist-exome/cmd-pipe/bam/twist-exome/$Tumor.bwa.sort.bam
#$Normal/trannel/proj/ctg/hedenfalk_twist-exome/cmd-pipe/bam/twist-exome/$Normal.bwa.sort.bam

/data/bnf/sw/bwa/0.7.15/bwa mem -t 6 -R @RG\tID:$Tumor\tSM:$Tumor\tPL:illumina\tLB:$Tumor\tPU:$Tumor /data/bnf/ref/b37/human_g1k_v37_decoy.fasta /trannel/proj/ctg//190920_A00681_0033_AHFYC3DRXX/Data/Intensities/BaseCalls/2019-72-IH/L32-T/L32-T_S24_L002_R2_001.fastq.gz /trannel/proj/ctg//190920_A00681_0033_AHFYC3DRXX/Data/Intensities/BaseCalls/2019-72-IH/L32-T/L32-T_S24_L002_R1_001.fastq.gz
# also for L32-N


#
# Realign and bam, QC tumor
#

$SENTIEON_INSTALL_DIR/bin/sentieon driver -t 18 -r /data/bnf/ref/b37/human_g1k_v37_decoy.fasta \
				   -i /trannel/proj/ctg/hedenfalk_twist-exome/cmd-pipe/bam/twist-exome/$Tumor.bwa.sort.bam --algo Realigner bam/$Tumor.realigned.bam
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t 18 -r /data/bnf/ref/b37/human_g1k_v37_decoy.fasta \
				   -i bam/$Tumor.realigned.bam --algo QualCal qual/$Tumor.RECAL_DATA.TABLE
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t 18 -r /data/bnf/ref/b37/human_g1k_v37_decoy.fasta -i bam/$Tumor.realigned.bam \
				   --algo GCBias --summary qual/$Tumor.GC_SUMMARY_TXT qual/$Tumor.GC_METRIC_TXT \
				   --algo MeanQualityByCycle qual/$Tumor.MQ_METRIC_TXT \
				   --algo QualDistribution qual/$Tumor.QD_METRIC_TXT \
				   --algo InsertSizeMetricAlgo qual/$Tumor.IS_METRIC_TXT  \
				   --algo AlignmentStat qual/$Tumor.ALN_METRIC_TXT

$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o qual/$Tumor.GC_METRIC.PDF qual/$Tumor.GC_METRIC_TXT
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o qual/$Tumor.MQ_METRIC.PDF qual/$Tumor.MQ_METRIC_TXT
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o qual/$Tumor.QD_METRIC.PDF qual/$Tumor.QD_METRIC_TXT
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o qual/$Tumor.IS_METRIC.PDF qual/$Tumor.IS_METRIC_TXT

#
# Realign and bam QC, Normal sample
#

$SENTIEON_INSTALL_DIR/bin/sentieon driver -t 18 -r /data/bnf/ref/b37/human_g1k_v37_decoy.fasta \
				   -i /trannel/proj/ctg/hedenfalk_twist-exome/cmd-pipe/bam/twist-exome/$Normal.bwa.sort.bam --algo Realigner bam/$Normal.realigned.bam
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t 18 -r /data/bnf/ref/b37/human_g1k_v37_decoy.fasta \
				   -i bam/$Normal.realigned.bam --algo QualCal qual/$Normal.RECAL_DATA.TABLE
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t 18 -r /data/bnf/ref/b37/human_g1k_v37_decoy.fasta -i bam/$Normal.realigned.bam \
				   --algo GCBias --summary qual/$Normal.GC_SUMMARY_TXT qual/$Normal.GC_METRIC_TXT \
				   --algo MeanQualityByCycle qual/$Normal.MQ_METRIC_TXT \
				   --algo QualDistribution qual/$Normal.QD_METRIC_TXT \
				   --algo InsertSizeMetricAlgo qual/$Normal.IS_METRIC_TXT  \
				   --algo AlignmentStat qual/$Normal.ALN_METRIC_TXT

$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o qual/$Normal.GC_METRIC.PDF qual/$Normal.GC_METRIC_TXT
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o qual/$Normal.MQ_METRIC.PDF qual/$Normal.MQ_METRIC_TXT
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o qual/$Normal.QD_METRIC.PDF qual/$Normal.QD_METRIC_TXT
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o qual/$Normal.IS_METRIC.PDF qual/$Normal.IS_METRIC_TXT


#
# Paired variant calling
#

$SENTIEON_INSTALL_DIR/bin/sentieon driver --interval /data/bnf/ref/twist/twexome/Exome_RefSeq_SpikeIn_chrM.bed -t 12 -r $REFERENCE -i bam/$Tumor.realigned.bam -q qual/$Tumor.RECAL_DATA.TABLE -i bam/$Normal.realigned.bam -q qual/$Normal.RECAL_DATA.TABLE  --algo TNscope --tumor_sample $Tumor --normal_sample $Normal --dbsnp /data/bnf/ref/b37/dbsnp_138.b37.vcf.gz vcf/CTG-L32.vcf


#
# Apply machine learning filter
#


$SENTIEON_INSTALL_DIR/bin/sentieon driver -t 12 -r $REFERENCE --algo TNModelApply -v vcf/CTG-L32.vcf -m $MODELFILE vcf/CTG-L32.ml.vcf

#
# Annotate vcf
#

/data/bnf/scripts/annotate_vcf_vep2.pl vcf/CTG-L32.ml.vcf vcf/CTG-L32.ml.vep.vcf
