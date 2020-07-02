#
# Create dirs and set paths/license
#

SENTIEON_INSTALL_DIR=/data/bnf/sw/sentieon/sentieon-genomics-201808.01
export SENTIEON_LICENSE=mtlucmds2.lund.skane.se:8990
REFERENCE="/data/bnf/ref/b37/human_g1k_v37_decoy.fasta"
MODELFILE="/data/bnf/ref/sentieon/Sentieon_GiAB_HighAF_LowFP_201711.05.model"

mkdir -p qual bam vcf


#
# Align using standard BWA
#

#CTG-L32-T/trannel/proj/ctg/hedenfalk_twist-exome/cmd-pipe/bam/twist-exome/CTG-L32-T.bwa.sort.bam
#CTG-L32-N/trannel/proj/ctg/hedenfalk_twist-exome/cmd-pipe/bam/twist-exome/CTG-L32-N.bwa.sort.bam

#/data/bnf/sw/bwa/0.7.15/bwa mem -t 6 -R @RG\tID:CTG-L32-T\tSM:CTG-L32-T\tPL:illumina\tLB:CTG-L32-T\tPU:CTG-L32-T /data/bnf/ref/b37/human_g1k_v37_decoy.fasta /trannel/proj/ctg//190920_A00681_0033_AHFYC3DRXX/Data/Intensities/BaseCalls/2019-72-IH/L32-T/L32-T_S24_L002_R2_001.fastq.gz /trannel/proj/ctg//190920_A00681_0033_AHFYC3DRXX/Data/Intensities/BaseCalls/2019-72-IH/L32-T/L32-T_S24_L002_R1_001.fastq.gz
# also for L32-N


#
# Realign and bam, QC tumor
#

$SENTIEON_INSTALL_DIR/bin/sentieon driver -t 18 -r /data/bnf/ref/b37/human_g1k_v37_decoy.fasta \
				   -i /trannel/proj/ctg/hedenfalk_twist-exome/cmd-pipe/bam/twist-exome/CTG-L32-T.bwa.sort.bam --algo Realigner bam/CTG-L32-T.realigned.bam
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t 18 -r /data/bnf/ref/b37/human_g1k_v37_decoy.fasta \
				   -i bam/CTG-L32-T.realigned.bam --algo QualCal qual/CTG-L32-T.RECAL_DATA.TABLE
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t 18 -r /data/bnf/ref/b37/human_g1k_v37_decoy.fasta -i bam/CTG-L32-T.realigned.bam \
				   --algo GCBias --summary qual/CTG-L32-T.GC_SUMMARY_TXT qual/CTG-L32-T.GC_METRIC_TXT \
				   --algo MeanQualityByCycle qual/CTG-L32-T.MQ_METRIC_TXT \
				   --algo QualDistribution qual/CTG-L32-T.QD_METRIC_TXT \
				   --algo InsertSizeMetricAlgo qual/CTG-L32-T.IS_METRIC_TXT  \
				   --algo AlignmentStat qual/CTG-L32-T.ALN_METRIC_TXT

$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o qual/CTG-L32-T.GC_METRIC.PDF qual/CTG-L32-T.GC_METRIC_TXT
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o qual/CTG-L32-T.MQ_METRIC.PDF qual/CTG-L32-T.MQ_METRIC_TXT
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o qual/CTG-L32-T.QD_METRIC.PDF qual/CTG-L32-T.QD_METRIC_TXT
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o qual/CTG-L32-T.IS_METRIC.PDF qual/CTG-L32-T.IS_METRIC_TXT

#
# Realign and bam QC, Normal sample
#

$SENTIEON_INSTALL_DIR/bin/sentieon driver -t 18 -r /data/bnf/ref/b37/human_g1k_v37_decoy.fasta \
				   -i /trannel/proj/ctg/hedenfalk_twist-exome/cmd-pipe/bam/twist-exome/CTG-L32-N.bwa.sort.bam --algo Realigner bam/CTG-L32-N.realigned.bam
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t 18 -r /data/bnf/ref/b37/human_g1k_v37_decoy.fasta \
				   -i bam/CTG-L32-N.realigned.bam --algo QualCal qual/CTG-L32-N.RECAL_DATA.TABLE
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t 18 -r /data/bnf/ref/b37/human_g1k_v37_decoy.fasta -i bam/CTG-L32-N.realigned.bam \
				   --algo GCBias --summary qual/CTG-L32-N.GC_SUMMARY_TXT qual/CTG-L32-N.GC_METRIC_TXT \
				   --algo MeanQualityByCycle qual/CTG-L32-N.MQ_METRIC_TXT \
				   --algo QualDistribution qual/CTG-L32-N.QD_METRIC_TXT \
				   --algo InsertSizeMetricAlgo qual/CTG-L32-N.IS_METRIC_TXT  \
				   --algo AlignmentStat qual/CTG-L32-N.ALN_METRIC_TXT

$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o qual/CTG-L32-N.GC_METRIC.PDF qual/CTG-L32-N.GC_METRIC_TXT
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o qual/CTG-L32-N.MQ_METRIC.PDF qual/CTG-L32-N.MQ_METRIC_TXT
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o qual/CTG-L32-N.QD_METRIC.PDF qual/CTG-L32-N.QD_METRIC_TXT
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o qual/CTG-L32-N.IS_METRIC.PDF qual/CTG-L32-N.IS_METRIC_TXT


#
# Paired variant calling
#

$SENTIEON_INSTALL_DIR/bin/sentieon driver --interval /data/bnf/ref/twist/twexome/Exome_RefSeq_SpikeIn_chrM.bed -t 12 -r $REFERENCE -i bam/CTG-L32-T.realigned.bam -q qual/CTG-L32-T.RECAL_DATA.TABLE -i bam/CTG-L32-N.realigned.bam -q qual/CTG-L32-N.RECAL_DATA.TABLE  --algo TNscope --tumor_sample CTG-L32-T --normal_sample CTG-L32-N --dbsnp /data/bnf/ref/b37/dbsnp_138.b37.vcf.gz vcf/CTG-L32.vcf


#
# Apply machine learning filter
#


$SENTIEON_INSTALL_DIR/bin/sentieon driver -t 12 -r $REFERENCE --algo TNModelApply -v vcf/CTG-L32.vcf -m $MODELFILE vcf/CTG-L32.ml.vcf

#
# Annotate vcf
#

/data/bnf/scripts/annotate_vcf_vep2.pl vcf/CTG-L32.ml.vcf vcf/CTG-L32.ml.vep.vcf
