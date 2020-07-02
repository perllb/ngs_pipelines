#!/bin/bash

rundir=/trannel/dev/per/proj/TWIST/2020_58_Hedenfalk/scr
bamdir=/data/bnf/dev/per/proj/TWIST/2020_58_Hedenfalk/bam
qualdir=/data/bnf/dev/per/proj/TWIST/2020_58_Hedenfalk/QC

cd $rundir
mkdir -p 3_scr_realign_QC
mkdir -p $qualdir

cd $bamdir

SENTIEON_INSTALL_DIR=/data/bnf/sw/sentieon/sentieon-genomics-201808.01
REFERENCE=/data/bnf/ref/b37/human_g1k_v37_decoy.fasta
MODELFILE=/data/bnf/ref/sentieon/Sentieon_GiAB_HighAF_LowFP_201711.05.model

for bam in *.markdup.bam

do echo $bam

   realigned=$(echo $bam | sed 's/.bam/.realigned.bam/g')
   id=$(echo $bam | cut -f1 -d".")
   echo "BAM: $bam"
   echo "ID: $id"
   echo "realigned: $realigned"

   echo "#!/bin/bash -l
#SBATCH -p trannel
#SBATCH -t 65:00:00
#SBATCH -J ${id}Sention_real_QC
#SBATCH -N 1
#SBATCH -c 18
#SBATCH -e ${id}Sention_real_QC.err
#SBATCH -o ${id}Sention_real_QC.out

export SENTIEON_LICENSE=localhost:8990

echo $SENTIEON_INSTALL_DIR
echo $REFERENCE
echo $MODELFILE

$SENTIEON_INSTALL_DIR/bin/sentieon driver -t 18 -r /data/bnf/ref/b37/human_g1k_v37_decoy.fasta \
				      -i $bamdir/$bam --algo Realigner $bamdir/$realigned
   
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t 18 -r /data/bnf/ref/b37/human_g1k_v37_decoy.fasta \
				      -i $bamdir/$realigned --algo QualCal ${qualdir}/$id.RECAL_DATA.TABLE

$SENTIEON_INSTALL_DIR/bin/sentieon driver -t 18 -r /data/bnf/ref/b37/human_g1k_v37_decoy.fasta -i $bamdir/$realigned \
				      --algo GCBias --summary ${qualdir}/$id.GC_SUMMARY_TXT ${qualdir}/$id.GC_METRIC_TXT \
				      --algo MeanQualityByCycle ${qualdir}/$id.MQ_METRIC_TXT \
				      --algo QualDistribution ${qualdir}/$id.QD_METRIC_TXT \
				      --algo InsertSizeMetricAlgo ${qualdir}/$id.IS_METRIC_TXT  \
				      --algo AlignmentStat ${qualdir}/$id.ALN_METRIC_TXT
   
$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o ${qualdir}/$id.GC_METRIC.PDF ${qualdir}/$id.GC_METRIC_TXT
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o ${qualdir}/$id.MQ_METRIC.PDF ${qualdir}/$id.MQ_METRIC_TXT
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o ${qualdir}/$id.QD_METRIC.PDF ${qualdir}/$id.QD_METRIC_TXT
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o ${qualdir}/$id.IS_METRIC.PDF ${qualdir}/$id.IS_METRIC_TXT

" > $rundir/3_scr_realign_QC/${id}_realign_QC.sh
done
