#!/bin/bash

SENTIEON_INSTALL_DIR=/data/bnf/sw/sentieon/sentieon-genomics-201911

taskcpus=20
genome_file="/fs1/resources/ref/hg19/fasta/human_g1k_v37_decoy.fasta"
regions="/trannel/dev/per/proj/TWIST/2020_025_JKarlsson/Design_JK_202002_customTwist.sort.bed"

subscrdir="/trannel/dev/per/proj/TWIST/2020_025_JKarlsson/scr/1.2_scr_sent_BQSR_QC"
QCoutdir="/data/bnf/dev/per/proj/TWIST/2020_025_JKarlsson/sentQC/"
BQSRoutdir="/data/bnf/dev/per/proj/TWIST/2020_025_JKarlsson/sentBQSR/"

mkdir -p $QCoutdir
mkdir -p $BQSRoutdir

mkdir -p $subscrdir

bamdir="/data/bnf/dev/per/proj/TWIST/2020_025_JKarlsson/bam"

all_bam=$(ls $bamdir/*umi.sort.bam)

for bam in $all_bam
do echo "> BAM: $bam"

    bam_file=$(basename $bam)
    echo ">Basename: $bam_file"

    # Get ID 
    id=$(basename $bam_file .bwa.umi.sort.bam)
    echo ">ID: $id"

    type=$(if [[ $id == *"normal"* ]]; then echo normal; else echo tumor; fi)
    echo "> Type: $type"


    echo "#!/bin/bash -l
#SBATCH -p trannel
#SBATCH -t 65:00:00
#SBATCH -J ${id}sentionBQSR
#SBATCH -N 1
#SBATCH -c $taskcpus
#SBATCH -e ${id}sentionBQSR.err
#SBATCH -o ${id}sentionBQSR.out

export SENTIEON_LICENSE=localhost:8990

# BQSR generation
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t ${taskcpus} -r $genome_file -i $bam --algo QualCal $BQSRoutdir/${id}.bqsr.table

# QC 
sentieon driver \\
  --interval $regions -r $genome_file -t ${taskcpus} -i ${bam} \\
  --algo MeanQualityByCycle ${QCoutdir}/${id}_mq_metrics.txt --algo QualDistribution ${QCoutdir}/${id}_qd_metrics.txt \\
  --algo GCBias --summary ${QCoutdir}/${id}_gc_summary.txt ${QCoutdir}/${id}_gc_metrics.txt --algo AlignmentStat ${QCoutdir}/${id}_aln_metrics.txt \\
  --algo InsertSizeMetricAlgo ${QCoutdir}/${id}_is_metrics.txt \\
  --algo CoverageMetrics --cov_thresh 1 --cov_thresh 10 --cov_thresh 30 --cov_thresh 100 --cov_thresh 250 --cov_thresh 500 ${QCoutdir}/${id}_cov_metrics.txt


" > $subscrdir/${id}_sentBQSR.sh

done
    


	     
