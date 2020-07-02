#!/bin/sh

rundir=/trannel/dev/per/proj/TWIST/2020_58_Hedenfalk/scr
bamdir=/data/bnf/dev/per/proj/TWIST/2020_58_Hedenfalk/bam
fastqdir=/data/NovaSeq/200512_A00681_0134_AHLT75DRXX/Data/Intensities/BaseCalls/2020_58
fastqdir2=/fs1/seqdata/NovaSeq/190920_A00681_0033_AHFYC3DRXX/Data/Intensities/BaseCalls/2019-72-IH

cd $fastqdir 

fqs=$(ls *R1*)

cd $fastqdir2
fqs2=$(ls *R1*)

cd $rundir
mkdir -p 1_scr_bwa

#
# BWA alignment raw fastq
#

for fq in $fqs 
do
    echo "Running $fq"
    sid=$(echo $fq | cut -f1 -d"R")
    id=$(echo $sid | cut -f1 -d"_")

    #
    # Rename files

    echo "ID: $id"             


    echo "New Name: $ctgid"

    fastqR1=${fastqdir}/${sid}R1_001.fastq.gz
    fastqR2=${fastqdir}/${sid}R2_001.fastq.gz

    echo "Fastq R1: $fastqR1"
    echo "Fastq R2: $fastqR2"

    echo "#!/bin/bash -l
#SBATCH -p normal
#SBATCH -t 65:00:00
#SBATCH -J ${sid}bwa
#SBATCH -N 1
#SBATCH -c 6
#SBATCH -e ${sid}bwa.err
#SBATCH -o ${sid}bwa.out

/data/bnf/sw/bwa/0.7.15/bwa mem -t 6 \
				-R \"@RG\tID:$id\tSM:$id\tPL:illumina\tLB:$id\tPU:$id\" \
				/data/bnf/ref/b37/human_g1k_v37_decoy.fasta \
				$fastqR1 \
				$fastqR2 > ${bamdir}/${id}.bwa_first.sam" > 1_scr_bwa/${id}_bwaRun.sh


done



for fq in $fqs2
do
    echo "Running $fq"
    sid=$(echo $fq | cut -f1 -d"R")
    id=$(echo $sid | cut -f1 -d"_")
    echo "ID: $id"

    fastqR1=${fastqdir2}/${sid}R1_001.fastq.gz
    fastqR2=${fastqdir2}/${sid}R2_001.fastq.gz

    echo "Fastq R1: $fastqR1"
    echo "Fastq R2: $fastqR2"

    echo "#!/bin/bash -l
#SBATCH -p normal
#SBATCH -t 65:00:00
#SBATCH -J ${sid}bwa
#SBATCH -N 1
#SBATCH -c 6
#SBATCH -e ${sid}bwa.err
#SBATCH -o ${sid}bwa.out

/data/bnf/sw/bwa/0.7.15/bwa mem -t 6 \
				-R \"@RG\tID:$id\tSM:$id\tPL:illumina\tLB:$id\tPU:$id\" \
				/data/bnf/ref/b37/human_g1k_v37_decoy.fasta \
				$fastqR1 \
				$fastqR2 > ${bamdir}/${id}.bwa_first.sam" > 1_scr_bwa/${id}_bwaRun.sh


done


