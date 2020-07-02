#!/bin/bash

taskcpus=20
genome_file="/fs1/resources/ref/hg19/fasta/human_g1k_v37_decoy.fasta"
bedfile="/trannel/dev/per/proj/TWIST/2020_025_JKarlsson/Design_JK_202002_customTwist.sort.bed"
subscrdir="/trannel/dev/per/proj/TWIST/2020_025_JKarlsson/scr/2_scr_freebayes"
outfile="/data/bnf/dev/per/proj/TWIST/2020_025_JKarlsson/vcf/freebayes_cohortall.vcf"

mkdir -p /data/bnf/dev/per/proj/TWIST/2020_025_JKarlsson/vcf/
mkdir -p $subscrdir

bamdir="/data/bnf/dev/per/proj/TWIST/2020_025_JKarlsson/bam"
bamlist=$bamdir/bamfiles_to_vc.txt

all_bam=$(ls $bamdir/*umi.sort.bam)
echo $all_bam | sed 's/ /\n/g' > $bamlist

for bam in $bamlist
do echo "> BAM: $bam"

    bam_file=$(basename $bam)
    echo ">Basename: $bam_file"

    # Get ID + Sxx
    id=$(basename $bam_file .umi.sort.bam)
    echo ">ID: $id"

    # Get ID only 
    id=$(echo $_file | cut -f1 -d"_")
    echo ">ID: $id"

    type=$(if [[ $id == *"normal"* ]]; then echo normal; else echo tumor; fi)
    echo "> Type: $type"

    r1="$fqdir/${sid}_R1_001.fastq.gz"
    r2="$fqdir/${sid}_R2_001.fastq.gz"

#    echo "> R1: $r1"
#    echo "> R2: $r2"

   



echo "#!/bin/bash -l
#SBATCH -p trannel
#SBATCH -t 65:00:00
#SBATCH -J ${id}sentionQC
#SBATCH -N 1
#SBATCH -c $taskcpus
#SBATCH -e ${id}sentionQC.err
#SBATCH -o ${id}sentionQC.out

sentieon driver -t ${taskcpus} -r $genome_file -i $bam --algo QualCal ${id}.bqsr.table

 > $outfile" > $subscrdir/run_freebayes.sh
    


	     
