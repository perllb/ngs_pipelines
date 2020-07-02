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

cat $bamlist

echo "#!/bin/bash -l
#SBATCH -p trannel
#SBATCH -t 65:00:00
#SBATCH -J ${id}freebayes
#SBATCH -N 1
#SBATCH -c $taskcpus
#SBATCH -e ${id}freebayes.err
#SBATCH -o ${id}freebayes.out

freebayes -f $genome_file -t $bedfile --pooled-continuous --pooled-discrete --min-repeat-entropy 1 -F 0.03 -L $bamlist > $outfile" > $subscrdir/run_freebayes.sh
    


	     
