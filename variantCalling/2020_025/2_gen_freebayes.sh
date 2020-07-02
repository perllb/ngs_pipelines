#!/bin/bash

taskcpus=20
genome_file="/fs1/resources/ref/hg19/fasta/human_g1k_v37_decoy.fasta"
bedfile="/trannel/dev/per/proj/TWIST/2020_025_JKarlsson/Design_JK_202002_customTwist.sort.bed"
subscrdir="/trannel/dev/per/proj/TWIST/2020_025_JKarlsson/scr/2_scr_freebayes"
outfile="/data/bnf/dev/per/proj/TWIST/2020_025_JKarlsson/vcf/freebayes_cohortall.vcf"
outfilePre="/data/bnf/dev/per/proj/TWIST/2020_025_JKarlsson/vcf/freebayes_cohortall"
scrdir="/trannel/dev/per/proj/TWIST/2020_025_JKarlsson/scr/"
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

freebayes -f $genome_file -t $bedfile --pooled-continuous --pooled-discrete --min-repeat-entropy 1 -F 0.03 -L $bamlist > ${outfile}.raw

/data/bnf/sw/vcflib/bin/vcffilter -F LowCov -f \"DP > 500\" -f \"QA > 1500\" ${outfile}.raw | vcffilter -F LowFrq -o -f \"AB > 0.05\" -f \"AB = 0\" | vcfglxgt > ${outfilePRE}.filt1.vcf

${scrdir}/filter_freebayes_unpaired.pl ${outfilePRE}.filt1.vcf > ${outfile}

" > $subscrdir/run_freebayes.sh
    


	     
