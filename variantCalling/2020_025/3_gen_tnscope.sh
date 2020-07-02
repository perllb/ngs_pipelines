#!/bin/bash

SENTIEON_INSTALL_DIR=/data/bnf/sw/sentieon/sentieon-genomics-201911

taskcpus=20
genome_file="/fs1/resources/ref/hg19/fasta/human_g1k_v37_decoy.fasta"
regions="/trannel/dev/per/proj/TWIST/2020_025_JKarlsson/Design_JK_202002_customTwist.sort.bed"

scrdir="/trannel/dev/per/proj/TWIST/2020_025_JKarlsson/scr/"
subscrdir="${scrdir}/3_scr_tnscope"

outdir="/data/bnf/dev/per/proj/TWIST/2020_025_JKarlsson/vcf_tnscope/"
BQSRoutdir="/data/bnf/dev/per/proj/TWIST/2020_025_JKarlsson/sentBQSR/"

mkdir -p $outdir
mkdir -p $subscrdir

## Get Sample IDs
bamdir="/data/bnf/dev/per/proj/TWIST/2020_025_JKarlsson/bam"
all_bam=$(ls $bamdir/*umi.sort.bam)

# get normal samples to get unique IDs
for bam in $all_bam
do 
    if [[ $bam == *"normal"* ]]; then
	echo
	echo 
	echo ">BAM: $bam"
	echo ">Type: normal"

	nid=$(basename $bam .normal.bwa.umi.sort.bam)
	echo "> Normal ID: $nid"
	# Patient ID
	pid=$(echo $nid | cut -f1 -d"-")
	echo "> PID: $pid"

	# get all tumor samples for corresponding normal id
	Tbams=$(echo $bamdir/${pid}*tumor*.umi.sort.bam)
	echo "> Tumor bams:"

	for Tbam in $Tbams
	do
	    echo ">"
	    echo $Tbam
	    tid=$(basename $Tbam .tumor.bwa.umi.sort.bam)
	    echo $tid

	    # Get bqsr files for tumor and normal sample
	    Tbqsr=$(echo $BQSRoutdir/${tid}*bqsr.table)
	    Nbqsr=$(echo $BQSRoutdir/${nid}*bqsr.table)

	  
	    
	    # Generate script!
	    echo "#!/bin/bash -l
#SBATCH -p trannel
#SBATCH -t 65:00:00
#SBATCH -J ${tid}tnscope
#SBATCH -N 1
#SBATCH -c $taskcpus
#SBATCH -e ${tid}tnscope.err
#SBATCH -o ${tid}tnscope.out

export SENTIEON_LICENSE=localhost:8990


$SENTIEON_INSTALL_DIR/bin/sentieon driver -t ${taskcpus} \\
  -r $genome_file \\
  -i $Tbam  -q ${Tbqsr} \\
  -i $bam -q ${Nbqsr} \\
  --interval $regions --algo TNscope \\
  --tumor_sample $tid --normal_sample $nid \\
  --clip_by_minbq 1 --max_error_per_read 3 --min_init_tumor_lod 2.0 \\
  --min_base_qual 10 --min_base_qual_asm 10 --min_tumor_allele_frac 0.0005 \\
  $outdir/tnscope_${tid}_${nid}.vcf.raw

echo
echo '>filter...'
echo 

perl ${scrdir}/filter_tnscope_somatic.pl ${outdir}/tnscope_${tid}_${nid}.vcf.raw ${tid} ${nid} > ${outdir}/tnscope_${tid}_${nid}.filter.vcf	    
	
" > $subscrdir/run_tnscope_${tid}_${nid}.sh



	done
    fi
done
