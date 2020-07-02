#!/bin/bash

taskcpus=20
genome_file="/fs1/resources/ref/hg19/fasta/human_g1k_v37_decoy.fasta"
regions="/trannel/dev/per/proj/TWIST/2020_025_JKarlsson/Design_JK_202002_customTwist.sort.bed"

scrdir="/trannel/dev/per/proj/TWIST/2020_025_JKarlsson/scr/"
subscrdir="${scrdir}/4_scr_vardict"

outdir="/data/bnf/dev/per/proj/TWIST/2020_025_JKarlsson/vcf_vardict/"

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
	    echo "> Tumor Bam: $Tbam"
	    tid=$(basename $Tbam .tumor.bwa.umi.sort.bam)
	    echo "> Tumor ID: $tid"

	    # Generate script!
	    echo "#!/bin/bash -l
#SBATCH -p trannel
#SBATCH -t 65:00:00
#SBATCH -J ${tid}vardict
#SBATCH -N 1
#SBATCH -c $taskcpus
#SBATCH -e ${tid}vardict.err
#SBATCH -o ${tid}vardict.out

vardict-java -G $genome_file -f 0.01 -N ${tid} -b \"${Tbam}|${bam}\" -c 1 -S 2 -E 3 -g 4 -U $regions \\
| /data/bnf/sw/VarDict/testsomatic.R | /data/bnf/sw/VarDict/var2vcf_paired.pl -N \"${tid}|${nid}\" -f 0.01 > ${outdir}/vardict_${tid}.vcf.raw

echo
echo '>filter...'
echo 

perl ${scrdir}/filter_vardict_somatic.pl ${outdir}/vardict_${tid}.vcf.raw ${tid} ${nid} > ${outdir}/vardict_${tid}.vcf 
" > ${subscrdir}/run_vardict_${tid}_${nid}.sh



	done
    fi
done
