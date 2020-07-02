#!/bin/bash

#
#
# NB! Run with "bash 1_gen_bwa.sh" NOT "sh .."
#
#

SENTIEON_INSTALL_DIR=/data/bnf/sw/sentieon/sentieon-genomics-201911


fqdir="/fs1/seqdata/NovaSeq/200616_A00681_0161_AHMLV5DRXX/Data/Intensities/BaseCalls/JKarlsson_2020_025/"

taskcpus=20
genome_file="/fs1/resources/ref/hg19/fasta/human_g1k_v37_decoy.fasta"

subscrdir="/trannel/dev/per/proj/TWIST/2020_025_JKarlsson/scr/1_scr_bwa"
mkdir $subscrdir

outdir="/data/bnf/dev/per/proj/TWIST/2020_025_JKarlsson/bam"

all_r1=$(ls $fqdir/*R1*gz)

for r1 in $all_r1
do
#    echo " ...."
#    echo "> $r1"
    r1_file=$(basename $r1)
#    echo ">Basename: $r1_file"

    # Get ID + Sxx
    sid=$(basename $r1_file _R1_001.fastq.gz)
    echo ">ID: $sid"

    # Get ID only 
    id=$(echo $r1_file | cut -f1 -d"_")
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
#SBATCH -J ${id}bwa
#SBATCH -N 1
#SBATCH -c $taskcpus
#SBATCH -e ${id}bwa.err
#SBATCH -o ${id}bwa.out

export SENTIEON_LICENSE=localhost:8990

export skip_coord_end=true

echo '>> Sentieon UMI extract, BWA mem, tee'

$SENTIEON_INSTALL_DIR/bin/sentieon umi extract -d 3M2S+T,3M2S+T $r1 $r2 \\
    | $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem \\
	    -R \"@RG\tID:$id\tSM:$id\tLB:$id\tPL:illumina\" \\
	    -t ${taskcpus} \\
	    -p -C $genome_file - \\
    | /data/bnf/sw/miniconda3/bin/tee -a $outdir/${id}_noumi.sam \\
    | $SENTIEON_INSTALL_DIR/bin/sentieon umi consensus -o $outdir/${id}_consensus.fastq.gz

echo '>> Sentieon BWA mem'


$SENTIEON_INSTALL_DIR/bin/sentieon bwa mem \\
    -R \"@RG\tID:$id\tSM:$id\tLB:$id\tPL:illumina\" \\
    -t ${taskcpus} \\
    -p -C $genome_file $outdir/${id}_consensus.fastq.gz \\
| $SENTIEON_INSTALL_DIR/bin/sentieon util sort -i - \\
    -o $outdir/${id}.${type}.bwa.umi.sort.bam \\
    --sam2bam

echo '>> Sentieon util sort'

$SENTIEON_INSTALL_DIR/bin/sentieon util sort -i $outdir/${id}_noumi.sam -o $outdir/${id}.${type}.bwa.sort.bam --sam2bam



rm $outdir/${id}_noumi.sam
	" > $subscrdir/${id}_run.bwa.sh
    
done

	     
