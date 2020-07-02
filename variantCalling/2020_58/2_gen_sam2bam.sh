#!/bin/bash

rundir=/trannel/dev/per/proj/TWIST/2020_58_Hedenfalk/scr
bamdir=/data/bnf/dev/per/proj/TWIST/2020_58_Hedenfalk/bam

cd $bamdir 

sams=$(ls *.sam)

cd $rundir
mkdir -p 2_scr_sam2bam_markdup

#
# sam to bam sorted
#

for sam in $sams

do

    echo "Running $sam"
    id=$(echo $sam | cut -f1 -d".")
    echo "ID: $id"
        
    echo "#!/bin/bash -l
#SBATCH -p normal
#SBATCH -t 65:00:00
#SBATCH -J ${id}sam2bam
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -e ${id}sam2bam.err
#SBATCH -o ${id}sam2bam.out

cd $bamdir

/data/bnf/sw/samtools/1.10/samtools view -S -b $bamdir/$sam > $bamdir/${id}.bam
/data/bnf/sw/samtools/1.10/samtools sort $bamdir/${id}.bam -o $bamdir/${id}.sort.bam 
/data/bnf/sw/samtools/1.10/samtools index $bamdir/${id}.sort.bam 

# MARKING DUPLICATES 
/data/bnf/sw/sambamba/0.6.7/sambamba markdup --tmpdir /data/tmp -t 6 $bamdir/${id}.sort.bam $bamdir/${id}.markdup.bam 
 " > 2_scr_sam2bam_markdup/${id}_sam2bam_markdup.sh

done
