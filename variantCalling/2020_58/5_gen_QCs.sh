#!/bin/sh

bamdir="/data/bnf/dev/per/proj/TWIST/2020_58_Hedenfalk/bam"
cd $bamdir
bamfiles=$(ls *sort.bam)

cd -

for file in $bamfiles
do echo $file
    echo "producing $file.."
    id=$(echo $file | cut -f1 -d".")

    # Post alignment QC

    echo "#!/bin/bash -l
#SBATCH -p normal
#SBATCH -t 65:00:00
#SBATCH -J ${sid}bwaQC
#SBATCH -N 1
#SBATCH -c 6
#SBATCH -e ${sid}bwaQC.err
#SBATCH -o ${sid}bwaQC.out

/data/bnf/scripts/postaln_qc.pl $bamdir/$file /data/bnf/ref/twist/twexome/Exome_RefSeq_SpikeIn_chrM.bed $id 6 /data/bnf/ref/twist/twexome/Exome_RefSeq_SpikeIn_chrM.bed /data/bnf/ref/b37/human_g1k_v37_decoy.fasta > /trannel/proj/ctg/hedenfalk_twist-exome/cmd-pipe/postmap/twist-exome/$id.2.bwa.QC " > 5_scr_QC/${id}_qc.sh

done
    

