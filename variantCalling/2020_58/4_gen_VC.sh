#!/bin/bash 

rundir=/trannel/dev/per/proj/TWIST/2020_58_Hedenfalk/scr
bamdir=/data/bnf/dev/per/proj/TWIST/2020_58_Hedenfalk/bam
homedir=/data/bnf/dev/per/proj/TWIST/2020_58_Hedenfalk

cd $rundir

mkdir -p 4_scr_VC
mkdir -p $homedir/vcf


SENTIEON_INSTALL_DIR=/data/bnf/sw/sentieon/sentieon-genomics-201808.01
REFERENCE=/data/bnf/ref/b37/human_g1k_v37_decoy.fasta
MODELFILE=/data/bnf/ref/sentieon/Sentieon_GiAB_HighAF_LowFP_201711.05.model


create_script(){
    ID=$1

    Tumor=$2
    Ttag=$3

    Normal=$4
    Ntag=$5

    pair=$6


    echo "ID: $ID"
    echo "Tumor: $Tumor"
    echo "Normal: $Normal"
    echo "pair: $pair"
    
    echo "#!/bin/bash -l
#SBATCH -p trannel
#SBATCH -t 65:00:00
#SBATCH -J ${ID}_${pair}_SenVC
#SBATCH -N 1
#SBATCH -c 18
#SBATCH -e ${ID}_${pair}_SenVC.err
#SBATCH -o ${ID}_${pair}_SenVC.out

export SENTIEON_LICENSE=localhost:8990

cd $homedir

#
# Paired variant calling
#

$SENTIEON_INSTALL_DIR/bin/sentieon driver --interval /data/bnf/ref/twist/twexome/Exome_RefSeq_SpikeIn_chrM.bed -t 12 -r $REFERENCE -i bam/$Tumor.markdup.realigned.bam -q QC/$Tumor.RECAL_DATA.TABLE -i bam/$Normal.markdup.realigned.bam -q QC/$Normal.RECAL_DATA.TABLE  --algo TNscope --tumor_sample $Ttag --normal_sample $Ntag --dbsnp /data/bnf/ref/b37/dbsnp_138.b37.vcf.gz vcf/${ID}_${pair}.vcf


#
# Apply machine learning filter
#


$SENTIEON_INSTALL_DIR/bin/sentieon driver -t 12 -r $REFERENCE --algo TNModelApply -v vcf/${ID}_${pair}.vcf -m $MODELFILE vcf/${ID}_${pair}.ml.vcf

#
# Annotate vcf
#

/data/bnf/scripts/annotate_vcf_vep2.pl vcf/${ID}_${pair}.ml.vcf vcf/${ID}_${pair}.ml.vep.vcf" > ./4_scr_VC/${ID}_${pair}_vc.sh

}


####### L13 ############

# T - N
metID=L13
ID=CTG-${metID}
Ttag=${metID}-T
Ntag=${metID}-N

Tumor=$ID-T
Normal=$ID-N

pair="T-N"

create_script $ID $Tumor $Ttag $Normal $Ntag $pair

# E - N 
metID=L13
ID=CTG-${metID}
Ttag=${metID}E
Ntag=${metID}-N

Tumor=$ID-E
Normal=$ID-N
pair="E-N"

create_script $ID $Tumor $Ttag $Normal $Ntag $pair





####### L15 ############

# T - N
metID=L15
ID=CTG-${metID}
Ttag=${metID}-T2
Ntag=${metID}-N

Tumor=$ID-T2
Normal=$ID-N
pair="T-N"

create_script $ID $Tumor $Ttag $Normal $Ntag $pair

# E - N 
metID=L15
ID=CTG-${metID}
Ttag=${metID}E
Ntag=${metID}-N

Tumor=$ID-E
Normal=$ID-N
pair="E-N"

create_script $ID $Tumor $Ttag $Normal $Ntag $pair




####### L18 ############

# T - N
metID=L18
ID=CTG-${metID}
Ttag=${metID}T
Ntag=${metID}N

Tumor=$ID-T
Normal=$ID-N
pair="T-N"

create_script $ID $Tumor $Ttag $Normal $Ntag $pair

# E - N 
metID=L18
ID=CTG-${metID}
Ttag=${metID}E
Ntag=${metID}N

Tumor=$ID-E
Normal=$ID-N
pair="E-N"

create_script $ID $Tumor $Ttag $Normal $Ntag $pair






####### L23 ############

# T - N
metID=L23
ID=CTG-${metID}
Ttag=${metID}-T
Ntag=${metID}-N

Tumor=$ID-T
Normal=$ID-N
pair="T-N"

create_script $ID $Tumor $Ttag $Normal $Ntag $pair


# E - N
metID=L23
ID=CTG-${metID}
Ttag=${metID}E
Ntag=${metID}-N

Tumor=$ID-E
Normal=$ID-N
pair="E-N"

create_script $ID $Tumor $Ttag $Normal $Ntag $pair






####### L32 ############

# T - N
metID=L32
ID=CTG-${metID}
Ttag=${metID}-T
Ntag=${metID}-N

Tumor=$ID-T
Normal=$ID-N
pair="T-N"

create_script $ID $Tumor $Ttag $Normal $Ntag $pair


# E - N 
metID=L32
ID=CTG-${metID}
Ttag=${metID}E
Ntag=${metID}-N

Tumor=$ID-E
Normal=$ID-N
pair="E-N"

create_script $ID $Tumor $Ttag $Normal $Ntag $pair

 
