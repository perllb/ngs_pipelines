#!/bin/bash

homedir="/fs1/per/proj/scRNAseq/iCell8/200423_els/"
subscr="$homedir/subscr"
outdir="$homedir/outdir"
rundir="/fs1/seqdata/NovaSeq/200423_A00681_0119_BHLVJ2DRXX"

contCR="/fs1/per/containers_singularity/cellranger-rna/cellranger-rna_3.1.0_pipe.sif"
contiCell="/fs1/per/containers_singularity/iCell8x/icell8x.sif"

mkdir -p $subscr
mkdir -p $outdir

bcl2fastq=0
mappa=0

usage="Usage:
> sh run_icell8x.sh -m <to run mappa_demuxer> -b <to run bcl2fastq> -h (help)"

while getopts ":b:mh" OPTION; do
    case "$OPTION" in
	b) bcl2fastq=1 ;;
	m) mappa=1 ;;
	h)
	    echo "$usage"
	    exit ;;
	:) echo "missing argument for -%s\n" "$OPTARG" >&2
	    echo "$usage" >&2
	    exit 1
	    ;;
	\?) echo "illegal option: -%s\n" "$OPTARG" >&2
	    echo "$usage" >&2
	    exit 1
       ;;
    esac
done

if [ $OPTIND -eq 1 ]; then echo "No options were passed"; echo $usage; exit; fi
shift $((OPTIND-1))
echo "$# non-option arguments"


echo ">> Set up for running: "
echo "> Bcl2fastq (-b):  $bcl2fastq"
echo "> Mappa (-m):      $mappa"

### BCL2FASTQ ###
if [ $bcl2fastq == 1 ]
then

    echo ">> Running BCL2FASTQ"

    jname="ic8_bcl2"
    echo "#!/bin/bash -l
#SBATCH -p normal
#SBATCH -t 5:00:00
#SBATCH -J $jname
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem 189840
#SBATCH -e $jname.err
#SBATCH -o $jname.out

ml singularity

cd $rundir 

singularity exec $cont bcl2fastq" > $subscr/$jname.run.sh

    sbatch $subscr/$jname.run.sh

    ## Move fastq files to curr dir
    mkdir -p $outdir/fastq
    mv /fs1/seqdata/NovaSeq/200423_A00681_0119_BHLVJ2DRXX/Data/Intensities/BaseCalls/Und*fastq.gz $outdir/fastq

fi

### MAPPA ####

if [ $mappa == 1 ]
then
    echo ">> Running Mappa"

    jname="ic8_map_demux"

    echo "#!/bin/bash -l
#SBATCH -p normal
#SBATCH -t 65:00:00
#SBATCH -J $jname
#SBATCH -N 1
#SBATCH -c 56
#SBATCH --mem 189840
#SBATCH -e $jname.err
#SBATCH -o $jname.out

ml singularity

cd $homedir

singularity exec $contiCell mappa_demuxer.py \
   -i outdir/fastq/Undetermined_S0_L002_R1_001.fastq.gz \
   -p outdir/fastq/Undetermined_S0_L002_R2_001.fastq.gz \
   -b 200423_Els_iCell8cx_test_2/117117-200402_Els_WellList.i5RevComp.TXT \
   -t Full_length \
   -o $outdir/mappa \
   -n 56  " > $subscr/$jname.run.sh

    sbatch $subscr/$jname.run.sh
    
fi

