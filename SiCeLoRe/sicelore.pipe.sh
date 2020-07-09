# BASH Pipeline to analyse single cell ONT + short reads 
# Using the Sicelore suite

ilpar=0
fstp=0
minim=0
addGNT=0
tagger=0
mergeBCUt=0
consensus=0

while getopts ":ifmatbc" OPTION; do
    case "$OPTION" in
	i) ilpar=1 ;;
	f) fstp=1 ;;
	m) minim=1 ;;
	a) addGNT=1 ;;
	t) tagger=1 ;;
	b) mergeBCUt=1 ;;
	c) consensus=1 ;;
    esac
done

if [ $OPTIND -eq 1 ]
then
    echo " "
    echo " "
    echo ">> ERROR: Pass at least one option for what module to run"
    echo "Sicelore pipeline usage:"
    echo "> sh sicelore.pipeline -i -f -m -a -t -c "
    echo " "
    echo "- Illumina Parser (-i) "
    echo "- Fastp (-f) "
    echo "- Minimap2 (-m) "
    echo "- AddGeneNameTag (-a) "
    echo "- BC+UMI tagging (-t) "
    echo "- Merge BC+UMI tag reads (-b) "
    echo "- Consensus gen. (-c) "
    exit 1;
fi

echo ">> Set up for running: "
echo "> Illumina Parser (-i):  $ilpar"
echo "> Fastp (-f):            $fstp"
echo "> Minimap2 (-m):         $minim"
echo "> AddGeneNameTag (-a):   $addGNT"
echo "> BC+UMI tagging (-t):   $tagger"
echo "> Merge BC+UMItag (-b):  $mergeBCUt"
echo "> Consensus gen (-c):    $consensus"
echo ". "
echo ".."
echo "..." 
 
# 0.1 - Define dirs
homedir="/fs1/per/proj/scRNAseq/Nanopore/Nano_10x/Run2_SampleA/"
scrdir="$homedir/sicelorePipe/scr"

# 10x home
x10home="$homedir/data/10x/"
# 10x bamfile
bam="${x10home}/possorted_genome_bam.bam"
# 10x barcode file
barcode="${x10home}/filtered_feature_bc_matrix/barcodes.tsv"

# Outdir for sicelore pipeline
outdir_sic="$homedir/outdir_fulldata/sicelore/"
mkdir -p $outdir_sic

# Input ONT fastq path
ONTpath="$homedir/data/Nanopore/"

#Minimap2 parameters
build='/opt/bin/sicelore/Gencode/gencode.v31.hg38'
genc="${build}.junctions.bed"
mmi='/opt/bin/sicelore/Data/hg38.mmi'
ref='/opt/bin/ref/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa'

# addGeneNameTag parameters
cellRef='/opt/bin/sicelore/Gencode/gencode.v31.hg38.refFlat.txt'

# Set singularity
ml singularity
singCnt="/fs1/per/containers_singularity/SiCeLoRe/SiCeLoRe.sif"
singCmd="singularity exec --bind /fs1/ $singCnt"

subscr="${scrdir}/shell_sicelorePipe/"
mkdir -p $subscr

# 1.1 - First sicelore step - illumina parser
if [ $ilpar == 1 ]
then

    jname="ill_parser"

    echo "> 1. Illumina Parser running"
    echo "> -jobname: $jname"

    if [ -f ${barcode}.gz ] 
    then
        zcat ${barcode}.gz > $barcode
    fi
    
    echo "> BAM: $bam"
    echo $bam

    echo "#!/bin/bash -l
#SBATCH -p normal
#SBATCH -t 40:00:00
#SBATCH -J $jname
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem 189840
#SBATCH -e $jname.err
#SBATCH -o $jname.out

$singCmd \
  java -Xmx15000m -jar /opt/bin/sicelore/Jar/IlluminaParser-1.0.jar \
              --inFileIllumina $bam \
              --tsv $barcode \
              --outFile $outdir_sic/parsedForNanopore_v1.0.obj --cellBCflag CB --umiFlag UB --geneFlag GN" > $subscr/$jname.run.sh
    
    sbatch $subscr/$jname.run.sh

fi


# 2 - Map Nanopore reads to reference genome (minimap)
# 2.1 Fastp 
# - split fastq into chunks for paralellization
fastpout=$outdir_sic/2_fastp

if [ $fstp == 1 ]
then

    jname="fastp"

    echo "> 2. Fastp running"
    echo "> -jobname: $jname"
    echo "> -output fir: $fastpout"

    mkdir -p $fastpout

    for file in $ONTpath/*fastq
    do

	baseFile=$(basename $file .fastq)
	echo "#!/bin/bash -l
#SBATCH -p normal
#SBATCH -t 5:00:00
#SBATCH -J $jname
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem 189840
#SBATCH -e $jname.err
#SBATCH -o $jname.out

$singCmd \
        fastp -i $file -Q -A --thread 20 --split_prefix_digits=3 --out1=$fastpout/$baseFile.fastq --split=24
    " > $subscr/$jname.$baseFile.run.sh

	sbatch $subscr/$jname.$baseFile.run.sh
	
    done

fi


# 2.2 Minimap2
# - Map the split fastq files to genome

# 2.2.2 Run minimap2
minimout=$outdir_sic/2_minimap2

if [ $minim == 1 ]
then

    jname="minimap"

    echo "> 2. Minimap2 running"
    echo "> -jobname: $jname"
    echo "> -output fir: $minimout"

    mkdir -p $minimout

    for file in $fastpout/*fastq
    do

        splt=$(basename $file .fastq)
        echo "#!/bin/bash -l
#SBATCH -p normal
#SBATCH -t 5:00:00
#SBATCH -J $jname.$splt
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem 189840
#SBATCH -e $jname.$splt.err
#SBATCH -o $jname.$splt.out


cd $fastpout

$singCmd minimap2 -ax splice -uf --MD -N 100 --sam-hit-only -t 20 --junc-bed ${genc}.junctions.bed $ref $splt.fastq > $minimout/${splt}.sam  

cd $minimout

$singCmd samtools view -Sb ${splt}.sam -o ${splt}.unsorted.bam

$singCmd samtools sort ${splt}.unsorted.bam -o ${splt}.bam

$singCmd samtools index ${splt}.bam ;

rm ${splt}.unsorted.bam

cd $srcdir 
    " > $subscr/$jname.$splt.run.sh

        sbatch $subscr/$jname.$splt.run.sh

    done

fi


# 3. Tag Nanopore SAM records 
# 3.1 With gene names
# 3.2 With QV and read sequences (from fastq splits fastp)
minimout=$outdir_sic/2_minimap2
tagout=$outdir_sic/3_tagBams

if [ $addGNT == 1 ]
then

    mkdir -p $tagout

    jname="addGeneNameTag"

    echo "> 3. AddGeneNameTag running"
    echo "> -jobname: $jname"
    echo "> -output dir: $tagout"

    for file in $minimout/*bam
    do

        splt=$(basename $file .bam)
	GEfile=$tagout/$splt.GE.bam
	geUSfile=$tagout/$splt.GEUS.bam

        echo "#!/bin/bash -l
#SBATCH -p normal
#SBATCH -t 5:00:00
#SBATCH -J $jname.$splt
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem 189840
#SBATCH -e $jname.$splt.err
#SBATCH -o $jname.$splt.out

$singCmd java -jar -Xmx12g  /opt/bin/sicelore/Jar/Sicelore-1.0.jar AddGeneNameTag I=$file O=$GEfile REFFLAT=$cellRef GENETAG=GE

$singCmd samtools index $GEfile

## TAG QV and readSequences (GEUS)

$singCmd java -jar -Xmx12g /opt/bin/sicelore/Jar/Sicelore-1.0.jar AddBamReadSequenceTag I=$GEfile O=$geUSfile FASTQ=$fastpout/$splt.fastq

$singCmd samtools index $geUSfile
    " > $subscr/$jname.$splt.run.sh

        sbatch $subscr/$jname.$splt.run.sh

    done

fi


### 4. BC and UMI finder
umiout=${outdir_sic}/4_bc_umi_bams
ill_pars=$outdir_sic/parsedForNanopore_v1.0.obj 

if [  $tagger == 1 ]
then

    mkdir -p $umiout

    jname="BC_UMI_tagger"

    echo "> 4. BC_UMI tagging running"
    echo "> -jobname: $jname"
    echo "> -output dir: $umiout"

    for file in $tagout/*.GEUS.bam
    do

        baseName=$(basename $file .GEUS.bam)
	outfile=$umiout/$baseName.bcUMI.bam
	echo "> Running BC+UMI taggin of $file!"
        echo "#!/bin/bash -l
#SBATCH -p normal
#SBATCH -t 5:00:00
#SBATCH -J $jname.$baseName
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem 189840
#SBATCH -e $jname.$baseName.err
#SBATCH -o $jname.$baseName.out

$singCmd \
    java -jar -Xmx30g /opt/bin/sicelore/Jar/NanoporeBC_UMI_finder-1.0.jar -i $file \\
      -o $outfile \\
      -k $ill_pars \\
      --maxUMIfalseMatchPercent 6 \\
      --maxBCfalseMatchPercent 6 \\
      --logFile $outfile.log
 
    " > $subscr/$jname.$splt.run.sh

        sbatch $subscr/$jname.$splt.run.sh

    done

fi








### 5. Merge data from BC_UMI_tag
mergout=${outdir_sic}/5_mergedBAM

if [  $mergeBCUt == 1 ]
then

    jname="merge_BCUtag"

    echo "> 5. Merge files from BC_UMI tagging"
    echo "> -jobname: $jname"
    echo "> -output dir: $mergout"

    ## Get different basefiles
    for file in $ONTpath/*fastq
    do
        baseName=$(basename $file .fastq)
#	echo $baseName
	outfile="$mergout/$baseName.bcUMI_umifound_.merged.bam"
	scrFile="$subscr/${jname}_${baseName}_umifound_run.sh"
	bamFiles="$mergout/$baseName.bamfiles_umifound.txt"

	touch $bamFiles
	## Add INPUT= for each split bam
	for base in $umiout/*$baseName.bcUMI_umifound_.bam
	do 
	    echo "$base" >> $bamFiles
	done

	echo "#!/bin/bash -l
#SBATCH -p normal
#SBATCH -t 5:00:00
#SBATCH -J $jname.$baseName
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem 189840
#SBATCH -e $jname.$baseName.err
#SBATCH -o $jname.$baseName.out

$singCmd samtools merge -b $bamFiles $outfile  

$singCmd samtools index $outfile" > $scrFile
	
    sbatch $scrFile

    done
    
fi # mergeBUMI tag done



### 6. Generate consensus seq
consout=${outdir_sic}/6_consensus

if [ $consensus == 1 ]
then
    jname="sic_cons"
    mkdir -p $consout
    echo "> 6. Create consensus from merge bams"
    echo "> -jobname: $jname"
    echo "> -output dir: $consout"

    declare -a jobs=('1' '2' 'X' 'MT' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' 'Y');

    #Testing
#    declare -a jobs=('1' '2') 

    for file in $ONTpath/*fastq
    do
        baseName=$(basename $file .fastq)
#	echo $baseName
	mergedfile="$mergout/$baseName.bcUMI_umifound_.merged.bam"

	for i in ${jobs[@]} 
	do
    	    chroutfile="$consout/$baseName.bcUMI_umifound_.chr.${i}.bam"
	    outfile="$consout/$baseName.bcUMI_umifound_.molecules.chr.${i}.bam"
	    scrFile="$subscr/${jname}_${baseName}.chr${i}.run.sh"
	    echo "#!/bin/bash -l
#SBATCH -p normal
#SBATCH -t 5:00:00
#SBATCH -J $jname.$baseName.chr${i}
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem 189840
#SBATCH -e $jname.$baseName.chr${i}.err
#SBATCH -o $jname.$baseName.chr${i}.out

$singCmd samtools view -Sb $mergedfile ${i} > $chroutfile

$singCmd samtools index $chroutfile

$singCmd java -jar -Xmx22g /opt/bin/sicelore/Jar/Sicelore-1.0.jar ComputeConsensus I=$chroutfile O=$outfile T=20 TMPDIR=$consout/tmp/ 

" > $scrFile

	    sbatch $scrFile
       done # finish iterating chr 

    done # finish basename itr

fi # consensus finish
