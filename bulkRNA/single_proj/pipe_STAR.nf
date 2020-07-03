#!/usr/bin/env nextFlow

// set variables
exp = params.experiment // param has to be set manually in nextflow.config file
metaID = params.metaid // param has to be set manually in nextflow.config file

// automatically assigned in config file
OUTDIR = params.outDir
FQDIR = params.fqDir
BAMDIR = params.bamDir
QUANTDIR = params.quantDir

lanes = params.lanes // should be set to lanes=0 if all lanes to be included. otherwise, set "1,3" etc.

// Read and process sample sheet
sheet = file(params.sheet)


// create new file for reading into channels that provide sample info!
newsheet = file("$exp/sample_sheet.nf.csv")

allLines = sheet.readLines()
writeB = false // if next lines has sample info
readS = false // if next line contain species info
newsheet.text=""

for ( line in allLines ) {

    // read if [Data]
    if ( writeB ) {
        newsheet.append(line + "\n")
    }
    if (line.contains("[Data]")) {
        writeB = true
    }	
}

println "======= Info =========="
println ">>> Demux + QC >>> "

println "> Experiment: $exp "
println "> Sample sheet: $sheet "
println "> Project ID: $metaID "

println "> output dir: $OUTDIR "
println "======================= "

// Samples
Channel
    .fromPath(newsheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.Sample_ID, row.Sample_Name, row.Sample_Project) }
    .unique()
    .tap{infoSamples}
    .set { star_ch }

// Projects
Channel
    .fromPath(newsheet)
    .splitCsv(header:true)
    .map { row -> row.Sample_Project }
    .unique()
    .tap{infoProject}
    .set { quant_proj }

infoProject.subscribe{ println "Info: $it" }

// Run STAR
process STAR  {

    publishDir "${OUTDIR}/${projid}/STAR/", mode: 'copy', overwrite: true

    cores=56

    input:
    set sid, sname, projid from star_ch

    output:
    set val(sname), file("${sname}*") into qualimap_map, quant_ch


    """

    # Get fastq of reads 
    read1=\$(echo ${OUTDIR}/${projid}/Fastq/$sid/${sname}*R1*fastq.gz)
    read2=\$(echo ${OUTDIR}/${projid}/Fastq/$sid/${sname}*R2*fastq.gz)

    echo ">READ1: \$read1"

    # Check if fastq is not containing wildcards (due to sample fastq are not put in sample id folder
    if [[ \$read1 == *"*R1*"* ]]; then
       read1=\$(echo ${OUTDIR}/${projid}/Fastq/${sname}*R1*fastq.gz)
       read2=\$(echo ${OUTDIR}/${projid}/Fastq/${sname}*R2*fastq.gz)
    fi

    echo ">READ1: \$read1"

    STAR --genomeDir ${params.ref_genome_dir} \\
    --readFilesIn \${read1} \${read2} \\
    --runThreadN ${task.cpus}  \\
    --outSAMtype BAM SortedByCoordinate \\
    --readFilesCommand zcat \\
    --limitBAMsortRAM 10000000000 \\
    --outFileNamePrefix ${sname}


    """

}


process FeatureCounts {


	input:
	val projid from quant_proj
	val sname from quant_ch.collect()

	output:
	val projid into postCount
	
	"""
	 
     	mkdir -p  ${OUTDIR}/${projid}/Quantification/
	
	# Get bam files for project
	bams=\$(echo ${OUTDIR}/${projid}/STAR/*.bam)
	
	outfile=${OUTDIR}/${projid}/Quantification/${projid}_genename.featureCounts.txt

     	echo "featureCounts -T ${task.cpus} -a ${params.gtf} -g gene_name -o \${outfile} -p -s ${params.stranded} \${bams}"

 	featureCounts -T ${task.cpus} -a ${params.gtf} -g gene_name -o \${outfile} -p -s ${params.stranded} \${bams}
		
	"""	

}


process editFeatureCount {


	input:
	val projid from postCount

	output:
	val "1" into edit_ch


	"""

	outfile=${OUTDIR}/${projid}/Quantification/${projid}_genename.featureCounts.txt
	
	bpath=${OUTDIR}/$projid/STAR/

	sedded=\$(echo \$bpath | sed "s/\\//\\\\//g')

	sed 's/'"$sedded"'//g' \$outfile | sed 's/Aligned.sortedByCoord.out.bam//g' > tmp.txt
	mv tmp.txt \$outfile

	sed 's/'"$sedded"'//g' \${outfile}.summary | sed 's/Aligned.sortedByCoord.out.bam//g' > tmp.txt
	mv tmp.txt \${outfile}.summary

	"""
	
}