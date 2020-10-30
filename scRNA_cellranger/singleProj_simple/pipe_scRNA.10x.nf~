#!/usr/bin/env nextFlow

// set variables
exp = params.experiment
metaID = params.metaid
OUTDIR = params.outDir
FQDIR = params.fqDir
CNTDIR = params.countDir
QCDIR = params.qcDir
SUMDIR = params.sumDir

// Read and process sample sheet
sheet = file(params.sheet)

// create new file for reading into channels that provide sample info!
newsheet = file("$exp/sample_sheet.nf.csv")

// Process samplesheet data
allLines = sheet.readLines()
writeB = false // if next lines has sample info
readS = false // if next line contain species info
newsheet.text=""     

for ( line in allLines ) {

    // Define species of samples
    if ( readS ) {
	species = line
	if ( species.contains("Human") || species.contains("human" )) {
	    genome = "/opt/refdata-cellranger-GRCh38-3.0.0/" 
	} else {
	    genome = "/opt/refdata-cellranger-mm10-3.0.0/"
	}
	readS = false
    }
    if ( writeB ) {
	newsheet.append(line + "\n")
    }
    if (line.contains("[Data]")) {
	writeB = true
    }
    if (line.contains("[Species]")) {
	readS = true
    }
}

// all samplesheet info
Channel
    .fromPath(newsheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.Sample_ID, row.Sample_Name, row.Sample_Project) }
    .unique()
    .tap{infoall}
    .into { cellrangerMKF; crCount_csv; infoch}

println " ========= INFO =========== "
println ">>> scRNAseq 10x Chromium >>>"
println "> Experiment: $exp "
println "> Sample sheet: $sheet "
println "> Project ID: $metaID "
println "> Species: $species "
println "> Reference data: $genome \n"
println "> output dir:  $OUTDIR "
println "> Fastq dir:   $FQDIR "
println "> Count dir:   $CNTDIR "
println "> QC dir:      $QCDIR "
println "> VDJ dir:     $VDJDIR "
println "> Summary dir: $SUMDIR "
println "============================="



infoall.subscribe{ println "Info: $it" }

// Run mkFastq
process mkfastq {
	
	input:
        val sheet 

	output:
	file("$metaID/outs/**/*.fastq.gz") into fqc_ch
	val 1 into crcount

	"""
	cellranger mkfastq \\
    		   --id=$metaID \\
        	   --run=$exp \\
		   --csv=$sheet \\
	           --jobmode=local \\
		   --localmem=120 \\
		   --qc \\
                   --output-dir $FQDIR 
		  
	"""
}


process count {

	publishDir "${CNTDIR}/", mode: "copy", overwrite: true

	input: 
	val x from crcount
        set sid, sname, projid from crCount_csv

	output:
	file "${sname}/*" into cellrangerPost

	"""
	cellranger count \\
	     --id=$sname \\
	     --fastqs=${FQDIR}/$projid/$sname/ \\
	     --sample=$sname \\
             --project=$projid \\
	     --transcriptome=$genome \\
             --localcores=56 --localmem=90 

        mkdir -p $QCDIR
        mkdir -p $SUMDIR

        cp ${sname}/outs/web_summary.html ${SUMDIR}/${sname}.web_summary.html
        cp ${sname}/outs/cloupe.cloupe ${SUMDIR}/${sname}.cloupe
	"""
}


process fastqc {

	publishDir "${QCDIR}/FastQC/", mode: 'copy', overwrite: true, pattern: "*_fastqc.*"

	input:
	file x from fqc_ch.flatten()

        output:
        file "*fastqc*" into qc_ch
    
	"""
        mkdir -p $QCDIR
        mkdir -p $QCDIR/FastQC/

	fastqc $x 
	"""
    
}

process multiqc {


    publishDir "${QCDIR}/", mode: 'copy', overwrite: true

    input:
    file x from qc_ch.collect()

    output:
    file "multiqc_report.html" into multiqc_outch
    file "multiqc_data"

    script:
    """
    multiqc $OUTDIR

    """
}
