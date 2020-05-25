#!/usr/bin/env nextFlow

// set variables
exp = params.experiment
metaID = params.metaid
OUTDIR = params.outDir
FQDIR = params.fqDir
CNTDIR = params.countDir
QCDIR = params.qcDir

sheet = file(params.sheet)

if ( params.sheetType == "simple" ) {
   
   Channel
	 .fromPath(params.sheet)
	 .splitCsv(header:false)
	 .map{ row -> val(row[1]) }
         .tap{info1}
	 .into { cellrangerMKF; cellrangerCount }
	 
   sheet = file(params.sheet)
}
else {

     newsheet = file("$exp/sample_sheet.nf.csv")

     allLines = sheet.readLines()
     
     writeB = false
     newsheet.text=""     

     for ( line in allLines ) {

	 if ( writeB ) {
	    newsheet.append(line + "\n")
	 }
	 if (line.contains("[Data]")) {
	    writeB = true
	 }
	
     }
     Channel
	.fromPath(newsheet)
	.splitCsv(header:true)
	.map{ row -> tuple(row.Sample_ID,row.Sample_Name,row.Sample_Project) }
        .tap{info1}
	.into { cellrangerMKF; crCount_csv; infoch}
     
}

println "================================="
println ">> Sample info from samplesheet: "
info1.subscribe{ "> Info: $it " }

/*
* Print Info on experiment and sample sheet
*/

process printInfo {

	input:
	val exp
	val sheet
	val metaID
	val OUTDIR

	output:
	stdout info

	"""
	printf "======= Info ==========\n"

	printf ">>> scATACseq 10x Chromium >>>\n"

	printf "> Experiment: $exp \n"
	printf "> Sample sheet: $sheet \n"
	printf "> Project ID: $metaID \n"
	printf "> output dir: $OUTDIR \n"

	printf "======================= \n"
	"""	
}

info.view {  it }


// Run mkFastq

process mkfastq {
	
	publishDir "${FQDIR}/", mode: 'copy', overwrite: true

	input:
        val sheet 

	output:
	file("$metaID/outs/**/*.fastq.gz") into (cellrangerCount, fqc_ch, multiqc_ch) mode flatten
	file("$metaID/outs/fastq_path/Stats/*") into bqcPaths

	"""
	cellranger-atac mkfastq \\
    		   --id=$metaID \\
        	   --run=$exp \\
		   --csv=$sheet \\
	           --jobmode=local \\
		   --localmem=120 \\
		   --qc
		  
	"""
}

/*
process count {
	publishDir "${CNTDIR}/", mode: "copy", overwrite: true

	input: 
	file x from cellrangerCount
	set sample_id, sample_name, proj from crCount_csv

	output:
	file "${proj}/${proj}.mri.tgz" into cellrangerPost

	
	"""
	cellranger-atac count \\
	     --id=$ID \\
	     --fastqs="${FQDIR}/${proj}/outs/fastq_path/${project}/" \\
	     --sample=$sample_name \\
	     --reference=/opt/refdata-cellranger-atac-GRCh38-1.2.0/

	"""
}

*/
process fastqc {

	publishDir "${QCDIR}/", mode: 'copy', overwrite: true

	input:
	file x from fqc_ch

        output:
        file "*fastqc.{zip,html}" into qc_ch mode flatten
    
	"""
        mkdir -p $QCDIR

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
