#!/usr/bin/env nextFlow

// set variables
exp = params.experiment
metaID = params.metaid
OUTDIR = params.outDir
FQDIR = params.fqDir
CNTDIR = params.countDir
QCDIR = params.qcDir
lanes = params.lanes // should be set to lanes=0 if all lanes to be included. otherwise, set "1,3" etc.

// Read and process sample sheet
sheet = file(params.sheet)

if ( params.sheetType == "simple" ) {
   
   Channel
	 .fromPath(params.sheet)
	 .splitCsv(header:false)
	 .map{ row -> tuple(row[0], row[1], row[2]) } // 0: samp.id; 1: samp.name; 2.samp.project
         .tap{info1}
	 .into { cellrangerMKF; crCount_csv }

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
// all samplesheet info
     Channel
	.fromPath(newsheet)
	.splitCsv(header:true)
        .map { row -> tuple( row.Sample_ID, row.Sample_Name, row.Sample_Project) }
        .unique()
        .tap{infoall}
	.into { cellrangerMKF; crCount_csv; infoch}
}

infoall.subscribe{ println "Info: $it" }
     
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

	printf ">>> Novaseq custom library BCL2FASTQ >>>\n"

	printf "> Experiment: $exp \n"
	printf "> Sample sheet: $sheet \n"
	printf "> Project ID: $metaID \n"
	printf "> output dir: $OUTDIR \n"

	printf "======================= \n"
	"""	
}

info.view {  it }

// Run bcl2fastq
process bcl2fastq {
	
	publishDir "${FQDIR}/", mode: 'copy', overwrite: true

	input:
        val sheet 

	output:
	file("$metaID/outs/**/*.fastq.gz") into (cellrangerCount, fqc_ch, multiqc_ch) mode flatten
	file("$metaID/outs/fastq_path/Stats/*") into bqcPaths

	"""
        bcl2fastq -R $exp \\
                  -o $OUTDIR \\
		  --sample-sheet $sheet \\
	         
        
                		  
	"""
}

/*
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
 //   file y from cellrangerPost

    output:
    file "multiqc_report.html" into multiqc_outch
    file "multiqc_data"

    script:
    """
    multiqc $OUTDIR

    """
}

*/
