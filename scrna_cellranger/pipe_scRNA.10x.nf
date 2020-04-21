#!/usr/bin/env nextFlow

// set variables
exp = params.experiment
metaID = params.metaid
OUTDIR = params.outDir
FQDIR = params.fqDir
CNTDIR = params.countDir
QCDIR = params.qcDir

// Read and process sample sheet
sheet = file(params.sheet)

// The "simple" type should be on the format: <sample id>, <sample name>, <sample project>
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
     // create new file for reading into channels that provide sample info!
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

	printf ">>> scRNAseq 10x Chromium >>>\n"

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
	
	input:
        val sheet 

	output:
	file("**.fastq.gz") into fqc_ch
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
	     --fastqs=${FQDIR} \\
	     --sample=$sname \\
	     --transcriptome=/opt/refdata-cellranger-GRCh38-3.0.0/ \\
             --localcores=56 --localmem=90

        mkdir -p $QCDIR
        cp ${sname}/outs/web_summary.html ${QCDIR}/${sname}.web_summary.html
	"""
}


process fastqc {

	publishDir "${QCDIR}/", mode: 'copy', overwrite: true, pattern: "*_fastqc.*"

	input:
	file x from fqc_ch 

        output:
        val "1" into qc_ch
    
	"""
        mkdir -p $QCDIR

	fastqc $x 
	"""
    
}

process multiqc {


    publishDir "${QCDIR}/", mode: 'copy', overwrite: true

    input:
    val x from qc_ch

    output:
    file "multiqc_report.html" into multiqc_outch
    file "multiqc_data"

    script:
    """
    multiqc $OUTDIR

    """
}
