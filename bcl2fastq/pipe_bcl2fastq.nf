#!/usr/bin/env nextFlow

// set variables
exp = params.experiment // param has to be set manually in nextflow.config file
metaID = params.metaid // param has to be set manually in nextflow.config file

// automatically assigned in config file
OUTDIR = params.outDir
FQDIR = params.fqDir
CNTDIR = params.countDir
QCDIR = params.qcDir
lanes = params.lanes // should be set to lanes=0 if all lanes to be included. otherwise, set "1,3" etc.

// Read and process sample sheet
sheet = file(params.sheet)

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


