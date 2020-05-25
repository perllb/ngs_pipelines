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

println "======= Info =========="
println ">>> Demux + QC >>> "
println "> Experiment: $exp "
println "> Sample sheet: $sheet "
println "> Project ID: $metaID "
println "> output dir: $OUTDIR "
println "======================= "


// Run bcl2fastq
process bcl2fastq {
	
	input:
        val sheet 

	output:
        val "1" into fqc_ch

	"""
        bcl2fastq -R $exp \\
                  --sample-sheet $sheet \\
                  --no-lane-splitting \\
                  --output-dir $OUTDIR        
   
                		  
	"""
}

process fastqc {

	publishDir "${QCDIR}/FastQC/", mode: 'copy', overwrite: true

	input:
	file x from fqc_ch

        output:
        file "*fastqc.{zip,html}" into qc_ch 
    
	"""
        mkdir -p $QCDIR
        mkdir -p $QCDIR/FastQC

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


