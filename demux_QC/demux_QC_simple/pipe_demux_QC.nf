#!/usr/bin/env nextFlow

// set variables
exp = params.experiment // param has to be set manually in nextflow.config file
metaID = params.metaid // param has to be set manually in nextflow.config file

// automatically assigned in config file
OUTDIR = params.outDir
FQDIR = params.fqDir
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
    val "$FQDIR/*fastq.gz" into fqc_ch
        
    """
    bcl2fastq -R $exp \\
              --sample-sheet $sheet \\
              --no-lane-splitting  \\
              -r 16  \\
              -p 16  \\
              -w 16  \\
              --output-dir $FQDIR
   
     """
}

fastqs = Channel.fromPath( "${FQDIR}/*.fastq.gz" )

// fastqc 
process fastqc {

	input:
	val bclfq from fqc_ch.collect()
        file fastq from fastqs.flatten()

        output:
        val "${QCDIR}/FastQC/*_fastqc*" into qc_ch

	"""
        fastqc ${FQDIR}/$fastq

        mkdir -p ${QCDIR}/
        mkdir -p ${QCDIR}/FastQC

        mv ${FQDIR}/*_fastqc* ${QCDIR}/FastQC/

	"""
    
}


process multiqc {


    publishDir "${QCDIR}/", mode: 'copy', overwrite: true

    input:
    val x from qc_ch.collect()

    output:
    file "multiqc_report.html" into multiqc_outch
    file "multiqc_data"

    script:
    """
    
    multiqc $OUTDIR

    """
}


