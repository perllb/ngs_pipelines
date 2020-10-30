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

// create new file for reading into channels that provide sample info!
newsheet = file("$exp/sample_sheet.nf.csv")

allLines = sheet.readLines()
writeB = false // if next lines has sample info
readS = false // if next line contain species info
newsheet.text=""
for ( line in allLines ) {

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
    .map { row -> tuple( row.Sample_ID, row.Sample_Name, row.Sample_Project ) }
    .unique()
    .tap{infoall}
    .into { fastQC_csv; multiQC_ch; }

// Projects
Channel
    .fromPath(newsheet)
    .splitCsv(header:true)
    .map { row -> row.Sample_Project }
    .unique()
    .tap{infoProject}
    .into { mvFastq_ch }


println "======= Info =========="
println ">>> Demux + QC >>> "
println "> Experiment: $exp "
println "> Sample sheet: $sheet "
println "> EXP ID: $metaID "
println "> output dir: $OUTDIR "
println "> Fastq dir: $FQDIR "
println "> QC dir: $QCDIR "
println "======================= "


//infoall.subscribe{ println "Info: $it" }

println " > Projects to process : "
infoProject.subscribe{ println "- Projects   : $it" }


println " > Starting pipeline ...."


// Run bcl2fastq
process bcl2fastq {

    input:
    val sheet

    output:
    val "x" into moveFastq

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

fastqs = Channel.fromPath( "${FQDIR}/**/*.fastq.gz" )

process moveFastq {

    input:
    val x from moveFastq
    val projid from mvFastq_ch

    output:
    val projid into fqc_ch

    """
    mkdir -p ${OUTDIR}/${projid}
    mkdir -p ${OUTDIR}/${projid}/Fastq_Raw

    mv ${FQDIR}/${projid} ${OUTDIR}/${projid}/Fastq_Raw

    """

}

// fastqc
process fastqc {

    input:
    val projid from fqc_ch

    output:
    val projid into multiqcProj


    """
    mkdir -p ${OUTDIR}/$projid/QC/

    cd $OUTDIR/$projid/QC
    for file in $OUTDIR/$projid/Fastq_Raw/*/*fastq.gz
         do fastqc \$file
    done
    """

}


process multiqc {


    publishDir "${OUTDIR}/$projid/QC/", mode: 'copy', overwrite: true

    input:
    val projid from multiqcProj

    output:
    file "multiqc_report.html" into multiqc_outch
    file "multiqc_data"

    script:
    """

    cd $OUTDIR/$projid

    multiqc -n $projid .

    """
}
