#!/usr/bin/env nextFlow

// set variables
exp = params.experiment
metaID = params.metaid
OUTDIR = params.outDir
FQDIR = params.fqDir
CNTDIR = params.countDir
QCDIR = params.qcDir
VDJDIR = params.vdjDir
SUMDIR = params.sumDir

// Read and process sample sheet
sheet = file(params.sheet)

// create new file for reading into channels that provide sample info!
newsheet = file("$exp/sample_sheet.nf.csv")


// read samplesheet. Get species, and sample info
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
	    vdj_genome = "/opt/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0"
	} else {
	    genome = "/opt/refdata-cellranger-mm10-3.0.0/"
	    vdj_genome = "/opt/refdata-cellranger-vdj-GRCm38-alts-ensembl-3.1.0"
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

// collect samplesheet info
Channel
    .fromPath(newsheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.Sample_ID, row.Sample_Name, row.Sample_Project, row.Library_Type ) }
    .unique()
    .tap{infoall}
    .into { crCount_csv; crVDJ_csv}

infoall.subscribe{ println "Info: $it" }

println "============================="
println ">>> scVDJseq 10x Chromium >>>"
println "> Experiment: $exp "
println "> Sample sheet: $sheet "
println "> Project ID: $metaID "
println "> Species: $species "
println "> Reference data: $genome \n"
println "> output dir: $OUTDIR "
println "> Fastq dir:   $FQDIR "
println "> Count dir:   $CNTDIR "
println "> QC dir:      $QCDIR "
println "> VDJ dir:     $VDJDIR "
println "> Summary dir: $SUMDIR "
println "============================="


// Run mkFastq
process mkfastq {
	
	input:
        val sheet 

	output:
	file("$metaID/outs/**/*.fastq.gz") into fqc_ch
	val 1 into crcount
        val 1 into crVdj


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

process crVDJ {

    publishDir "${VDJDIR}/", mode: "copy", overwrite: true
    
    input:
    val x from crVdj
    set sid, sname, projid, libtype from crVDJ_csv

    output: 
    file "${metaID}/*" into vdjPost

    when:
    libtype == "VDJ"

    """
    cellranger vdj \\
    --id=$metaID \\
    --fastqs=${FQDIR} \\
    --reference=$vdj_genome \\
    --sample=${sname} \\
    --localcores=24 \\
    --localmem=64 

    mkdir -p $SUMDIR
    cp ${metaID}/outs/web_summary.html ${SUMDIR}/${sname}.VDJ_web_summary.html
    cp ${metaID}/outs/vloupe.vloupe ${SUMDIR}/${sname}.VDJ_vloupe.vloupe


    """
}

process count {

	publishDir "${CNTDIR}/", mode: "copy", overwrite: true

	input: 
	val x from crcount
        set sid, sname, projid, libtype from crCount_csv

	output:
	file "${sname}/*" into cellrangerPost

        when:
        libtype == "mRNA"

	"""
	cellranger count \\
	     --id=$sname \\
	     --fastqs=${FQDIR} \\
	     --sample=$sname \\
             --project=$projid \\
	     --transcriptome=$genome \\
             --localcores=56 --localmem=90 

        mkdir -p $SUMDIR
        cp ${sname}/outs/web_summary.html ${SUMDIR}/${sname}.mRNA_web_summary.html
        cp ${sname}/outs/cloupe.cloupe ${SUMDIR}/${sname}.mRNA_cloupe.cloupe
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
