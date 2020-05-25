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

// create new file for reading into channels that provide sample info!
newsheet = file("$exp/sample_sheet.nf.csv")

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
    .map { row -> tuple( row.Sample_ID, row.Sample_Name, row.Sample_Project, row.Sample_Species) }
    .unique()
    .tap{infoall}
    .into { cellrangerMKF; crCount_csv; infoch; merge_csv }

// Projects
Channel
    .fromPath(newsheet)
    .splitCsv(header:true)
    .map { row -> row.Sample_Project }
    .unique()
    .tap{infoProject}
    .into { mvFastq_csv }




println " ========== INFO ============= "
println ">>> scRNAseq 10x Chromium >>>"
println "> Experiment: $exp "
println "> Sample sheet: $sheet "
println "> Project ID: $metaID "
println "> output dir: $OUTDIR "
println "> Species: $species "
println "> Reference data: $genome \n"
println " ========================== "


infoall.subscribe{ println "Info: $it" }

println " > Projects to process : "
infoProject.subscribe{ println "Info Projects: $it" }
println " > Starting pipeline ...." 

// Run mkFastq
process mkfastq {
	
	input:
        val sheet 

	output:
	val 1 into moveFastq

	"""
	cellranger mkfastq \\
    		   --id=$metaID \\
        	   --run=$exp \\
		   --csv=$sheet \\
	           --jobmode=local \\
		   --localmem=120 \\
		   --qc \\
                   --barcode-mismatches=0 \\
                   --output-dir $FQDIR 
		  
	"""
}

process moveFastq {

    input:
    val x from moveFastq
    val projid from mvFastq_csv

    output:
    val "1" into crCount
    val projid into fqc_ch

    """
    mkdir -p ${OUTDIR}/${projid}
    mkdir -p ${OUTDIR}/${projid}/Fastq_Raw
    mv ${FQDIR}/${projid} ${OUTDIR}/${projid}/Fastq_Raw/
    """

}

process count {

	publishDir "${OUTDIR}/${projid}/Count_CR/", mode: "copy", overwrite: true

	input: 
	val x from crCount
        set sid, sname, projid, species from crCount_csv

	output:
        file "${sname}" into mergecount

	"""
        if [ $species == "Human" ] || [ $species == "human" ]
        then
            genome="/opt/refdata-cellranger-GRCh38-3.0.0/"
        elif [ $species == "Murine"  ] || [ $species == "murine" ] || [ $species == "mouse" ] || [ $species == "Mouse" ]
        then
            genome="/opt/refdata-cellranger-mm10-3.0.0/"
        else
            echo ">SPECIES NOT RECOGNIZED!"
            genome="ERR"
        fi

        mkdir -p ${OUTDIR}/${projid}/Count_CR/

	cellranger count \\
	     --id=$sname \\
	     --fastqs=${OUTDIR}/$projid/Fastq_Raw/$sname \\
	     --sample=$sname \\
             --project=$projid \\
	     --transcriptome=\$genome \\
             --localcores=56 --localmem=90 


        mkdir -p "${OUTDIR}/${projid}/QC"
        cp ${sname}/outs/web_summary.html "${OUTDIR}/${projid}/QC/${sname}.web_summary.html
	"""
}

process fastqc {

	publishDir "${OUTDIR}/$projid/QC/", mode: 'copy', overwrite: true, pattern: "*_fastqc.*"

	input:
	val projid from fqc_ch

        output:
        val "1" into mergeqc
        val "2" into multiqc
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
    val x from multiqc.collect()
    val projid from multiqcProj

    output:
    file "multiqc_report.html" into multiqc_outch
    file "multiqc_data"

    script:
    """
    cd $OUTDIR/$projid
    multiqc .

    """
}



