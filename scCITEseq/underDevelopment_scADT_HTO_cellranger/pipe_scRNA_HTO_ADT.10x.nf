#!/usr/bin/env nextFlow

// set variables
exp = params.experiment
metaID = params.metaid
OUTDIR = params.outDir
FQDIR = params.fqDir
CNTDIR = params.countDir
QCDIR = params.qcDir
CSDIR = params.citeseqDir

// Read and process sample sheet
sheet = file(params.sheet)
adtsheet = file(params.adt_ssheet)
htosheet = file(params.hto_ssheet)


// create new file for reading into channels that provide sample info!
newsheet = file("$exp/sample_sheet.nf.htoadt.csv")

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
    .map { row -> tuple( row.Sample_ID, row.Sample_Name, row.Sample_Project ) }
    .unique()
    .tap{infoall}
    .into { cellrangerMKF; crCount_csv; infoch; }

println " ========= INFO =========== "
println ">>> scRNAseq 10x Chromium >>>"
println "> Experiment: $exp "
println "> Sample sheet mRNA: $sheet "
println "> Sample sheet ADT: $adtsheet "
println "> Sample sheet HTO: $htosheet "
println "> Project ID: $metaID "
println "> output dir: $OUTDIR "
println "> Species: $species "
println "> Reference data: $genome \n"
println " ========================== "

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


process bcl2fastq_adt {

    input:
    val adtsheet

    output:
    val "$FQDIR/ADT/${metaid}/*.fastq.gz" into fqc_adt
    val "x" into citeseqcount_adt

    when:
    params.adt
    
    """

    bcl2fastq -R $exp \\
            --sample-sheet $adtsheet \\
            --no-lane-splitting \\
            --output-dir $FQDIR/ADT/
    """

}


process bcl2fastq_hto {

    input:
    val htosheet

    output:
    val "$FQDIR/HTO/${metaid}/*.fastq.gz" into fqc_hto
    val "x" into citeseqcount_hto

    when:
    params.hto

    """
    bcl2fastq -R $exp \\
            --sample-sheet $htosheet \\
            --no-lane-splitting \\
            --output-dir $FQDIR/HTO/
    """

}


process count {

	publishDir "${CNTDIR}/", mode: "copy", overwrite: true

	input: 
	val x from crcount
        set sid, sname, projid from crCount_csv

	output:
	file "${sname}/*" into cellrangerPost
        val "${CNTDIR}/${sname}/outs/filtered_feature_bc_matrix/barcodes.tsv" into citeseq_bc


	"""
	cellranger count \\
	     --id=$sname \\
	     --fastqs=${FQDIR}/$projid/$sname/ \\
	     --sample=$sname \\
             --project=$projid \\
	     --transcriptome=$genome \\
             --localcores=56 --localmem=90 

        mkdir -p $QCDIR
        cp ${sname}/outs/web_summary.html ${QCDIR}/${sname}.web_summary.html
	"""
}


process fastqc {

	publishDir "${QCDIR}/FastQC/", mode: 'copy', overwrite: true, pattern: "*_fastqc.*"

	input:
        file x from fqc_ch.flatten()
        
        output:
        val "*fastqc*" into qc_ch

    
	"""
        mkdir -p $QCDIR
        mkdir -p $QCDIR/FastQC/

	fastqc $x 
	"""
    
}

process fastqc_adt {

    input:
    file x from fqc_adt

    output:
    val "x" into adt_multiqc

    when:
    params.adt && $x =~ /^Undetermined*/ 

    """
    mkdir -p $QCDIR
    mkdir -p $QCDIR/FastQC/

    fastqc $x
    """

}

process fastqc_hto {

    input:
    file x from fqc_hto

    output:
    val "x" into hto_multiqc

    when:
    params.hto && $x =~ /^Undetermined*/ 
    
    """
    mkdir -p $QCDIR
    mkdir -p $QCDIR/FastQC/

    fastqc $x
    """

}

process multiqc {


    publishDir "${QCDIR}/", mode: 'copy', overwrite: true

    input:
    val x from qc_ch.collect()
    val y from adt_multiqc.collect()
    val z from hto_multiqc.collect()

    output:
    file "multiqc_report.html" into multiqc_outch
    file "multiqc_data"

    script:
    """
    multiqc $OUTDIR

    """
}

/*
process citeseqCount {

    input:
    val x from citeseqcount_adt


    when:
    params.citeseqCount


    """
    source activate CITEseq

    if [ -f ${bcs}.gz ]
    then
        gunzip ${bcs}.gz
    fi

    CITE-seq-Count \
       -R1 $r1 \
       -R2 $r2 \
       -t params.citeseq_tags \
       -cbf 1 \
       -cbl 16 \
       -umif 17 \
       -umil 26 \
       --whitelist $bcs \
       -o $CSDIR

    """
}

*/
