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
sheet_RNA_in = file(params.sheet)
sheet_ADT_in = file(params.citeseq_ssheet)

// create new file for reading into channels that provide sample info!
RNAsheet = file("$exp/sample_sheet.nf.csv")
ADTsheet = file("$exp/sample_sheet.nf_adt.csv")

// Process samplesheet data
//
// RNA
//
allLines = sheet_RNA_in.readLines()
writeB = false // if next lines has sample info
readS = false // if next line contain species info
RNAsheet.text=""     

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
	RNAsheet.append(line + "\n")
    }
    if (line.contains("[Data]")) {
	writeB = true
    }
    if (line.contains("[Species]")) {
	readS = true
    }
}

// Process samplesheet data
//
// ADT
//
allLines = sheet_ADT_in.readLines()
writeB = false // if next lines has sample info
ADTsheet.text=""     

for ( line in allLines ) {

    if ( writeB ) {
	ADTsheet.append(line + "\n")
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
    .fromPath(RNAsheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.Sample_ID, row.Sample_Name, row.Sample_Project ) }
    .unique()
    .tap{infoall}
    .into { cellrangerMKF; crCount_csv; infoch; mRNA_sample_CSC }


// ADT samplesheet 
Channel
    .fromPath(ADTsheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.Sample_ID, row.Sample_Name, row.Sample_Project ) }
    .unique()
    .tap{infoADT}
    .into { ADT_sample_CSC }

println " ========= INFO =========== "
println ">>> scRNAseq 10x Chromium >>>"
println "> Experiment: $exp "
println "> Sample sheet mRNA: $sheet_RNA_in "
println "> Sample sheet ADT: $sheet_ADT_in "
println "> Project ID: $metaID "
println "> output dir: $OUTDIR "
println "> Species: $species "
println "> Reference data: $genome \n"
println " ========================== "


infoall.subscribe{ println "Info RNA Samples\n: $it" }

infoADT.subscribe{ println "Info ADT Samples\n: $it" }

// Run mkFastq
process mkfastq {
	
	input:
        val sheet_RNA_in

	output:
	file("$metaID/outs/**/*.fastq.gz") into fqc_ch
	val 1 into crcount

	"""
	cellranger mkfastq \\
    		   --id=$metaID \\
        	   --run=$exp \\
		   --csv=$sheet_RNA_in \\
	           --jobmode=local \\
		   --localmem=120 \\
		   --qc \\
                   --output-dir $FQDIR 
		  
	"""
}

process bcl2fastq_adt {

    input:
    val sheet_ADT_in

    output:
    val "$FQDIR/ADT/*/*.fastq.gz" into fqc_adt
    val "x" into citeseqcount_adt

    """
    mkdir -p $FQDIR
    mkdir -p $FQDIR/ADT

    bcl2fastq -R $exp \\
            --sample-sheet $sheet_ADT_in \\
            --no-lane-splitting \\
            --output-dir $FQDIR/ADT/

    rm -r -I $FQDIR/ADT/Undetermined*

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
        file "*fastqc*" into qc_ch
    
	"""
        mkdir -p $QCDIR
        mkdir -p $QCDIR/FastQC/

	fastqc $x 
	"""
    
}

process fastqc_adt {

    publishDir "${QCDIR}/FastQC/", mode: 'copy', overwrite: true, pattern: "*fastqc.*"

    input:
    val x from fqc_adt

    output:
    val "*fastqc*" into adt_multiqc

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
    file y from adt_multiqc.collect()

    output:
    file "multiqc_report.html" into multiqc_outch
    file "multiqc_data"

    script:
    """
    multiqc $OUTDIR

    """
}


process citeseqCount {

    input:
    val x from citeseqcount_adt.collect()
    set sid, sname, projid from ADT_sample_CSC
    set msid, msname, mprojid from mRNA_sample_CSC    
    val y from citeseq_bc

    when:
    params.citeseqCount

    bcs="${CNTDIR}/${msname}/outs/filtered_feature_bc_matrix/barcodes.tsv"
    citeseqtags=params.citeseq_tags

    shell:
    """
    source activate CITEseq

    echo $bcs
    
    if [ -f ${bcs}.gz ]
    then
        zcat ${bcs}.gz > $bcs
    fi

    r1=\$(echo ${FQDIR}/ADT/$projid/${sname}*R1*fastq.gz)
    r2=\$(echo ${FQDIR}/ADT/$projid/${sname}*R2*fastq.gz)

    CITE-seq-Count \
       -R1 \$r1 \
       -R2 \$r2 \
       -t $citeseqtags \
       -cbf 1 \
       -cbl 16 \
       -umif 17 \
       -umil 28 \
       --whitelist $bcs \
       -cells \$(cat $bcs | wc -l) \
       -T 56 \
       -o $CSDIR/${sname}_${msname}/

    """
}

