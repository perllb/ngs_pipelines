#!/usr/bin/env nextFlow

// set variables
exp = params.experiment // param has to be set manually in nextflow.config file
metaID = params.metaid // param has to be set manually in nextflow.config file

// automatically assigned in config file
OUTDIR = params.outDir
FQDIR = params.fqDir
QCDIR = params.qcDir
BAMDIR = params.bamDir
QDIR = params.quantDir

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


println "======= Info =========="
println ">>> Demux + QC >>> "
println "> Experiment: $exp "
println "> Sample sheet: $sheet "
println "> Project ID: $metaID "
println "> output dir: $OUTDIR "
println "> Fastq dir:  $FQDIR "
println "> QC dir:     $QCDIR "
println "> BAM dir:    $BAMDIR "
println "> Quant dir:  $QDIR "
println "======================= "

// Projects
Channel
    .fromPath(newsheet)
    .splitCsv(header:true)
    .map { row -> row.Sample_Project }
    .unique()
    .tap{infoProject}
    .set{quant_proj}

infoProject.subscribe{ println "Projects: $it" }

// sample info
Channel
    .fromPath(newsheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.Sample_ID, row.Sample_Name, row.Sample_Project) }
    .unique()
    .tap{infoSamples}
    .set{star_ch }

infoSamples.subscribe{ println "Samples: $it" }

// Run bcl2fastq
process bcl2fastq {

    input:
    val sheet 

    output:
    val "$FQDIR/*/*/*fastq.gz" into fqc_ch
        
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

// fastqc 
process fastqc {

	input:
	val fq from fqc_ch

        output:
        val "${QCDIR}/FastQC/*_fastqc*" into starGO

	"""
        fastqc $fq

        mkdir -p ${QCDIR}/
        mkdir -p ${QCDIR}/FastQC

        mv ${FQDIR}/*/*/*_fastqc* ${QCDIR}/FastQC/

       
	"""
    
}


// Run STAR
process STAR  {

    publishDir "${BAMDIR}/", mode: 'copy', overwrite: true

    input:
    file a from starGO.collect()
    set sid, sname, projid from star_ch

    output:
    set val(sname), file("${sname}*") into postStar
    file "${sname}Aligned.sortedByCoord.out.bam" into count

    
    """

    read1=\$(echo ${FQDIR}/$projid/$sid/${sname}*R1*fastq.gz)
    read2=\$(echo ${FQDIR}/$projid/$sid/${sname}*R2*fastq.gz)

    # Check if fastq is not containing wildcards (due to sample fastq are not put in sample id folder
    if [[ \$read1 == *"*R1*"* ]]; then
       read1=\$(echo ${FQDIR}/${projid}/${sname}*R1*fastq.gz)
       read2=\$(echo ${FQDIR}/${projid}/${sname}*R2*fastq.gz)
    fi

    STAR --genomeDir ${params.ref_genome_dir} \\
    --readFilesIn \${read1} \${read2} \\
    --runThreadN ${task.cpus}  \\
    --outSAMtype BAM SortedByCoordinate \\
    --readFilesCommand zcat \\
    --limitBAMsortRAM 10000000000 \\
    --outFileNamePrefix ${sname}


    """ 

}



process featureCounts {


	input:
	val projid from quant_proj
	file bams from count.collect()

	output:
	val projid into postCount
	
	"""
	 
     	mkdir -p  ${QDIR}
	
        
     	featureCounts -T ${task.cpus} -t ${params.feature} -a ${params.gtf} -g gene_name -o ${QDIR}/${projid}_genename.featureCounts.txt ${params.paired} -s ${params.stranded} ${bams}

        
        python /fs1/per/github/ngs_pipelines/bulkRNA/bin/fCounts2fpkm.py -r -i ${QDIR}/${projid}_genename.featureCounts.txt

	"""	

}

process multiqc {



    input:
    val projid from postCount.collect()

    output:
    val "multiqc_report.html" into multiqc_outch
    
    script:
    """
    cd $OUTDIR
    multiqc .

    mv multiqc* ${QCDIR}
    """
}



