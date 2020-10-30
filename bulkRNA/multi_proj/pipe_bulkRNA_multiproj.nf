#!/usr/bin/env nextFlow

// set variables
exp = params.experiment // param has to be set manually in nextflow.config file
launchdir = params.workdir
metaID = params.metaid // param has to be set manually in nextflow.config file

// automatically assigned in config file
OUTDIR = params.outDir // top outdir (output from individual projects will be added in separate folders herein)
FQDIR = params.fqDir

lanes = params.lanes // should be set to lanes=0 if all lanes to be included. otherwise, set "1,3" etc.

// Read and process sample sheet
sheet = file(params.sheet)

// create new file for reading into channels that provide sample info!
newsheet = file("${launchdir}/sample_sheet.nf.csv")

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
println ">>> Bulk RNA pipeline - multiple projects >>> "
println "> Experiment 	       : $exp "
println "> Sample sheet	       : $sheet "
println "> Experiment ID       : $metaID "
println "> Output dir 	       : $OUTDIR "
println "> Common Fastq dir    : $FQDIR "

println "======================= "

// Projects
Channel
    .fromPath(newsheet)
    .splitCsv(header:true)
    .map { row -> row.Sample_Project }
    .unique()
    .tap{infoProject}
    .into{mvFastq_ch; quant_proj}

infoProject.subscribe{ println "Projects: $it" }

// sample info
Channel
    .fromPath(newsheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.Sample_ID, row.Sample_Name, row.Sample_Project ) }
    .unique()
    .tap{infoSamples}
    .into{star_ch; fcount_ch}

infoSamples.subscribe{ println "Samples: $it" }

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
              -r $task.cpus \\
              -p $task.cpus  \\
              -w $task.cpus  \\
              --output-dir $FQDIR
   
     """
}

process moveFastq {

    input:
    val x from moveFastq
    val projid from mvFastq_ch

    output:
    val projid into fqc_ch, star_go

    """
    mkdir -p ${OUTDIR}/${projid}
    mkdir -p ${OUTDIR}/${projid}/Fastq_Raw

    mv ${FQDIR}/${projid}/* ${OUTDIR}/${projid}/Fastq_Raw/

    """

}



// fastqc 
process fastqc {

	input:
	val projid from fqc_ch

        output:
        val projid into multiqc_fastqc
	
	"""

        mkdir -p ${OUTDIR}/$projid/QC/
        mkdir -p ${OUTDIR}/$projid/QC/FastQC

    	cd $OUTDIR/$projid/QC
        for file in $OUTDIR/$projid/Fastq_Raw/**/*fastq.gz
             do fastqc -t $task.cpus \$file
  	done

        mv ${OUTDIR}/${projid}/Fastq_Raw/**/*fastqc* $OUTDIR/$projid/QC/FastQC 

           
	"""
    
}


// Run STAR
process STAR  {

    publishDir "${OUTDIR}/$projid/STAR/", mode: 'copy', overwrite: true

    input:
    val projDone from star_go.collect()
    set sid, sname, projid, species from star_ch

    output:
    set val(sname), file("${sname}*") into postStar
    file "${sname}Aligned.sortedByCoord.out.bam" into count

    when:
    params.align
   
    
    """

    read1=\$(echo ${OUTDIR}/$projid/Fastq_Raw/$sid/${sname}*R1*fastq.gz)
    read2=\$(echo ${OUTDIR}/$projid/Fastq_Raw/$sid/${sname}*R2*fastq.gz)

    # Check if fastq is not containing wildcards (due to sample fastq are not put in sample id folder
    if [[ \$read1 == *"*R1*"* ]]; then
       read1=\$(echo ${OUTDIR}/${projid}/Fastq_Raw/${sname}*R1*fastq.gz)
       read2=\$(echo ${OUTDIR}/${projid}/Fastq_Raw/${sname}*R2*fastq.gz)
    fi

    # Species 
    if [ $species == "Human" ] || [ $species == "human" ]; then
       genome="/opt/star_ref_index/"

    elif [ $species == "Murine"  ] || [ $species == "murine" ] || [ $species == "mouse" ] || [ $species == "Mous\
e" ]; then
	# Mouse reference not yet in container
        genome=params.ref_genome_dir_mouse
    else
	echo ">SPECIES NOT RECOGNIZED - using nextflow.config params!"
        genome=params.ref_genome_dir
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

	when:
	params.align
	params.quant
	
	"""
 	mkdir -p  ${OUTDIR}/$projid/Quantification/

        # Get bam files for project
        cd ${OUTDIR}/$projid/STAR 
	bams=\$(echo *.bam)
	echo "BAMS: $bams"

        outfile="${OUTDIR}/$projid/Quantification/${projid}_genename.featureCounts.txt"

	# Species 
	if [ $species == "Human" ] || [ $species == "human" ]; then
           genome="/opt/gencode.v33.annotation.gtf.gz"

	elif [ $species == "Murine"  ] || [ $species == "murine" ] || [ $species == "mouse" ] || [ $species == "Mouse" ]; then
	   # Mouse reference not yet in container
           genome=params.ref_genome_dir_mouse
	else 
      	   echo ">SPECIES NOT RECOGNIZED - using nextflow.config params!"
	   genome=params.ref_genome_dir
      	fi
		
        
	featureCounts -T ${task.cpus} -t ${params.feature} -a ${params.gtf} -g gene_name -o \${outfile} -p -s ${params.stranded} \${bams}

	
   	
	"""

}

process multiqc {



    input:
    val projid from postCount

    output:
    val "multiqc_report.html" into multiqc_outch
    
    script:
    """
    cd $OUTDIR/$projid
    multiqc -n $projid -o ${OUTDIR}/$projid/QC/ .

    
    """
}



