#!/usr/bin/env nextFlow

// set variables
exp = params.experiment
metaID = params.metaid
OUTDIR = params.outDir
FQDIR = params.fqDir
CNTDIR = params.countDir
QCDIR = params.qcDir
AGGDIR = params.aggDir
SUMDIR = params.sumDir


// Read and process sample sheet
sheet = file(params.sheet)

// create new file for reading into channels that provide sample info!
newsheet = file("$exp/sample_sheet.nf.csv")

// Read and process sample sheet
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
    .tap{infoall}
    .into { crCount_csv; crAgg_ch }

// Projects
Channel
    .fromPath(newsheet)
    .splitCsv(header:true)
    .map { row -> row.Sample_Project }
    .unique()
    .tap{infoProject}
    .into { mvFastq_csv ; agg_ch }

println "============================="
println ">>> scRNAseq 10x Chromium >>>"
println "> Experiment: $exp "
println "> Sample sheet: $sheet "
println "> Project ID: $metaID "
println "> Species: $species "
println "> Reference data: $genome \n"
println "> output dir:  $OUTDIR "
println "> Fastq dir:   $FQDIR "
println "> Count dir:   $CNTDIR "
println "> QC dir:      $QCDIR "
println "> Aggregate dir: $AGGDIR "
println "> Summary dir: $SUMDIR "
println "============================="


println " > Samples to process: "
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
    val "y" into crCount
    val projid into fqc_ch

    """
    mkdir -p ${OUTDIR}/${projid}
    mkdir -p ${OUTDIR}/${projid}/Fastq_Raw
    mv ${FQDIR}/${projid}/* ${OUTDIR}/${projid}/Fastq_Raw/
    
    """

}

process count {

	publishDir "${OUTDIR}/${projid}/Count_CR/", mode: "copy", overwrite: true

	input: 
	val y from crCount.collect()
        set sid, sname, projid, species from crCount_csv

	output:
        file "${sname}" into samplename
        val "x" into count_agg

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

        mkdir -p ${OUTDIR}/${projid}/Summaries
        cp ${sname}/outs/web_summary.html ${OUTDIR}/${projid}/Summaries/${sname}.web_summary.html
        cp ${sname}/outs/cloupe.cloupe ${OUTDIR}/${projid}/Summaries/${sname}_cloupe.cloupe

	"""

}

process fastqc {

	input:
        val projid from fqc_ch

        output:
        val projid into mqc_cha

	"""

        mkdir -p ${OUTDIR}/${projid}/QC
        mkdir -p ${OUTDIR}/${projid}/QC/FastQC

        for file in ${OUTDIR}/${projid}/Fastq_Raw/*/*fastq.gz
            do fastqc \$file --outdir=${OUTDIR}/${projid}/QC/FastQC
        done
	"""
    
}

process multiqc {

    input:
    val projid from mqc_cha

    output:
    val "x" into multiqc_outch

    script:
    """
    cd $OUTDIR/$projid
    multiqc . --outdir ${OUTDIR}/$projid/QC/

    """
}


// Aggregation

process gen_aggCSV {

    input:
    set sid, sname, projid from crAgg_ch

    output:
    val projid into crAggregate

    when:
    params.aggregate

    """
    
    aggdir=$OUTDIR/$projid/Aggregate

    mkdir -p \$aggdir

    aggcsv=\$aggdir/${projid}_libraries.csv

    if [ -f \$aggcsv ]
    then
        if grep -q $sname \$aggcsv
        then
             echo ""
        else
             echo "${sname},${OUTDIR}/${projid}/Count_CR/${sname}/outs/molecule_info.h5" >> \$aggcsv
        fi
    else
        echo "library_id,molecule_h5" > \$aggcsv
        echo "${sname},${OUTDIR}/${projid}/Count_CR/${sname}/outs/molecule_info.h5" >> \$aggcsv
    fi


    """
}



process aggregate {

    publishDir "${OUTDIR}/${projid}/Aggregate/", mode: 'move', overwrite: true
  
    input:
    val y from crAggregate.collect()
    val x from count_agg.collect()
    val projid from agg_ch

    output:
    file "${projid}_agg" 

    when:
    params.aggregate

    """

    aggdir="$OUTDIR/$projid/Aggregate"

    cellranger aggr \
       --id=${projid}_agg \
       --csv=\${aggdir}/${projid}_libraries.csv \
       --normalize=mapped

    cp \${aggdir}/${projid}_agg/outs/cloupe.cloupe ${OUTDIR}/${projid}/Summaries/${projid}_Aggregate_cloupe.cloupe
    cp \${aggdir}/${projid}_agg/outs/web_summary.html ${OUTDIR}/${projid}/Summaries/${projid}_Aggregate_web_summary.html

  
    """

}

