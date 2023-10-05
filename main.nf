#!/user/bin/nextflow
nextflow.enable.dsl=2

params.db = "gorg_dark_v1"

/// INPUTS
DIR_input = "/mnt/scgc_nfs/lab/ggavelis/notebooks/classifier/gom/gom_3/input/contigs/"
SUFFIX_of_input_files = ".fasta"

// OUTPUT
params.outdir = "results"

/// DATABASES
params.prokka = "/mnt/scgc_nfs/ref/uniprot_swissprot_prokka.fasta"
params.gtdb = "/mnt/scgc_nfs/ref/gtdb/release207" //release207 = db version for gtdbtk_2_0_0
DIR_taxdump = "/mnt/scgc_nfs/ref/taxdump/release207/" // these files were downloaded from
// nodes=file("${DIR_taxdump}/nodes.dmp")
// names=file("${DIR_taxdump}/names.dmp")

/// Default settings 
params.publishmode = 'Copy'
params.cpus = 8
params.mismatches = 3  // kaiju
params.minlength = 11  // kaiju


// DEV MODE STUFF
params.dev = false // Lets user testrun nextflow command (by adding flag '--dev') which will have this pipeline run on JUST ONE SAG (again, as a test)
params.num_inputs = 10
params.well = false // To run pipeline on a specific sample e.g 'AM-306-M06' use like: "nextflow run [this_script].nf --well='AM-306-M06' "
params.classif_by_ani = false

workflow {
	// Look for the reference database to classify the reads against
	if (params.dbdir == "/mnt/scgc_nfs/ref/gorg_classifier_2022/gorg_dark_v1") {
    	VAL_gtdb_release = 'R207'
    	annotations=file("${dbdir}/annotations.tsv")
    	lineages=file("${dbdir}/gtdb_R207_lineages.tsv")
    	clusters=file("${dbdir}/gorg_dark_60minid_m80.tsv")
    	nodes=file("${dbdir}/nodes.dmp")
    	names=file("${dbdir}/names.dmp")
    	fmi=file("${dbdir}/index.fmi") }

    // Print settings to screen
	log.info """

    =======================================================================
    GORG Classifier V2, Single Cell Genome Center, Bigelow Laboratory
    =======================================================================

    #### Authors
    Greg Gavelis <zoark0@gmail.com>
    -----------------------------------------------------------------------

    Sequences          (*.fq/*.fna)   : ${params.seqs}
    Mode                              : ${params.db}
    Nodes              (nodes.dmp)    : ${nodes}
    Names              (names.dmp)    : ${names}
    Kaiju Index        (.fmi)         : ${fmi}
    GORG Annotations   (.tsv)         : ${annotations}
    Output directory                  : ${params.outdir}
    Kaiju mismatches                  : ${params.mismatches}
    Kaiju minimum alignment length    : ${params.minlength}
    Kaiju CPUs                        : ${params.cpus}
    -----------------------------------------------------------------------

    """.stripIndent()

	if (!params.seqs) { exit 1, "--seqs is not defined" }
    seqs = Channel
        .fromFilePairs(params.seqs, size: -1, checkIfExists: true, flat: true)
        .map { it ->
            if (it.size == 2) {
                it.add([])
            }
            return it
        }

    KAIJU_v1_9_2(seqs, nodes, fmi)
    CLASSIFY_AND_ANNOTATE_READS(KAIJU_v1_9_2.out, annotations, clusters, lineages)
    CLUSTERWISE_ANNOTATION(CLASSIFY_AND_ANNOTATE_READS.out, annotations, clusters, lineages)
}

workflow make_new_db {

	if (!params.dbdir) { exit 1, "--dbdir is not defined"}
	CH_taxdump_config_files = Channel.fromPath("${DIR_taxdump}/{names,nodes,merged,delnodes}.dmp").collect() // Looks in DIR_taxdump for files called names.dmp, nodes.dmp, merged.dmp, and delnodes.dmp
	GTDB_LINEAGE_TABLE(CH_taxdump_config_files)
	CH_lineage_AND_taxid = GTDB_LINEAGE_TABLE.out.splitCsv(header: ['taxid', 'lineage'], skip: 1 ).map{ it -> tuple(it.lineage, it.taxid)}
	// CH_lineage_AND_taxid.view()
	
  CH_sample_AND_fasta = channel
    .fromPath("$DIR_input/*$SUFFIX_of_input_files", checkIfExists:true) // --> path(fasta)
    .filter({it.isEmpty() != true})            // ignore empty files
    .map { it -> tuple(it.getBaseName(), it) } // --> tuple( val(sample), path(fasta))    #  1. Extracts sampleID using .getBaseName(), then stores it in a tuple (with the inputfile path)
    .filter({it[0].length() < 11})             // input filename  must be under 11 characters (not counting file suffix, e.g. ".fasta")
    .filter({it[0].contains('_') != true})     // input filename must not contain underscore
    .filter({it[0].contains(';') != true})     // input filename must not contain semicolon
    //.view()

  // #### TROUBLESHOOT (1) Use '--dev' flag to run pipeline ON 1st SAMPLE ONLY ########
  // 'nextflow run this_pipeline.nf --dev'

  // #### TROUBLESHOOT (2) Use '--well=' flag to run pipeline ON 1 SPECIFIC SAMPLE #######
  // E.g. if target sample is "AG-510-B11.fasta", run pipeline like:
  // 'nextflow run this_pipeline.nf --well="AG-510-B11"'
  if( params.well != false )
		CH_sample_AND_fasta = CH_sample_AND_fasta .filter( { it[0] == params.well} )

	// ##### CLASSIFY SAGS ######### ________________________________________________________________________

	RENAME_CONTIGS(CH_sample_AND_fasta.take( params.dev ? params.num_inputs : -1))
	GTDBTK_v2_0_0(RENAME_CONTIGS.out)
	PARSE_GTDBTK(GTDBTK_v2_0_0.out)
	// PARSE_GTDBTK.out.collectFile(name: 'sample.txt', newLine: false)
  //   .subscribe {
  //       println "Entries are saved to file: $it"
  //       println "File content is: ${it.text}"
  //   }

	LOG_GTDBTK(PARSE_GTDBTK.out.collect())
	ADD_TAXIDS(LOG_GTDBTK.out, GTDB_LINEAGE_TABLE.out)

	CH_sample_AND_classif_AND_taxid = ADD_TAXIDS.out.splitCsv(header: ['sag', 'classification', 'taxid'], skip: 1).map { it -> tuple(it.sag, it.classification, it.taxid)}
	//     Emits (e.g. [AH-707-C14, d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Halieaceae;g__Congregibacter;s__])

	// Move forward only with samples that were classified by GTDB-Tk  --> CH_ok_sample_AND_fasta
	CH_ok_sample_AND_classif_AND_taxid = CH_sample_AND_classif_AND_taxid.filter( { it[1] != '0'})  // Unclassifieds have '0' in their gtdb 'classification' field

	CH_ok_sample_AND_classif_AND_taxid.count().subscribe {  println "$it samples were classified by GTDB-Tk, then had taxids assigned by Taxdump"  }

	CH_ok_sample_AND_classif_AND_taxid_AND_fasta = CH_ok_sample_AND_classif_AND_taxid.join(CH_sample_AND_fasta)
	CH_ok_sample_AND_fasta = CH_ok_sample_AND_classif_AND_taxid_AND_fasta.map{it -> tuple(it[0],it[3])}

	//CLASSIF_TO_TAXID(CH_id_AND_classif_AND_fasta.take( params.dev ? params.num_inputs : -1), CH_taxdump_config_files)

	// ##### USE '--classif_by_ani' flag to TRY TO CLASSIFY SAGS ON WHICH GTDB-TK FAILED ##### 
	// Run pipeline like:
	//  'nextflow run this_pipeline.nf --classif_by_ani'
	if( params.classif_by_ani != false )
	{
		// Which samples weren't classified?
		CH_unclass_sample = CH_sample_AND_classif.filter( { it[1] == '0'}) //.view()
		// Find those fastas
		CH_unclass_sample_AND_fasta = CH_unclass_sample.join(CH_sample_AND_fasta).map( it -> tuple(it[0], it[2]))//.view()
		
		// Now for *classified* samples, find those that GTDBTk classified to genus or species (will use as references to map unclassifieds against)
		CH_classif_sample = CH_sample_AND_classif.filter( { it[1] != '0'})
		CH_samples_known_genus = CH_classif_sample
			.filter({it[1].split(';')[5] != 'g__'}) // genus is known if it's more than just "g__" 

		// Get those genus-classified fastas (as reference genomes)
		CH_fastas_known_genus = CH_samples_known_genus.join(CH_sample_AND_fasta).map( it -> it[2]).collect()

		// FastANI needs paths to those reference genomes written as a list
		LIST_GENUS_FASTAS(CH_fastas_known_genus)
		
		// Run FastANI. To ask: For each *un*classified fasta, what amount of ANI does it shared with the classified ones?
		FASTANI_AGAINST_KNOWN_GENERA_v1_3_2(CH_unclass_id_AND_fasta, LIST_GENUS_FASTAS.out)
	}

	// ######## RENAME_CONTIGS ############ ____________________________________________________________
	

	PROKKA_v1_14_6(CH_ok_sample_AND_fasta)

	// combine prokka seqs and annots with gtdb_classif and taxid
	PROKKA_2_TABLE(PROKKA_v1_14_6.out.join(CH_sample_AND_classif_AND_taxid))
}


process ADD_TAXIDS{
		container = 'brwnj/kmernorm:v1.0.0'
		publishDir "${params.outdir}/", mode: params.publishmode
		errorStrategy 'finish'
		input:
			path(classifications)
			path(lineages)
		output:
			path("sag_2_taxid.csv")
		script: template "add_taxids.py" }

process CLASSIFY_AND_ANNOTATE_READS {
    tag "$sample"
    publishDir "$params.outdir", mode: 'copy'
    conda '/mnt/scgc/scgc_nfs/opt/common/anaconda3a'
    memory = { 8.GB * task.attempt }

    input:
    tuple val(sample), path(kaiju_hits)
    path(annotations)
    path(clusters)
    path(lineages)
    output:
    tuple val(sample), path("${sample}_reads_classified_v_${params.mode}_GTDB_${VAL_gtdb_release}.tsv.gz")

    script:
    template 'classify_and_annotate_reads.py' }

process CLUSTERWISE_ANNOTATION {
    tag "$sample"
    publishDir "$params.outdir", mode: 'copy'
    errorStrategy = 'finish'
    conda '/mnt/scgc/scgc_nfs/opt/common/anaconda3a'
    memory = { 8.GB * task.attempt }

    input:
    tuple val(sample), path(classified_reads)
    path(annotations)
    path(clusters)
    path(lineages)
    output:
    tuple val(sample), path("${sample}_clusters_v_${params.db}_GTDB_${VAL_gtdb_release}.tsv.gz")

    script:
    template 'clusterwise_annotation.py' }

process GTDBTK_v2_0_0  {
	errorStrategy = 'ignore'
	tag "${sample}"
	memory='10 GB'
  	cpus=6
  	container='docker://quay.io/biocontainers/gtdbtk:2.0.0--pyhdfd78af_1'
	//publishDir "${params.outdir}/${sample}/annotation_${sample}", mode: params.publishmode
	input: tuple val(sample), path(contigs)
	output: tuple val(sample), path("gtdbtk_classification_${sample}")
	"""
	mkdir tmp_dir; cp ${contigs} ./tmp_dir/final_contigs_${sample}.fasta; export GTDBTK_DATA_PATH=${params.gtdb}
	gtdbtk classify_wf --genome_dir tmp_dir --out_dir gtdbtk_classification_${sample} --cpus ${task.cpus} -x fasta
	rm -r tmp_dir
	""" }

process GTDB_LINEAGE_TABLE{
	publishDir "${params.outdir}/", mode: params.publishmode
	conda '/mnt/scgc/scgc_nfs/opt/common/anaconda3a/envs/taxonkit'
	input: path(taxdump_config_files)
	output: path("gtdb_lineages.tsv")
	script:
	"""
	# Converts names.dmp (via taxdump) to a table and fills in skipped GTDB lineages.
	taxonkit list --data-dir . --ids 1 -I "" \
    | taxonkit lineage  --data-dir . \
    | taxonkit reformat --data-dir . -I 1 -F -P --prefix-k d__ -p "" \
    | sed 's/ phylum//' | sed 's/ class//' | sed 's/ order//' | sed 's/ family//' | sed 's/ genus//' | sed 's/ species//' \
    | cut -f 1,3 \
    | csvtk add-header -t -n taxid,lineage \
 	>  gtdb_lineages.tsv
	""" }

process KAIJU_v1_9_2 {
    tag "$sample"
    // publishDir "$params.outdir"
    cpus params.cpus
    memory = { 16.GB * task.attempt }
    //container = 'brwnj/classifier-nf:1.1.0'
    container='docker://quay.io/biocontainers/kaiju:1.9.2--h43eeafb_3'

    input:
    tuple val(sample), path(r1), path(r2)
    path(nodes)
    path(fmi)

    output:
    tuple val(sample), path("${sample}_hits.txt")

    script:
    def r2path = r2 ? "-j ${r2}" : ""
    """
    kaiju -z ${task.cpus} -v -m ${params.minlength} \
        -e ${params.mismatches} -t $nodes -f $fmi \
        -i ${r1} ${r2path} -o ${sample}_hits.txt
    """ }

process PROKKA_2_TABLE{
	tag "${sample}"
	container='brwnj/kmernorm:v1.0.0'
	publishDir "${params.outdir}/test/", mode: params.publishmode, pattern: "*tsv"
	memory = { 1.GB * task.attempt }
  errorStrategy = 'finish'
  maxRetries = 3
	input: tuple val(sample), path(faa), path(gff)
	output:
		tuple val(sample), path("${sample}_contigs_proteins_annot.tsv"), emit: tsv
	script: template "prokka_2_table.py" }

process PROKKA_v1_14_6{
	//publishDir "${params.outdir}/test/", mode: params.publishmode
	tag "${sample}"
	memory = { 10.GB * task.attempt }
  errorStrategy = 'retry'
  cpus=6
  maxRetries = 3
  container='docker://quay.io/biocontainers/prokka:1.14.6--pl5262hdfd78af_1'
	input: tuple val(sample), path(contigs)
	output: tuple val(sample), path("prokka_${sample}/${sample}.faa"), path("prokka_${sample}/${sample}.gff")
	script: "prokka --outdir prokka_${sample} --prefix ${sample} --locustag ${sample} --quiet --compliant --force --proteins ${params.prokka} --cpus ${task.cpus} ${contigs}" }

// process TABULATE_CONTIG_PROTEIN_ANNOT{
//	tag "${sample}"
// }

process RENAME_CONTIGS {
	tag "${sample}"
	//publishDir "${params.outdir}/renamed/", mode: params.publishmode
	errorStrategy = 'finish'
	container = 'brwnj/kmernorm:v1.0.0'
	input: tuple val(sample), path(contigs)
	output: tuple val(sample), path("${sample}_tmp_renamed.fasta")
	script:
	"""
	#!/usr/bin/env python
	import Bio
	from Bio import SeqIO
	NUM_contig = 0
	with open("${sample}_tmp_renamed.fasta", "w") as handle:
		for record in SeqIO.parse("${contigs}", "fasta"):
			NUM_contig += 1
			STR_short_sample_id = "${sample}".replace('-','')                # remove dashes
			STR_new_seqID = ">" + STR_short_sample_id + "_C" + str(NUM_contig)  # make new SeqID
			handle.write(STR_new_seqID+'\\n')
			handle.write(str(record.seq+'\\n'))
	""" }

process PARSE_GTDBTK{
	tag "${sample}"
	container = 'brwnj/kmernorm:v1.0.0'
	// publishDir "${params.outdir}/${sample}/logs_${sample}", mode: params.publishmode
	input: tuple val(sample), path(gtdbtk_dir)
	output: path("gtdbtk_stats_${sample}.csv")
	script:
	"""
	#!/usr/bin/env python
	from os.path import exists; import pandas as pd
	classification_via_GTDBTk=0

	TSV_gtdb_bac = "${gtdbtk_dir}/gtdbtk.bac120.summary.tsv"
	TSV_gtdb_ar  = "${gtdbtk_dir}/gtdbtk.ar53.summary.tsv"

	if exists(TSV_gtdb_bac) and not exists(TSV_gtdb_ar): #If its classified as bacterial
		DF_gtdb = pd.read_csv(TSV_gtdb_bac, sep='\\t')
		classification_via_GTDBTk = DF_gtdb.loc[0,'classification']
	elif exists(TSV_gtdb_ar) and not exists(TSV_gtdb_bac): # If archaeal
		DF_gtdb = pd.read_csv(TSV_gtdb_ar, sep='\\t')
		classification_via_GTDBTk = DF_gtdb.loc[0,'classification']
	elif exists(TSV_gtdb_ar) and exists(TSV_gtdb_bac): # If both
		print('warning: It may be ambiguous as to whether this is archaeal or bacterial!')
	else:
		print('warning: No classification files were found')

	with open("gtdbtk_stats_${sample}.csv", "w") as out:
		print("${sample}",classification_via_GTDBTk,sep=",",file=out)
	""" }

process LOG_GTDBTK {
	publishDir "${params.outdir}/", mode: "copy"
	container = 'brwnj/kmernorm:v1.0.0'
	input: path(countfiles)
	output: path("gtdbtk_classified.csv"), emit: csv
	shell: ''' echo "sag,classification" > gtdbtk_classified.csv; for LINE in !{countfiles}; do cat ${LINE} >> gtdbtk_classified.csv; done ''' }

process LIST_GENUS_FASTAS{
	errorStrategy = 'finish'
	tag "writing list of reference sags for FastANI"
	//container = 'brwnj/kmernorm:v1.0.0'
	publishDir "${params.outdir}/fastani", mode: "copy"
	input: path(all_fastas)
	output: path("list_all_fastas.txt")
	shell:
	'''
	touch list_all_fastas.txt
	for PATH in !{all_fastas}; do /usr/bin/readlink -f ${PATH} >> list_all_fastas.txt; done
	'''}

process FASTANI_AGAINST_KNOWN_GENERA_v1_3_2{
	tag "${sample}"
	errorStrategy = 'ignore'
	publishDir "${params.outdir}/fastani", mode: "copy"
	conda='/mnt/scgc/scgc_nfs/opt/common/anaconda3a/envs/fastani_1.32'
	input:
	tuple val(sample), path(contigs)
	path(list_fastas)
	output: path("${sample}_fastani.txt")
	script: "fastANI -q ${contigs} --rl ${list_fastas} -o ${sample}_fastani.txt"}

