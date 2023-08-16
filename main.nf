nextflow.enable.dsl=2

params.db = "gorg_dark_v1"

/// Default settings 

params.outdir = "results"
params.cpus = 8
params.mismatches = 3
params.minlength = 11

/// For the selected mode (e.g. 'gorg_dark_v1') find the required accessory files

params.seqs = false
if (!params.seqs) { exit 1, "--seqs is not defined" }
if (params.db == "gorg_dark_v1") {
    VAL_gtdb_release = 'R207'
    refdir = "/mnt/scgc_nfs/ref/gorg_classifier_2022/gorg_dark_v1"
    annotations=file("${refdir}/annotations.tsv")
    lineages=file("${refdir}/gtdb_R207_lineages.tsv")
    clusters=file("${refdir}/gorg_dark_60minid_m80.tsv")
    nodes=file("${refdir}/nodes.dmp")
    names=file("${refdir}/names.dmp")
    fmi=file("${refdir}/index.fmi") }

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

workflow {
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
