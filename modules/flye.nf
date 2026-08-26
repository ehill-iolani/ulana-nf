process FLYE_ASSEMBLY {
    tag "$sample"
    label 'process_high'
    container 'quay.io/biocontainers/flye:2.9.6--py312h734f728_1'
    publishDir(path: { "${params.outdir}/${sample}/assembly/flye" }, mode: 'copy')

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("assembly.fasta"), emit: assembly
    tuple val(sample), path("assembly_info.txt"), emit: info
    tuple val(sample), path("assembly_graph.gfa"), emit: gfa

    script:
    """
    flye ${params.flye_mode} ${reads} --out-dir . --threads ${task.cpus}
    """
}
