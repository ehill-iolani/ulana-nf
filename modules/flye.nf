process FLYE_ASSEMBLY {
    tag "$sample"
    label 'process_high'
    container "${params.container_registry}/biocontainers/flye:2.9.6--py312h734f728_1"
    // the .version is for the GenBank assembly comment, not a result
    publishDir(path: { "${params.outdir}/${sample}/assembly/flye" }, mode: 'copy', saveAs: { it.endsWith('.version') ? null : it })

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("assembly.fasta"), emit: assembly
    tuple val(sample), path("assembly_info.txt"), emit: info
    tuple val(sample), path("assembly_graph.gfa"), emit: gfa
    tuple val(sample), path("flye.version"), emit: version

    script:
    """
    flye ${params.flye_mode} ${reads} --out-dir . --threads ${task.cpus}
    # "2.9.6-b1802" -> "2.9.6"
    flye --version | cut -d- -f1 > flye.version
    """
}
