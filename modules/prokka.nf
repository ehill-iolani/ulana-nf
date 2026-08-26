process PROKKA {
    tag "$sample"
    label 'process_medium'
    container 'quay.io/biocontainers/prokka:1.14.6--pl5321hdfd78af_4'
    publishDir(path: { "${params.outdir}/${sample}/annotation/prokka" }, mode: 'copy')

    input:
    tuple val(sample), path(assembly)

    output:
    tuple val(sample), path("prokka_annotation.tsv"), emit: tsv
    tuple val(sample), path("prokka_annotation.faa"), emit: faa
    tuple val(sample), path("prokka_annotation.ffn"), emit: ffn
    tuple val(sample), path("prokka_annotation.gff"), emit: gff

    script:
    """
    prokka --outdir . ${assembly} --prefix prokka_annotation --cpus ${task.cpus} --force
    """
}
