process CHOPPER {
    tag "$sample"
    label 'process_low'
    container 'quay.io/biocontainers/chopper:0.7.0--hdcf5f25_0'
    publishDir(path: { "${params.outdir}/${sample}/reads" }, mode: 'copy')

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.filtered.fastq.gz"), emit: filtered

    script:
    """
    zcat -f ${fastq} \\
      | chopper -q ${params.chopper_q} --minlength ${params.chopper_minlength} --threads ${task.cpus} \\
      | gzip > ${sample}.filtered.fastq.gz
    """
}
