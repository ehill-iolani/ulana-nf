process MERGE_FASTQ {
    tag "$sample"
    label 'process_low'
    container 'ubuntu:22.04'
    publishDir(path: { "${params.outdir}/${sample}/reads" }, mode: 'copy')

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.merged.fastq.gz"), emit: merged

    script:
    // a samplesheet row's fastq field may be a single file or a glob matching
    // several ONT part-files; cat them into one file per sample so downstream
    // steps don't need to care how many parts came in
    def infiles = fastq instanceof List ? fastq : [fastq]
    def cmds = infiles.collect { f -> f.name.endsWith('.gz') ? "zcat -f ${f}" : "cat ${f}" }.join(' ; ')
    """
    { ${cmds} ; } | gzip > ${sample}.merged.fastq.gz
    """
}
