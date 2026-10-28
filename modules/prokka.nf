process PROKKA {
    tag "$sample"
    label 'process_medium'
    container "${params.container_registry}/biocontainers/prokka:1.14.6--pl5321hdfd78af_4"
    // the .fsa/.tbl only feed GENBANK_PREP, which publishes its own fixed
    // copies; the .version is for the GenBank assembly comment
    publishDir(path: { "${params.outdir}/${sample}/annotation/prokka" }, mode: 'copy', pattern: '*.{tsv,faa,ffn,gff}')

    input:
    tuple val(sample), path(assembly), val(meta)

    output:
    tuple val(sample), path("prokka_annotation.tsv"), emit: tsv
    tuple val(sample), path("prokka_annotation.faa"), emit: faa
    tuple val(sample), path("prokka_annotation.ffn"), emit: ffn
    tuple val(sample), path("prokka_annotation.gff"), emit: gff
    tuple val(sample), path("prokka_annotation.fsa"), path("prokka_annotation.tbl"), emit: genbank
    tuple val(sample), path("prokka.version"), emit: version

    script:
    // for a GenBank submission the locus tags take the sample's
    // NCBI-assigned prefix instead of Prokka's made-up one, and the organism
    // goes into the .gbk; the samplesheet values were checked in main.nf
    def genbank = ''
    if (params.enable_genbank) {
        def words = meta.organism.split(/\s+/)
        genbank = "--locustag ${meta.locus_tag} --genus '${words[0]}' --species '${words[1..-1].join(' ')}'"
        if (meta.strain) {
            genbank += " --strain '${meta.strain}'"
        }
    }
    """
    prokka --outdir . ${assembly} --prefix prokka_annotation --cpus ${task.cpus} --force ${genbank}
    prokka --version 2>&1 | awk '{ print \$NF }' > prokka.version
    """
}
