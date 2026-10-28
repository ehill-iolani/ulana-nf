// GenBank submission files from the Prokka annotation: GENBANK_PREP fixes
// what NCBI's validator rejects in Prokka's .fsa/.tbl (see
// bin/genbank_prep.py), then TABLE2ASN builds the .sqn and runs NCBI's own
// validation + discrepancy checks on it, so problems show up in the results
// instead of in an email from GenBank weeks later

process GENBANK_PREP {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/pandas:2.2.1"
    publishDir(path: { "${params.outdir}/${sample}/genbank" }, mode: 'copy')

    input:
    tuple val(sample), path(fsa), path(tbl), path(assembly_info), val(meta), val(assembly_method), val(annotation)

    output:
    tuple val(sample), path("${sample}.fsa"), path("${sample}.tbl"), path("${sample}.cmt"), emit: files
    tuple val(sample), path("submission_info.txt"), emit: info

    script:
    """
    genbank_prep.py \\
        --fsa ${fsa} \\
        --tbl ${tbl} \\
        --assembly-info ${assembly_info} \\
        --organism '${meta.organism}' \\
        --strain '${meta.strain}' \\
        --locus-tag '${meta.locus_tag}' \\
        --assembly-method '${assembly_method}' \\
        --annotation '${annotation}' \\
        --out-prefix ${sample}
    """
}

process TABLE2ASN {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/table2asn:1.28.1179--he45da00_1"
    publishDir(path: { "${params.outdir}/${sample}/genbank" }, mode: 'copy')

    input:
    tuple val(sample), path(fsa), path(tbl), path(cmt), val(annotation)
    path(template)   // the submitter template (.sbt), or [] for none

    output:
    tuple val(sample), path("${sample}.sqn"), emit: sqn
    tuple val(sample), path("${sample}.val"), path("${sample}.stats"), path("${sample}.dr"), emit: reports
    tuple val(sample), path("${sample}.gbf"), emit: gbf

    script:
    // -M n: genome submission checks (validation + discrepancy report, as
    // NCBI runs them), -J: proteins as a delayed genomic product set, -c f:
    // NCBI's own product-name cleanup, -V b: also write the GenBank flatfile
    // so the record can be read before it's sent. Without a template the
    // .sqn has no submitter or citation, and the .val says so (NoPubFound);
    // table2asn exits 0 either way, its findings are in the .val/.dr
    def template_arg = template ? "-t ${template}" : ''
    """
    table2asn -i ${fsa} -f ${tbl} -w ${cmt} -o ${sample}.sqn \\
        -M n -J -c f -Z -V b ${template_arg} \\
        -y 'Annotated using ${annotation} from https://github.com/tseemann/prokka'
    """
}
