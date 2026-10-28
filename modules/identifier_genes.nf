process IDENTIFIER_GENES {
    tag "$sample"
    label 'process_low'
    container "${params.dockerhub_registry}/library/ubuntu:22.04"
    publishDir(path: { "${params.outdir}/${sample}/id_genes" }, mode: 'copy')

    input:
    tuple val(sample), path(ffn)

    output:
    tuple val(sample), path("16S_rRNA.fasta"), path("dnaA.fasta"), path("rpoB.fasta"), emit: genes

    script:
    // one record per gene copy (see bin/pick_genes.sh). The rpoB pattern is
    // anchored at the end so it doesn't also match rpoC ("... subunit beta'");
    // the 16S one also takes Prokka's "16S ribosomal RNA (partial)"
    """
    pick_genes.sh ${ffn} '^16S ribosomal RNA' '16S ribosomal RNA' 16S_rRNA.fasta
    pick_genes.sh ${ffn} '^Chromosomal replication initiator protein DnaA[[:space:]]*\$' 'Chromosomal replication initiator protein DnaA' dnaA.fasta
    pick_genes.sh ${ffn} '^DNA-directed RNA polymerase subunit beta[[:space:]]*\$' 'DNA-directed RNA polymerase subunit beta' rpoB.fasta
    """
}
