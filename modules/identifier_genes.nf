process IDENTIFIER_GENES {
    tag "$sample"
    label 'process_low'
    container 'ubuntu:22.04'
    publishDir(path: { "${params.outdir}/${sample}/id_genes" }, mode: 'copy')

    input:
    tuple val(sample), path(ffn)

    output:
    tuple val(sample), path("16S_rRNA.fasta"), path("dnaA.fasta"), path("rpoB.fasta"), emit: genes

    script:
    """
    echo '>16S ribosomal RNA' > 16S_rRNA.fasta
    awk '/16S ribosomal RNA/{flag=1;next}/^>/{flag=0}flag' ${ffn} >> 16S_rRNA.fasta

    echo '>Chromosomal replication initiator protein DnaA' > dnaA.fasta
    awk '/Chromosomal replication initiator protein DnaA/{flag=1;next}/^>/{flag=0}flag' ${ffn} >> dnaA.fasta

    echo '>DNA-directed RNA polymerase subunit beta RpoB' > rpoB.fasta
    awk '/DNA-directed RNA polymerase subunit beta RpoB/{flag=1;next}/^>/{flag=0}flag' ${ffn} >> rpoB.fasta
    """
}
