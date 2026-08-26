process AMRFINDER {
    tag "$sample"
    label 'process_medium'
    container 'quay.io/biocontainers/ncbi-amrfinderplus:4.2.7--hf69ffd2_0'
    publishDir(path: { "${params.outdir}/${sample}/amr" }, mode: 'copy')

    input:
    tuple val(sample), path(faa), path(gff)
    path amrfinder_db

    output:
    tuple val(sample), path("amrfinder_pro_results.tsv"), emit: tsv
    tuple val(sample), path("amrfinder_pro.fasta"), emit: proteins

    script:
    """
    amrfinder \\
      -a prokka \\
      --threads ${task.cpus} \\
      -p ${faa} \\
      -g ${gff} \\
      -d ${amrfinder_db} \\
      --print_node --plus \\
      --protein_output amrfinder_pro.fasta \\
      > amrfinder_pro_results.tsv
    """
}
