process AMRFINDER {
    tag "$sample"
    label 'process_medium'
    container "${params.container_registry}/biocontainers/ncbi-amrfinderplus:4.2.7--hf69ffd2_0"
    publishDir(path: { "${params.outdir}/${sample}/amr" }, mode: 'copy')

    input:
    tuple val(sample), path(faa), path(gff), path(fna)
    path amrfinder_db

    output:
    tuple val(sample), path("amrfinder_pro_results.tsv"), emit: tsv
    tuple val(sample), path("amrfinder_pro.fasta"), emit: proteins

    script:
    // Combined mode: -n adds a translated search of the contigs themselves
    // (Prokka's .fna, so the names match the GFF), on top of the search of
    // Prokka's proteins. A protein-only search misses a gene that is in the
    // assembly but that Prokka called short (a start at an internal codon)
    // or that runs off a contig end; those now come back with a method
    // ending in X (EXACTX, ALLELEX, PARTIAL_CONTIG_ENDX...) and Protein id NA.
    """
    amrfinder \\
      -a prokka \\
      --threads ${task.cpus} \\
      -p ${faa} \\
      -g ${gff} \\
      -n ${fna} \\
      -d ${amrfinder_db} \\
      --print_node --plus \\
      --protein_output amrfinder_pro.fasta \\
      > amrfinder_pro_results.tsv
    """
}
