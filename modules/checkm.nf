process CHECKM {
    tag "$sample"
    // not process_medium/process_high -- lineage_wf's pplacer tree-placement
    // step is known to need well above the 8GB process_medium gives it, even
    // with --reduced_tree (the Snakemake original never hit this because it
    // had no per-process memory cap; Nextflow's `memory` directive becomes a
    // hard `docker --memory` limit, so undersizing it here means a silent
    // OOM-kill deep inside pplacer, not a normal "out of memory" error).
    // Override per-environment with `-c` / `process { withName: CHECKM {...} }`.
    // cpus   = 4
    // memory = '16.GB'
    // time   = '4.h'
    container 'quay.io/biocontainers/checkm-genome:1.2.2--pyhdfd78af_1'
    publishDir(path: { "${params.outdir}/${sample}/qc/checkm" }, mode: 'copy')

    input:
    tuple val(sample), path(assembly)
    path checkm_db

    output:
    tuple val(sample), path("summary.tsv"), emit: summary

    script:
    // checkm insists on scanning a whole directory of genomes (-x fasta), not
    // a single file path, so it has to be staged into its own bin/ dir first.
    // `data setRoot` has to be re-run against CHECKM_DB's published directory
    // in every task -- each task starts from a fresh container, so nothing
    // from a prior `setRoot` call persists between tasks
    """
    mkdir -p bin
    cp ${assembly} bin/${sample}.fasta

    checkm data setRoot ${checkm_db}
    checkm lineage_wf --reduced_tree -t ${task.cpus} -x fasta bin .
    checkm qa -o 2 --tab_table -t ${task.cpus} -f summary.tsv lineage.ms .
    """
}
