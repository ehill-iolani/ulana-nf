process AMRFINDER_DB {
    label 'process_low'
    container 'quay.io/biocontainers/ncbi-amrfinderplus:4.2.7--hf69ffd2_0'
    publishDir "${params.outdir}/db", mode: 'copy'

    output:
    path "amrfinder_db", emit: db

    script:
    // amrfinder_update -d writes a dated subdir plus an "amrfinder_db/latest"
    // symlink; Nextflow's own containerized process boundary is fresh for
    // every task (unlike the shared conda env the Snakemake rule relied on),
    // so the resolved db has to be published as a plain directory here and
    // passed explicitly to AMRFINDER below, rather than relying on amrfinder's
    // internal default db path
    """
    amrfinder_update -d raw_db
    mkdir -p amrfinder_db
    cp -rL raw_db/latest/. amrfinder_db/
    """
}
