process CHECKM_DB {
    label 'process_low'
    container 'quay.io/biocontainers/checkm-genome:1.2.2--pyhdfd78af_1'
    publishDir "${params.outdir}/db", mode: 'copy'

    output:
    path "checkm_data", emit: db

    script:
    // CheckM's reference tree/marker-set data (~275MB) isn't bundled in the
    // biocontainer and normally requires an interactive `checkm data setRoot`
    // one-time setup -- but each CHECKM task below gets a fresh container, so
    // that has to be re-run inside every task against this published
    // directory instead (see CHECKM's script for why)
    """
    mkdir -p checkm_data
    curl -sSL https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz \\
      | tar -xz -C checkm_data
    """
}
