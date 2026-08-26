process MEDAKA_POLISH {
    tag "$sample"
    label 'process_high'
    // ONT's own multi-arch image, not biocontainers -- see edna-ont-nf's
    // MEDAKA module for why (amd64-only biocontainers build SIGILLs under
    // Docker emulation on Apple Silicon). Pinned to a digest tag (not the
    // floating `latest`) rather than a numbered release: as of this writing
    // the newest numbered tag is v1.11.3 (medaka 1.x), whose bundled model
    // list tops out at v4.3.0 and rejects the v5.0.0 default model this
    // pipeline uses (matching ulana-ht's config.yaml); this digest is
    // medaka 2.2.1, which knows v5.0.0.
    container 'ontresearch/medaka:shaf39439188c323053e66cf79618fae1ab33d1b38a'
    publishDir(path: { "${params.outdir}/${sample}/polish/medaka" }, mode: 'copy')

    input:
    tuple val(sample), path(reads), path(draft)

    output:
    tuple val(sample), path("polished_consensus.fasta"), emit: polished

    script:
    """
    medaka_consensus -i ${reads} -d ${draft} -o medaka_out -t ${task.cpus} -m ${params.medaka_model}
    medaka sequence medaka_out/*.hdf ${draft} polished_consensus.fasta
    """
}
