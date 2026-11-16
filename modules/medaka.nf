process MEDAKA_POLISH {
    tag "$sample"
    label 'process_high'
    // ONT's own multi-arch image, not biocontainers -- see edna-ont-nf's
    // MEDAKA module for why (amd64-only biocontainers build SIGILLs under
    // Docker emulation on Apple Silicon). Pinned to a digest tag (not the
    // floating `latest`) rather than a numbered release: as of this writing
    // the newest numbered tag is v1.11.3 (medaka 1.x), whose bundled model
    // list tops out at v4.3.0 and lacks both v5.0.0 models and the default
    // r1041_e82_400bps_bacterial_methylation; this digest is medaka 2.2.1,
    // which knows both.
    container "${params.dockerhub_registry}/ontresearch/medaka:shaf39439188c323053e66cf79618fae1ab33d1b38a"
    // the .version is for the GenBank assembly comment, not a result
    publishDir(path: { "${params.outdir}/${sample}/polish/medaka" }, mode: 'copy', saveAs: { it.endsWith('.version') ? null : it })

    input:
    tuple val(sample), path(reads), path(draft)

    output:
    tuple val(sample), path("polished_consensus.fasta"), emit: polished
    tuple val(sample), path("medaka.version"), emit: version

    script:
    """
    medaka_consensus -i ${reads} -d ${draft} -o medaka_out -t ${task.cpus} -m ${params.medaka_model}
    medaka sequence medaka_out/*.hdf ${draft} polished_consensus.fasta
    # "medaka 2.2.1" -> "2.2.1"
    medaka --version 2>&1 | awk '{ print \$NF }' > medaka.version
    """
}
