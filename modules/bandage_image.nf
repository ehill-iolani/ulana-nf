process BANDAGE_IMAGE {
    tag "$sample"
    label 'process_low'
    container 'quay.io/biocontainers/bandage:0.9.0--h9948957_0'
    publishDir(path: { "${params.outdir}/${sample}/assembly/flye" }, mode: 'copy')

    input:
    tuple val(sample), path(gfa)

    output:
    tuple val(sample), path("assembly_graph.png"), emit: image

    script:
    // Bandage is a Qt GUI app -- QT_QPA_PLATFORM=offscreen is required for it
    // to render without a display (there's no windowing system inside a
    // container); without it, `Bandage image` fails outright rather than
    // just looking worse
    """
    QT_QPA_PLATFORM=offscreen Bandage image ${gfa} assembly_graph.png --height ${params.bandage_height} --lengths --depth --fontsize 20
    """
}
