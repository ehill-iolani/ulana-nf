// Caps each sample's read depth before assembly. Very deep runs make Flye fail:
// a ~91 kb phage at ~5,000x stopped with "No disjointigs were assembled", while
// a 1/20 subsample of the same reads assembled fine.
//
// DOWNSAMPLE_PLAN counts the filtered bases and gets the genome size (from
// --genome_size, or else estimated by LRGE from the reads), then decides:
//   - reads above --target_coverage x genome size -> DOWNSAMPLE (Rasusa) takes
//     a random subset of about that many bases;
//   - otherwise, or if the size can't be estimated -> every filtered read is
//     used, as before this step existed.
// Rasusa picks reads at random rather than by quality or length, so small
// plasmids keep their share of the reads (length-weighted tools like Filtlong
// can drop them). The plan is published as {sample}.downsampling.tsv either way.

process DOWNSAMPLE_PLAN {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/lrge:1.0.0--hf1b6044_0"
    publishDir(path: { "${params.outdir}/${sample}/reads" }, mode: 'copy')

    input:
    tuple val(sample), path(reads)
    val genome_size   // bp from --genome_size, or 0 to estimate it

    output:
    tuple val(sample), env('DECISION'), env('TARGET_BASES'), emit: plan
    tuple val(sample), path("${sample}.downsampling.tsv"), emit: summary

    script:
    """
    # chopper writes 4-line records, so every 4th line from the 2nd is sequence
    READ_BASES=\$(zcat -f ${reads} | awk 'NR % 4 == 2 { b += length(\$0) } END { printf "%d", b }')

    if [ "${genome_size}" -gt 0 ]; then
        GENOME_SIZE=${genome_size}
        SOURCE=given
    else
        # seeded so a rerun makes the same estimate; a failed or empty estimate
        # (e.g. too few reads) means "unknown", and then every read is kept
        GENOME_SIZE=\$(lrge -t ${task.cpus} -s 1 ${reads} 2> lrge.log || true)
        SOURCE=estimated
        case "\$GENOME_SIZE" in
            ''|*[!0-9]*|0) GENOME_SIZE=''; SOURCE=unknown ;;
        esac
    fi

    if [ -n "\$GENOME_SIZE" ]; then
        TARGET_BASES=\$(( ${params.target_coverage} * GENOME_SIZE ))
    else
        TARGET_BASES=0
    fi
    if [ "\$TARGET_BASES" -gt 0 ] && [ "\$READ_BASES" -gt "\$TARGET_BASES" ]; then
        DECISION=downsample
        ASSEMBLY_BASES=\$TARGET_BASES
    else
        DECISION=all
        ASSEMBLY_BASES=\$READ_BASES
    fi

    cov() { [ -n "\$GENOME_SIZE" ] && awk -v b="\$1" -v g="\$GENOME_SIZE" 'BEGIN { printf "%.1f", b / g }'; }
    {
        printf 'sample\\tread_bases\\tgenome_size\\tgenome_size_source\\tread_coverage\\ttarget_coverage\\treads_used\\tassembly_bases\\tassembly_coverage\\n'
        printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
            "${sample}" "\$READ_BASES" "\$GENOME_SIZE" "\$SOURCE" "\$(cov \$READ_BASES)" \\
            "${params.target_coverage}" "\$([ \$DECISION = downsample ] && echo subsample || echo all)" \\
            "\$ASSEMBLY_BASES" "\$(cov \$ASSEMBLY_BASES)"
    } > ${sample}.downsampling.tsv
    """
}

process DOWNSAMPLE {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/rasusa:5.1.0--hfa8f182_0"

    input:
    tuple val(sample), path(reads), val(target_bases)

    output:
    tuple val(sample), path("${sample}.downsampled.fastq.gz"), emit: reads

    script:
    // fixed seed: the same reads give the same subset on a rerun
    """
    rasusa reads --bases ${target_bases} --seed 1 -o ${sample}.downsampled.fastq.gz ${reads}
    """
}
