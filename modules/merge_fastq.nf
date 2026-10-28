process MERGE_FASTQ {
    tag "$sample"
    label 'process_low'
    container "${params.dockerhub_registry}/library/ubuntu:22.04"
    publishDir(path: { "${params.outdir}/${sample}/reads" }, mode: 'copy')

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.merged.fastq.gz"), emit: merged

    script:
    // a samplesheet row's fastq field may be a single file or a glob matching
    // several ONT part-files; cat them into one file per sample so downstream
    // steps don't need to care how many parts came in.
    //
    // Each gzipped part is checked first (this image has GNU gzip). Without
    // the check a damaged part only truncates the sample: zcat stops at the
    // bad part, but the pipe's status is the final gzip's, so the task still
    // succeeds with every read after that part silently lost. Plain FASTQ
    // parts (no gzip magic) are passed through unchecked, as zcat -f does.
    def infiles = fastq instanceof List ? fastq : [fastq]
    def cmds = infiles.collect { f -> f.name.endsWith('.gz') ? "zcat -f ${f}" : "cat ${f}" }.join(' ; ')
    """
    for f in ${infiles.join(' ')}; do
        if [ "\$(head -c 2 "\$f" | od -An -tx1 | tr -d ' \\n')" = "1f8b" ]; then
            if ! msg=\$(gzip -t "\$f" 2>&1); then
                echo "ERROR: sample ${sample}: input file \$f is damaged and can't be read (gzip -t: \${msg##*: })." >&2
                echo "Copy it again from the sequencer, or leave it out of the samplesheet." >&2
                exit 1
            fi
        fi
    done

    set -o pipefail
    { ${cmds} ; } | gzip > ${sample}.merged.fastq.gz
    """
}
