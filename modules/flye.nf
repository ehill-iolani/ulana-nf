// Assembles each sample with Flye.
//
// Second try on very deep reads: DOWNSAMPLE_PLAN caps the depth before this
// step, but only as well as it knows the genome size, and a failed LRGE
// estimate keeps every read. ulana-nf v0.1.8 (before downsampling) stopped
// with "No disjointigs were assembled" at ~4,300x; ulana-phage-nf hit the same
// on a host-heavy sample LRGE sized wrong (583x, 2026-11-10). Flye measures
// that coverage itself before it gives up, so when it fails that way at more
// than --target_coverage, this takes a seeded random subset of target/measured
// of the reads and runs Flye once more. Any other failure, or a second
// failure, still fails the task. Only with --enable_downsampling, since
// turning that off asks for every read.
//
// Medaka still polishes with every read this step was given. flye_reads.tsv
// records which reads were assembled.
//
// One copy per genome: Flye sometimes reports a small circular plasmid as one
// contig holding two or more back-to-back copies. K. pneumoniae MGH 78578's
// 3,478 bp pKPN7 came out at 6,956 bp (validation, 2026-11-10); the cut copy
// matches RefSeq base for base. bin/collapse_copies.py cuts such a contig to
// one copy before Medaka, so polishing, Prokka and AMRFinder see each plasmid
// once. On the 20 validation drafts it cut only plasmid-size contigs (1.8-46 kb
// copies, 2-5 per contig) and never a chromosome. copies.tsv lists what was
// cut; assembly_info.txt and the graph stay as Flye made them. Every contig is
// checked: a 6.3 Mb chromosome takes under a second.
//
// The same module as ulana-phage-nf's apart from the --max-length limit
// (there, contigs over twice the phage size cap are host DNA and skipped).
// Keep fixes in step.

process FLYE_ASSEMBLY {
    tag "$sample"
    label 'process_high'
    container "${params.container_registry}/biocontainers/flye:2.9.6--py312h734f728_1"
    // the .version is for the GenBank assembly comment, not a result
    publishDir(path: { "${params.outdir}/${sample}/assembly/flye" }, mode: 'copy', saveAs: { it.endsWith('.version') ? null : it })

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("assembly.fasta"), emit: assembly
    tuple val(sample), path("assembly_info.txt"), emit: info
    tuple val(sample), path("assembly_graph.gfa"), emit: gfa
    tuple val(sample), path("flye_reads.tsv"), emit: reads_used
    tuple val(sample), path("copies.tsv"), emit: copies
    tuple val(sample), path("flye.version"), emit: version

    script:
    """
    run_flye() {
        rm -rf flye_run
        flye ${params.flye_mode} "\$1" --out-dir flye_run --threads ${task.cpus}
    }
    # Flye's own depth estimate from read overlaps ("Overlap-based coverage: 583")
    overlap_cov() {
        grep -o 'Overlap-based coverage: [0-9]*' flye_run/flye.log | tail -n 1 | grep -o '[0-9]*\$' || true
    }

    READS_USED=all
    FRACTION=1
    FIRST_COV=''
    if ! run_flye ${reads}; then
        COV=\$(overlap_cov)
        if [ "${params.enable_downsampling}" = true ] \\
            && grep -q 'No disjointigs were assembled' flye_run/flye.log \\
            && [ -n "\$COV" ] && [ "\$COV" -gt ${params.target_coverage} ]; then
            FRACTION=\$(awk -v t=${params.target_coverage} -v c="\$COV" 'BEGIN { printf "%.4f", t / c }')
            echo "Flye assembled no disjointigs at \${COV}x; retrying on a seeded \$FRACTION subset of the reads" >&2
            FIRST_COV=\$COV
            # 4-line fastq records (chopper/Rasusa output); seed 1 so a rerun
            # picks the same reads
            python3 - ${reads} "\$FRACTION" <<'PY'
import gzip, random, sys
src, frac = sys.argv[1], float(sys.argv[2])
rng = random.Random(1)
opener = gzip.open if src.endswith('.gz') else open
with opener(src, 'rb') as fin, gzip.open('flye_subset.fastq.gz', 'wb', compresslevel=1) as fout:
    while True:
        rec = [fin.readline() for _ in range(4)]
        if not rec[0]:
            break
        if rng.random() < frac:
            fout.writelines(rec)
PY
            READS_USED=subset
            run_flye flye_subset.fastq.gz
        else
            exit 1
        fi
    fi
    mv flye_run/assembly_info.txt flye_run/assembly_graph.gfa .

    # a contig that is two or more copies of a plasmid -> one copy (see above);
    # 20 Mb is above any bacterial genome, so no contig is skipped
    collapse_copies.py --sample "${sample}" --in flye_run/assembly.fasta \\
        --out assembly.fasta --table copies.tsv --max-length 20000000

    # first_try_coverage = Flye's overlap coverage on the reads it was given,
    # assembled_coverage = the same on the reads it assembled
    {
        printf 'sample\\treads_used\\tsubset_fraction\\tfirst_try_coverage\\tassembled_coverage\\n'
        printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "${sample}" "\$READS_USED" "\$FRACTION" \\
            "\${FIRST_COV:-\$(overlap_cov)}" "\$(overlap_cov)"
    } > flye_reads.tsv

    # "2.9.6-b1802" -> "2.9.6"
    flye --version | cut -d- -f1 > flye.version
    """
}
