#!/usr/bin/env nextflow
/*
 * ULANA-HT-NF -- bacterial whole-genome assembly and characterization from
 * Oxford Nanopore long reads. Nextflow DSL2 reimplementation of ulana-ht
 * (https://github.com/ehill-iolani/ulana-ht), keeping the same steps/toggles
 * and following edna-ont-nf's conventions so both pipelines share one
 * samplesheet schema and one launch pattern for the shared submission hub.
 */

nextflow.enable.dsl = 2

include { ULANA_WGS   } from './workflows/ulana_wgs.nf'
include { MERGE_FASTQ } from './modules/merge_fastq.nf'
include { CHOPPER     } from './modules/chopper.nf'
include { READ_STATS; READ_STATS_REPORT } from './modules/read_stats.nf'

// ---- top-level params (override via -params-file or --flag) ----
// tool params (chopper_q, flye_mode, medaka_model, ...) and the enable_*
// toggles all live in nextflow.config, not here -- anything read directly
// inside an included module or workflow script must be set there to be
// reliably visible by the time that module/workflow script runs
params.input  = null   // path to samplesheet.csv (sample,fastq[,organism,strain,locus_tag]); fastq globs resolve relative to the launch dir, not the CSV's location
params.outdir = "results"
params.help   = false

def helpMessage() {
    log.info """
    ULANA-HT-NF -- bacterial WGS assembly and characterization (ONT)
    ------------------------------------------------------------------
    Usage:
      nextflow run main.nf --input samplesheet.csv --outdir results

    Required:
      --input       CSV: sample,fastq_path (fastq may be a single file or a glob), plus
                    organism,strain,locus_tag for --enable_genbank

    Key optional toggles (all default to the ulana-ht Snakemake defaults):
      --enable_read_stats read length/Q-score summary before vs. after filtering (default ${params.enable_read_stats})
      --enable_medaka     polish the Flye assembly with Medaka (default ${params.enable_medaka})
      --enable_bandage    render a PNG of the Flye assembly graph (default ${params.enable_bandage})
      --enable_prokka     annotate with Prokka (default ${params.enable_prokka})
      --enable_checkm     QC with CheckM (default ${params.enable_checkm})
      --enable_id_genes   extract 16S/dnaA/rpoB from Prokka output, requires --enable_prokka (default ${params.enable_id_genes})
      --enable_amrfinder  AMR detection with AMRFinderPlus, requires --enable_prokka (default ${params.enable_amrfinder})
      --enable_genbank    GenBank submission files (.sqn, feature table, validation report), requires --enable_prokka
                          and organism + locus_tag samplesheet columns (default ${params.enable_genbank})
      --genbank_template  NCBI submission template (.sbt) to build into the .sqn (optional)

    Key tool params:
      --chopper_q / --chopper_minlength   chopper quality/length filtering (default ${params.chopper_q} / ${params.chopper_minlength})
      --flye_mode                         Flye read-type flag (default ${params.flye_mode})
      --medaka_model                      Medaka model name (default ${params.medaka_model})

    See nextflow_schema.json for the full parameter list.
    """.stripIndent()
}

def samplesheetToChannel(path) {
    Channel
        .fromPath(path, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            // checkIfExists doesn't catch a glob that matches zero files, so check explicitly --
            // otherwise this silently produces an empty file list downstream instead of failing
            def fq = file(row.fastq, checkIfExists: true)
            def files = fq instanceof List ? fq : [fq]
            if (files.isEmpty()) {
                error "Samplesheet row '${row.sample}': no files matched '${row.fastq}' (glob paths resolve relative to the launch directory, not the samplesheet's location)"
            }
            tuple(row.sample, fq)
        }
}

// The optional organism,strain,locus_tag columns, as tuple(sample, meta map).
// They're only read with --enable_genbank, and then organism and locus_tag
// are required: NCBI rejects a record without them, so fail here rather than
// after hours of assembly
def samplesheetMeta(path) {
    // allowlist, since the values go into Prokka's command line and into the
    // [key=value] modifiers of the GenBank FASTA headers
    def safe = /^[A-Za-z0-9 .\-_\/():,+#]+$/
    Channel
        .fromPath(path, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                organism : (row.organism ?: '').trim(),
                strain   : (row.strain ?: '').trim(),
                locus_tag: (row.locus_tag ?: '').trim(),
            ]
            if (params.enable_genbank) {
                def where = "Samplesheet row '${row.sample}'"
                if (!meta.organism) {
                    error "${where}: --enable_genbank needs an organism column (e.g. 'Pseudomonas aeruginosa', or 'Pseudomonas sp.' if the species isn't known)"
                }
                if (meta.organism.split(/\s+/).size() < 2) {
                    error "${where}: organism '${meta.organism}' needs a genus and a species, e.g. '${meta.organism} sp.' if the species isn't known"
                }
                if (!(meta.organism ==~ safe) || (meta.strain && !(meta.strain ==~ safe))) {
                    error "${where}: organism and strain may only use letters, digits, spaces and . - _ / ( ) : , + #"
                }
                // NCBI names an unnamed species by its strain ("Pseudomonas sp. KEA-6")
                if (meta.organism.endsWith(' sp.') && !meta.strain) {
                    error "${where}: organism '${meta.organism}' needs a strain column, since NCBI names an unidentified species by its strain"
                }
                if (!meta.locus_tag) {
                    error "${where}: --enable_genbank needs a locus_tag column, the prefix NCBI assigns when you register the BioSample"
                }
                // NCBI's rule for locus tag prefixes
                if (!(meta.locus_tag ==~ /^[A-Za-z][A-Za-z0-9]{2,11}$/)) {
                    error "${where}: locus_tag '${meta.locus_tag}' isn't a valid NCBI prefix (3-12 letters and digits, starting with a letter, no underscore)"
                }
            }
            tuple(row.sample, meta)
        }
}

workflow {
    if (params.help || !params.input) {
        helpMessage()
        exit 0
    }

    ULANA_WGS(samplesheetToChannel(params.input), samplesheetMeta(params.input))
}

// dev/debug entry point -- run just the merge step in isolation, e.g.:
//   nextflow run main.nf -entry MERGE_ONLY --input samplesheet.csv -profile docker
workflow MERGE_ONLY {
    if (!params.input) {
        log.error "MERGE_ONLY requires --input samplesheet.csv"
        exit 1
    }
    MERGE_FASTQ(samplesheetToChannel(params.input))
}

// Read QC only: merge -> chopper -> before/after length & Q-score summary, no
// assembly. Handy for picking --chopper_q/--chopper_minlength before
// committing to a full run:
//   nextflow run main.nf -entry READ_QC_ONLY --input samplesheet.csv --chopper_q 12 -profile docker
workflow READ_QC_ONLY {
    if (!params.input) {
        log.error "READ_QC_ONLY requires --input samplesheet.csv"
        exit 1
    }
    MERGE_FASTQ(samplesheetToChannel(params.input))
    CHOPPER(MERGE_FASTQ.out.merged)
    READ_STATS(MERGE_FASTQ.out.merged.join(CHOPPER.out.filtered))
    READ_STATS_REPORT(
        READ_STATS.out.stats.map { sample, stats, hist -> stats }.collect(),
        READ_STATS.out.stats.map { sample, stats, hist -> hist }.collect()
    )
}
