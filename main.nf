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
params.input  = null   // path to samplesheet.csv (sample,fastq); fastq globs resolve relative to the launch dir, not the CSV's location
params.outdir = "results"
params.help   = false

def helpMessage() {
    log.info """
    ULANA-HT-NF -- bacterial WGS assembly and characterization (ONT)
    ------------------------------------------------------------------
    Usage:
      nextflow run main.nf --input samplesheet.csv --outdir results

    Required:
      --input       CSV: sample,fastq_path (fastq may be a single file or a glob)

    Key optional toggles (all default to the ulana-ht Snakemake defaults):
      --enable_read_stats read length/Q-score summary before vs. after filtering (default ${params.enable_read_stats})
      --enable_medaka     polish the Flye assembly with Medaka (default ${params.enable_medaka})
      --enable_bandage    render a PNG of the Flye assembly graph (default ${params.enable_bandage})
      --enable_prokka     annotate with Prokka (default ${params.enable_prokka})
      --enable_checkm     QC with CheckM (default ${params.enable_checkm})
      --enable_id_genes   extract 16S/dnaA/rpoB from Prokka output, requires --enable_prokka (default ${params.enable_id_genes})
      --enable_amrfinder  AMR detection with AMRFinderPlus, requires --enable_prokka (default ${params.enable_amrfinder})

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

workflow {
    if (params.help || !params.input) {
        helpMessage()
        exit 0
    }

    ULANA_WGS(samplesheetToChannel(params.input))
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
