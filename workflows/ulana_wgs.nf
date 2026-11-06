include { MERGE_FASTQ      } from '../modules/merge_fastq.nf'
include { CHOPPER          } from '../modules/chopper.nf'
include { DOWNSAMPLE_PLAN; DOWNSAMPLE } from '../modules/downsample.nf'
include { READ_STATS; READ_STATS_REPORT } from '../modules/read_stats.nf'
include { FLYE_ASSEMBLY    } from '../modules/flye.nf'
include { BANDAGE_IMAGE    } from '../modules/bandage_image.nf'
include { MEDAKA_POLISH    } from '../modules/medaka.nf'
include { PROKKA           } from '../modules/prokka.nf'
include { CHECKM_DB        } from '../modules/checkm_db.nf'
include { CHECKM           } from '../modules/checkm.nf'
include { IDENTIFIER_GENES } from '../modules/identifier_genes.nf'
include { AMRFINDER_DB     } from '../modules/amrfinder_db.nf'
include { AMRFINDER        } from '../modules/amrfinder.nf'
include { GENBANK_PREP; TABLE2ASN } from '../modules/genbank.nf'

// --genome_size as a whole number of bases: "92kb", "5.5 Mb", "4600000" ->
// 92000 / 5500000 / 4600000; blank -> 0 (estimate it). Stops the run at launch
// on anything else, rather than after merging and filtering every sample
def genomeSizeBp(value) {
    def text = (value == null || value == true) ? '' : value.toString().trim()
    if (!text) {
        return 0
    }
    def m = text =~ /^(\d+(?:\.\d+)?)\s*([kKmMgG]?)(?:[bB][pP]?)?$/
    if (!m.matches()) {
        error "--genome_size '${text}' isn't a size: give a number of bases, e.g. 92kb, 5.5mb or 4600000"
    }
    def scale = [k: 1e3, m: 1e6, g: 1e9][m.group(2).toLowerCase()] ?: 1
    def bp = (new BigDecimal(m.group(1)) * scale).toLong()
    if (bp < 1000 || bp > 50_000_000) {
        error "--genome_size '${text}' is ${bp} bases; expected 1 kb to 50 Mb${bp < 1000 ? ' (did you leave off the kb or mb?)' : ''}"
    }
    return bp
}

workflow ULANA_WGS {

    take:
    reads_ch   // tuple(sample, fastq)
    meta_ch    // tuple(sample, [organism, strain, locus_tag]), read only with --enable_genbank

    main:
    def genome_size = genomeSizeBp(params.genome_size)
    if (!(params.target_coverage instanceof Integer) || params.target_coverage < 1) {
        error "--target_coverage must be a whole number of 1 or more, not '${params.target_coverage}'"
    }

    // 1. merge multi-part fastq(.gz) files per sample
    MERGE_FASTQ(reads_ch)

    // 2. quality + length filter
    CHOPPER(MERGE_FASTQ.out.merged)

    // 2b. read length / Q-score summary, before vs. after filtering -- a
    // side branch, nothing downstream depends on it
    if (params.enable_read_stats) {
        READ_STATS(MERGE_FASTQ.out.merged.join(CHOPPER.out.filtered))
        READ_STATS_REPORT(
            READ_STATS.out.stats.map { sample, stats, hist -> stats }.collect(),
            READ_STATS.out.stats.map { sample, stats, hist -> hist }.collect()
        )
    }

    // 2c. cap the depth at --target_coverage before assembly (modules/downsample.nf);
    // a sample that doesn't reach it keeps every filtered read. Flye and Medaka
    // both get these reads -- polishing at thousands of x is slow and adds nothing
    if (params.enable_downsampling) {
        DOWNSAMPLE_PLAN(CHOPPER.out.filtered, genome_size)
        def planned = CHOPPER.out.filtered
            .join(DOWNSAMPLE_PLAN.out.plan)
            .branch { sample, reads, decision, target_bases ->
                cut: decision == 'downsample'
                keep: true
            }
        DOWNSAMPLE(planned.cut.map { sample, reads, decision, target_bases -> tuple(sample, reads, target_bases as long) })
        assembly_reads = planned.keep
            .map { sample, reads, decision, target_bases -> tuple(sample, reads) }
            .mix(DOWNSAMPLE.out.reads)
    } else {
        assembly_reads = CHOPPER.out.filtered
    }

    // 3. assemble
    FLYE_ASSEMBLY(assembly_reads)

    // 3b. render the assembly graph -- optional, visualizes the Flye GFA
    // directly (unaffected by Medaka, which only polishes the sequence, not
    // the graph, so this always runs off FLYE_ASSEMBLY's output regardless
    // of --enable_medaka)
    if (params.enable_bandage) {
        BANDAGE_IMAGE(FLYE_ASSEMBLY.out.gfa)
    }

    // 4. polish -- optional; downstream steps consume whichever of these is current
    if (params.enable_medaka) {
        MEDAKA_POLISH(assembly_reads.join(FLYE_ASSEMBLY.out.assembly))
        assembly_ch = MEDAKA_POLISH.out.polished
    } else {
        assembly_ch = FLYE_ASSEMBLY.out.assembly
    }

    // 5. annotate -- optional; id_genes and amrfinder both depend on this
    if (params.enable_prokka) {
        PROKKA(assembly_ch.join(meta_ch))
    }

    // 6. QC -- optional
    if (params.enable_checkm) {
        CHECKM_DB()
        CHECKM(assembly_ch, CHECKM_DB.out.db.first())
    }

    // 7. marker gene extraction -- optional, requires prokka
    if (params.enable_id_genes) {
        if (!params.enable_prokka) {
            error "--enable_id_genes requires --enable_prokka (identifier genes are extracted from the Prokka .ffn output)"
        }
        IDENTIFIER_GENES(PROKKA.out.ffn)
    }

    // 8. AMR detection -- optional, requires prokka
    if (params.enable_amrfinder) {
        if (!params.enable_prokka) {
            error "--enable_amrfinder requires --enable_prokka (AMRFinder needs Prokka's .faa/.gff output)"
        }
        AMRFINDER_DB()
        AMRFINDER(
            PROKKA.out.faa.join(PROKKA.out.gff),
            AMRFINDER_DB.out.db.first()
        )
    }

    // 9. GenBank submission files -- optional, requires prokka
    if (params.enable_genbank) {
        if (!params.enable_prokka) {
            error "--enable_genbank requires --enable_prokka (the GenBank feature table is Prokka's annotation)"
        }
        // tool versions for the assembly comment, e.g. "Flye v. 2.9.6; Medaka v. 2.2.1"
        def method_ch = FLYE_ASSEMBLY.out.version.map { sample, f -> tuple(sample, "Flye v. ${f.text.trim()}") }
        if (params.enable_medaka) {
            method_ch = method_ch
                .join(MEDAKA_POLISH.out.version)
                .map { sample, method, f -> tuple(sample, "${method}; Medaka v. ${f.text.trim()}") }
        }
        def annotation_ch = PROKKA.out.version.map { sample, f -> tuple(sample, "Prokka ${f.text.trim()}") }

        GENBANK_PREP(
            PROKKA.out.genbank
                .join(FLYE_ASSEMBLY.out.info)
                .join(meta_ch)
                .join(method_ch)
                .join(annotation_ch)
        )
        TABLE2ASN(
            GENBANK_PREP.out.files.join(annotation_ch),
            params.genbank_template ? file(params.genbank_template, checkIfExists: true) : []
        )
    }

    emit:
    assembly = assembly_ch
}
