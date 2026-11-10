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
include { SPECIES_CHECK_ANI; SPECIES_CHECK_16S } from '../modules/species_check.nf'
include { SPECIES_TREE_PICK; SPECIES_TREE; SPECIES_TREE_ANI } from '../modules/species_tree.nf'

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

    // 10. species check -- optional, needs a type-strain snapshot
    // (--species_db, built by scripts/build_typestrain_db.py; the google_batch
    // profile points it at the project's copy). Without one the step is
    // skipped with a warning rather than failing, so a local run with no
    // 50 GB snapshot still works. The 16S half needs IDENTIFIER_GENES' 16S copies
    if (params.enable_species_check) {
        if (!params.species_db) {
            log.warn "Species check skipped: no --species_db given (a type-strain snapshot from scripts/build_typestrain_db.py)"
        } else {
            def db = params.species_db.toString().replaceAll('/+$', '')
            SPECIES_CHECK_ANI(
                assembly_ch,
                file("${db}/sketch", type: 'dir', checkIfExists: true),
                file("${db}/typestrains.tsv", checkIfExists: true),
                file("${db}/build_info.json", checkIfExists: true)
            )
            if (params.enable_id_genes) {
                SPECIES_CHECK_16S(
                    IDENTIFIER_GENES.out.genes.map { sample, sixteen_s, dnaA, rpoB -> tuple(sample, sixteen_s) },
                    files("${db}/16S/16S_ribosomal_RNA.*", checkIfExists: true)
                )
            }

            // 10b. species tree (modules/species_tree.nf) -- optional; needs the
            // snapshot's lineage.tsv + gtotree_hmms/ (scripts/build_typestrain_lineage.py),
            // so an older snapshot without them skips it with a warning. Without
            // --enable_id_genes the neighbours come from skani alone
            if (params.enable_species_tree) {
                def lineage = file("${db}/lineage.tsv")
                if (!lineage.exists()) {
                    log.warn "Species tree skipped: ${db} has no lineage.tsv (add it with scripts/build_typestrain_lineage.py)"
                } else {
                    def wide_ch = params.enable_id_genes
                        ? SPECIES_CHECK_16S.out.wide
                        : SPECIES_CHECK_ANI.out.hits.map { sample, hits -> tuple(sample, []) }
                    SPECIES_TREE_PICK(
                        SPECIES_CHECK_ANI.out.hits.join(wide_ch),
                        file("${db}/typestrains.tsv", checkIfExists: true),
                        lineage
                    )
                    // a sample with too few neighbours stops here (its
                    // tree_info.json says why); the rest stage their
                    // neighbours' genomes from the snapshot's genomes/
                    def tree_plan = SPECIES_TREE_PICK.out.plan
                        .map { sample, pick, genomes ->
                            tuple(sample, pick, new groovy.json.JsonSlurper().parseText(pick.text), genomes)
                        }
                        .filter { sample, pick, plan, genomes -> plan.status == 'build' }
                        .map { sample, pick, plan, genomes ->
                            def staged = genomes.readLines().findAll { acc -> acc.trim() }.collect { acc ->
                                file("${db}/genomes/${acc.trim()}.fna.gz", checkIfExists: true)
                            }
                            tuple(sample, staged, file("${db}/gtotree_hmms/${plan.gtotree_set}.hmm", checkIfExists: true), pick)
                        }
                    tree_input = assembly_ch.join(tree_plan)
                    SPECIES_TREE(tree_input)
                    SPECIES_TREE_ANI(tree_input.map { sample, contigs, staged, hmm, pick -> tuple(sample, contigs, staged) })
                }
            }
        }
    }

    emit:
    assembly = assembly_ch
}
