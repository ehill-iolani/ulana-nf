include { MERGE_FASTQ      } from '../modules/merge_fastq.nf'
include { CHOPPER          } from '../modules/chopper.nf'
include { FLYE_ASSEMBLY    } from '../modules/flye.nf'
include { BANDAGE_IMAGE    } from '../modules/bandage_image.nf'
include { MEDAKA_POLISH    } from '../modules/medaka.nf'
include { PROKKA           } from '../modules/prokka.nf'
include { CHECKM_DB        } from '../modules/checkm_db.nf'
include { CHECKM           } from '../modules/checkm.nf'
include { IDENTIFIER_GENES } from '../modules/identifier_genes.nf'
include { AMRFINDER_DB     } from '../modules/amrfinder_db.nf'
include { AMRFINDER        } from '../modules/amrfinder.nf'

workflow ULANA_WGS {

    take:
    reads_ch   // tuple(sample, fastq)

    main:
    // 1. merge multi-part fastq(.gz) files per sample
    MERGE_FASTQ(reads_ch)

    // 2. quality + length filter
    CHOPPER(MERGE_FASTQ.out.merged)

    // 3. assemble
    FLYE_ASSEMBLY(CHOPPER.out.filtered)

    // 3b. render the assembly graph -- optional, visualizes the Flye GFA
    // directly (unaffected by Medaka, which only polishes the sequence, not
    // the graph, so this always runs off FLYE_ASSEMBLY's output regardless
    // of --enable_medaka)
    if (params.enable_bandage) {
        BANDAGE_IMAGE(FLYE_ASSEMBLY.out.gfa)
    }

    // 4. polish -- optional; downstream steps consume whichever of these is current
    if (params.enable_medaka) {
        MEDAKA_POLISH(CHOPPER.out.filtered.join(FLYE_ASSEMBLY.out.assembly))
        assembly_ch = MEDAKA_POLISH.out.polished
    } else {
        assembly_ch = FLYE_ASSEMBLY.out.assembly
    }

    // 5. annotate -- optional; id_genes and amrfinder both depend on this
    if (params.enable_prokka) {
        PROKKA(assembly_ch)
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

    emit:
    assembly = assembly_ch
}
