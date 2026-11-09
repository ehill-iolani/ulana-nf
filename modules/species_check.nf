// Species check: how close the assembly is to the nearest named species.
// The assembly is compared with every bacterial type-strain genome in NCBI
// RefSeq (skani ANI, the main evidence), and each full-length 16S copy is
// BLASTed against NCBI's 16S type-strain database (supporting evidence).
//
// Both read a type-strain snapshot built once by scripts/build_typestrain_db.py
// (--species_db: sketch/, typestrains.tsv, build_info.json, 16S/). Nothing is
// downloaded per run. Neither step decides anything: the results page turns
// the hits into "same species" / "borderline" / "possible new species" / "no
// close type strain", gated on CheckM, so the cutoffs can change without a
// pipeline release.
//
// Outputs, in {sample}/species_check/:
//   genome_hits.tsv      top type-strain hits by ANI, with names (bin/species_hits.sh)
//   16S_hits.tsv         top 16S hits per full-length copy (bin/sixteen_s_hits.sh)
//   typestrain_db.json   which snapshot was searched (its build_info.json)

process SPECIES_CHECK_ANI {
    tag "$sample"
    // process_low (4 vCPU / 16 GB on Batch) is plenty: in the Batch test a
    // 6 Mb genome took ~18 s and peaked at 1.5 GB, the marker index included
    label 'process_low'
    container "${params.container_registry}/biocontainers/skani:0.4.0--hec9b1f2_0"
    publishDir(path: { "${params.outdir}/${sample}/species_check" }, mode: 'copy')

    input:
    tuple val(sample), path(assembly)
    path sketch       // the snapshot's consolidated skani database (sketch/)
    path typestrains  // typestrains.tsv
    path build_info   // build_info.json

    output:
    tuple val(sample), path("genome_hits.tsv"), emit: hits
    path "typestrain_db.json"

    script:
    // skani's defaults, as in the leave-one-out validation the verdict cutoffs
    // came from. -n 20 rather than 10: a genome can sit at >= 96% to several
    // type strains at once (E. coli and the four Shigella), and the page lists
    // every one of them. The sketch is read in place, not copied: skani only
    // reads it, and only the candidate genomes' records past the marker screen
    """
    skani search -d ${sketch} -q ${assembly} -n 20 -t ${task.cpus} -o skani.tsv
    species_hits.sh skani.tsv ${typestrains} genome_hits.tsv
    cp ${build_info} typestrain_db.json
    """
}

process SPECIES_CHECK_16S {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/blast:2.15.0--pl5321h6f7f691_1"
    publishDir(path: { "${params.outdir}/${sample}/species_check" }, mode: 'copy')

    input:
    tuple val(sample), path(sixteen_s)
    path db_files     // the snapshot's 16S/16S_ribosomal_RNA.* (not its 300 MB taxdb: titles carry the names)

    output:
    tuple val(sample), path("16S_hits.tsv"), emit: hits

    script:
    """
    sixteen_s_hits.sh ${sixteen_s} 16S_ribosomal_RNA ${task.cpus} 16S_hits.tsv
    """
}
