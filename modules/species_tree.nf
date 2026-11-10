// Species tree: where the assembly sits among its closest named relatives.
// Runs after the species check and reads the same type-strain snapshot
// (--species_db), plus two files the tree adds to it: lineage.tsv (NCBI
// lineage and GToTree marker set per taxid, scripts/build_typestrain_lineage.py)
// and gtotree_hmms/ (the marker sets themselves, so a run never downloads them).
//
//   SPECIES_TREE_PICK  which type genomes to draw it with, and which marker set
//                      (bin/species_tree_pick.py: skani hits, then 16S hits)
//   SPECIES_TREE       GToTree marker-gene alignment + IQ-TREE with UFBoot
//                      (bin/species_tree.py)
//   SPECIES_TREE_ANI   pairwise ANI between the same genomes (skani triangle),
//                      for the page's ANI view and 95% species groups
//
// Outputs, in {sample}/species_check/tree/:
//   neighbours.tsv   the query and each neighbour, with where it came from
//   tree.nwk         marker-gene tree, UFBoot at forks ("query" = the assembly)
//   tree_report.txt  IQ-TREE's report
//   ani.tsv          pairwise ANI (pairs under ~80% are missing)
//   tree_info.json   status (built / no_neighbours / too_few_neighbours /
//                    too_few_marker_genes), marker set, model, sizes, versions
//
// Like the species check, no verdict here: the page draws the tree and says
// what it shows.

process SPECIES_TREE_PICK {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/pandas:2.2.1"
    publishDir(path: { "${params.outdir}/${sample}/species_check/tree" }, mode: 'copy',
               pattern: '{neighbours.tsv,tree_info.json}')

    input:
    tuple val(sample), path(genome_hits), path(sixteen_s)  // sixteen_s may be [] (no --enable_id_genes)
    path typestrains
    path lineage

    output:
    tuple val(sample), path("pick.json"), path("genomes.txt"), emit: plan
    tuple val(sample), path("neighbours.tsv"), emit: neighbours
    path "tree_info.json", optional: true   // only when no tree will be built

    script:
    def wide = sixteen_s ? "--sixteen-s ${sixteen_s}" : ''
    """
    species_tree_pick.py --sample '${sample}' --genome-hits ${genome_hits} ${wide} \\
        --typestrains ${typestrains} --lineage ${lineage}
    """
}

process SPECIES_TREE {
    tag "$sample"
    // GToTree runs Prodigal + HMMER per genome in parallel and IQ-TREE uses
    // every core. Hand tests (amd64 emulation on a laptop, 4 threads): 2-6 min
    // for GToTree, 1-9 min for IQ-TREE, the longest for a distant
    // neighbourhood (more variable sites)
    label 'process_medium'
    // GToTree's image carries IQ-TREE (3.1.2) and Python too
    container "${params.container_registry}/biocontainers/gtotree:1.8.19--h9ee0642_0"
    publishDir(path: { "${params.outdir}/${sample}/species_check/tree" }, mode: 'copy',
               pattern: '{tree.nwk,tree_report.txt,tree_info.json}')

    input:
    tuple val(sample), path(assembly), path(genomes, stageAs: 'genomes/*'), path(hmm), path(pick)

    output:
    tuple val(sample), path("tree.nwk"), emit: tree, optional: true
    path "tree_report.txt", optional: true
    path "tree_info.json"
    // not published; declared so they come back from a Batch VM to the work dir
    path "{gtotree,iqtree}.log", optional: true

    script:
    """
    species_tree.py --query ${assembly} --genomes ${genomes} --hmm ${hmm} \\
        --pick ${pick} --cpus ${task.cpus}
    """

    stub:
    // -stub-run (CI): the test genomes are random sequence with no marker
    // genes, so a real run can only end in too_few_marker_genes. This writes
    // a star tree over the real labels instead, so the wiring and the output
    // contract are what CI checks; run -profile test without -stub-run for the
    // real step
    def labels = ['query'] + (genomes instanceof List ? genomes : [genomes]).collect { g -> g.name.replaceAll(/^.*\//, '').replaceAll(/\.(fna|fa|fasta)(\.gz)?$/, '') }
    """
    echo '(${labels.collect { l -> "${l}:0.1" }.join(',')});' > tree.nwk
    echo 'stub run: no IQ-TREE report' > tree_report.txt
    python3 - <<'EOF'
import json
info = json.load(open('${pick}'))
info.update(status='built', stub=True, genomes_in_tree=${labels.size()}, model=None)
json.dump(info, open('tree_info.json', 'w'), indent=2)
EOF
    """
}

process SPECIES_TREE_ANI {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/skani:0.4.0--hec9b1f2_0"
    publishDir(path: { "${params.outdir}/${sample}/species_check/tree" }, mode: 'copy', pattern: 'ani.tsv')

    input:
    tuple val(sample), path(assembly), path(genomes, stageAs: 'genomes/*')

    output:
    tuple val(sample), path("ani.tsv"), emit: ani

    script:
    """
    species_tree_ani.sh ${assembly} ${task.cpus} ani.tsv ${genomes}
    """
}
