#!/usr/bin/env bash
#
# Pairwise ANI between the assembly and its species-tree neighbours (skani
# triangle), for the tree's second view and its 95% species groups. skani
# reports nothing below about 80% ANI, so for distant neighbourhoods (the 16S
# fallback) most pairs are simply missing: the page treats a missing pair as
# "under ~80%".
#
# Writes ani.tsv: a, b, ani, af_a, af_b (af = % of that genome aligned), one
# row per pair skani reported, labelled as in neighbours.tsv ("query" = the
# assembly).
#
#   species_tree_ani.sh <assembly> <threads> <out.tsv> <genome>...
set -euo pipefail

query=$1
threads=$2
out=$3
shift 3

printf '%s\n' "$query" "$@" > genomes.list
skani triangle -E -l genomes.list -t "$threads" -o skani.tsv

printf 'a\tb\tani\taf_a\taf_b\n' > "$out"
# skani -E: Ref_file Query_file ANI Align_fraction_ref Align_fraction_query ...
awk -F '\t' -v OFS='\t' -v query="$query" '
    function label(path) {
        if (path == query) return "query"
        sub(/^.*\//, "", path)
        sub(/\.(fna|fa|fasta)(\.gz)?$/, "", path)
        return path
    }
    NR > 1 { print label($1), label($2), $3, $4, $5 }
' skani.tsv >> "$out"
