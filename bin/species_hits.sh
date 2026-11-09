#!/usr/bin/env bash
#
# Turns skani's search output into the species check's genome_hits.tsv: one row
# per type-strain hit, named from the type-strain set's typestrains.tsv, best
# ANI first. skani names each hit by its file (genomes/<accession>.fna.gz),
# which is the join key. Header only when skani found nothing (it reports no
# hit below about 80% ANI).
#
# No verdict here: the results page decides "same species" / "possible new
# species" from these numbers, so the cutoffs can change without a release.
#
#   species_hits.sh <skani search .tsv> <typestrains.tsv> <out.tsv>
set -euo pipefail

skani=$1
meta=$2
out=$3

printf 'accession\tname\tspecies\tsubspecies\tstrain\ttaxid\tspecies_taxid\tani\taf_query\taf_ref\tassembly_level\trelation_to_type_material\tgenome_size\n' > "$out"

awk -F '\t' -v OFS='\t' '
    # columns found by name, so a reordered table still joins
    FNR == 1 { for (i = 1; i <= NF; i++) col[FILENAME, $i] = i; next }
    FILENAME == ARGV[1] {
        acc = $col[FILENAME, "accession"]
        info[acc] = $col[FILENAME, "name"] OFS $col[FILENAME, "species"] OFS $col[FILENAME, "subspecies"] OFS \
            $col[FILENAME, "strain"] OFS $col[FILENAME, "taxid"] OFS $col[FILENAME, "species_taxid"]
        extra[acc] = $col[FILENAME, "assembly_level"] OFS $col[FILENAME, "relation_to_type_material"] OFS \
            $col[FILENAME, "genome_size"]
        next
    }
    {
        acc = $col[FILENAME, "Ref_file"]
        sub(/^.*\//, "", acc)
        sub(/\.(fna|fa|fasta)(\.gz)?$/, "", acc)
        # a hit missing from the table means the sketch and the table are from
        # different builds; keep it, named as skani saw it
        named = (acc in info) ? info[acc] : $col[FILENAME, "Ref_name"] OFS "" OFS "" OFS "" OFS "" OFS ""
        more = (acc in extra) ? extra[acc] : "" OFS "" OFS ""
        print acc, named, $col[FILENAME, "ANI"], $col[FILENAME, "Align_fraction_query"], \
            $col[FILENAME, "Align_fraction_ref"], more
    }
' "$meta" "$skani" \
    | sort -t "$(printf '\t')" -k8,8gr -k9,9gr \
    >> "$out"
