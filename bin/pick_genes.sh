#!/usr/bin/env bash
#
# Pulls every Prokka gene whose product matches a pattern out of a .ffn file,
# one FASTA record per copy (header kept as Prokka wrote it: ">LOCUS_TAG
# product"), longest copy first so a partial copy at a contig end never comes
# first. A genome with seven 16S copies gives seven records, not one record
# with all seven run together.
#
# When nothing matches, the output is the placeholder header alone, which the
# results page reads as "the pipeline found no such gene".
#
#   pick_genes.sh <ffn> <product regex, matched against the whole product> <placeholder header> <out.fasta>
set -euo pipefail

ffn=$1
pattern=$2
placeholder=$3
out=$4

awk -v pat="$pattern" '
    function flush() { if (header != "") print length(seq) "\t" header "\t" seq }
    /^>/ {
        flush()
        header = ""; seq = ""
        product = index($0, " ") ? substr($0, index($0, " ") + 1) : ""
        if (product ~ pat) header = $0
        next
    }
    header != "" { seq = seq $0 }
    END { flush() }
' "$ffn" \
    | sort -t "$(printf '\t')" -k1,1nr -s \
    | awk -F '\t' '{ print $2; for (i = 1; i <= length($3); i += 60) print substr($3, i, 60) }' \
    > "$out"

if [ ! -s "$out" ]; then
    echo ">$placeholder" > "$out"
fi
