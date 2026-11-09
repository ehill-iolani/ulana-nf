#!/usr/bin/env bash
#
# BLASTs each full-length 16S copy (IDENTIFIER_GENES' 16S_rRNA.fasta, one
# record per copy) against NCBI's 16S type-strain database and writes the
# species check's 16S_hits.tsv: the top hits per copy with identity and how
# much of the copy they cover. Supporting evidence only; the genome ANI is
# the verdict.
#
# Copies under 1,200 bp are skipped: a partial 16S at a contig end can't
# separate species, and NCBI's own entries are close to full length. Header
# only when no copy is long enough (or Prokka found none: the placeholder
# record has no sequence).
#
#   sixteen_s_hits.sh <16S_rRNA.fasta> <blast db prefix> <threads> <out.tsv>
set -euo pipefail

genes=$1
db=$2
threads=$3
out=$4
min_length=1200

awk -v min="$min_length" '
    function flush() { if (id != "" && length(seq) >= min) print ">" id "\n" seq }
    /^>/ { flush(); id = substr($1, 2); seq = ""; next }
    { seq = seq $0 }
    END { flush() }
' "$genes" > full_length.fasta

printf 'copy\tcopy_length\taccession\ttaxid\torganism\tidentity\talign_length\tquery_cover\tbitscore\ttitle\n' > "$out"

if [ -s full_length.fasta ]; then
    blastn -task megablast -query full_length.fasta -db "$db" \
        -max_target_seqs 10 -max_hsps 1 -num_threads "$threads" \
        -outfmt '6 qseqid qlen sacc staxid pident length qcovs bitscore stitle' \
        > blast.tsv
    # "Escherichia coli strain U 5/41 16S ribosomal RNA, complete sequence"
    # -> organism "Escherichia coli strain U 5/41"; the full title stays too
    awk -F '\t' -v OFS='\t' '{
        organism = $9
        sub(/ 16S ribosomal RNA.*$/, "", organism)
        print $1, $2, $3, $4, organism, $5, $6, $7, $8, $9
    }' blast.tsv >> "$out"
fi
