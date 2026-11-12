#!/usr/bin/env python3
"""
Cuts a contig that is two or more back-to-back copies of the same genome
down to one copy.

Flye sometimes reports a phage as a circular contig about twice its length:
8NZ0, EGVZ and ZTYS (PRJNA1502461, validation 2026-11-11) came out at
~85 kb for ~43 kb genomes. The results then show a phage of the wrong size,
Pharokka annotates every gene twice, and Medaka polishes badly because each
read fits either copy equally well (ZTYS: 444 differences per 100 kb after
polishing). This runs between Flye and Medaka, so Medaka, CheckV and Pharokka
all see one copy.

How it decides: k-mers sampled every STEP bp along the contig are looked up
everywhere else in it, and each later match votes for an offset. A contig
made of copies has one offset (the genome length) that most k-mers agree on.
It's cut to that length only when
  - at least MIN_SHARED of the k-mers that could recur do recur at that
    offset (copies agree; the dimers above: 0.89-1.00, other contigs <= 0.3),
  - the repeat spans at least MIN_SPAN of a copy (contig >= 1.9 copies).
    A phage's terminal repeat (DTR, a few hundred bp at both ends, which is
    what lets CheckV call a genome complete) fails this and is left alone.
  - the copy is at least MIN_COPY bp.
The first copy is kept, cut where the contig's start comes round again: on a
circular contig that's the whole genome, and Medaka then polishes it with
every read. Where the copies disagree, the kept copy keeps its own version
(ZTYS: a ~170 bp stretch near 26.8 kb stays wrong after polishing, against
368 differences when the dimer was polished as it was).

Pure standard library. Usage:

  collapse_copies.py --sample S --in flye_raw.fasta --out assembly.fasta \\
      --table copies.tsv [--max-length 1000000]

copies.tsv has one row per contig that was cut (only the header when none):
sample, contig, flye_length, copy_length, copies, shared_kmers.
"""
import argparse
import collections

K = 21
STEP = 50
MIN_SHARED = 0.8
MIN_SPAN = 0.9
MIN_COPY = 1000
OFFSET_SLACK = 20  # bp; the copies' errors shift k-mers by a few bases
START_WINDOW = 2000  # bp at the contig's start that place the cut
START_SLACK = 0.05  # how far (fraction of a copy) the cut may sit from the majority offset


def read_fasta(path):
    records, name, header, parts = [], None, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    records.append((name, header, ''.join(parts)))
                header = line[1:]
                name = header.split()[0] if header.split() else ''
                parts = []
            elif line:
                parts.append(line.strip())
    if name is not None:
        records.append((name, header, ''.join(parts)))
    return records


def find_copy_length(seq):
    """(copy_length, shared) if seq is >= 2 copies of one sequence, else None."""
    seq = seq.upper()
    n = len(seq)
    if n < 2 * MIN_COPY:
        return None
    # index only the sampled k-mers, then scan every position against them
    sampled = collections.defaultdict(list)
    for i in range(0, n - K + 1, STEP):
        sampled[seq[i:i + K]].append(i)
    votes = collections.Counter()
    start_hits = []
    for j in range(n - K + 1):
        starts = sampled.get(seq[j:j + K])
        if starts:
            for i in starts:
                if j > i:
                    votes[j - i] += 1
                    if i < START_WINDOW:
                        start_hits.append((i, j))
    if not votes:
        return None
    best, support = None, 0
    for offset, _ in votes.most_common(20):
        near = sum(c for o, c in votes.items() if abs(o - offset) <= OFFSET_SLACK)
        if near > support:
            best, support = offset, near
    if best < MIN_COPY:
        return None
    eligible = len(range(0, n - best - K + 1, STEP))
    shared = support / eligible if eligible else 0.0
    if shared < MIN_SHARED or (n - best) < MIN_SPAN * best:
        return None
    # The cut goes where the contig's start comes round again. Copies that
    # differ by an indel give two offsets (ZTYS: one copy lacks ~355 bp at
    # 26.8 kb), and the majority one would leave the start repeated at the
    # end; the k-mers in the first START_WINDOW bp say where the start recurs.
    near_start = collections.Counter()
    for (i, j) in start_hits:
        offset = j - i
        if abs(offset - best) <= START_SLACK * best:
            near_start[offset] += 1
    if near_start:
        cut, cut_support = best, 0
        for offset, _ in near_start.most_common(20):
            near = sum(c for o, c in near_start.items() if abs(o - offset) <= OFFSET_SLACK)
            if near > cut_support:
                cut, cut_support = offset, near
        best = cut
    return best, shared


def wrap(seq, width=60):
    return '\n'.join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--sample', required=True)
    ap.add_argument('--in', dest='src', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--table', required=True)
    ap.add_argument('--max-length', type=int, default=1_000_000,
                    help='skip longer contigs (host DNA); default 1 Mb = two copies of the 500 kb size cap')
    args = ap.parse_args()

    rows = []
    with open(args.out, 'w') as out:
        for name, header, seq in read_fasta(args.src):
            found = find_copy_length(seq) if len(seq) <= args.max_length else None
            if found:
                length, shared = found
                rows.append((args.sample, name, len(seq), length,
                             f'{len(seq) / length:.2f}', f'{shared:.2f}'))
                seq = seq[:length]
            out.write(f'>{header}\n{wrap(seq)}\n')

    with open(args.table, 'w') as tsv:
        tsv.write('sample\tcontig\tflye_length\tcopy_length\tcopies\tshared_kmers\n')
        for row in rows:
            tsv.write('\t'.join(map(str, row)) + '\n')


if __name__ == '__main__':
    main()
