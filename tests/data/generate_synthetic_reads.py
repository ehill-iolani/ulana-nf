#!/usr/bin/env python3
"""Regenerate the tiny synthetic ONT-style dataset used by -profile test.

Not run as part of the pipeline -- a one-off generator, re-run only if the
checked-in fixtures under tests/data/reads/ need to change. Produces two
samples, each simulated as noisy long-read coverage of its own small random
"genome" so Flye has something trivial but real to assemble.
"""
import gzip
import random

random.seed(0)

BASES = "ACGT"


def random_genome(length):
    return "".join(random.choice(BASES) for _ in range(length))


def mutate(seq, sub_rate=0.03, indel_rate=0.01):
    out = []
    for base in seq:
        r = random.random()
        if r < indel_rate:
            continue  # deletion
        if r < indel_rate * 2:
            out.append(random.choice(BASES))  # insertion
        out.append(random.choice(BASES) if random.random() < sub_rate else base)
    return "".join(out)


def simulate_reads(genome, coverage, read_len_range=(2000, 4000)):
    reads = []
    total_len = 0
    target = len(genome) * coverage
    while total_len < target:
        rl = random.randint(*read_len_range)
        start = random.randint(0, max(0, len(genome) - 1))
        # wrap around so reads can start near the end (circular-ish genome)
        raw = (genome + genome)[start:start + rl]
        reads.append(mutate(raw))
        total_len += rl
    return reads


def write_fastq_gz(path, reads, ids):
    with gzip.open(path, "wt") as fh:
        for read_id, seq in zip(ids, reads):
            qual = "I" * len(seq)
            fh.write(f"@{read_id}\n{seq}\n+\n{qual}\n")


def main():
    # genome sizes/read lengths tuned so Flye reliably assembles a single
    # contig in CI: Flye's --min-overlap floor is 1000bp (not configurable
    # below that), so reads need to be long enough (multi-kb) relative to a
    # tiny synthetic genome for neighboring reads to clear it
    for sample, genome_len in [("sample_a", 30000), ("sample_b", 25000)]:
        genome = random_genome(genome_len)
        reads = simulate_reads(genome, coverage=30)
        # IDs numbered across the whole sample (not restarted per part-file)
        # so they stay unique after MERGE_FASTQ concatenates part_0 + part_1
        # -- Flye rejects duplicate read IDs outright
        ids = [f"{sample}_{i}" for i in range(len(reads))]
        # split into two part-files per sample, like real ONT barcode output,
        # to exercise MERGE_FASTQ's multi-file-per-sample glob handling
        half = len(reads) // 2
        write_fastq_gz(f"reads/{sample}/part_0.fastq.gz", reads[:half], ids[:half])
        write_fastq_gz(f"reads/{sample}/part_1.fastq.gz", reads[half:], ids[half:])


if __name__ == "__main__":
    main()
