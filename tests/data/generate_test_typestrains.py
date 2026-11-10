#!/usr/bin/env python3
"""Build the tiny type-strain snapshot used by -profile test (tests/data/typestrains/).

Not run as part of the pipeline -- a one-off generator, like
generate_synthetic_reads.py, re-run only if the fixture needs to change. The
real snapshot (scripts/build_typestrain_db.py) is ~51 GB, so the test one
holds four made-up "type strains": sample_a's synthetic genome, rebuilt from
the same seed, and three relatives made from it by random substitution
(about 94, 90 and 86% ANI), each its own species. That gives the species
check and the species tree their outcomes in CI:

  sample_a  assembles to the first genome: ~99-100% ANI hit plus the three
            relatives, so the species tree has neighbours to build with
  sample_b  a different random genome: no hit at all, so no tree

The genomes are kept (genomes/), since the species tree stages its
neighbours' sequences, not just their sketches. They're random sequence with
no marker genes, so a real tree run ends in "too_few_marker_genes"; CI runs
-stub-run, where the tree step writes a star tree instead.

Also written: lineage.tsv (every test taxid -> the Bacteria marker set) and
gtotree_hmms/Bacteria.hmm, three profiles cut from GToTree's real Bacteria
set (the full one is 5 MB) so GToTree has a real HMM file to search with.

The 16S BLAST db is a single made-up record so the files exist and the
process runs; the synthetic genomes have no 16S gene, so nothing is BLASTed.

Run from tests/data/ (needs docker for skani, makeblastdb and hmmfetch, and
network for GToTree's Bacteria.hmm):

    python3 generate_test_typestrains.py
"""
import csv
import gzip
import hashlib
import json
import os
import random
import shutil
import subprocess
import urllib.request

SKANI_IMAGE = 'quay.io/biocontainers/skani:0.4.0--hec9b1f2_0'
BLAST_IMAGE = 'quay.io/biocontainers/blast:2.15.0--pl5321h6f7f691_1'
GTOTREE_IMAGE = 'quay.io/biocontainers/gtotree:1.8.19--h9ee0642_0'
OUT = 'typestrains'
NAME = 'Testibacter syntheticus'
STRAIN = 'T1'
# accession, name, strain, taxid (= species_taxid), substitution rate from sample_a's genome
GENOMES = [
    ('TEST_000001.1', NAME, STRAIN, '1', 0.0),
    ('TEST_000002.1', 'Testibacter proximus', 'T2', '900002', 0.06),
    ('TEST_000003.1', 'Testibacter medius', 'T3', '900003', 0.10),
    ('TEST_000004.1', 'Testibacter remotus', 'T4', '900004', 0.14),
]
HMM_PROFILES = ['Ribosomal_L1', 'Ribosomal_L13', 'Ribosomal_S2']


def docker(image, *args, stdout=None):
    subprocess.run(['docker', 'run', '--rm', '--platform', 'linux/amd64',
                    '-v', f'{os.path.abspath(OUT)}:/w', '-w', '/w', image, *args], check=True, stdout=stdout)


def mutate(genome, rate, seed):
    rnd = random.Random(seed)
    return ''.join(rnd.choice([b for b in 'ACGT' if b != base]) if rnd.random() < rate else base
                   for base in genome)


def bacteria_hmm():
    """GToTree's Bacteria.hmm, by the link + md5 in its own hmm-sources-and-info.tsv."""
    info = subprocess.run(['docker', 'run', '--rm', '--platform', 'linux/amd64', GTOTREE_IMAGE, 'cat',
                           '/usr/local/share/gtotree/hmm_sets/hmm-sources-and-info.tsv'],
                          check=True, capture_output=True, text=True).stdout
    row = next(r for r in csv.DictReader(info.splitlines(), delimiter='\t') if r['file'] == 'Bacteria.hmm')
    with urllib.request.urlopen(row['link'], timeout=300) as r:
        data = r.read()
    if hashlib.md5(data).hexdigest() != row['md5']:
        raise SystemExit('Bacteria.hmm: md5 mismatch')
    return data


def main():
    # generate_synthetic_reads.py's first draws after seed 0 are sample_a's genome
    random.seed(0)
    genome = ''.join(random.choice('ACGT') for _ in range(30000))

    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(f'{OUT}/genomes')
    os.makedirs(f'{OUT}/16S')
    os.makedirs(f'{OUT}/gtotree_hmms')
    seqs = {}
    for i, (acc, name, strain, taxid, rate) in enumerate(GENOMES):
        seqs[acc] = mutate(genome, rate, seed=100 + i) if rate else genome
        # mtime 0 so a rerun gives byte-identical files
        with open(f'{OUT}/genomes/{acc}.fna.gz', 'wb') as raw, \
                gzip.GzipFile(fileobj=raw, mode='wb', mtime=0) as fh:
            fh.write(f'>{acc} {name} strain {strain}\n{seqs[acc]}\n'.encode())
    with open(f'{OUT}/genomes.txt', 'w') as fh:
        fh.write(''.join(f'genomes/{g[0]}.fna.gz\n' for g in GENOMES))
    docker(SKANI_IMAGE, 'skani', 'sketch', '-l', 'genomes.txt', '-o', 'sketch')
    os.remove(f'{OUT}/genomes.txt')

    rnd = random.Random(16)
    with open(f'{OUT}/16S/test_16S.fasta', 'w') as fh:
        fh.write(f'>NR_000001.1 {NAME} strain {STRAIN} 16S ribosomal RNA, complete sequence\n')
        fh.write(''.join(rnd.choice('ACGT') for _ in range(1500)) + '\n')
    docker(BLAST_IMAGE, 'makeblastdb', '-in', '16S/test_16S.fasta', '-dbtype', 'nucl',
           '-parse_seqids', '-taxid', '1', '-out', '16S/16S_ribosomal_RNA', '-title', 'test 16S')
    os.remove(f'{OUT}/16S/test_16S.fasta')

    columns = ['accession', 'name', 'species', 'subspecies', 'strain', 'organism_name',
               'taxid', 'species_taxid', 'relation_to_type_material', 'assembly_level',
               'refseq_category', 'genome_size', 'gc_percent', 'seq_rel_date', 'file']
    with open(f'{OUT}/typestrains.tsv', 'w') as fh:
        fh.write('\t'.join(columns) + '\n')
        for acc, name, strain, taxid, _ in GENOMES:
            s = seqs[acc]
            gc = round(100 * (s.count('G') + s.count('C')) / len(s), 1)
            fh.write('\t'.join([acc, name, name, '', strain, f'{name} {strain}', taxid, taxid,
                                'assembly from type material', 'Complete Genome', 'na', str(len(s)),
                                str(gc), '2026-11-01', f'genomes/{acc}.fna.gz']) + '\n')

    # same columns as scripts/build_typestrain_lineage.py
    with open(f'{OUT}/lineage.tsv', 'w') as fh:
        fh.write('taxid\tcurrent_taxid\tspecies_taxid\tspecies\tgenus\tfamily\torder\tclass\tphylum\tdomain\tgtotree_sets\n')
        for acc, name, strain, taxid, _ in GENOMES:
            fh.write('\t'.join([taxid, taxid, taxid, name, 'Testibacter', 'Testibacteraceae', '', '', '',
                                'Bacteria', 'Bacteria']) + '\n')

    with open(f'{OUT}/gtotree_hmms/Bacteria_full.hmm', 'wb') as fh:
        fh.write(bacteria_hmm())
    with open(f'{OUT}/gtotree_hmms/Bacteria.hmm', 'w') as fh:
        for profile in HMM_PROFILES:
            docker(GTOTREE_IMAGE, 'hmmfetch', 'gtotree_hmms/Bacteria_full.hmm', profile, stdout=fh)
    os.remove(f'{OUT}/gtotree_hmms/Bacteria_full.hmm')

    with open(f'{OUT}/build_info.json', 'w') as fh:
        json.dump({
            'snapshot': 'test',
            'built': None,
            'source': 'tests/data/generate_test_typestrains.py',
            'genomes': len(GENOMES),
            'species': len(GENOMES),
            'failed_downloads': 0,
            'subset': True,
            'skani_image': SKANI_IMAGE,
        }, fh, indent=2)
        fh.write('\n')


if __name__ == '__main__':
    main()
