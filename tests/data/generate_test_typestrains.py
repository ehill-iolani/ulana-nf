#!/usr/bin/env python3
"""Build the tiny type-strain snapshot used by -profile test (tests/data/typestrains/).

Not run as part of the pipeline -- a one-off generator, like
generate_synthetic_reads.py, re-run only if the fixture needs to change. The
real snapshot (scripts/build_typestrain_db.py) is ~51 GB, so the test one
holds a single made-up "type strain": sample_a's synthetic genome, rebuilt
from the same seed. That gives both species-check outcomes in CI:

  sample_a  assembles to that genome, so skani reports a ~99-100% ANI hit
  sample_b  a different random genome, so no hit at all

The 16S BLAST db is a single made-up record so the files exist and the
process runs; the synthetic genomes have no 16S gene, so nothing is BLASTed.

Run from tests/data/ (needs docker for skani and makeblastdb):

    python3 generate_test_typestrains.py
"""
import gzip
import json
import os
import random
import shutil
import subprocess

SKANI_IMAGE = 'quay.io/biocontainers/skani:0.4.0--hec9b1f2_0'
BLAST_IMAGE = 'quay.io/biocontainers/blast:2.15.0--pl5321h6f7f691_1'
OUT = 'typestrains'
ACC = 'TEST_000001.1'
NAME = 'Testibacter syntheticus'
STRAIN = 'T1'


def docker(image, *args):
    subprocess.run(['docker', 'run', '--rm', '--platform', 'linux/amd64',
                    '-v', f'{os.path.abspath(OUT)}:/w', '-w', '/w', image, *args], check=True)


def main():
    # generate_synthetic_reads.py's first draws after seed 0 are sample_a's genome
    random.seed(0)
    genome = ''.join(random.choice('ACGT') for _ in range(30000))

    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(f'{OUT}/genomes')
    os.makedirs(f'{OUT}/16S')
    with gzip.open(f'{OUT}/genomes/{ACC}.fna.gz', 'wt') as fh:
        fh.write(f'>{ACC} {NAME} strain {STRAIN}\n{genome}\n')
    with open(f'{OUT}/genomes.txt', 'w') as fh:
        fh.write(f'genomes/{ACC}.fna.gz\n')
    docker(SKANI_IMAGE, 'skani', 'sketch', '-l', 'genomes.txt', '-o', 'sketch')
    # the sketch holds everything skani needs; the genome itself isn't read at run time
    shutil.rmtree(f'{OUT}/genomes')
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
    gc = round(100 * (genome.count('G') + genome.count('C')) / len(genome), 1)
    row = [ACC, NAME, NAME, '', STRAIN, f'{NAME} {STRAIN}', '1', '1',
           'assembly from type material', 'Complete Genome', 'na', str(len(genome)),
           str(gc), '2026-11-01', f'genomes/{ACC}.fna.gz']
    with open(f'{OUT}/typestrains.tsv', 'w') as fh:
        fh.write('\t'.join(columns) + '\n' + '\t'.join(row) + '\n')

    with open(f'{OUT}/build_info.json', 'w') as fh:
        json.dump({
            'snapshot': 'test',
            'built': None,
            'source': 'tests/data/generate_test_typestrains.py',
            'genomes': 1,
            'species': 1,
            'failed_downloads': 0,
            'subset': True,
            'skani_image': SKANI_IMAGE,
        }, fh, indent=2)
        fh.write('\n')


if __name__ == '__main__':
    main()
