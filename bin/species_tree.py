#!/usr/bin/env python3
"""
Builds the species tree: the assembly and its type-strain neighbours
(species_tree_pick.py) on a marker-gene protein tree with support values.

  GToTree    finds the marker set's single-copy genes in every genome
             (Prodigal + HMMER), aligns each gene and joins them into one
             protein alignment. A genome with under half the markers is
             dropped (GToTree's default), which is how an assembly too
             fragmented or too unusual for the set shows up
  IQ-TREE    maximum-likelihood tree, ultrafast bootstrap (1,000, with -bnni,
             which keeps UFBoot from overstating support when the model fits
             imperfectly). The model is chosen by BIC from a short list
             (MODEL_ARGS): in the hand tests (2026-11-10) it picked the same
             model and tree as the full ModelFinder on the distant case in
             one seventh of the time, and the model never changed a tree's
             shape, only support at weak forks

Writes:
  tree.nwk         the tree, leaves labelled as in neighbours.tsv ("query" =
                   the assembly), UFBoot (0-100) at each fork; unrooted
  tree_report.txt  IQ-TREE's report (model, likelihood, alignment summary)
  tree_info.json   status built / too_few_marker_genes, marker set, genes in
                   the set and found in the assembly, alignment size, model,
                   genomes GToTree dropped, tool versions

  species_tree.py --query assembly.fasta --genomes genomes/*.fna.gz \\
      --hmm Bacteria.hmm --pick pick.json --cpus 8
"""
import argparse
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys

MODEL_ARGS = ['-m', 'MFP', '-mset', 'LG,WAG,JTT,Q.pfam', '-mfreq', 'FU,F', '-mrate', 'G4,R4']
QUERY = 'query'


def run(cmd, log):
    with open(log, 'w') as out:
        return subprocess.run(cmd, stdout=out, stderr=subprocess.STDOUT).returncode


def version(cmd, pattern):
    out = subprocess.run(cmd, capture_output=True, text=True)
    m = re.search(pattern, out.stdout + out.stderr)
    return m.group(1) if m else ''


def label_of(path):
    return re.sub(r'\.(fna|fa|fasta)(\.gz)?$', '', os.path.basename(path))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--query', required=True)
    ap.add_argument('--genomes', nargs='+', required=True)
    ap.add_argument('--hmm', required=True)
    ap.add_argument('--pick', required=True)
    ap.add_argument('--cpus', type=int, default=2)
    args = ap.parse_args()

    with open(args.pick) as f:
        info = json.load(f)
    info.pop('status', None)

    # Every genome goes to GToTree as a plain inputs/<label>.fna: GToTree
    # names a genome after its file, so the leaves come out as "query" and the
    # accessions with no -m label map. It also decides "gzipped?" with
    # `file`, which on Google Batch sees the staged inputs as symlinks into the
    # bucket mount, not gzip: it then skipped its own gunzip and named every
    # neighbour "<accession>.fna", so the map no longer matched (2026-11-10)
    os.makedirs('inputs', exist_ok=True)
    plain = []
    for src, label in [(args.query, QUERY)] + [(g, label_of(g)) for g in args.genomes]:
        dest = os.path.join('inputs', f'{label}.fna')
        with open(src, 'rb') as f:
            gz = f.read(2) == b'\x1f\x8b'
        with (gzip.open(src, 'rb') if gz else open(src, 'rb')) as fin, open(dest, 'wb') as fout:
            shutil.copyfileobj(fin, fout)
        plain.append(dest)
    with open('genomes.list', 'w') as f:
        f.write(''.join(p + '\n' for p in plain))
    with open(args.hmm) as f:
        genes_in_set = sum(1 for line in f if line.startswith('NAME '))

    run(['GToTree', '-f', 'genomes.list', '-H', args.hmm, '-N',
         '-j', str(args.cpus), '-n', '1', '-o', 'gtotree'], 'gtotree.log')

    # GToTree calls it Aligned_SCGs.faa, or Aligned_SCGs_mod_names.faa when
    # given a label map (-m, not used any more)
    aln = next((p for p in ('gtotree/Aligned_SCGs.faa', 'gtotree/Aligned_SCGs_mod_names.faa')
                if os.path.exists(p)), 'gtotree/Aligned_SCGs.faa')
    labels = []
    if os.path.exists(aln):
        with open(aln) as f:
            labels = [line[1:].strip() for line in f if line.startswith('>')]

    # markers found in the assembly (SCG_hit_counts.tsv: one row per genome,
    # a hit count per gene; GToTree uses only genes hit exactly once)
    found = None
    hits = 'gtotree/SCG_hit_counts.tsv'
    if os.path.exists(hits):
        with open(hits) as f:
            for r in csv.reader(f, delimiter='\t'):
                if r and r[0] == QUERY:
                    found = sum(1 for x in r[1:] if x.strip() == '1')

    info.update({
        'genes_in_set': genes_in_set,
        'genes_found_in_query': found,
        'dropped': sorted({QUERY, *map(label_of, args.genomes)} - set(labels)),
        'gtotree_version': version(['GToTree', '-v'], r'v?(\d+\.\d+\.\d+)'),
    })

    if QUERY not in labels or len(labels) < 4:
        info['status'] = 'too_few_marker_genes'
        with open('tree_info.json', 'w') as f:
            json.dump(info, f, indent=2)
        print(f"no tree: {'the assembly' if QUERY not in labels else 'too few genomes'} kept "
              f"after the marker search ({len(labels)} genomes in the alignment)", file=sys.stderr)
        return

    code = run(['iqtree3', '-s', aln, *MODEL_ARGS, '-B', '1000', '-bnni',
                '-T', str(args.cpus), '--seed', '1', '--prefix', 'iqtree'], 'iqtree.log')
    if code != 0:
        sys.exit(f'IQ-TREE failed (exit {code}); see iqtree.log')

    with open('iqtree.iqtree') as f:
        report = f.read()

    def grab(pattern, cast=str):
        m = re.search(pattern, report)
        return cast(m.group(1)) if m else None

    os.replace('iqtree.treefile', 'tree.nwk')
    os.replace('iqtree.iqtree', 'tree_report.txt')
    info.update({
        'status': 'built',
        'genomes_in_tree': len(labels),
        'alignment_columns': grab(r'(\d+) amino-acid sites', int),
        'informative_sites': grab(r'Number of parsimony informative sites: (\d+)', int),
        'model': grab(r'Best-fit model according to BIC: (\S+)'),
        'support': 'UFBoot, 1000 replicates',
        'iqtree_version': version(['iqtree3', '--version'], r'version (\d+\.\d+\.\d+)'),
    })
    with open('tree_info.json', 'w') as f:
        json.dump(info, f, indent=2)


if __name__ == '__main__':
    main()
