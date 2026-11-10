#!/usr/bin/env python3
"""
Picks the type-strain genomes the species tree is drawn with, and the GToTree
marker-gene set to draw it from. Standard library only.

Neighbours, up to --max (10):
  1. skani hits (genome_hits.tsv), best ANI first, one genome per species,
     except the closest species keeps every type strain it has (its
     subspecies), since a split there is what the tree is for
     (KEA-6: chelonae subsp. bovistauri/gwanakae vs the chelonae type);
  2. then, to fill up the rest, 16S hits (the wide search, 16S_wide.tsv)
     with query cover >= 90, best identity first, one genome per species,
     taken through lineage.tsv (16S taxids are strain-level) to a type genome
     of the same species; species with no genome, and archaea, are skipped.
     With no skani hit at all this is where every neighbour comes from (a
     possible new genus); otherwise it adds a few more distant relatives.

Marker set: the most specific GToTree set (lineage.tsv's gtotree_sets) that
the closest neighbour and every other neighbour share, so a neighbourhood
that spans classes falls back to the phylum set, then Bacteria.

Fewer than --min (3) neighbours: no tree (IQ-TREE needs four sequences).

Writes:
  neighbours.tsv  the query and each neighbour: label (accession; "query" for
                  the assembly), names, species_taxid, genus, family, where it
                  came from (skani / 16S) with its ANI or 16S identity
  genomes.txt     the neighbours' accessions, one per line (staged from the
                  snapshot's genomes/)
  pick.json       status (build / no_neighbours / too_few_neighbours),
                  marker set, counts
  tree_info.json  only when no tree will be built: the same status, so the
                  results page can say why

  species_tree_pick.py --sample S --genome-hits genome_hits.tsv \\
      [--sixteen-s 16S_wide.tsv] --typestrains typestrains.tsv --lineage lineage.tsv
"""
import argparse
import csv
import json
import sys

MIN_QUERY_COVER = 90
COLUMNS = ['label', 'accession', 'name', 'species', 'species_taxid', 'genus', 'family',
           'source', 'ani', 'af_query', 'identity_16S', 'query_cover_16S']


def read_tsv(path):
    with open(path, encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--sample', required=True)
    ap.add_argument('--genome-hits', required=True)
    ap.add_argument('--sixteen-s', help="wide 16S BLAST output (outfmt 6: qseqid qlen sacc staxid pident length qcovs bitscore stitle); absent or empty is fine")
    ap.add_argument('--typestrains', required=True)
    ap.add_argument('--lineage', required=True)
    ap.add_argument('--max', type=int, default=10)
    ap.add_argument('--min', type=int, default=3)
    args = ap.parse_args()

    lineage = {r['taxid']: r for r in read_tsv(args.lineage)}
    genomes_of = {}
    for r in read_tsv(args.typestrains):
        genomes_of.setdefault(r['species_taxid'], []).append(r)

    def lin(species_taxid):
        return lineage.get(species_taxid) or {}

    picked = []          # dicts in COLUMNS order
    species_seen = set()
    closest = None

    def add(genome, source, **numbers):
        sp = genome['species_taxid']
        L = lin(sp)
        picked.append({
            'label': genome['accession'], 'accession': genome['accession'], 'name': genome['name'],
            'species': L.get('species') or genome.get('species', ''), 'species_taxid': sp,
            'genus': L.get('genus', ''), 'family': L.get('family', ''), 'source': source,
            'ani': numbers.get('ani', ''), 'af_query': numbers.get('af_query', ''),
            'identity_16S': numbers.get('identity', ''), 'query_cover_16S': numbers.get('qcov', ''),
        })

    # 1. skani hits, best ANI first (genome_hits.tsv is sorted that way)
    typestrain = {g['accession']: g for gs in genomes_of.values() for g in gs}
    for h in read_tsv(args.genome_hits):
        if len(picked) >= args.max:
            break
        g = typestrain.get(h['accession'])
        if not g:   # sketch and table from different builds; nothing to stage
            continue
        sp = g['species_taxid']
        if closest is None:
            closest = sp
        if sp in species_seen and sp != closest:
            continue
        species_seen.add(sp)
        add(g, 'skani', ani=h['ani'], af_query=h['af_query'])

    # 2. 16S hits to fill up the rest
    rows = []
    if args.sixteen_s:
        try:
            with open(args.sixteen_s, encoding='utf-8') as f:
                for r in csv.reader(f, delimiter='\t'):
                    if len(r) >= 7 and float(r[6]) >= MIN_QUERY_COVER:
                        rows.append((float(r[4]), float(r[6]), r[3]))
        except FileNotFoundError:
            pass
    rows.sort(key=lambda r: -r[0])
    for identity, qcov, taxid in rows:
        if len(picked) >= args.max:
            break
        L = lineage.get(taxid)
        sp = L and L['species_taxid']
        if not sp or sp in species_seen:
            continue
        species_seen.add(sp)
        if L['domain'] != 'Bacteria' or sp not in genomes_of:
            continue
        candidates = genomes_of[sp]
        if closest is None:
            closest = sp      # the closest species keeps all its type strains
        else:
            candidates = [next((g for g in candidates if g['taxid'] == taxid), candidates[0])]
        for g in candidates:
            if len(picked) < args.max:
                add(g, '16S', identity=f'{identity:g}', qcov=f'{qcov:g}')

    # marker set shared by all
    gene_set = 'Bacteria'
    if picked:
        lists = [lin(p['species_taxid']).get('gtotree_sets', 'Bacteria').split(',') for p in picked]
        gene_set = next((s for s in lists[0] if all(s in other for other in lists[1:])), 'Bacteria')

    status = ('build' if len(picked) >= args.min
              else 'no_neighbours' if not picked else 'too_few_neighbours')
    info = {
        'status': status,
        'neighbours': len(picked),
        'from_skani': sum(p['source'] == 'skani' for p in picked),
        'from_16S': sum(p['source'] == '16S' for p in picked),
        'gtotree_set': gene_set,
        'min_neighbours': args.min,
        'max_neighbours': args.max,
    }

    with open('neighbours.tsv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, COLUMNS, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerow({c: '' for c in COLUMNS} | {'label': 'query', 'name': args.sample})
        w.writerows(picked)
    with open('genomes.txt', 'w', encoding='utf-8') as f:
        f.write(''.join(p['accession'] + '\n' for p in picked))
    with open('pick.json', 'w', encoding='utf-8') as f:
        json.dump(info, f, indent=2)
    if status != 'build':
        with open('tree_info.json', 'w', encoding='utf-8') as f:
            json.dump(info, f, indent=2)
    print(f"{args.sample}: {status}, {len(picked)} neighbours ({info['from_skani']} skani, "
          f"{info['from_16S']} 16S), marker set {gene_set}", file=sys.stderr)


if __name__ == '__main__':
    main()
