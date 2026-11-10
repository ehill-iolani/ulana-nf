#!/usr/bin/env python3
"""
Adds what the species tree needs to a type-strain snapshot (see
build_typestrain_db.py), without a rebuild:

  lineage.tsv     the NCBI lineage of every taxid the species check and the
                  tree can meet, and which GToTree marker-gene sets fit it,
                  from NCBI's taxdump
  gtotree_hmms/   those marker sets (<set>.hmm, ~160 MB), so a run never
                  fetches them (GToTree would download one from Zenodo per run)

Rows: every taxid and species_taxid in typestrains.tsv, plus every taxid in
the snapshot's 16S BLAST database (strain-level taxids; the tree's fallback
turns 16S hits into neighbour genomes through species_taxid).

Columns:
  taxid           as found in typestrains.tsv / the 16S db
  current_taxid   the same taxon now (differs when NCBI merged it)
  species_taxid   current species-rank ancestor ('' when there is none)
  species         its current scientific name, e.g. "Mycobacteroides chelonae"
  genus family order class phylum domain   current names ('' when missing)
  gtotree_sets    GToTree sets that fit, most specific first, comma-separated,
                  always ending in Bacteria (Archaea for archaea)

GToTree's set names predate NCBI's 2021 phylum renames, so they are matched by
taxid (Firmicutes = Bacillota 1239, Proteobacteria = Pseudomonadota 1224, ...).
Epsilonproteobacteria is now the phylum Campylobacterota, outside
Pseudomonadota, so it gets no Proteobacteria fallback.

The tree step picks the first set shared by the query's closest hit and all of
its neighbours (a neighbourhood that spans classes falls back to the phylum,
then Bacteria).

Usage (a few minutes; the taxdump is ~80 MB, the marker sets ~160 MB):

  build_typestrain_lineage.py --snapshot /data/typestrains --out lineage.tsv --hmms gtotree_hmms
  gcloud storage cp -r lineage.tsv gtotree_hmms gs://nextflow-jobs/ulana/db_cache/typestrains/2026-11/

--snapshot needs typestrains.tsv and 16S/ (or pass --typestrains and
--sixteen-s-taxids yourself). The 16S taxids come from blastdbcmd, on PATH or
through docker (BLAST_IMAGE, keep in step with modules/species_check.nf).
--taxdump takes an unpacked taxdump directory instead of downloading one.
--hmms fetches the marker sets by the links and md5s in GToTree's own
hmm-sources-and-info.tsv (read from GTOTREE_IMAGE, keep in step with
modules/species_tree.nf). Pure standard library otherwise.
"""
import argparse
import collections
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request

TAXDUMP_URL = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
BLAST_IMAGE = 'quay.io/biocontainers/blast:2.15.0--pl5321h6f7f691_1'
GTOTREE_IMAGE = 'quay.io/biocontainers/gtotree:1.8.19--h9ee0642_0'
USER_AGENT = 'ulana-nf-typestrain-build'

# GToTree 1.8.19 HMM set -> NCBI taxid it was built for (checked against the
# 2026-11-10 taxdump). Archaea and Universal sets are left out on purpose.
GTOTREE_SETS = {
    'Actinobacteria': '201174',         # phylum Actinomycetota
    'Alphaproteobacteria': '28211',
    'Bacteroidetes': '976',             # phylum Bacteroidota
    'Betaproteobacteria': '28216',
    'Chlamydiae': '204428',             # phylum Chlamydiota
    'Cyanobacteria': '1117',            # phylum Cyanobacteriota
    'Epsilonproteobacteria': '29547',   # now phylum Campylobacterota
    'Firmicutes': '1239',               # phylum Bacillota
    'Gammaproteobacteria': '1236',
    'Proteobacteria': '1224',           # phylum Pseudomonadota
    'Tenericutes': '544448',            # phylum Mycoplasmatota
    'Bacteria': '2',
    'Archaea': '2157',
}
RANKS = ['species', 'genus', 'family', 'order', 'class', 'phylum']
COLUMNS = ['taxid', 'current_taxid', 'species_taxid', 'species', *RANKS[1:], 'domain', 'gtotree_sets']


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def dmp(path):
    with open(path, encoding='utf-8') as f:
        for line in f:
            yield line.rstrip('\t|\n').split('\t|\t')


def load_taxdump(d):
    parent, rank, name = {}, {}, {}
    for r in dmp(os.path.join(d, 'nodes.dmp')):
        parent[r[0]] = r[1]
        rank[r[0]] = r[2]
    for r in dmp(os.path.join(d, 'names.dmp')):
        if r[3] == 'scientific name':
            name[r[0]] = r[1]
    merged = {r[0]: r[1] for r in dmp(os.path.join(d, 'merged.dmp'))}
    return parent, rank, name, merged


def fetch_taxdump(dest):
    tgz = os.path.join(dest, 'taxdump.tar.gz')
    log(f'downloading {TAXDUMP_URL}')
    req = urllib.request.Request(TAXDUMP_URL, headers={'User-Agent': USER_AGENT})
    with urllib.request.urlopen(req, timeout=300) as r, open(tgz, 'wb') as f:
        shutil.copyfileobj(r, f)
    with tarfile.open(tgz) as t:
        for m in ('nodes.dmp', 'names.dmp', 'merged.dmp'):
            t.extract(m, dest, filter='data')
    return dest


def sixteen_s_taxids(db_dir):
    """Taxids in the 16S BLAST db, via blastdbcmd on PATH or the pipeline's image."""
    db_dir = os.path.abspath(db_dir)
    args = ['blastdbcmd', '-db', '16S_ribosomal_RNA', '-entry', 'all', '-outfmt', '%T']
    if not shutil.which('blastdbcmd'):
        args = ['docker', 'run', '--rm', '--platform', 'linux/amd64', '-v', f'{db_dir}:/db', '-w', '/db',
                BLAST_IMAGE, *args]
    out = subprocess.run(args, cwd=db_dir, check=True, capture_output=True, text=True).stdout
    return {t.strip() for t in out.split() if t.strip()}


def fetch_hmms(dest):
    """Every GToTree set lineage.tsv can name (Archaea aside: the tree skips
    archaea), checked against GToTree's md5s. Ones already there and right are kept."""
    info = subprocess.run(['docker', 'run', '--rm', '--platform', 'linux/amd64', GTOTREE_IMAGE, 'cat',
                           '/usr/local/share/gtotree/hmm_sets/hmm-sources-and-info.tsv'],
                          check=True, capture_output=True, text=True).stdout
    rows = {r['file']: r for r in csv.DictReader(info.splitlines(), delimiter='\t')}
    os.makedirs(dest, exist_ok=True)
    for name in GTOTREE_SETS:
        if name == 'Archaea':
            continue
        row, path = rows[f'{name}.hmm'], os.path.join(dest, f'{name}.hmm')
        if os.path.exists(path) and md5(path) == row['md5']:
            continue
        log(f'downloading {name}.hmm')
        req = urllib.request.Request(row['link'], headers={'User-Agent': USER_AGENT})
        with urllib.request.urlopen(req, timeout=600) as r, open(path, 'wb') as f:
            shutil.copyfileobj(r, f)
        if md5(path) != row['md5']:
            os.remove(path)
            sys.exit(f'{name}.hmm: md5 mismatch, removed; rerun to retry')


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument('--snapshot', help='snapshot folder with typestrains.tsv and 16S/')
    ap.add_argument('--typestrains', help='typestrains.tsv (default: <snapshot>/typestrains.tsv)')
    ap.add_argument('--sixteen-s-taxids', help='file of 16S db taxids, one per line (default: read <snapshot>/16S)')
    ap.add_argument('--taxdump', help='unpacked taxdump directory (default: download it)')
    ap.add_argument('--out', required=True, help='lineage.tsv to write')
    ap.add_argument('--hmms', help='also fetch the GToTree marker sets into this folder (gtotree_hmms/)')
    args = ap.parse_args()
    if args.hmms:
        fetch_hmms(args.hmms)

    ts_path = args.typestrains or (args.snapshot and os.path.join(args.snapshot, 'typestrains.tsv'))
    if not ts_path:
        ap.error('give --snapshot or --typestrains')
    with open(ts_path, encoding='utf-8') as f:
        ts = list(csv.DictReader(f, delimiter='\t'))
    wanted = {r[k] for r in ts for k in ('taxid', 'species_taxid')}
    n_ts = len(wanted)

    if args.sixteen_s_taxids:
        with open(args.sixteen_s_taxids, encoding='utf-8') as f:
            s16 = {line.strip() for line in f if line.strip()}
    elif args.snapshot and os.path.isdir(os.path.join(args.snapshot, '16S')):
        s16 = sixteen_s_taxids(os.path.join(args.snapshot, '16S'))
    else:
        ap.error('no 16S taxids: give --snapshot with a 16S/ folder, or --sixteen-s-taxids')
    wanted |= s16

    with tempfile.TemporaryDirectory() as tmp:
        parent, rank, name, merged = load_taxdump(args.taxdump or fetch_taxdump(tmp))
    set_of = {taxid: s for s, taxid in GTOTREE_SETS.items()}

    rows, missing = [], []
    for taxid in sorted(wanted, key=int):
        cur = taxid if taxid in parent else merged.get(taxid)
        if not cur:
            missing.append(taxid)
            continue
        lineage, sets, t = {}, [], cur
        for _ in range(100):
            if rank.get(t) in RANKS and rank[t] not in lineage:
                lineage[rank[t]] = t
            if rank.get(t) in ('domain', 'superkingdom'):
                lineage['domain'] = t
            if t in set_of:
                sets.append(set_of[t])
            if t == '1' or t not in parent:
                break
            t = parent[t]
        rows.append({
            'taxid': taxid, 'current_taxid': cur,
            'species_taxid': lineage.get('species', ''),
            'species': name.get(lineage.get('species'), ''),
            **{r: name.get(lineage.get(r), '') for r in [*RANKS[1:], 'domain']},
            'gtotree_sets': ','.join(sets),
        })

    with open(args.out, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, COLUMNS, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)

    # summary: per type genome, the most specific set its species gets
    by_taxid = {r['taxid']: r for r in rows}
    first = collections.Counter(
        (by_taxid.get(r['species_taxid']) or by_taxid.get(r['taxid']) or {}).get('gtotree_sets', '?').split(',')[0]
        for r in ts)
    log(f'{len(rows)} taxids written to {args.out} ({n_ts} from typestrains.tsv, {len(s16)} from 16S)')
    log(f'merged taxids followed: {sum(r["taxid"] != r["current_taxid"] for r in rows)}; '
        f'not in the taxdump: {len(missing)} {missing[:10]}')
    log('type genomes by most specific GToTree set:')
    for s, n in first.most_common():
        log(f'  {s:24} {n:6}  {n / len(ts):.1%}')
    log(json.dumps({'no_species_rank': sum(not r['species_taxid'] for r in rows)}))


if __name__ == '__main__':
    main()
