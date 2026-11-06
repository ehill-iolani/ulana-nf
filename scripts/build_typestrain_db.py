#!/usr/bin/env python3
"""
Builds the type-strain reference set for the species check: every bacterial
type-strain genome in NCBI RefSeq (one per species or subspecies), sketched
with skani, plus NCBI's 16S type-strain BLAST database. Run once per snapshot
(about yearly); runs read the result from GCS and never touch NCBI.

Output (in --out, and under --dest/<snapshot>/ when given):

  typestrains.tsv     one row per genome kept: accession, names, strain, taxids,
                      type relation, assembly level, size, GC. The pipeline joins
                      skani hits to this by accession (skani reports the file path)
  excluded.tsv        every type-material assembly left out, and why
  genomes/            <accession>.fna.gz, the genomes themselves (the tree step
                      needs the neighbours' sequences, not just their sketches)
  sketch/             skani consolidated database, for `skani search`
  sketches/           the same as one .sketch per genome, for `skani triangle`
                      on a handful of them (triangle can't read sketch/)
  16S/                NCBI 16S_ribosomal_RNA BLAST database (type strains)
  assembly_summary.txt  NCBI's table as downloaded, for provenance
  build_info.json     snapshot name, source dates, counts, tool versions

Selection: rows of NCBI's RefSeq bacteria assembly_summary.txt with a
relation_to_type_material and version_status "latest". NCBI often lists the
same type strain more than once (several assemblies, or the strain under a
second name like "Pseudomonas protegens CHA0"), so rows are grouped by name
cut to "Genus species [subsp. x]" (strain, serovar and pathovar words dropped)
and one is kept per name. A species-level row joins its nominate subspecies'
group (P. chlororaphis with P. chlororaphis subsp. chlororaphis), since they
share one type strain. Preference within a group: not an engineered derivative
("..._operon_deletion"); assembly from type material > neotype or reftype >
heterotypic synonym > other; NCBI's reference genome; Complete > Chromosome >
Scaffold > Contig; then the newest.

Downloads are resumable: a genome already in genomes/ was checked when it
landed (full gzip read, FASTA header) and is skipped; a rerun fetches only
what's missing or failed. Sketches are rebuilt from scratch each run (minutes).

Where to run the full build: about 40 GB of genomes plus a few GB of
sketches, so not a laptop. A us-central1 VM next to the bucket, e.g.

  gcloud compute instances create typestrain-build --zone us-central1-a \\
      --machine-type e2-standard-8 --boot-disk-size 120GB \\
      --image-family debian-12 --image-project debian-cloud --scopes cloud-platform
  # on the VM: install docker (or skani), copy this script over, then
  python3 build_typestrain_db.py --out /data/typestrains \\
      --dest gs://nextflow-jobs/ulana/db_cache/typestrains
  # delete the VM afterwards

Smaller builds, for testing or the -profile test fixture:

  build_typestrain_db.py --out ts --select-only           # just typestrains.tsv + excluded.tsv
  build_typestrain_db.py --out ts --accessions ids.txt    # only these accessions
  build_typestrain_db.py --out ts --limit 50              # first 50 selected

skani: the `skani` on PATH if there is one, else the pipeline's own container
through docker (SKANI_IMAGE below; keep it in step with the skani module).
Pure standard library otherwise.
"""
import argparse
import concurrent.futures as cf
import csv
import datetime as dt
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.error
import urllib.request

SUMMARY_URL = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt'
SIXTEEN_S_URL = 'https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz'
SKANI_IMAGE = 'quay.io/biocontainers/skani:0.4.0--hec9b1f2_0'
USER_AGENT = 'ulana-nf-typestrain-build'

RELATION_RANK = {
    'assembly from type material': 0,
    'assembly designated as neotype': 1,
    'assembly designated as reftype': 1,
    'assembly from heterotypic synonym type material': 2,
}
LEVEL_RANK = {'Complete Genome': 0, 'Chromosome': 1, 'Scaffold': 2, 'Contig': 3}

# "Genus species" with an optional "subsp. x"; genus may be bracketed
# ("[Clostridium] scindens") or prefixed with Candidatus
NAME_RE = re.compile(
    r'^((?:Candidatus )?\[?[A-Z][a-z]+\]? [a-z][a-z-]+)(?: subsp\. ([a-z][a-z-]+))?')

ALTERED_RE = re.compile(r'delet|knock|mutant|::|Δ|derivative|cured', re.IGNORECASE)

TSV_COLUMNS = [
    'accession', 'name', 'species', 'subspecies', 'strain', 'organism_name',
    'taxid', 'species_taxid', 'relation_to_type_material', 'assembly_level',
    'refseq_category', 'genome_size', 'gc_percent', 'seq_rel_date', 'file',
]


def log(msg):
    print(f'[{dt.datetime.now():%H:%M:%S}] {msg}', file=sys.stderr, flush=True)


def fetch(url, dest, tries=4):
    """Download url to dest via a .part file; returns the Last-Modified header."""
    part = dest + '.part'
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': USER_AGENT})
            with urllib.request.urlopen(req, timeout=120) as res, open(part, 'wb') as out:
                shutil.copyfileobj(res, out, 1 << 20)
                modified = res.headers.get('Last-Modified')
            os.replace(part, dest)
            return modified
        except urllib.error.HTTPError as e:
            if e.code == 404 or attempt == tries - 1:
                raise
        except (urllib.error.URLError, OSError):
            if attempt == tries - 1:
                raise
        time.sleep(2 ** attempt)


def read_summary(path):
    """Rows of assembly_summary.txt as dicts keyed by its own header."""
    with open(path, newline='') as fh:
        header = None
        for line in fh:
            if line.startswith('#assembly_accession') or line.startswith('# assembly_accession'):
                header = line.lstrip('# ').rstrip('\n').split('\t')
                break
        if header is None:
            sys.exit(f'{path}: no "#assembly_accession" header line, is this assembly_summary.txt?')
        for line in fh:
            yield dict(zip(header, line.rstrip('\n').split('\t')))


def name_key(organism):
    m = NAME_RE.match(organism)
    if not m:
        return organism, organism, ''
    species, subsp = m.group(1), m.group(2) or ''
    return (f'{species} subsp. {subsp}' if subsp else species), species, subsp


def select(rows):
    """(kept, excluded): one row per species/subspecies name, best first."""
    by_name = {}
    excluded = []
    for r in rows:
        rel = r.get('relation_to_type_material', '')
        if not rel or rel == 'na':
            continue
        if r.get('version_status') != 'latest':
            excluded.append((r, 'not the latest version'))
            continue
        if not r.get('ftp_path') or r['ftp_path'] == 'na':
            excluded.append((r, 'no ftp_path'))
            continue
        name, species, subsp = name_key(r['organism_name'])
        r['_name'], r['_species'], r['_subsp'] = name, species, subsp
        by_name.setdefault(name, []).append(r)

    # The subspecies that repeats the species epithet (P. chlororaphis subsp.
    # chlororaphis) has the species' own type strain, so a species-level row
    # and that subspecies' rows are one strain: pool them under the subspecies
    for name in list(by_name):
        if name not in by_name:  # a species-level group already pooled into its subspecies
            continue
        rows_ = by_name[name]
        species, subsp = rows_[0]['_species'], rows_[0]['_subsp']
        if subsp and subsp == species.split()[-1] and species in by_name:
            for r in by_name.pop(species):
                r['_name'], r['_subsp'] = name, subsp
                rows_.append(r)

    def rank(r):
        return (
            # an engineered derivative of the type strain ("ATCC 13985_hcn_operon_deletion") goes last
            1 if ALTERED_RE.search(r.get('infraspecific_name', '')) else 0,
            RELATION_RANK.get(r['relation_to_type_material'], 3),
            0 if r.get('refseq_category') == 'reference genome' else 1,
            LEVEL_RANK.get(r.get('assembly_level'), 4),
            # newest first: negate YYYY/MM/DD as a number
            -int(re.sub(r'\D', '', r.get('seq_rel_date', '')) or 0),
            r['assembly_accession'],
        )

    kept = []
    for name in sorted(by_name):
        group = sorted(by_name[name], key=rank)
        kept.append(group[0])
        for r in group[1:]:
            excluded.append((r, f'same name as {group[0]["assembly_accession"]} ({name}), ranked lower'))
    return kept, excluded


def genome_url(r):
    path = r['ftp_path'].rstrip('/').replace('ftp://', 'https://')
    return f'{path}/{path.rsplit("/", 1)[1]}_genomic.fna.gz'


def check_fasta_gz(path):
    """Read the whole file: a truncated or non-FASTA download raises."""
    with gzip.open(path, 'rb') as fh:
        if fh.read(1) != b'>':
            raise ValueError('not a FASTA file')
        while fh.read(1 << 20):
            pass


def download_one(r, genomes_dir):
    acc = r['assembly_accession']
    dest = os.path.join(genomes_dir, f'{acc}.fna.gz')
    if os.path.exists(dest):
        return acc, None
    tmp = dest + '.download'
    try:
        fetch(genome_url(r), tmp)
        check_fasta_gz(tmp)
        os.replace(tmp, dest)
        return acc, None
    except Exception as e:  # noqa: BLE001 -- recorded per genome, run carries on
        for p in (tmp, tmp + '.part'):
            if os.path.exists(p):
                os.remove(p)
        return acc, f'{type(e).__name__}: {e}'


def skani_cmd(out_dir):
    """argv prefix that runs skani with out_dir as the working directory."""
    if shutil.which('skani'):
        return ['skani'], out_dir
    if not shutil.which('docker'):
        sys.exit('needs skani on PATH, or docker to run ' + SKANI_IMAGE)
    return ['docker', 'run', '--rm', '--platform', 'linux/amd64',
            '-u', f'{os.getuid()}:{os.getgid()}',
            '-v', f'{os.path.abspath(out_dir)}:/w', '-w', '/w', SKANI_IMAGE, 'skani'], None


def run_skani(out_dir, args, threads):
    prefix, cwd = skani_cmd(out_dir)
    cmd = prefix + args + (['-t', str(threads)] if threads else [])
    log('$ ' + ' '.join(cmd))
    res = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if res.returncode != 0:
        sys.exit(f'skani failed ({res.returncode}):\n{res.stderr[-2000:]}')
    return res.stdout


def sketch(out_dir, accessions, threads):
    listing = os.path.join(out_dir, 'genomes.txt')
    with open(listing, 'w') as fh:
        for acc in accessions:
            fh.write(f'genomes/{acc}.fna.gz\n')
    for target, extra in (('sketch', []), ('sketches', ['--separate-sketches'])):
        final = os.path.join(out_dir, target)
        tmp = final + '.new'
        shutil.rmtree(tmp, ignore_errors=True)
        run_skani(out_dir, ['sketch', '-l', 'genomes.txt', '-o', target + '.new'] + extra, threads)
        shutil.rmtree(final, ignore_errors=True)
        os.replace(tmp, final)
    return run_skani(out_dir, ['--version'], 0).strip()


def fetch_16s(out_dir):
    dest_dir = os.path.join(out_dir, '16S')
    tarball = os.path.join(out_dir, '16S_ribosomal_RNA.tar.gz')
    modified = fetch(SIXTEEN_S_URL, tarball)
    shutil.rmtree(dest_dir, ignore_errors=True)
    os.makedirs(dest_dir)
    with tarfile.open(tarball) as tf:
        try:
            tf.extractall(dest_dir, filter='data')
        except TypeError:  # Python before 3.11.4 / 3.10.12 (Debian 12 ships 3.11.2): no filter
            for m in tf.getmembers():
                if m.name.startswith('/') or '..' in m.name.split('/') or not (m.isfile() or m.isdir()):
                    sys.exit(f'unexpected entry in {SIXTEEN_S_URL}: {m.name}')
            tf.extractall(dest_dir)
    os.remove(tarball)
    return modified


def write_tsv(path, header, rows):
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(header)
        w.writerows(rows)


def kept_row(r):
    strain = r.get('infraspecific_name', '')
    strain = strain.split('=', 1)[1] if strain.startswith('strain=') else (strain if strain != 'na' else '')
    return [
        r['assembly_accession'], r['_name'], r['_species'], r['_subsp'], strain,
        r['organism_name'], r['taxid'], r['species_taxid'],
        r['relation_to_type_material'], r['assembly_level'], r.get('refseq_category', ''),
        r.get('genome_size', ''), r.get('gc_percent', ''), r.get('seq_rel_date', ''),
        f'genomes/{r["assembly_accession"]}.fna.gz',
    ]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog='See the top of this file for the layout and where to run it.')
    ap.add_argument('--out', required=True, help='working + output directory (reused by reruns)')
    ap.add_argument('--snapshot', default=dt.date.today().strftime('%Y-%m'),
                    help='snapshot name, the GCS subfolder under --dest (default: this month, %(default)s)')
    ap.add_argument('--dest', help='gs:// prefix to upload to as <dest>/<snapshot>/ (needs gcloud)')
    ap.add_argument('--summary', help='use this local assembly_summary.txt instead of downloading it')
    ap.add_argument('--accessions', help='file of accessions to keep (one per line), for small builds')
    ap.add_argument('--limit', type=int, help='keep only the first N selected genomes')
    ap.add_argument('--select-only', action='store_true', help='write typestrains.tsv + excluded.tsv and stop')
    ap.add_argument('--jobs', type=int, default=8, help='parallel downloads (default %(default)s; be kind to NCBI)')
    ap.add_argument('--threads', type=int, default=os.cpu_count() or 4, help='skani threads (default %(default)s)')
    ap.add_argument('--skip-16s', action='store_true', help="don't fetch the 16S BLAST database")
    a = ap.parse_args()

    os.makedirs(os.path.join(a.out, 'genomes'), exist_ok=True)
    summary = os.path.join(a.out, 'assembly_summary.txt')
    summary_modified = None
    if a.summary:
        if os.path.abspath(a.summary) != os.path.abspath(summary):
            shutil.copyfile(a.summary, summary)
    else:
        log(f'downloading {SUMMARY_URL}')
        summary_modified = fetch(SUMMARY_URL, summary)

    kept, excluded = select(read_summary(summary))
    log(f'{len(kept)} type-strain genomes selected, {len(excluded)} left out')
    if a.accessions:
        wanted = {l.strip() for l in open(a.accessions) if l.strip() and not l.startswith('#')}
        missing = wanted - {r['assembly_accession'] for r in kept}
        if missing:
            log(f'not among the selected genomes (see excluded.tsv): {", ".join(sorted(missing))}')
        kept = [r for r in kept if r['assembly_accession'] in wanted]
    if a.limit:
        kept = kept[:a.limit]

    write_tsv(os.path.join(a.out, 'excluded.tsv'),
              ['accession', 'organism_name', 'relation_to_type_material', 'assembly_level', 'reason'],
              [[r['assembly_accession'], r['organism_name'], r['relation_to_type_material'],
                r.get('assembly_level', ''), why] for r, why in excluded])

    if a.select_only:
        write_tsv(os.path.join(a.out, 'typestrains.tsv'), TSV_COLUMNS, [kept_row(r) for r in kept])
        log(f'wrote {len(kept)} rows to {a.out}/typestrains.tsv (--select-only, nothing downloaded)')
        return

    genomes_dir = os.path.join(a.out, 'genomes')
    todo = [r for r in kept if not os.path.exists(os.path.join(genomes_dir, f'{r["assembly_accession"]}.fna.gz'))]
    log(f'{len(kept) - len(todo)} genomes already here, downloading {len(todo)} with {a.jobs} at a time')
    failed = {}
    done = 0
    started = time.time()
    with cf.ThreadPoolExecutor(a.jobs) as pool:
        for acc, err in pool.map(lambda r: download_one(r, genomes_dir), todo):
            done += 1
            if err:
                failed[acc] = err
            if done % 500 == 0 or done == len(todo):
                rate = done / max(time.time() - started, 1e-9)
                log(f'  {done}/{len(todo)} ({len(failed)} failed, {rate:.1f}/s)')

    by_acc = {r['assembly_accession']: r for r in kept}
    got = [acc for acc in by_acc if acc not in failed]
    write_tsv(os.path.join(a.out, 'failed.tsv'), ['accession', 'error'], sorted(failed.items()))
    write_tsv(os.path.join(a.out, 'typestrains.tsv'), TSV_COLUMNS, [kept_row(by_acc[acc]) for acc in got])
    if failed:
        log(f'{len(failed)} downloads failed (failed.tsv); they are left out of this build, rerun to retry them')

    # a genome from an earlier, bigger build in the same --out isn't part of this one
    stray = sorted(set(f[:-len('.fna.gz')] for f in os.listdir(genomes_dir) if f.endswith('.fna.gz')) - set(got))
    if stray:
        log(f'{len(stray)} genomes in genomes/ aren\'t in this build and are left out of the sketches')

    log(f'sketching {len(got)} genomes')
    skani_version = sketch(a.out, got, a.threads)

    sixteen_s_modified = None
    if not a.skip_16s:
        log(f'downloading {SIXTEEN_S_URL}')
        sixteen_s_modified = fetch_16s(a.out)

    info = {
        'snapshot': a.snapshot,
        'built': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'),
        'source': SUMMARY_URL if not a.summary else f'local file {a.summary}',
        'assembly_summary_last_modified': summary_modified,
        'genomes': len(got),
        'species': len({by_acc[acc]['_species'] for acc in got}),
        'failed_downloads': len(failed),
        'subset': bool(a.accessions or a.limit),
        'skani': skani_version,
        'skani_image': None if shutil.which('skani') else SKANI_IMAGE,
        'sixteen_s_db': None if a.skip_16s else SIXTEEN_S_URL,
        'sixteen_s_last_modified': sixteen_s_modified,
    }
    with open(os.path.join(a.out, 'build_info.json'), 'w') as fh:
        json.dump(info, fh, indent=2)
        fh.write('\n')
    log(f'built {a.out}: {info["genomes"]} genomes, {info["species"]} species')

    if a.dest:
        dest = f'{a.dest.rstrip("/")}/{a.snapshot}'
        if not shutil.which('gcloud'):
            sys.exit(f'--dest needs gcloud; upload by hand with:\n  gcloud storage rsync -r {a.out} {dest}')
        # stray genomes from an earlier build stay local: only what this build used goes up
        for sub in ('sketch', 'sketches', '16S'):
            if os.path.isdir(os.path.join(a.out, sub)):
                subprocess.run(['gcloud', 'storage', 'rsync', '-r', '--delete-unmatched-destination-objects',
                                os.path.join(a.out, sub), f'{dest}/{sub}'], check=True)
        listing = os.path.join(a.out, 'upload_genomes.txt')
        with open(listing, 'w') as fh:
            for acc in got:
                fh.write(os.path.join(genomes_dir, f'{acc}.fna.gz') + '\n')
        with open(listing) as fh:
            subprocess.run(['gcloud', 'storage', 'cp', '-n', '-I', f'{dest}/genomes/'], stdin=fh, check=True)
        for f in ('typestrains.tsv', 'excluded.tsv', 'failed.tsv', 'build_info.json', 'assembly_summary.txt'):
            subprocess.run(['gcloud', 'storage', 'cp', os.path.join(a.out, f), f'{dest}/{f}'], check=True)
        log(f'uploaded to {dest}/')


if __name__ == '__main__':
    main()
