#!/usr/bin/env python3
"""
Turns Prokka's .fsa + .tbl into the files table2asn needs for a GenBank
genome submission, fixing what NCBI's validator rejects in Prokka's output:

  - every CDS/RNA gets a `gene` feature carrying its locus_tag (Prokka's
    --addgenes would do this, but it also adds gene rows to the .tsv and .gff
    the rest of the pipeline publishes, so it's done here instead)
  - every CDS gets a protein_id (gnl|ncbi|<locus_tag>)
  - "hypothetical protein" CDSs (incl. "putative protein" and the like,
    renamed to it) lose their EC number and gene name
  - gene names lose Prokka's _1, _2 ... duplicate suffixes
  - barrnap's "(partial)" rRNA fragments are dropped (NCBI flags them as too
    short with a non-standard name), as are Prokka's non-standard `score`
    qualifiers
  - contigs under 200 bp, which GenBank doesn't accept, are dropped with
    their features

and writes the FASTA headers' source modifiers (organism, strain, genetic
code, topology=circular for contigs Flye closed) plus the
Genome-Assembly-Data structured comment (assembly method, coverage,
sequencing technology) NCBI asks every genome submission for.

Pure standard library. Outputs <prefix>.fsa, <prefix>.tbl, <prefix>.cmt and
submission_info.txt (what was changed, and what to enter in NCBI's portal):

  genbank_prep.py --fsa prokka.fsa --tbl prokka.tbl --assembly-info assembly_info.txt \\
      --organism "Pseudomonas sp." --strain KEA-6 --locus-tag KEAT \\
      --assembly-method "Flye v. 2.9.6; Medaka v. 2.2.1" --annotation "Prokka 1.14.6" \\
      --out-prefix KEA-6
"""
import argparse
import re
import sys

MIN_CONTIG = 200
# features that get a `gene` feature of their own
GENE_PARENTS = {'CDS', 'rRNA', 'tRNA', 'tmRNA', 'ncRNA', 'misc_RNA'}
HYPOTHETICAL = 'hypothetical protein'
# products table2asn's own cleanup (-c f) renames to "hypothetical protein",
# after which an EC number or gene name on them fails validation; renamed
# here first so those get removed too
UNKNOWN_PRODUCTS = {'putative protein', 'protein', 'uncharacterized protein', 'unknown protein', HYPOTHETICAL}


def read_fasta(path):
    """[(id, defline, sequence)] in file order."""
    records = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip('\n')
            if line.startswith('>'):
                defline = line[1:]
                records.append([defline.split()[0], defline, []])
            elif records:
                records[-1][2].append(line.strip())
    return [(seq_id, defline, ''.join(parts)) for seq_id, defline, parts in records]


def read_tbl(path):
    """{contig: [feature]} where a feature is {'locs': [[start, stop]], 'key', 'quals': [[name, value]]}."""
    contigs = {}
    features = None
    with open(path) as handle:
        for line in handle:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>Feature'):
                features = contigs.setdefault(line.split(None, 1)[1].strip(), [])
                continue
            cols = line.split('\t')
            if cols[0]:
                # a location line: "start stop key" opens a feature, "start
                # stop" adds an interval to the current one
                if len(cols) >= 3 and cols[2]:
                    features.append({'locs': [cols[:2]], 'key': cols[2], 'quals': []})
                else:
                    features[-1]['locs'].append(cols[:2])
            else:
                name = cols[3] if len(cols) > 3 else ''
                value = cols[4] if len(cols) > 4 else ''
                features[-1]['quals'].append([name, value])
    return contigs


def read_assembly_info(path):
    """{contig: (length, coverage, circular)} from Flye's assembly_info.txt."""
    info = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.rstrip('\n').split('\t')
            info[cols[0]] = (int(cols[1]), float(cols[2]), cols[3] == 'Y')
    return info


def qual(feature, name):
    return next((value for key, value in feature['quals'] if key == name), '')


def clean_features(features, log):
    """The fixed feature list for one contig, gene features added."""
    out = []
    for feature in features:
        key = feature['key']
        if key == 'gene':
            # regenerated below, from the feature each one belongs to
            continue
        product = qual(feature, 'product')
        locus_tag = qual(feature, 'locus_tag')
        if key == 'rRNA' and '(partial)' in product:
            log['partial_rrna'].append(f'{locus_tag} {product}')
            continue

        if key == 'CDS' and product.lower() in UNKNOWN_PRODUCTS and product != HYPOTHETICAL:
            log['renamed'] += 1
            product = HYPOTHETICAL
        quals = []
        for name, value in feature['quals']:
            if name == 'product':
                value = product
            if name == 'score':
                continue
            if product == HYPOTHETICAL and name in ('EC_number', 'gene'):
                log['hypothetical_' + name] += 1
                continue
            if name == 'gene':
                stripped = re.sub(r'_\d+$', '', value)
                if stripped != value:
                    log['gene_suffix'] += 1
                value = stripped
            quals.append([name, value])
        if key == 'CDS' and locus_tag and not any(name == 'protein_id' for name, _ in quals):
            quals.append(['protein_id', f'gnl|ncbi|{locus_tag}'])

        if key in GENE_PARENTS and locus_tag:
            gene_quals = [[name, value] for name, value in quals if name == 'gene']
            gene_quals.append(['locus_tag', locus_tag])
            # a gene spans its feature: first interval's start to last's stop
            locs = feature['locs']
            out.append({'locs': [[locs[0][0], locs[-1][1]]], 'key': 'gene', 'quals': gene_quals})
        out.append({'locs': feature['locs'], 'key': key, 'quals': quals})
    return out


def write_tbl(path, contigs):
    with open(path, 'w') as handle:
        for contig, features in contigs:
            handle.write(f'>Feature {contig}\n')
            for feature in features:
                first, *rest = feature['locs']
                handle.write(f'{first[0]}\t{first[1]}\t{feature["key"]}\n')
                for start, stop in rest:
                    handle.write(f'{start}\t{stop}\n')
                for name, value in feature['quals']:
                    handle.write(f'\t\t\t{name}\t{value}\n')


def organism_name(organism, strain):
    # NCBI wants an unnamed species' strain in the name itself
    # ("Pseudomonas sp. KEA-6"); a bare "Pseudomonas sp." fails validation
    if strain and organism.endswith(' sp.'):
        return f'{organism} {strain}'
    return organism


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--fsa', required=True)
    parser.add_argument('--tbl', required=True)
    parser.add_argument('--assembly-info', required=True)
    parser.add_argument('--organism', required=True)
    parser.add_argument('--strain', default='')
    parser.add_argument('--locus-tag', required=True)
    parser.add_argument('--assembly-method', required=True)
    parser.add_argument('--annotation', required=True)
    parser.add_argument('--out-prefix', required=True)
    args = parser.parse_args(argv)

    organism = organism_name(args.organism, args.strain)
    info = read_assembly_info(args.assembly_info)
    tbl = read_tbl(args.tbl)
    log = {'short_contigs': [], 'partial_rrna': [], 'renamed': 0, 'hypothetical_EC_number': 0, 'hypothetical_gene': 0, 'gene_suffix': 0}

    kept = []
    for seq_id, _, seq in read_fasta(args.fsa):
        if len(seq) < MIN_CONTIG:
            log['short_contigs'].append(f'{seq_id} ({len(seq)} bp)')
            continue
        kept.append((seq_id, seq))
    if not kept:
        sys.exit(f'no contig of {MIN_CONTIG} bp or more in {args.fsa}')

    circular = [seq_id for seq_id, _ in kept if info.get(seq_id, (0, 0, False))[2]]
    with open(f'{args.out_prefix}.fsa', 'w') as handle:
        for seq_id, seq in kept:
            mods = [f'[organism={organism}]']
            if args.strain:
                mods.append(f'[strain={args.strain}]')
            mods.append('[gcode=11]')
            if seq_id in circular:
                # topology only: [completeness=complete] fails validation
                # (UnwantedCompleteFlag) without the chromosome/plasmid
                # labels the portal adds, so the validator's warning that a
                # circle isn't flagged complete (CompleteCircleProblem) stays
                mods.append('[topology=circular]')
            handle.write(f'>{seq_id} {" ".join(mods)}\n')
            for i in range(0, len(seq), 60):
                handle.write(seq[i:i + 60] + '\n')

    cleaned = [(seq_id, clean_features(tbl.get(seq_id, []), log)) for seq_id, _ in kept]
    write_tbl(f'{args.out_prefix}.tbl', cleaned)

    # length-weighted mean depth over the contigs kept; Flye's own number
    # for each contig, so it's the depth of the reads Flye assembled from
    sizes = [(len(seq), info[seq_id][1]) for seq_id, seq in kept if seq_id in info]
    total = sum(length for length, _ in sizes)
    coverage = f'{round(sum(length * cov for length, cov in sizes) / total)}x' if total else ''
    comment = [
        ('StructuredCommentPrefix', '##Genome-Assembly-Data-START##'),
        ('Assembly Method', args.assembly_method),
        ('Genome Coverage', coverage),
        ('Sequencing Technology', 'Oxford Nanopore'),
        ('StructuredCommentSuffix', '##Genome-Assembly-Data-END##'),
    ]
    with open(f'{args.out_prefix}.cmt', 'w') as handle:
        handle.writelines(f'{key}\t{value}\n' for key, value in comment if value)

    genome_size = sum(len(seq) for _, seq in kept)
    n_cds = sum(1 for _, features in cleaned for f in features if f['key'] == 'CDS')
    lines = [
        f'GenBank submission files for {args.out_prefix}',
        '',
        'What NCBI\'s Genome Submission Portal asks for:',
        f'  Organism                {organism}',
        f'  Strain                  {args.strain or "(none given)"}',
        f'  Locus tag prefix        {args.locus_tag}',
        f'  Assembly method         {args.assembly_method}',
        f'  Genome coverage         {coverage or "(unknown)"}',
        '  Sequencing technology   Oxford Nanopore',
        f'  Annotation              {args.annotation}',
        f'  Genome size             {genome_size:,} bp in {len(kept)} contig{"s" * (len(kept) != 1)}',
        f'  Circular contigs        {", ".join(circular) or "none"}',
        f'  CDS                     {n_cds:,}',
        '',
        'Changed from Prokka\'s annotation, so NCBI\'s validator accepts it:',
        '  - gene features added (one per CDS/RNA), and a protein_id on every CDS',
        f'  - {log["renamed"]} "putative protein" or similar product(s) renamed "hypothetical protein"',
        f'  - EC number removed from {log["hypothetical_EC_number"]} hypothetical protein(s)',
        f'  - gene name removed from {log["hypothetical_gene"]} hypothetical protein(s)',
        f'  - _1/_2 duplicate suffix removed from {log["gene_suffix"]} gene name(s)',
    ]
    if log['partial_rrna']:
        lines.append(f'  - {len(log["partial_rrna"])} partial rRNA fragment(s) left out: {"; ".join(log["partial_rrna"])}')
    if log['short_contigs']:
        lines.append(f'  - {len(log["short_contigs"])} contig(s) under {MIN_CONTIG} bp left out: {", ".join(log["short_contigs"])}')
    lines += [
        '',
        f'Upload {args.out_prefix}.sqn. Check {args.out_prefix}.val first: every ERROR line',
        'there has to be fixed before NCBI will take the file. A circular contig is',
        'marked circular only; whether it is the chromosome or a plasmid is entered in',
        'the portal.',
        '',
    ]
    with open('submission_info.txt', 'w') as handle:
        handle.write('\n'.join(lines))


if __name__ == '__main__':
    main()
