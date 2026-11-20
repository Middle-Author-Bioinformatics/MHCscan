"""Explicit annotation selection and conservative protein/CDS mapping (1-based inclusive)."""
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import unquote
import re
from Bio import SeqIO
from Bio.SeqFeature import ExactPosition


@dataclass
class Gene:
    id: str
    contig: str
    start: int
    end: int
    strand: str
    protein_id: str = ''
    partial: str = 'unknown'
    product: str = ''
    compound: bool = False


@dataclass
class Genome:
    id: str
    proteins: dict = field(default_factory=dict)
    genes: list = field(default_factory=list)
    contigs: dict = field(default_factory=dict)
    sources: list = field(default_factory=list)
    annotation: str = ''
    warnings: list = field(default_factory=list)


def fasta(path):
    records = list(SeqIO.parse(path, 'fasta'))
    if not records or len({r.id for r in records}) != len(records):
        raise ValueError(f'{path}: empty FASTA or duplicate sequence identifiers')
    return records


def proteins(path):
    result = {}
    for r in fasta(path):
        seq = str(r.seq).upper().rstrip('*')
        if not seq or '*' in seq or not re.fullmatch('[ACDEFGHIKLMNPQRSTVWYBXZJUO]+', seq):
            raise ValueError(f'{path}: invalid protein or internal stop in {r.id}')
        result[r.id] = seq
    return result


def attributes(text):
    return {unquote(k): unquote(v) for k, v in (part.split('=', 1) for part in text.split(';') if '=' in part)}


def gff(genome, path, id_attribute=None):
    rows = []
    parents = {}
    with open(path) as f:
        for line in f:
            if line.startswith('##FASTA'):
                break
            if line.startswith('##sequence-region'):
                _, c, lo, hi = line.split()
                genome.contigs[c] = dict(length=int(hi), circular=False)
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 9:
                raise ValueError(f'{path}: GFF3 row must have 9 columns')
            a = attributes(parts[8])
            if a.get('ID'):
                parents[a['ID']] = a
            if parts[2] == 'region' and a.get('Is_circular', '').lower() == 'true':
                genome.contigs[parts[0]] = dict(length=int(parts[4]), circular=True)
            if parts[2] == 'CDS':
                rows.append((parts, a))
    grouped = {}
    used_proteins = set()
    group_identity = {}
    for n, (p, a) in enumerate(rows, 1):
        # Prefer explicit protein_id, then user-selected field. Resolve Parent
        # only via an exact identifier, never by substring or stripping versions.
        keys = [id_attribute] if id_attribute else ['protein_id', 'ID', 'locus_tag', 'Parent']
        candidates = []
        for key in keys:
            value = a.get(key, '')
            for alias in value.split(','):
                if alias in genome.proteins:
                    candidates.append(alias)
                pa = parents.get(alias, {})
                if pa.get('protein_id') in genome.proteins:
                    candidates.append(pa['protein_id'])
        candidates = list(dict.fromkeys(candidates))
        if len(candidates) > 1:
            raise ValueError(f'{path}: ambiguous CDS mapping {a}; use --gff-id-attribute')
        pid = candidates[0] if candidates else ''
        gid = a.get('ID') or a.get('Parent') or a.get('locus_tag') or f'CDS_{n}'
        groupkey = (p[0], pid or gid)
        start, end = int(p[3]), int(p[4])
        if start < 1 or end < start or p[6] not in ['+', '-', '.','?']:
            raise ValueError(f'{path}: invalid CDS coordinates/strand: {gid}')
        partial = a.get('partial', 'unknown')
        identity = a.get('Parent') or gid
        if groupkey in grouped:
            if group_identity[groupkey] != identity:
                raise ValueError(f'{path}: protein {pid} assigned to different CDS IDs without a common Parent')
            gene = grouped[groupkey]
            if gene.strand != p[6]:
                raise ValueError(f'{path}: inconsistent split CDS strands for {gid}')
            gene.start, gene.end = min(gene.start, start), max(gene.end, end)
            gene.compound = True
        else:
            if pid and pid in used_proteins:
                raise ValueError(f'{path}: protein {pid} maps to multiple contigs/CDSs')
            group_identity[groupkey] = identity
            grouped[groupkey] = Gene(gid, p[0], start, end, p[6], pid, partial, a.get('product', ''))
            if pid:
                used_proteins.add(pid)
    genome.genes = list(grouped.values())
    missing = set(genome.proteins) - used_proteins
    if missing:
        genome.warnings.append(f'{len(missing)} proteins have no exact CDS mapping; their neighborhoods are unavailable')
    if not rows:
        genome.warnings.append('No CDS features; sequence screening only')
    genome.sources.append(str(path.resolve()))


def genbank(path):
    genome = Genome(path.stem, sources=[str(path.resolve())], annotation='GenBank; no gene prediction')
    for rec in SeqIO.parse(path, 'genbank'):
        if rec.id in genome.contigs:
            raise ValueError(f'{path}: duplicate GenBank record {rec.id}')
        genome.contigs[rec.id] = dict(length=len(rec), circular=rec.annotations.get('topology') == 'circular')
        for n, feat in enumerate(rec.features):
            if feat.type != 'CDS':
                continue
            q = feat.qualifiers
            gid = q.get('locus_tag', q.get('protein_id', [f'{rec.id}_CDS{n}']))[0]
            pid = q.get('protein_id', [gid])[0] if 'translation' in q else ''
            if pid:
                if pid in genome.proteins:
                    raise ValueError(f'{path}: duplicate protein {pid}')
                genome.proteins[pid] = q['translation'][0].upper().rstrip('*')
            partial = str(any(not isinstance(x, ExactPosition) for part in feat.location.parts for x in [part.start, part.end]) or 'partial' in q).lower()
            genome.genes.append(Gene(gid, rec.id, int(feat.location.start)+1, int(feat.location.end),
                                     '+' if feat.location.strand == 1 else '-' if feat.location.strand == -1 else '.',
                                     pid, partial, q.get('product', [''])[0], len(feat.location.parts) > 1))
    if not genome.contigs:
        raise ValueError(f'{path}: no GenBank records')
    return genome


def predict(path, out, translation_table=11):
    try:
        import pyrodigal
    except ImportError as e:
        raise ValueError('Genome prediction needs pyrodigal: pip install "mhcscan-eet[analysis]"') from e
    genome = Genome(path.stem, sources=[str(path.resolve())], annotation=f'Pyrodigal meta; translation table {translation_table}')
    bins = pyrodigal.MetagenomicBins([b for b in pyrodigal.METAGENOMIC_BINS if b.training_info.translation_table == translation_table])
    if not len(bins):
        raise ValueError(f'No Pyrodigal meta training bins for translation table {translation_table}')
    finder = pyrodigal.GeneFinder(meta=True, metagenomic_bins=bins)
    with open(out / 'predicted.gff3', 'w') as gf:
        gf.write('##gff-version 3\n')
        for rec in fasta(path):
            genome.contigs[rec.id] = dict(length=len(rec), circular=False)
            genes = finder.find_genes(str(rec.seq).encode())
            gf.write(f'##sequence-region {rec.id} 1 {len(rec)}\n')
            for i, gene in enumerate(genes, 1):
                pid = f'{rec.id}_{i}'
                genome.proteins[pid] = gene.translate(translation_table=translation_table).rstrip('*')
                partial = str(gene.partial_begin or gene.partial_end).lower()
                strand = '+' if gene.strand == 1 else '-'
                genome.genes.append(Gene(pid, rec.id, gene.begin, gene.end, strand, pid, partial))
                gf.write(f'{rec.id}\tPyrodigal\tCDS\t{gene.begin}\t{gene.end}\t.\t{strand}\t0\tID={pid};partial={partial}\n')
    return genome


def load(path, kind, out, id_attribute=None, translation_table=11):
    if kind == 'auto':
        if path.suffix.lower() in ['.gb', '.gbk', '.gbff', '.genbank']:
            kind = 'genbank'
        elif path.suffix.lower() in ['.faa', '.pep']:
            kind = 'proteins'
        else:
            kind = 'genomes'
    if kind == 'genbank':
        genome = genbank(path)
    elif kind == 'proteins':
        genome = Genome(path.stem, proteins=proteins(path), sources=[str(path.resolve())], annotation='Existing protein FASTA')
        matches = [p for ext in ['.gff', '.gff3'] if (p := path.with_suffix(ext)).exists()]
        if len(matches) > 1:
            raise ValueError(f'{path}: both .gff and .gff3 found; keep one matching annotation')
        if matches:
            gff(genome, matches[0], id_attribute)
        else:
            genome.warnings.append('No matching GFF3; sequence screening only; neighborhoods unavailable')
    else:
        # Never silently ignore adjacent annotations and predict a second set.
        if any(path.with_suffix(e).exists() for e in ['.gff', '.gff3', '.gbk', '.gbff', '.faa']):
            raise ValueError(f'{path}: adjacent annotations/proteins found; use --input-type proteins or genbank in a separate input folder')
        genome = predict(path, out, translation_table)
    seen = set()
    for gene in genome.genes:
        if (gene.contig, gene.id) in seen:
            raise ValueError(f'{path}: duplicate CDS ID {gene.id}')
        seen.add((gene.contig, gene.id))
        if gene.contig in genome.contigs and gene.end > genome.contigs[gene.contig]['length']:
            raise ValueError(f'{path}: CDS outside contig: {gene.id}')
    return genome
