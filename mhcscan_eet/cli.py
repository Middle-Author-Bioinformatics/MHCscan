"""CLI orchestration, source-specific results and content-addressed resume."""
import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import csv
import hashlib
import importlib.metadata
import json
import logging
from pathlib import Path
import platform
import subprocess
import sys
from collections import defaultdict
from . import __version__
from .domains import FAMILIES, digest, setup_db, search, fegenie_models, hmm_metadata
from .inputs import load
from .motifs import characterize, sorting_architecture
from .neighborhoods import contexts
from .evidence import FEATURES, import_evidence, run_feet, table, check_ids
from .report import report, svg

LOG = logging.getLogger('mhc-eet')


def j(value):
    return json.dumps(value, separators=(',', ':'), ensure_ascii=True)


def write_tsv(path, rows, fields=None):
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row))
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for row in rows:
            w.writerow({k: j(v) if isinstance(v, (dict, list)) else v for k,v in row.items()})


def revision(path):
    if not path:
        return None
    try:
        return dict(path=str(path.resolve()), commit=subprocess.check_output(['git','-C',str(path),'rev-parse','HEAD'],stderr=subprocess.DEVNULL,text=True).strip(),
                    dirty=bool(subprocess.check_output(['git','-C',str(path),'status','--porcelain'],stderr=subprocess.DEVNULL,text=True)))
    except (OSError, subprocess.CalledProcessError):
        return dict(path=str(path.resolve()), commit='unavailable (not a git checkout)')


def packages():
    result = dict(python=platform.python_version(), mhcscan_eet=__version__)
    for name in ['biopython','pyhmmer','pyrodigal']:
        try:
            result[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            result[name] = 'not_installed'
    return result


def choose_inputs(directory, kind):
    extensions = dict(auto=['.fna','.fa','.fasta','.faa','.pep','.gb','.gbk','.gbff','.genbank'],
                      genomes=['.fna','.fa','.fasta'], proteins=['.faa','.pep','.fa','.fasta'], genbank=['.gb','.gbk','.gbff','.genbank'])
    paths = sorted(p for p in directory.iterdir() if p.is_file() and p.suffix.lower() in extensions[kind])
    if not paths:
        raise ValueError(f'{directory}: no {kind} input files found (uncompressed files required)')
    ids = [p.stem for p in paths]
    if len(set(ids)) != len(ids):
        raise ValueError('Duplicate genome stems; use an explicit --input-type and one annotation source per genome')
    if kind == 'auto':
        modes = {('proteins' if p.suffix in ['.faa','.pep'] else 'genbank' if p.suffix in extensions['genbank'] else 'genomes') for p in paths}
        if len(modes)>1:
            raise ValueError('Mixed annotation types; choose --input-type proteins/genbank/genomes explicitly or use separate folders')
    return paths


def circular_overrides(path, genome):
    if not path:
        return
    for row in table(path):
        if row.get('genome') == genome.id:
            if row.get('contig') not in {g.contig for g in genome.genes} | set(genome.contigs):
                raise ValueError(f'Unknown circular-contig ID {row.get("contig")} in {genome.id}')
            length = int(row['length'])
            existing = genome.contigs.get(row['contig'], {})
            if length < 1 or (existing.get('length') and existing['length'] != length):
                raise ValueError(f'Circular length mismatch: {genome.id}/{row["contig"]}')
            genome.contigs[row['contig']] = dict(length=length, circular=True, evidence='explicit user circular-contigs TSV')


def run(args):
    if args.genes == 0 and args.bp is None:
        raise ValueError('--genes 0 requires --bp')
    args.output.mkdir(parents=True, exist_ok=True)
    handler = logging.FileHandler(args.output/'run.log')
    handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
    LOG.addHandler(handler)
    try:
        return execute(args)
    finally:
        LOG.removeHandler(handler)
        handler.close()


def execute(args):
    paths = choose_inputs(args.input, args.input_type)
    parameters = {k: str(v.resolve()) if isinstance(v, Path) else v for k,v in vars(args).items() if k not in ['resume','force','func']}
    # Hash all inputs, matching annotations, databases, imports, and wrapper source.
    files = set(paths)
    for p in paths:
        files.update(p.with_suffix(ext) for ext in ['.gff','.gff3'] if p.with_suffix(ext).exists())
    if args.surface_db:
        if not args.surface_db.exists() and not args.skip_surface:
            raise ValueError(f'{args.surface_db} not found. Run mhc-eet setup-db --output databases; or explicitly --skip-surface')
        if not args.skip_surface:
            files.add(args.surface_db)
            meta = args.surface_db.parent/'database_manifest.json'
            if meta.exists():
                files.add(meta)
    if args.evidence_dir:
        files.update(p for p in args.evidence_dir.rglob('*') if p.is_file())
    if args.circular_contigs:
        files.add(args.circular_contigs)
    fg_models = fegenie_models(args.fegenie_repo, args.fegenie_categories.split(',')) if args.fegenie_repo else []
    files.update(p for p,cut,cat in fg_models)
    if args.fegenie_repo:
        files.update(args.fegenie_repo.glob('**/HMM-bitcutoffs.txt'))
    if args.feet_repo:
        files.add(args.feet_repo/'FEET-pipeline-scripts/Find_multihemes_CNO.py')
    files.update(Path(__file__).parent.glob('*.py'))
    hashes = {str(p.resolve()):digest(p) for p in sorted(files)}
    versions = packages()
    signature = hashlib.sha256(j(dict(parameters=parameters, hashes=hashes, versions=versions)).encode()).hexdigest()
    completed = args.output/'manifest.json'
    if completed.exists():
        old = json.loads(completed.read_text())
        if args.resume and not args.force and old.get('signature') == signature:
            expected = old.get('output_sha256', {})
            if expected and all((args.output/p).exists() and digest(args.output/p)==sha for p,sha in expected.items()) and (args.output/'report.html').exists() and (args.output/'report.sha256').exists() and (args.output/'report.sha256').read_text().strip()==digest(args.output/'report.html'):
                LOG.info('Resume: verified inputs, parameters, versions and output checksums; reusing completed run')
                return
        if not (args.force or args.resume):
            raise ValueError('Output already contains a completed run; use --resume or a new directory')
        LOG.info('Inputs/parameters/outputs changed; recomputing analysis with fresh source tables')
    # Remove old completion marker before any mutation, so a failed run cannot be resumed as complete.
    completed.unlink(missing_ok=True)
    (args.output/'report.html').unlink(missing_ok=True)
    manifest = dict(schema_version=1, wrapper_version=__version__, signature=signature,
                    started_utc=datetime.now(timezone.utc).isoformat(), parameters=parameters,
                    software=versions, input_sha256=hashes,
                    repositories=dict(MHCscan=revision(Path(__file__).parent.parent), FeGenie=revision(args.fegenie_repo), FEET=revision(args.feet_repo)),
                    commands=[dict(argv=sys.argv, cwd=str(Path.cwd()))], models=[], genomes=[], skipped_analyses=[], warnings=[],
                    scientific_scope='Candidate discovery only; proximity and motifs do not establish EET, physical interaction, localization, or an operon')
    if args.surface_db.exists() and not args.skip_surface:
        manifest['surface_database'] = dict(path=str(args.surface_db.resolve()), sha256=digest(args.surface_db), models=hmm_metadata(args.surface_db.read_text()))
        provenance = args.surface_db.parent/'database_manifest.json'
        if provenance.exists():
            audit = json.loads(provenance.read_text())
            if audit.get('surface_sha256') == digest(args.surface_db):
                manifest['surface_database']['setup_provenance'] = audit
            else:
                manifest['warnings'].append('Adjacent database_manifest.json hash mismatches surface HMM; provenance not adopted')
    progress = args.output/'manifest.inprogress.json'
    progress.write_text(json.dumps(manifest, indent=2))
    all_candidates, all_neighbors, all_domains, all_localizations, summaries, diagrams = [], [], [], [], [], []
    candidate_sequences = []
    for path in paths:
        LOG.info('Analyzing %s', path.name)
        out = args.output/'genomes'/path.stem
        out.mkdir(parents=True, exist_ok=True)
        genome = load(path, args.input_type, out, args.gff_id_attribute, args.translation_table)
        if genome.annotation.startswith('Pyrodigal'):
            manifest['commands'].append(dict(api='pyrodigal.GeneFinder.find_genes', genome=genome.id, meta=True, translation_table=args.translation_table, bins='restricted_to_requested_translation_table'))
        circular_overrides(args.circular_contigs, genome)
        genome.genes.sort(key=lambda g:(g.contig,g.start,g.end,g.id))
        write_tsv(out/'cds_mapping.tsv', [asdict(g) for g in genome.genes],
                  ['id','contig','start','end','strand','protein_id','partial','product','compound'])
        protein_path = out/'proteins.faa'
        protein_path.write_text(''.join(f'>{p}\n{s}\n' for p,s in genome.proteins.items()))
        evidence_path = args.evidence_dir/genome.id if args.evidence_dir else None
        localization, external, raw_imports = import_evidence(evidence_path, genome.proteins)
        if args.feet_repo:
            external['feet_run'] = run_feet(args.feet_repo, protein_path, out, manifest['commands'])
        (out/'imported_source_rows.json').write_text(json.dumps(raw_imports, indent=2))
        domain_hits = []
        if genome.proteins and not args.skip_surface:
            hits, records = search(args.surface_db, protein_path, out/'surface.domtblout', args.threads, evalue=args.domain_evalue, resume=args.resume and not args.force)
            domain_hits.extend(hits)
            manifest['models'].extend(dict(genome=genome.id, **r) for r in records)
            manifest['commands'].append(dict(api='pyhmmer.hmmer.hmmsearch', genome=genome.id, database=str(args.surface_db), cpus=args.threads, thresholds='GA or explicit fallback; see models'))
        if args.skip_surface:
            manifest['skipped_analyses'].append(dict(genome=genome.id, analysis='surface HMM search', reason='--skip-surface; imported domains may still be used'))
        imported_path = evidence_path/'domains.tsv' if evidence_path else None
        if imported_path and imported_path.exists():
            hits = table(imported_path)
            check_ids(hits, genome.proteins)
            for h in hits:
                required = ['source','accession','name','target_from','target_to','domain_bitscore','independent_evalue','model_coverage','threshold']
                if any(not h.get(k) for k in required):
                    raise ValueError(f'{imported_path}: missing required domain-hit columns: {required}')
                h['target_from'],h['target_to'] = int(h['target_from']),int(h['target_to'])
                if not 1 <= h['target_from'] <= h['target_to'] <= len(genome.proteins[h['protein_id']]):
                    raise ValueError(f'{imported_path}: invalid domain coordinates')
                h['role'] = FAMILIES.get(h['accession'].split('.')[0], ('','user_domain',''))[1]
            domain_hits.extend(hits)
        fg_hits = []
        if genome.proteins:
            for model,cutoff,category in fg_models:
                hits, records = search(model, protein_path, out/(category+'_'+model.stem+'.domtblout'), args.threads, source='FeGenie:'+category, threshold=cutoff, resume=args.resume and not args.force)
                manifest['commands'].append(dict(api='pyhmmer.hmmer.hmmsearch', genome=genome.id, source='FeGenie:'+category, database=str(model.resolve()), cpus=args.threads, sequence_bitscore_cutoff=cutoff))
                fg_hits.extend(hits)
                manifest['models'].extend(dict(genome=genome.id, **r) for r in records)
        if not fg_models:
            manifest['skipped_analyses'].append(dict(genome=genome.id, analysis='FeGenie HMM search', reason='--fegenie-repo not supplied; imported results retained if supplied'))
        write_tsv(out/'fegenie_hits.tsv', fg_hits, ['protein_id','source','accession','name','sequence_bitscore','sequence_evalue','domain_bitscore','independent_evalue','target_from','target_to','model_coverage','threshold'])
        by_domain, by_fg, by_loc = defaultdict(list), defaultdict(list), defaultdict(list)
        for h in domain_hits:
            by_domain[h['protein_id']].append(h)
        for h in fg_hits:
            by_fg[h['protein_id']].append(h)
        for h in localization:
            by_loc[h['protein_id']].append(h)
        for h in domain_hits+fg_hits:
            all_domains.append(dict(genome=genome.id, **h))
        for h in localization:
            all_localizations.append(dict(genome=genome.id, **h))
        mapping = {g.protein_id:g for g in genome.genes if g.protein_id}
        genome_candidates = []
        motif_map = {pid:characterize(sequence) for pid,sequence in genome.proteins.items()}
        for pid,sequence in genome.proteins.items():
            motifs = motif_map[pid]
            if motifs['canonical_count'] < args.min_motifs:
                continue
            neighbors, limitations = contexts(genome, pid, args.genes, args.bp)
            same = by_domain[pid]
            association = []
            for n in neighbors:
                nhits = by_domain[n['neighbor_id']] if n['neighbor_has_protein'] else []
                row = dict(genome=genome.id, **n, neighbor_domains=nhits,
                           neighbor_fegenie=by_fg[n['neighbor_id']] if n['neighbor_has_protein'] else [],
                           neighbor_canonical_motifs=motif_map.get(n['neighbor_id'], {}).get('canonical_count', 'unavailable'),
                           neighbor_localization=by_loc[n['neighbor_id']] if n['neighbor_has_protein'] else 'unavailable',
                           qualifies_surface_association=bool(nhits), evidence='separate_neighbor_hypothesis' if nhits else 'context_only',
                           limitations=limitations)
                association.append(row)
            all_neighbors.extend(association)
            categories = ['canonical_motif_candidate']
            for role in sorted({d['role'] for d in same}):
                categories.append('same_protein:'+role)
            for role in sorted({d['role'] for n in association for d in n['neighbor_domains']}):
                categories.append('neighbor:'+role)
            fg = by_fg[pid]
            if fg or external['fegenie'].get(pid):
                categories.append('FeGenie_family_evidence')
            architecture = sorting_architecture(sequence)
            if any(s['architecture_pass'] for s in architecture):
                categories.append('LPXTG_C_terminal_architecture_heuristic')
            gene = mapping.get(pid)
            loc = by_loc[pid]
            loc_values = {f:[r for r in loc if r['feature']==f] or 'unavailable' for f in FEATURES}
            if not loc:
                limitations.append('localization_not_assessed')
            if 'coordinates_unavailable' in limitations or any('neighborhood_unavailable' in x for x in limitations):
                neighborhood_status = 'unavailable'
            else:
                neighborhood_status = 'assessed'
            row = dict(genome=genome.id, protein_id=pid, length_aa=len(sequence),
                       motif_density_per_100aa=round(100*motifs['canonical_count']/len(sequence), 4),
                       surface_domain_status='searched_all_proteins' if not args.skip_surface else 'imported_hits_scope_unspecified' if imported_path and imported_path.exists() else 'not_assessed',
                       fegenie_status='searched_all_proteins' if fg_models else 'imported_results_only' if external['fegenie'] else 'not_assessed',
                       contig=gene.contig if gene else '', start=gene.start if gene else '', end=gene.end if gene else '',
                       strand=gene.strand if gene else '', cds_id=gene.id if gene else '',
                       **motifs, feet_full_pipeline_status='not_assessed_Gram_positive_adaptation',
                       feet_original_script_rows=external.get('feet_run',{}).get(pid, []), feet_imported_rows=external['feet'].get(pid, []),
                       mhcscan_imported_rows=external['mhcscan'].get(pid, []),
                       fegenie_hits=fg, fegenie_imported_rows=external['fegenie'].get(pid, []), same_protein_domains=same,
                       **loc_values, sorting_signal_architecture=architecture,
                       neighborhood_status=neighborhood_status, neighbor_count=len(neighbors),
                       surface_neighbor_ids=[n['neighbor_id'] for n in association if n['qualifies_surface_association']],
                       nearby_mhc_ids=[n['neighbor_id'] for n in association if isinstance(n['neighbor_canonical_motifs'], int) and n['neighbor_canonical_motifs'] >= args.min_motifs],
                       evidence_categories=categories, limitations=limitations)
            all_candidates.append(row)
            genome_candidates.append(row)
            candidate_sequences.append(f'>{genome.id}|{pid}\n{sequence}\n')
            if gene and neighborhood_status=='assessed':
                anchor = dict(protein_id=pid, contig=gene.contig, start=gene.start, end=gene.end, strand=gene.strand,
                              contig_length=genome.contigs.get(gene.contig,{}).get('length', 0))
                diagrams.append((f'{genome.id} / {pid}', svg(anchor, association, {n['neighbor_id'] for n in association if n['qualifies_surface_association']})))
        skipped_features = [f for f in FEATURES if not any(r['feature']==f for r in localization)]
        manifest['skipped_analyses'].append(dict(genome=genome.id, analysis='localization features without predictions', features=skipped_features, reason='no imported predictions; built-in architecture is only a heuristic'))
        if not args.feet_repo:
            manifest['skipped_analyses'].append(dict(genome=genome.id, analysis='original FEET executable', reason='--feet-repo absent; independently calculated code-compatible and intended motif screens still reported'))
        summaries.append(dict(genome=genome.id, proteins=len(genome.proteins), total_cds=len(genome.genes),
                              mhc_candidates=len(genome_candidates), feet_code_motif_pass=sum(r['feet_code_motif_pass'] for r in genome_candidates),
                              same_protein_surface_domains=sum(bool(r['same_protein_domains']) for r in genome_candidates),
                              mhcs_with_surface_neighbors=sum(bool(r['surface_neighbor_ids']) for r in genome_candidates),
                              mhcs_with_localization=sum(r['localization']!='unavailable' for r in genome_candidates),
                              annotation=genome.annotation, warnings=genome.warnings))
        manifest['genomes'].append(dict(id=genome.id, contigs=genome.contigs, annotation=genome.annotation, warnings=genome.warnings))
        progress.write_text(json.dumps(manifest, indent=2))
        for warning in genome.warnings:
            LOG.warning('%s: %s', genome.id, warning)
        LOG.info('%s: %d candidate MHCs; %d CDSs; %d surface domains', genome.id, len(genome_candidates), len(genome.genes), len(domain_hits))
    write_tsv(args.output/'mhcs.tsv', all_candidates, list(all_candidates[0]) if all_candidates else ['genome','protein_id','canonical_count','evidence_categories','limitations'])
    write_tsv(args.output/'neighbors.tsv', all_neighbors, list(all_neighbors[0]) if all_neighbors else ['genome','mhc_id','neighbor_id','intervening_cds','bp_separation','neighbor_domains'])
    write_tsv(args.output/'domain_hits.tsv', all_domains, list(dict.fromkeys(k for h in all_domains for k in h)) if all_domains else ['genome','protein_id','source','accession','target_from','target_to','domain_bitscore','independent_evalue','model_coverage','threshold'])
    write_tsv(args.output/'localization_evidence.tsv', all_localizations, list(dict.fromkeys(k for h in all_localizations for k in h)) if all_localizations else ['genome','protein_id','source','version','feature','value'])
    write_tsv(args.output/'genomes.tsv', summaries)
    (args.output/'candidates.faa').write_text(''.join(candidate_sequences))
    manifest['completed_utc'] = datetime.now(timezone.utc).isoformat()
    manifest['output_sha256'] = {str(p.relative_to(args.output)):digest(p) for p in args.output.rglob('*') if p.is_file() and p.name not in ['manifest.json','manifest.inprogress.json','report.html','report.sha256','run.log']}
    completed.write_text(json.dumps(manifest, indent=2))
    report(args.output, all_candidates, all_neighbors, diagrams)
    (args.output/'report.sha256').write_text(digest(args.output/'report.html')+'\n')
    progress.unlink(missing_ok=True)
    LOG.info('Completed: %s', args.output/'report.html')


def parser():
    p = argparse.ArgumentParser(description='MHCscan Gram-positive candidate EET evidence wrapper (no confirmation or confidence score)')
    p.add_argument('--version', action='version', version=__version__)
    sub = p.add_subparsers(dest='command', required=True)
    db = sub.add_parser('setup-db', help='Download current official Pfam HMMs with version and hash audit, or extract a local release')
    db.add_argument('--output', type=Path, default=Path('databases'))
    db.add_argument('--accessions', default=','.join(FAMILIES))
    db.add_argument('--pfam-hmm', type=Path, help='Local Pfam-A.hmm[.gz]; no network needed')
    r = sub.add_parser('run', help='Analyze a folder of genomes, GenBank annotations, or matched proteins/GFF3')
    r.add_argument('--input', type=Path, required=True)
    r.add_argument('--input-type', choices=['auto','genomes','genbank','proteins'], default='auto')
    r.add_argument('--output', '-o', type=Path, required=True)
    r.add_argument('--surface-db', type=Path, default=Path('databases/surface.hmm'))
    r.add_argument('--skip-surface', action='store_true')
    r.add_argument('--domain-evalue', type=float, help='Explicit sequence/domain E-value fallback only for models lacking GA')
    r.add_argument('--fegenie-repo', type=Path)
    r.add_argument('--fegenie-categories', default='iron_reduction,iron_oxidation')
    r.add_argument('--feet-repo', type=Path, help='Run original FEET motif script; full CELLO/PCC pipeline is not required')
    r.add_argument('--evidence-dir', type=Path, help='Per-genome imports; see docs/GRAM_POSITIVE_EET.md')
    r.add_argument('--gff-id-attribute', help='Exact GFF3 attribute mapping to protein FASTA IDs')
    r.add_argument('--circular-contigs', type=Path, help='TSV: genome, contig, length; explicit circularity only')
    r.add_argument('--translation-table', type=int, choices=[4,11], default=11)
    r.add_argument('--min-motifs', type=int, default=2)
    r.add_argument('--genes', type=int, default=8, help='CDSs on EACH genomic side; 0 for bp-only')
    r.add_argument('--bp', type=int, help='Maximum inter-CDS gap in bp; intersects gene window when both set')
    r.add_argument('--threads','-t', type=int, default=4)
    r.add_argument('--resume', action='store_true', help='Reuse identical completed analysis after checking content hashes')
    r.add_argument('--force', action='store_true', help='Recompute a completed output directory')
    return p


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
    p = parser()
    args = p.parse_args(argv)
    try:
        if args.command == 'setup-db':
            audit = setup_db(args.output, args.accessions.split(','), args.pfam_hmm)
            LOG.info('Saved %d models to %s', len(audit['models']), args.output/'surface.hmm')
        else:
            if args.threads < 1 or args.min_motifs < 1 or args.genes < 0 or (args.bp is not None and args.bp < 0):
                raise ValueError('threads/min-motifs must be positive; genes/bp must be nonnegative')
            if args.domain_evalue is not None and args.domain_evalue <= 0:
                raise ValueError('--domain-evalue must be positive')
            run(args)
    except (ValueError, OSError, subprocess.CalledProcessError) as e:
        LOG.error('%s', e)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main())
