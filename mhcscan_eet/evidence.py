"""Strict source-preserving imports; absent predictions remain unavailable."""
import csv
import subprocess
import sys
from Bio import SeqIO

FEATURES = ['localization', 'secretion_signal', 'transmembrane_topology', 'lipoprotein_signal',
            'membrane_attachment', 'cell_wall_association', 'extracellular_exposure']


def table(path, delimiter='\t'):
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def check_ids(rows, known, column='protein_id'):
    for row in rows:
        if row.get(column) not in known:
            raise ValueError(f'Imported {column}={row.get(column)!r} has no exact protein FASTA match')


def import_evidence(directory, known):
    evidence = []
    external = {'fegenie': {}, 'feet': {}, 'mhcscan': {}}
    raw = {}
    if not directory or not directory.exists():
        return evidence, external, raw
    path = directory/'localization.tsv'
    if path.exists():
        rows = table(path)
        check_ids(rows, known)
        for r in rows:
            if not all(r.get(c) for c in ['source', 'version', 'feature', 'value']):
                raise ValueError(f'{path}: requires protein_id, source, version, feature, value')
            if r['feature'] not in FEATURES:
                raise ValueError(f'{path}: unsupported feature {r["feature"]}')
        evidence.extend(rows)
        raw[path.name] = rows
    path = directory/'phobius.short'
    if path.exists():
        rows = []
        for line in path.read_text().splitlines():
            cols = line.split()
            if not cols or cols[0] == 'SEQENCE' or cols[0] == 'SEQUENCE' or line.startswith('#'):
                continue
            if len(cols) < 3 or not cols[1].isdigit() or cols[2] not in ['Y','0','N']:
                raise ValueError(f'{path}: expected native Phobius short output ID TM SP PREDICTION')
            pid, tm, sp = cols[:3]
            rows.append(dict(protein_id=pid, tm_count=int(tm), sp=sp, prediction=' '.join(cols[3:])))
            evidence.extend([dict(protein_id=pid, source='Phobius', version='user_import_unspecified', feature='transmembrane_topology', value=' '.join(cols[3:]), tm_count=int(tm)),
                             dict(protein_id=pid, source='Phobius', version='user_import_unspecified', feature='secretion_signal', value='predicted' if sp=='Y' else 'not_predicted')])
        check_ids(rows, known)
        raw[path.name] = rows
    path = directory/'cello.tsv'
    if path.exists():
        rows = table(path)
        for r in rows:
            pid = r.get('#SeqName')
            location = r.get('#Most-likely-Location')
            if not pid or not location:
                raise ValueError(f'{path}: expected native CELLO Results headers #SeqName and #Most-likely-Location; use normalized localization.tsv otherwise')
            evidence.append(dict(protein_id=pid, source='CELLO_user_import_model_organism_unspecified', version='2.5_import',
                                 feature='localization', value=location,
                                 limitation='Confirm organism model: legacy FEET used Gram-negative categories; not a Gram-positive exposure inference'))
        check_ids(evidence, known)
        raw[path.name] = rows
    path = directory/'psortb.tsv'
    if path.exists():
        rows = table(path)
        # PSORTb native terse output uses SeqID, Final_Localization, Final_Score.
        for r in rows:
            if not r.get('SeqID') or 'Final_Localization' not in r:
                raise ValueError(f'{path}: expected PSORTb terse header SeqID and Final_Localization')
            evidence.append(dict(protein_id=r['SeqID'], source='PSORTb_Gram_positive_user_import', version='3.x_import',
                                 feature='localization', value=r['Final_Localization'], score=r.get('Final_Score','')))
        check_ids(evidence, known)
        raw[path.name] = rows
    path = directory/'signalp.tsv'
    if path.exists():
        # SignalP 6 prediction_results.txt is tabular with '# ID' header.
        lines = [l for l in path.read_text().splitlines() if l.startswith('# ID\t') or not l.startswith('#')]
        rows = list(csv.DictReader([l.removeprefix('# ') for l in lines], delimiter='\t'))
        for r in rows:
            pid, prediction = r.get('ID'), r.get('Prediction')
            if not pid or not prediction:
                raise ValueError(f'{path}: expected SignalP 6 # ID / Prediction columns')
            evidence.append(dict(protein_id=pid, source='SignalP', version='6_import', feature='secretion_signal', value=prediction, cleavage_site=r.get('CS Position','')))
            evidence.append(dict(protein_id=pid, source='SignalP', version='6_import', feature='lipoprotein_signal', value='predicted' if 'LIPO' in prediction else 'not_predicted'))
        check_ids(evidence, known)
        raw[path.name] = rows
    for filename, source, delimiter, key in [('fegenie.csv','fegenie',',','orf'),
                                              ('feet_heme_count.txt','feet','\t','gene_oid'),
                                              ('mhcscan.csv','mhcscan',',','locus')]:
        path = directory/filename
        if not path.exists():
            continue
        # Legacy MHCscan separates clusters with #### lines.
        lines = [l for l in path.read_text().splitlines() if l and not l.startswith('#')]
        rows = list(csv.DictReader(lines, delimiter=delimiter))
        check_ids(rows, known, key)
        genome_names = {r.get('genome/assembly', r.get('cell', '')) for r in rows} - {''}
        if len(genome_names) > 1:
            raise ValueError(f'{path}: contains multiple genomes; split per genome before importing')
        for r in rows:
            external[source].setdefault(r[key], []).append(r)
        raw[filename] = rows
    return evidence, external, raw


def run_feet(repo, protein_path, directory, commands):
    script = repo/'FEET-pipeline-scripts/Find_multihemes_CNO.py'
    if not script.exists():
        raise ValueError(f'{repo}: FEET-pipeline-scripts/Find_multihemes_CNO.py not found')
    stage = directory/'feet_external/ORF_calls'
    stage.mkdir(parents=True, exist_ok=True)
    target = stage/'proteins.faa'
    target.write_bytes(protein_path.read_bytes())
    command = [sys.executable, str(script.resolve()), str(target.resolve()), '3', '5']
    commands.append(dict(argv=command, cwd=str(directory.resolve()), source='FEET original script'))
    with open(directory/'feet_external.log', 'w') as log:
        subprocess.run(command, cwd=directory, stdout=log, stderr=subprocess.STDOUT, check=True)
    path = directory/'feet_external/MHCs/proteins_3_heme_count.txt'
    rows = table(path)
    check_ids(rows, {r.id for r in SeqIO.parse(protein_path, 'fasta')}, 'gene_oid')
    return {r['gene_oid']: [r] for r in rows}
