"""Model setup and HMMER searches. No external models are redistributed."""
from pathlib import Path
from datetime import datetime, timezone
import gzip
import hashlib
import json
import urllib.request

FAMILIES = {
 'PF01473': ('Choline_bind_1', 'binding_repeat', 'Putative choline/cell-wall binding repeat; ligand and boundaries can vary.'),
 'PF04122': ('CW_binding_2', 'binding_domain', 'CWB2; tandem architecture matters; PSII binding demonstrated in specific clostridial proteins.'),
 'PF00395': ('SLH', 'binding_domain', 'S-layer homology; supports secondary cell-wall polymer association, not EET function.'),
 'PF04203': ('Sortase', 'anchoring_machinery', 'Sortase transpeptidase machinery; does not identify a neighboring substrate.'),
 'PF01476': ('LysM', 'binding_domain', 'Peptidoglycan-associated binding module; also common in cell-wall hydrolases.'),
 'PF07501': ('G5', 'surface_context', 'G5; extracellular/adhesive or cell-wall metabolism context; not a universal anchor.'),
}


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def fetch(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'MHCscan-EET/0.1 (research model setup)'})
    with urllib.request.urlopen(req, timeout=90) as r:
        return r.read()


def hmm_metadata(text):
    models = []
    for block in text.split('//'):
        fields = {}
        for line in block.splitlines():
            if line.startswith(('NAME ', 'ACC ', 'DESC ', 'LENG ', 'GA ', 'DATE ')):
                key, value = line.split(None, 1)
                fields[key] = value
        if 'NAME' in fields:
            models.append(fields)
    return models


def setup_db(output, accessions, local_hmm=None):
    output.mkdir(parents=True, exist_ok=True)
    accessions = [a.split('.')[0] for a in accessions]
    content = ''
    audit = dict(retrieved_utc=datetime.now(timezone.utc).isoformat(), models=[], sources=[])
    local = {}
    if local_hmm:
        raw = local_hmm.read_bytes()
        if raw[:2] == b'\x1f\x8b':
            raw = gzip.decompress(raw)
        for block in raw.decode().split('//'):
            meta = hmm_metadata(block)
            if meta:
                key = meta[0].get('ACC', meta[0]['NAME']).split('.')[0]
                local[key] = block.strip() + '\n//\n'
        audit['local_hmm'] = dict(path=str(local_hmm.resolve()), sha256=digest(local_hmm))
    else:
        url = 'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz'
        audit['pfam_release'] = gzip.decompress(fetch(url)).decode().strip()
        audit['sources'].append(url)
    for acc in accessions:
        url = f'https://www.ebi.ac.uk/interpro/api/entry/pfam/{acc}/'
        if local_hmm:
            if acc not in local:
                raise ValueError(f'{acc} not found in {local_hmm}; no replacement model selected silently')
            block = local[acc]
            metadata = {'status': 'user-supplied; live InterPro status not checked'}
        else:
            metadata = json.loads(fetch(url))['metadata']
            if metadata.get('accession') != acc:
                raise ValueError(f'InterPro returned a different accession for {acc}')
            data = fetch(url + '?annotation=hmm')
            if data[:2] == b'\x1f\x8b':
                data = gzip.decompress(data)
            block = data.decode()
            audit['sources'].extend([url, url+'?annotation=hmm'])
        header = hmm_metadata(block)
        if len(header) != 1 or header[0].get('ACC', '').split('.')[0] != acc:
            raise ValueError(f'{acc}: response was not its HMM; use --pfam-hmm with an official Pfam-A.hmm release')
        content += block.strip() + '\n'
        audit['models'].append(dict(accession=acc, header=header[0], interpro=metadata,
                                    hmm_sha256=hashlib.sha256(block.strip().encode()).hexdigest()))
    tmp = output/'surface.hmm.tmp'
    tmp.write_text(content)
    tmp.replace(output/'surface.hmm')
    audit['surface_sha256'] = digest(output/'surface.hmm')
    (output/'database_manifest.json').write_text(json.dumps(audit, indent=2))
    return audit


def require_pyhmmer():
    try:
        import pyhmmer
        return pyhmmer
    except ImportError as e:
        raise ValueError('HMM searches need pyhmmer: pip install "mhcscan-eet[analysis]"; alternatively import --domain-hits') from e


def decode(value):
    return value.decode() if isinstance(value, bytes) else str(value or '')


def search(hmm_path, protein_path, raw_path, threads, source='Pfam', threshold=None, evalue=None, resume=False):
    ph = require_pyhmmer()
    cache_path = Path(str(raw_path)+'.cache.json')
    cache_key = hashlib.sha256(json.dumps(dict(hmm=digest(hmm_path), proteins=digest(protein_path), source=source, threshold=threshold, evalue=evalue, pyhmmer=ph.__version__), sort_keys=True).encode()).hexdigest()
    if resume and cache_path.exists() and Path(raw_path).exists():
        cached = json.loads(cache_path.read_text())
        if cached.get('key') == cache_key and cached.get('raw_sha256') == digest(raw_path):
            return cached['hits'], cached['models']
    with ph.plan7.HMMFile(str(hmm_path)) as f:
        models = list(f)
    if not models:
        raise ValueError(f'{hmm_path}: empty HMM collection')
    hits = []
    model_records = []
    with ph.easel.SequenceFile(str(protein_path), digital=True, alphabet=ph.easel.Alphabet.amino()) as f:
        sequences = f.read_block()
    with open(raw_path, 'wb') as raw:
        for model in models:
            ga = model.cutoffs.gathering
            if threshold is not None:
                mode = f'FeGenie sequence bitscore >= {threshold}'
                options = dict(T=threshold, incT=threshold, domE=10.0, incdomE=10.0)
            elif ga is not None:
                mode = f'Pfam GA sequence={ga[0]} domain={ga[1]}'
                options = dict(bit_cutoffs='gathering')
            elif evalue is not None:
                mode = f'user E <= {evalue} (sequence and independent domain)'
                options = dict(E=evalue, domE=evalue, incE=evalue, incdomE=evalue)
            else:
                raise ValueError(f'{hmm_path}: {decode(model.name)} lacks GA; supply --domain-evalue explicitly')
            accession = decode(model.accession) or decode(model.name)
            record = dict(source=source, accession=accession, name=decode(model.name), model_path=str(hmm_path.resolve()),
                          length=model.M, threshold=mode, sha256=digest(hmm_path))
            model_records.append(record)
            top = next(ph.hmmer.hmmsearch(model, sequences, cpus=threads, **options))
            top.write(raw, format='domains', header=True)
            for hit in top:
                if not hit.included:
                    continue
                for domain in hit.domains:
                    if not domain.included:
                        continue
                    aln = domain.alignment
                    hits.append(dict(protein_id=decode(hit.name), source=source, accession=accession,
                                     name=decode(model.name), role=FAMILIES.get(accession.split('.')[0], ('', 'user_domain', ''))[1],
                                     sequence_bitscore=round(hit.score, 4), sequence_evalue=hit.evalue,
                                     domain_bitscore=round(domain.score, 4), independent_evalue=domain.i_evalue,
                                     conditional_evalue=domain.c_evalue, hmm_from=aln.hmm_from, hmm_to=aln.hmm_to,
                                     target_from=aln.target_from, target_to=aln.target_to,
                                     envelope_from=domain.env_from, envelope_to=domain.env_to,
                                     model_length=model.M, protein_length=hit.length,
                                     model_coverage=round((aln.hmm_to-aln.hmm_from+1)/model.M, 5),
                                     protein_coverage=round((aln.target_to-aln.target_from+1)/hit.length, 5), threshold=mode))
    cache_path.write_text(json.dumps(dict(key=cache_key, raw_sha256=digest(raw_path), hits=hits, models=model_records)))
    return hits, model_records


def fegenie_models(repo, categories):
    roots = [repo/'hmms/iron', repo/'iron']
    root = next((r for r in roots if (r/'HMM-bitcutoffs.txt').exists()), None)
    if root is None:
        raise ValueError(f'{repo}: FeGenie iron/HMM-bitcutoffs.txt not found')
    cutoffs = {}
    for line in (root/'HMM-bitcutoffs.txt').read_text().splitlines()[1:]:
        cols = line.split()
        if len(cols) >= 2:
            cutoffs[cols[0]] = float(cols[1])
    models = []
    for category in categories:
        if not (root/category).is_dir():
            raise ValueError(f'{root}: category {category} not found')
        for path in sorted((root/category).glob('*.hmm')):
            if path.stem not in cutoffs:
                raise ValueError(f'{path.name}: missing FeGenie model-specific cutoff')
            models.append((path, cutoffs[path.stem], category))
    if not models:
        raise ValueError('No FeGenie HMMs selected')
    return models
