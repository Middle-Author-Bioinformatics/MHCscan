"""Optional integration tests; external repositories are never required for unit tests."""
import csv
import json
import os
from pathlib import Path
import pytest
from Bio import SeqIO
from mhcscan_eet.cli import main
from mhcscan_eet.evidence import run_feet
from mhcscan_eet.motifs import characterize
from mhcscan_eet.domains import search, setup_db


def repo(name):
    value=os.environ.get('MHC_EET_'+name.upper()+'_REPO')
    if not value:
        pytest.skip(f'MHC_EET_{name.upper()}_REPO not set')
    return Path(value).resolve()


def test_original_feet_motif_compatibility(tmp_path):
    feet=repo('feet')
    seqs={'two':'CAACHGG'*2,'five':'CAACHGG'*5,'mixed':'CAACHGG'*3+'CACHGG'*2,
          'literal_long':'CAACHGG'*3+'CACHGG'+'C'+'A'*11+'C'+'AA'+'CH',
          'intended_long':'CAACHGG'*3+'CACHGG'+'C'+'A'*14+'CH'}
    p=tmp_path/'proteins.faa';p.write_text(''.join(f'>{k}\n{v}\n' for k,v in seqs.items()))
    rows=run_feet(feet,p,tmp_path,[])
    assert set(rows)=={pid for pid,seq in seqs.items() if characterize(seq)['feet_code_motif_pass']}
    for pid,values in rows.items():
        assert int(values[0]['Total_Num_of_heme'])==characterize(seqs[pid])['feet_code_total']


def test_thermincola_references_and_all_tools(tmp_path):
    fg=repo('fegenie');feet=repo('feet')
    db=os.environ.get('MHC_EET_SURFACE_DB')
    if not db:pytest.skip('MHC_EET_SURFACE_DB not set')
    path=fg/'hmms/iron/iron_reduction/non-aligned/TherJR_SLCs.faa'
    inp=tmp_path/'input';inp.mkdir()
    inp.joinpath('Thermincola.faa').write_bytes(path.read_bytes())
    expected={'646799199':9,'646797728':6,'646796949':10}
    for r in SeqIO.parse(path,'fasta'):
        assert characterize(str(r.seq))['canonical_count']==expected[r.id]
    out=tmp_path/'output'
    assert main(['run','--input',str(inp),'--input-type','proteins','--output',str(out),'--surface-db',db,'--fegenie-repo',str(fg),'--feet-repo',str(feet),'--threads','2'])==0
    rows=list(csv.DictReader(open(out/'mhcs.tsv'),delimiter='\t'))
    assert len(rows)==3
    for r in rows:
        assert int(r['canonical_count'])==expected[r['protein_id']]
        assert json.loads(r['feet_original_script_rows'])
        assert r['localization']=='unavailable' and r['neighborhood_status']=='unavailable'
    assert any(json.loads(r['fegenie_hits']) for r in rows)
    manifest=json.load(open(out/'manifest.json'))
    assert len(manifest['models'])==49  # 6 Pfam + 43 FeGenie category-model entries
    assert manifest['repositories']['FEET']['commit']


def test_official_models_local_setup_and_consensus_repeats(tmp_path):
    ph=pytest.importorskip('pyhmmer')
    db=os.environ.get('MHC_EET_SURFACE_DB')
    if not db:pytest.skip('MHC_EET_SURFACE_DB not set')
    # Consensuses are model-derived positive controls, NOT characterized MHCs.
    audit=setup_db(tmp_path/'db',['PF01473','PF04122','PF00395','PF04203','PF01476','PF07501'],Path(db))
    assert len(audit['models'])==6
    with ph.plan7.HMMFile(db) as f:
        models=list(f)
    seqs={}
    for model in models:
        seqs[model.accession]=model.consensus.upper()+'GGGGGGGGGG'+model.consensus.upper()
    faa=tmp_path/'consensuses.faa'
    faa.write_text(''.join(f'>{k}\n{v}\n' for k,v in seqs.items()))
    hits,_=search(Path(db),faa,tmp_path/'raw.domtblout',2)
    for model in models:
        own=[h for h in hits if h['accession']==model.accession and h['protein_id']==model.accession]
        assert len(own)>=2
        assert all(h['threshold'].startswith('Pfam GA') for h in own)
