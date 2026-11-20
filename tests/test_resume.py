import json
from pathlib import Path
import pytest
from mhcscan_eet.domains import search


def test_search_cache_reuses_and_invalidates(tmp_path,monkeypatch):
    ph=pytest.importorskip('pyhmmer')
    alphabet=ph.easel.Alphabet.amino()
    seq='ACDEFGHIKLMNPQRSTVWY'*3
    hmm,_,_=ph.plan7.Builder(alphabet).build(ph.easel.TextSequence(name=b'test',sequence=seq).digitize(alphabet),ph.plan7.Background(alphabet))
    hmm.cutoffs.gathering=(10.,10.)
    db=tmp_path/'m.hmm'
    with open(db,'wb') as f:hmm.write(f)
    faa=tmp_path/'p.faa';faa.write_text('>p\n'+seq+'\n')
    raw=tmp_path/'m.domtblout'
    expected=search(db,faa,raw,1)
    def unavailable(*a,**kw):raise RuntimeError('search deliberately disabled')
    monkeypatch.setattr(ph.hmmer,'hmmsearch',unavailable)
    assert search(db,faa,raw,1,resume=True)==expected
    raw.write_text('corrupted raw output')
    with pytest.raises(RuntimeError,match='deliberately'):search(db,faa,raw,1,resume=True)
