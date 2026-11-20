import csv
import json
from pathlib import Path
import subprocess
import sys
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, BeforePosition
from mhcscan_eet.motifs import characterize, sorting_architecture, legacy_mhc_counts
from mhcscan_eet.inputs import Gene, Genome, gff, genbank, predict, load
from mhcscan_eet.neighborhoods import contexts
from mhcscan_eet.domains import search, setup_db
from mhcscan_eet.evidence import import_evidence
from mhcscan_eet.cli import main


def test_motifs_thresholds_and_legacy():
    assert characterize('CAACH') ['canonical_count'] == 1
    broad = characterize('CAACHGGCAACH')
    assert broad['canonical_count'] == 2
    assert not broad['feet_code_motif_pass']
    assert characterize('CAACHGG'*5)['feet_code_motif_pass']
    mixed = characterize('CAACHGG'*3 + 'CACHGGCACH')
    assert mixed['feet_code_total'] == 5
    assert mixed['feet_code_motif_pass']
    assert legacy_mhc_counts('CAACHGGCAAACH') == 2


def test_feet_long_regex_discrepancy():
    # Published notation's != is not Python negation. Record both outcomes.
    literal = characterize('C'+'A'*11+'C'+'AA'+'CH')
    intended = characterize('C'+'A'*11+'A'+'AA'+'CH')
    assert literal['feet_code_counts'][3] == 1
    assert literal['feet_intended_long_count'] == 0
    assert intended['feet_code_counts'][3] == 0
    assert intended['feet_intended_long_count'] == 1


def test_sorting_architecture_not_isolated_motif():
    assert not sorting_architecture('M'+'A'*20+'LPATG'+'E'*20)[0]['architecture_pass']
    assert sorting_architecture('M'+'A'*20+'LPATG'+'E'*3+'L'*15+'RKRK')[0]['architecture_pass']
    assert sorting_architecture('LPATG'+'A'*100) == []


def fixture_genome(n=20):
    genes = [Gene(f'g{i}', 'c', i*100+1, i*100+60, '+' if i%2==0 else '-', f'p{i}', 'false') for i in range(n)]
    return Genome('test', genes=genes, proteins={g.protein_id:'MCAACHCAACH' for g in genes}, contigs={'c':dict(length=n*100,circular=False)})


def test_gene_order_all_cds_contig_and_gene_window():
    g = fixture_genome()
    # A CDS without a translated protein is still counted.
    g.genes[10].protein_id = ''
    g.genes.append(Gene('other','other_contig',1,300,'+','unrelated'))
    g.genes.reverse()
    rows, limits = contexts(g,'p9')
    assert len(rows) == 16
    assert {r['signed_cds_step'] for r in rows} == set(range(-8,0))|set(range(1,9))
    assert all(r['contig']=='c' for r in rows)
    assert next(r for r in rows if r['neighbor_id']=='p11')['intervening_cds']==1
    assert any(r['neighbor_id']=='g10' and not r['neighbor_has_protein'] for r in rows)
    assert not limits


def test_spacing_strand_and_overlap():
    g = fixture_genome(3)
    rows,_ = contexts(g,'p1',8)
    assert rows[0]['bp_separation'] == 40
    assert rows[0]['orientation']=='convergent'
    assert rows[0]['relative_to_mhc_transcription']=='downstream'
    assert rows[1]['orientation']=='divergent'
    g.genes[2].start=150
    rows,_ = contexts(g,'p1',8)
    r=next(r for r in rows if r['neighbor_id']=='p2')
    assert r['bp_separation']==-11
    assert r['overlap_bp']==11


def test_bp_window_intersection_and_bp_only():
    g=fixture_genome()
    rows,_=contexts(g,'p9',8,40)
    assert len(rows)==2
    rows,_=contexts(g,'p9',0,140)
    assert len(rows)==4


def test_contig_edges_and_explicit_circularity():
    g=fixture_genome()
    rows,limits=contexts(g,'p0')
    assert len(rows)==8
    assert 'contig_edge_gene_window_truncated' in limits
    assert not any(r['circular_wrap'] for r in rows)
    g.contigs['c']['circular']=True
    rows,limits=contexts(g,'p0')
    assert len(rows)==16 and not limits
    prev=next(r for r in rows if r['neighbor_id']=='p19')
    assert prev['circular_wrap'] and prev['bp_separation']==40 and prev['intervening_cds']==0


def test_missing_coordinates_partial_compound():
    g=fixture_genome(3)
    assert contexts(g,'unmapped')[1]==['coordinates_unavailable']
    g.genes[0].partial='true'
    assert 'anchor_partial_or_completeness_unknown' in contexts(g,'p0')[1]
    g.genes[1].compound=True
    assert contexts(g,'p0')[1]==['compound_CDS_on_contig_neighborhood_unavailable']


def test_gff_exact_mapping_parent_and_unmapped_cds(tmp_path):
    p=tmp_path/'g.gff3'
    p.write_text('##gff-version 3\n##sequence-region c 1 1000\nc\tx\tmRNA\t100\t300\t.\t+\t.\tID=t1;protein_id=p1\nc\tx\tCDS\t100\t300\t.\t+\t0\tID=cds1;Parent=t1\nc\tx\tCDS\t400\t600\t.\t-\t0\tID=noncoding_pseudo\n')
    g=Genome('g', proteins={'p1':'MCAACHCAACH','p10':'MA'})
    gff(g,p)
    assert len(g.genes)==2
    assert g.genes[0].protein_id=='p1'
    assert g.genes[1].protein_id==''
    assert g.warnings and not g.contigs['c']['circular']


def test_gff_ambiguous_mapping_fails(tmp_path):
    p=tmp_path/'g.gff3';p.write_text('c\tx\tCDS\t1\t30\t.\t+\t0\tID=p1;protein_id=p2\n')
    with pytest.raises(ValueError,match='ambiguous'):
        gff(Genome('g',proteins={'p1':'M','p2':'M'}),p)
    g=Genome('g',proteins={'p1':'M','p2':'M'})
    gff(g,p,'protein_id')
    assert g.genes[0].protein_id=='p2'


def test_genbank_translation_and_partial(tmp_path):
    rec=SeqRecord(Seq('ATG'*200),id='c',annotations={'molecule_type':'DNA','topology':'circular'})
    rec.features=[SeqFeature(FeatureLocation(BeforePosition(0),90,strand=-1),type='CDS',qualifiers={'protein_id':['p'],'locus_tag':['tag'],'translation':['MCAACHCAACH']}),SeqFeature(FeatureLocation(100,190,strand=1),type='CDS',qualifiers={'locus_tag':['pseudo'],'pseudo':['']})]
    p=tmp_path/'g.gbk';SeqIO.write(rec,p,'genbank')
    g=genbank(p)
    assert g.proteins['p']=='MCAACHCAACH'
    assert len(g.genes)==2 and g.genes[0].start==1 and g.genes[0].strand=='-'
    assert g.genes[0].partial=='true' and g.contigs['c']['circular']


def test_annotation_mixing_rejected(tmp_path):
    p=tmp_path/'g.fna';p.write_text('>c\nATGATGATG\n')
    p.with_suffix('.faa').write_text('>p\nMCAACHCAACH\n')
    with pytest.raises(ValueError,match='adjacent annotations'):
        load(p,'genomes',tmp_path)


def test_prediction_is_reproducible(tmp_path):
    pytest.importorskip('pyrodigal')
    p=tmp_path/'g.fna'
    # Deterministic coding sequence embedded in noncoding flanks.
    p.write_text('>c\n'+'TAA'*30+'ATG'+'GCTGATGGT'*80+'TAA'+'TAA'*30+'\n')
    a=predict(p,tmp_path)
    b=predict(p,tmp_path)
    assert a.proteins==b.proteins and [vars(x) for x in a.genes]==[vars(x) for x in b.genes]
    assert len(a.genes)>0 and all(0<x.start<=x.end<=len(next(SeqIO.parse(p,'fasta')) ) for x in a.genes)


def test_hmmer_repeated_domains_and_ga(tmp_path):
    ph=pytest.importorskip('pyhmmer')
    alphabet=ph.easel.Alphabet.amino()
    # Build a reproducible test model, not an alleged biological domain model.
    seq='ACDEFGHIKLMNPQRSTVWY'*3
    msa=ph.easel.TextMSA(name=b'test_repeat', sequences=[ph.easel.TextSequence(name=b'a',sequence=seq),ph.easel.TextSequence(name=b'b',sequence=seq)])
    hmm,_,_=ph.plan7.Builder(alphabet).build_msa(msa.digitize(alphabet), ph.plan7.Background(alphabet))
    hmm.accession=b'PF99999.1';hmm.cutoffs.gathering=(10.,10.)
    db=tmp_path/'test.hmm'
    with open(db,'wb') as f:hmm.write(f)
    faa=tmp_path/'p.faa';faa.write_text('>p\n'+seq+'GGGGGGGGGG'+seq+'\n')
    hits,models=search(db,faa,tmp_path/'raw',1)
    assert len(hits)==2 and hits[0]['target_to']<hits[1]['target_from']
    assert all(h['model_coverage']>0.9 for h in hits)
    assert models[0]['threshold'].startswith('Pfam GA')
    hmm.cutoffs.gathering=None
    with open(db,'wb') as f:hmm.write(f)
    with pytest.raises(ValueError,match='lacks GA'):search(db,faa,tmp_path/'raw',1)
    assert search(db,faa,tmp_path/'raw',1,evalue=1e-4)[0]


def write_fixture(tmp_path):
    inp=tmp_path/'input';inp.mkdir()
    inp.joinpath('sample.faa').write_text('>mhc\nM'+'CAACHGG'*5+'\n>anchor\n'+'A'*60+'\n>none\nMAAAAA\n>lonely\nMCAACHGGCAACH\n')
    inp.joinpath('sample.gff3').write_text('##gff-version 3\n##sequence-region c 1 3000\n##sequence-region other 1 1000\nc\tx\tCDS\t501\t610\t.\t+\t0\tID=mhc\nc\tx\tCDS\t701\t880\t.\t-\t0\tID=anchor\nc\tx\tCDS\t901\t920\t.\t+\t0\tID=none\nother\tx\tCDS\t101\t150\t.\t+\t0\tID=lonely\n')
    ev=tmp_path/'evidence'/'sample';ev.mkdir(parents=True)
    ev.joinpath('domains.tsv').write_text('protein_id\tsource\taccession\tname\ttarget_from\ttarget_to\tdomain_bitscore\tindependent_evalue\tmodel_coverage\tthreshold\n'+'mhc\tPfam\tPF00395.1\tSLH\t1\t10\t40\t1e-10\t0.9\tGA\n'+'mhc\tPfam\tPF00395.1\tSLH\t11\t20\t42\t1e-11\t0.95\tGA\n'+'anchor\tPfam\tPF04203.1\tSortase\t1\t40\t50\t1e-15\t0.9\tGA\n')
    ev.joinpath('fegenie.csv').write_text('category,genome/assembly,orf,HMM,bitscore,bitscore_cutoff\niron_reduction,sample,mhc,example_family,100,80\n')
    ev.joinpath('feet_heme_count.txt').write_text('gene_oid\tNum_of_C..CH\tTotal_Num_of_heme\nmhc\t5\t5\n')
    ev.joinpath('mhcscan.csv').write_text('locus,periplasmic_or_outer_membrane,hemes,seq\nmhc,Y,5,MCAACH\n####\n')
    return inp,ev.parent


def test_cli_evidence_separation_no_localization_and_retained_lonely(tmp_path):
    inp,ev=write_fixture(tmp_path);out=tmp_path/'out'
    args=['run','--input',str(inp),'--input-type','proteins','--output',str(out),'--skip-surface','--evidence-dir',str(ev)]
    assert main(args)==0
    rows=list(csv.DictReader(open(out/'mhcs.tsv'),delimiter='\t'))
    assert len(rows)==2
    a,b=rows
    assert len(json.loads(a['same_protein_domains']))==2
    assert json.loads(a['surface_neighbor_ids'])==['anchor']
    assert 'same_protein:binding_domain' in a['evidence_categories']
    assert 'neighbor:anchoring_machinery' in a['evidence_categories']
    assert 'same_protein:anchoring_machinery' not in a['evidence_categories']
    assert a['localization']=='unavailable' and a['extracellular_exposure']=='unavailable'
    assert json.loads(a['fegenie_imported_rows'])[0]['HMM']=='example_family'
    assert json.loads(a['feet_imported_rows'])[0]['gene_oid']=='mhc'
    assert json.loads(a['mhcscan_imported_rows'])[0]['locus']=='mhc'
    assert json.loads(b['surface_neighbor_ids'])==[] and b['neighbor_count']=='0'
    assert '<svg' in (out/'report.html').read_text() and 'data:application/octet-stream;base64' in (out/'report.html').read_text()
    stamp=(out/'mhcs.tsv').stat().st_mtime_ns
    assert main(args+['--resume'])==0 and (out/'mhcs.tsv').stat().st_mtime_ns==stamp
    (out/'mhcs.tsv').write_text('corrupt')
    assert main(args+['--resume'])==0
    assert 'protein_id' in (out/'mhcs.tsv').read_text()


def test_native_localization_imports_and_unknown_ids(tmp_path):
    ev=tmp_path
    (ev/'phobius.short').write_text('SEQENCE ID TM SP PREDICTION\np 1 Y n1-3c4/5o6-50i\n')
    (ev/'psortb.tsv').write_text('SeqID\tFinal_Localization\tFinal_Score\np\tCellwall\t9.8\n')
    (ev/'signalp.tsv').write_text('# SignalP-6.0\n# ID\tPrediction\tCS Position\np\tLIPO(Sec/SPII)\tCS pos: 21-22\n')
    evidence,_,_=import_evidence(ev,{'p'})
    assert any(x['feature']=='localization' and x['value']=='Cellwall' for x in evidence)
    assert any(x['feature']=='lipoprotein_signal' and x['value']=='predicted' for x in evidence)
    assert not any(x['feature']=='extracellular_exposure' for x in evidence)
    with pytest.raises(ValueError,match='no exact'):import_evidence(ev,{'other'})


def test_cello_native_results_exact_ids(tmp_path):
    (tmp_path/'cello.tsv').write_text('#SeqNO.\textra\touter\tperi\tinner\tcyto\t#Most-likely-Location\t#SeqName\n#1\t1.9\t0.1\t0.1\t0.1\t0.1\tExtracellular\tp\n')
    evidence,_,raw=import_evidence(tmp_path,{'p'})
    assert evidence[0]['value']=='Extracellular'
    assert evidence[0]['feature']=='localization'
    assert 'Gram-negative' in evidence[0]['limitation']
    assert raw['cello.tsv'][0]['#SeqName']=='p'


def test_sequence_only_empty_and_threshold_cli(tmp_path):
    inp=tmp_path/'in';inp.mkdir()
    (inp/'g.faa').write_text('>one\nMCAACH\n>two\nMCAACHGGCAACH\n')
    out=tmp_path/'out'
    assert main(['run','--input',str(inp),'--output',str(out),'--skip-surface'])==0
    rows=list(csv.DictReader(open(out/'mhcs.tsv'),delimiter='\t'))
    assert len(rows)==1 and rows[0]['protein_id']=='two' and rows[0]['neighborhood_status']=='unavailable'
    out2=tmp_path/'none'
    assert main(['run','--input',str(inp),'--output',str(out2),'--skip-surface','--min-motifs','3'])==0
    assert len(list(csv.DictReader(open(out2/'mhcs.tsv'),delimiter='\t')))==0
    assert (out2/'report.html').exists()
    assert main(['run','--input',str(inp),'--output',str(tmp_path/'bad'),'--skip-surface','--genes','0'])==2


def test_negative_unknown_incomplete_neighbor_and_bp_boundary():
    g=fixture_genome(3)
    g.genes[0].end=180
    g.genes[1].start=100;g.genes[1].end=140
    rows,_=contexts(g,'p0',8)
    r=next(r for r in rows if r['neighbor_id']=='p1')
    assert r['overlap_bp']==41 and r['bp_separation']==-41
    assert len(contexts(g,'p0',8,0)[0])==1


def test_split_gff_mapping_preserved_conservatively(tmp_path):
    p=tmp_path/'g.gff3';p.write_text('c\tx\tCDS\t1\t30\t.\t+\t0\tID=cds;protein_id=p\nc\tx\tCDS\t50\t70\t.\t+\t0\tID=cds;protein_id=p\n')
    g=Genome('g',proteins={'p':'MCAACHCAACH'})
    gff(g,p)
    assert len(g.genes)==1 and g.genes[0].compound
    assert contexts(g,'p')[1]==['compound_CDS_on_contig_neighborhood_unavailable']
