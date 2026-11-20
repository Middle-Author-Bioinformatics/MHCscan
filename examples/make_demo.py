#!/usr/bin/env python3
"""Generate explicitly SYNTHETIC annotation and imported-hit fixtures for the report.

Usage: python examples/make_demo.py /tmp/mhc_demo
Do not use these synthetic hits as biological evidence or a sensitivity benchmark.
"""
import csv
from pathlib import Path
import sys

root=Path(sys.argv[1]);inp=root/'proteins';ev=root/'imports'/'synthetic'
inp.mkdir(parents=True,exist_ok=True);ev.mkdir(parents=True,exist_ok=True)
with open(inp/'synthetic.faa','w') as faa, open(inp/'synthetic.gff3','w') as gff:
    gff.write('##gff-version 3\n##sequence-region chromosome 1 40000\n##sequence-region fragment 1 3000\n')
    for i in range(1,30):
        pid=f'CDS_{i}'
        seq='M'+'ACDEFGHIKLMNPQRSTVWY'*8
        if i in [15,27]:seq='M'+('CAACHGG'*5)+'A'*120
        faa.write(f'>{pid}\n{seq}\n')
        gff.write(f'chromosome\tsynthetic\tCDS\t{i*1000}\t{i*1000+3*len(seq)-1}\t.\t{"+" if i%3 else "-"}\t0\tID={pid};partial=false\n')
    faa.write('>fragment_MHC\nMCAACHGGCAACH\n')
    gff.write('fragment\tsynthetic\tCDS\t1\t39\t.\t-\t0\tID=fragment_MHC;partial=true\n')
fields=['protein_id','source','accession','name','target_from','target_to','domain_bitscore','independent_evalue','model_coverage','threshold']
with open(ev/'domains.tsv','w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(fields)
    w.writerow(['CDS_15','synthetic_fixture_NOT_real_search','PF00395.27','SLH',45,84,40,'1e-10',0.9,'synthetic_fixture'])
    w.writerow(['CDS_15','synthetic_fixture_NOT_real_search','PF00395.27','SLH',90,129,42,'1e-11',0.9,'synthetic_fixture'])
    w.writerow(['CDS_17','synthetic_fixture_NOT_real_search','PF04203.20','Sortase',1,120,80,'1e-20',0.96,'synthetic_fixture'])
    w.writerow(['CDS_13','synthetic_fixture_NOT_real_search','PF07501.19','G5',1,70,38,'1e-9',0.93,'synthetic_fixture'])
print(f'Synthetic input: {inp}; synthetic imports: {root/"imports"}')
