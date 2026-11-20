"""Offline report with embedded downloadable tables and coordinate-scaled SVGs."""
import base64
import html


def esc(x):
    return html.escape(str(x), quote=True)


def svg(anchor, neighbors, hit_ids):
    rows = [dict(neighbor_id=anchor['protein_id'], neighbor_start=anchor['start'], neighbor_end=anchor['end'], neighbor_strand=anchor['strand'], signed_cds_step=0)] + neighbors
    # Wrapped features are shifted into the same unrolled coordinate interval.
    normalized = []
    for r in rows:
        start, end = r['neighbor_start'], r['neighbor_end']
        if r.get('circular_wrap'):
            shift = anchor['contig_length'] * (1 if r['signed_cds_step'] > 0 else -1)
            start, end = start+shift, end+shift
        normalized.append((r, start, end))
    low = min(s for _,s,e in normalized)
    high = max(e for _,s,e in normalized)
    span = max(high-low+1, 1)
    scale = 960/span
    shapes = []
    for r,s,e in normalized:
        x, width = 20+(s-low)*scale, max(2,(e-s+1)*scale)
        arrow = min(8, width/2)
        color = '#b3432e' if r['signed_cds_step']==0 else '#26706e' if r['neighbor_id'] in hit_ids else '#bfc5c8'
        y = 42 if r['signed_cds_step'] % 2 == 0 else 75
        if r['neighbor_strand']=='+':
            points = [(x,y-9),(x+width-arrow,y-9),(x+width,y),(x+width-arrow,y+9),(x,y+9)]
        elif r['neighbor_strand']=='-':
            points = [(x+arrow,y-9),(x+width,y-9),(x+width,y+9),(x+arrow,y+9),(x,y)]
        else:
            points = [(x,y-9),(x+width,y-9),(x+width,y+9),(x,y+9)]
        title = f'{r["neighbor_id"]}: {s}..{e}; strand {r["neighbor_strand"]}; gap {r.get("bp_separation",0)} bp; domains {r.get("neighbor_domains", "[]")}'
        pts = ' '.join(f'{a:.2f},{b}' for a,b in points)
        shapes.append(f'<polygon points="{pts}" fill="{color}"><title>{esc(title)}</title></polygon><text x="{x:.2f}" y="{y+22}" font-size="9">{esc(r["neighbor_id"])}</text>')
    return f'<svg viewBox="0 0 1000 125" role="img" aria-label="Gene neighborhood for {esc(anchor["protein_id"])}"><text x="20" y="15" font-size="12">{esc(anchor["contig"])}: {low}–{high} bp (unrolled if circular)</text>'+''.join(shapes)+'</svg>'


def report(output, candidates, associations, diagrams):
    downloads = []
    for name in ['mhcs.tsv', 'neighbors.tsv', 'domain_hits.tsv', 'genomes.tsv', 'candidates.faa', 'manifest.json', 'localization_evidence.tsv']:
        path = output/name
        if path.exists():
            data = base64.b64encode(path.read_bytes()).decode()
            downloads.append(f'<a download="{name}" href="data:application/octet-stream;base64,{data}">{name}</a>')
    def table(title, rows):
        fields = list(rows[0]) if rows else []
        head = ''.join(f'<th>{esc(f)}</th>' for f in fields)
        body = ''.join('<tr>'+''.join(f'<td>{esc(r.get(f,""))}</td>' for f in fields)+'</tr>' for r in rows)
        return f'<h2>{title} ({len(rows)})</h2><input class="filter" aria-label="Search {title}" placeholder="Search this table"><div class="scroll"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'
    text = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>MHCscan candidate EET report</title><style>
body{font:15px system-ui;margin:30px;color:#203437;background:#f8faf9}h1,h2{color:#164d50}a{margin-right:16px;display:inline-block;margin-bottom:10px}table{border-collapse:collapse;background:white;font-size:12px}td,th{padding:8px;border:1px solid #dce4e3;max-width:460px;overflow-wrap:anywhere}th{position:sticky;top:0;background:#dbe9e6}input{padding:10px;width:320px;margin:10px 0}.scroll{overflow:auto;max-height:550px}svg{width:100%;min-width:700px}details{background:white;padding:12px;margin:12px 0}.note{padding:14px;background:#e6eeeb;border-left:4px solid #26706e}</style>
<h1>Candidate EET proteins in Gram-positive bacteria</h1><p class="note">Motifs identify candidate multiheme proteins, not confirmed EET proteins. Same-protein binding domains, nearby binding proteins, sortase machinery, secretion, membrane attachment, and extracellular exposure are distinct evidence. Proximity does not establish localization, physical interaction, or an operon. Missing localization is unavailable, never a negative prediction. G5 is contextual evidence; sortase is machinery.</p>'''
    text += '<nav>'+''.join(downloads)+'</nav>'+table('MHC candidates', candidates)+table('MHC–neighbor associations', associations)
    text += '<h2>Neighborhoods</h2><p>Red: MHC. Teal: surface-domain hit. Grey: other CDS. Arrow direction shows strand; horizontal spacing uses nucleotide coordinates. Hover for distances and domains. Labels may overlap in dense regions.</p>'
    text += ''.join(f'<details><summary>{esc(label)}</summary><div class="scroll">{diagram}</div></details>' for label, diagram in diagrams)
    text += '''<script>document.querySelectorAll('.filter').forEach(input=>input.addEventListener('input',()=>{const q=input.value.toLowerCase();input.nextElementSibling.querySelectorAll('tbody tr').forEach(r=>{r.hidden=!r.textContent.toLowerCase().includes(q)})}));</script></html>'''
    (output/'report.html').write_text(text)
