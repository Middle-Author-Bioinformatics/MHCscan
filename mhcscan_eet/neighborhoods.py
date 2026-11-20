"""Coordinate-ordered neighborhoods; all CDSs counted, no operon inference."""


def contexts(genome, candidate, gene_window=8, bp_window=None):
    mapping = {g.protein_id: g for g in genome.genes if g.protein_id}
    anchor = mapping.get(candidate)
    if anchor is None:
        return [], ['coordinates_unavailable']
    # Compound CDS bounding boxes can cross the origin or span intervening CDSs.
    # Disable affected neighborhoods rather than manufacture a gene order.
    group = sorted([g for g in genome.genes if g.contig == anchor.contig], key=lambda g: (g.start, g.end, g.id))
    if any(g.compound for g in group):
        return [], ['compound_CDS_on_contig_neighborhood_unavailable']
    ci = genome.contigs.get(anchor.contig, {})
    circular = bool(ci.get('circular'))
    n = len(group)
    i = group.index(anchor)
    limits = []
    if anchor.partial not in ['false', '00', '0', 'False']:
        limits.append('anchor_partial_or_completeness_unknown')
    if not circular:
        if gene_window and (i < gene_window or n-1-i < gene_window):
            limits.append('contig_edge_gene_window_truncated')
        if bp_window is not None and ('length' not in ci or anchor.start-1 < bp_window or ci['length']-anchor.end < bp_window):
            limits.append('contig_edge_bp_window_truncated_or_length_unknown')
    result = []
    for j, neighbor in enumerate(group):
        if i == j:
            continue
        delta = j-i
        if circular:
            options = [delta, delta-n, delta+n]
            delta = min(options, key=lambda x: (abs(x), -x))
        if gene_window and abs(delta) > gene_window:
            continue
        wrapped = circular and ((delta > 0 and j < i) or (delta < 0 and j > i))
        if delta > 0:
            gap = neighbor.start - anchor.end - 1 + (ci['length'] if wrapped else 0)
        else:
            gap = anchor.start - neighbor.end - 1 + (ci['length'] if wrapped else 0)
        overlap = max(0, min(anchor.end, neighbor.end)-max(anchor.start, neighbor.start)+1) if not wrapped else 0
        if overlap:
            gap = -overlap
        if bp_window is not None and max(0, gap) > bp_window:
            continue
        if anchor.strand == neighbor.strand and anchor.strand in ['+', '-']:
            orientation = 'co-directional'
        elif anchor.strand not in ['+', '-'] or neighbor.strand not in ['+', '-']:
            orientation = 'unknown'
        else:
            left = anchor if delta > 0 else neighbor
            orientation = 'convergent' if left.strand == '+' else 'divergent'
        rel = ('downstream' if delta > 0 else 'upstream') if anchor.strand != '-' else ('upstream' if delta > 0 else 'downstream')
        if anchor.strand not in ['+', '-']:
            rel = 'unknown'
        result.append(dict(mhc_id=candidate, mhc_cds_id=anchor.id, contig=anchor.contig,
                           mhc_start=anchor.start, mhc_end=anchor.end, mhc_strand=anchor.strand,
                           neighbor_id=neighbor.protein_id or neighbor.id, neighbor_cds_id=neighbor.id,
                           neighbor_has_protein=bool(neighbor.protein_id), neighbor_start=neighbor.start,
                           neighbor_end=neighbor.end, neighbor_strand=neighbor.strand,
                           signed_cds_step=delta, intervening_cds=abs(delta)-1, bp_separation=gap,
                           overlap_bp=overlap, genomic_side='increasing_coordinate_direction' if delta > 0 else 'decreasing_coordinate_direction',
                           relative_to_mhc_transcription=rel, orientation=orientation, circular_wrap=wrapped,
                           neighbor_partial=neighbor.partial, neighbor_product=neighbor.product))
    return sorted(result, key=lambda x:x['signed_cds_step']), limits
