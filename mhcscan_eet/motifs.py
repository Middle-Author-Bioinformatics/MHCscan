"""Sequence characterization shared with legacy MHCscan.

FEET code-compatible counts are independent regex evaluations, not heme occupancy.
Its [!=C] class matches C, ! or =; it is NOT Python's negated character class.
"""
import re


def legacy_mhc_counts(sequence):
    """Original cytoscan.py non-overlapping CXXCH plus CXXXCH counts."""
    return len(re.findall(r'C(..)CH', sequence)) + len(re.findall(r'C(...)CH', sequence))


def characterize(sequence):
    seq = sequence.upper()
    starts = [m.start() + 1 for m in re.finditer(r'(?=C..CH)', seq)]
    patterns = [r'C..CH', r'C.CH', r'C...CH', r'C...........[!=C]..CH']
    counts = [len(re.findall(p, seq)) for p in patterns]
    intended_long = len(re.findall(r'C.{11}[^C].{2}CH', seq))
    return dict(canonical_count=len(starts), canonical_positions=starts,
                mhcscan_legacy_motif_count=legacy_mhc_counts(seq),
                feet_code_counts=counts, feet_code_total=sum(counts),
                feet_code_motif_pass=counts[0] >= 3 and sum(counts) >= 5,
                feet_intended_total=sum(counts[:3]) + intended_long,
                feet_intended_motif_pass=counts[0] >= 3 and sum(counts[:3]) + intended_long >= 5,
                feet_intended_long_count=intended_long)


def sorting_architecture(seq):
    """Conservative LPXTG heuristic; separate from secretion/localization prediction.

    Require motif in last 80 aa, a downstream 15-aa window with >=11 hydrophobic
    residues, and >=2 K/R in a terminal tail of 3..15 aa. These are heuristics,
    not SignalP/TMHMM predictions or proof of sortase processing.
    """
    hydrophobic = set('AILMFWVY')
    result = []
    for m in re.finditer(r'LP.TG', seq):
        if len(seq) - m.start() > 80:
            continue
        membrane = None
        for i in range(m.end(), len(seq) - 14):
            if sum(a in hydrophobic for a in seq[i:i+15]) >= 11:
                tail = seq[i+15:]
                if 3 <= len(tail) <= 15 and sum(a in 'KR' for a in tail) >= 2:
                    membrane = [i+1, i+15]
                    break
        result.append(dict(motif=m.group(), start=m.start()+1, end=m.end(),
                           downstream_hydrophobic_window=membrane,
                           architecture_pass=membrane is not None))
    return result
