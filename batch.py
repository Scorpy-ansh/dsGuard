"""Scan many target genes in one pass per organism.

All genes are joined into one query separated by an ambiguity code (so no match can
span two genes) and the reference is streamed past it ONCE. Matches are then split
back onto their genes. 30 genes cost about the same number of reference scans as 1.
"""
import numpy as np

from .seq import SEP


def concat_genes(arrays):
    gap = np.array([SEP], dtype=np.uint8)
    starts, parts, pos = [], [], 0
    for a in arrays:
        starts.append(pos)
        parts += [a, gap]
        pos += a.size + 1
    return np.concatenate(parts), np.array(starts, dtype=np.int64)


def split_matches(matches, starts):
    """Per-gene lists of Match with coordinates local to that gene."""
    out = [[] for _ in starts]
    for m in matches:
        g = int(np.searchsorted(starts, m.q_start, side="right")) - 1
        off = int(starts[g])
        out[g].append(m._replace(q_start=m.q_start - off, q_end=m.q_end - off))
    return out
