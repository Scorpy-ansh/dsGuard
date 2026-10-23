"""Continuous-match engine.

Strided seeding: put a seed on every k-th reference position. A perfect match of
length L >= 2k-1 always contains a grid point (the interval of valid seed starts
has L-k+1 >= k integers, so it contains a multiple of k). With k = 11 every match
of >= 21 nt is therefore found for certain, scanning 1/11 of the positions.
Matches shorter than 2k-1 may or may not be reported (best effort only).

The query (the short target gene) is indexed; the reference streams past in
chunks, so memory does not depend on how many organisms are in the panel.
"""
from collections import namedtuple

import numpy as np

from .seq import revcomp

K = 11
Match = namedtuple("Match", "q_start q_end r_start strand")  # q coords on the sense query


def threshold(k=K):
    if not isinstance(k, (int, np.integer)) or isinstance(k, bool) or not 1 <= k <= 16:
        raise ValueError("k must be an integer between 1 and 16")
    return 2 * k - 1

def pack_kmers(seq, k=K):
    """Rolling 2-bit packing -> (uint32 k-mer per start, validity mask)."""
    n = seq.size - k + 1
    if n <= 0:
        return np.zeros(0, np.uint32), np.zeros(0, bool)
    out = np.zeros(n, np.uint32)
    bad = np.zeros(n, bool)
    for j in range(k):
        s = seq[j:j + n]
        bad |= s > 3
        out = (out << np.uint32(2)) | (s & 3).astype(np.uint32)
    return out, ~bad


def _extend(query, ref, qi, r, k):
    a = 0
    while qi - a - 1 >= 0 and r - a - 1 >= 0 and query[qi - a - 1] < 4 \
            and query[qi - a - 1] == ref[r - a - 1]:
        a += 1
    b = k
    while qi + b < query.size and r + b < ref.size and query[qi + b] < 4 \
            and query[qi + b] == ref[r + b]:
        b += 1
    return (qi - a, qi + b, r - a)


def find_matches(query, ref, k=K, chunk=None):
    """Set of maximal perfect matches (q_start, q_end, r_start), one orientation."""
    if not isinstance(k, (int, np.integer)) or isinstance(k, bool) or not 1 <= k <= 16:
        raise ValueError("k must be an integer between 1 and 16")
    chunk = chunk or k * 1_000_000
    chunk -= chunk % k  # grid points must line up across chunks
    chunk = max(chunk, k)
    qk, qok = pack_kmers(query, k)
    qpos = np.flatnonzero(qok)
    if qpos.size == 0 or ref.size < k:
        return set()
    order = np.argsort(qk[qpos], kind="stable")
    q_sort = qk[qpos][order]
    q_pos = qpos[order]
    found = set()
    for start in range(0, ref.size, chunk):
        blk = ref[start:min(start + chunk + k - 1, ref.size)]
        if blk.size < k:
            break
        rk, rok = pack_kmers(blk, k)
        rk, rok = rk[::k], rok[::k]  # <-- the stride
        lo = np.searchsorted(q_sort, rk, "left")
        hi = np.searchsorted(q_sort, rk, "right")
        for gi in np.flatnonzero((hi > lo) & rok):  # real hits are rare
            r = start + int(gi) * k
            for qi in q_pos[lo[gi]:hi[gi]]:
                found.add(_extend(query, ref, int(qi), r, k))
    return found


def scan_both_strands(query, ref, k=K, chunk=None):
    """Screen both strands of the dsRNA: sense (+) and antisense (-).
    Antisense hits are mapped back onto sense coordinates."""
    n = query.size
    out = [Match(a, b, r, "+") for a, b, r in find_matches(query, ref, k, chunk)]
    for a, b, r in find_matches(revcomp(query), ref, k, chunk):
        out.append(Match(n - b, n - a, r, "-"))
    return sorted(out)


def coverage_profile(intervals, n):
    """Per position: length of the longest match covering it (for plots)."""
    cov = np.zeros(n, np.int32)
    for a, b in intervals:
        np.maximum(cov[a:b], b - a, out=cov[a:b])
    return cov


def window_longest(intervals, n, L):
    """For every window [w, w+L): the longest match lying inside it (clipped to
    the window). A window holds a >=21-nt hazard iff this is >= 21."""
    if not isinstance(L, (int, np.integer)) or isinstance(L, bool) or L <= 0:
        raise ValueError("L must be a positive integer")
    nw = n - L + 1
    ws = np.arange(nw)
    best = np.zeros(nw, np.int32)
    for a, b in intervals:
        lo, hi = max(0, a - L + 1), min(nw, b)  # windows that overlap [a, b)
        if hi <= lo:
            continue
        w = ws[lo:hi]
        ov = np.minimum(b, w + L) - np.maximum(a, w)
        np.maximum(best[lo:hi], ov, out=best[lo:hi])
    return best
