"""Slow, obviously-correct reference implementation. Used only to test engine.py."""
import numpy as np


def brute_matches(query, ref, min_len):
    out = set()
    n, m = query.size, ref.size
    for d in range(-(n - 1), m):  # diagonal: r = q + d
        qs = max(0, -d)
        rs = qs + d
        L = min(n - qs, m - rs)
        if L < min_len:
            continue
        qa, ra = query[qs:qs + L], ref[rs:rs + L]
        eq = ((qa == ra) & (qa < 4)).astype(np.int8)
        diff = np.diff(np.concatenate(([0], eq, [0])))
        for s, e in zip(np.flatnonzero(diff == 1), np.flatnonzero(diff == -1)):
            if e - s >= min_len:
                out.add((qs + int(s), qs + int(e), rs + int(s)))
    return out
