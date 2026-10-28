from collections import namedtuple



import numpy as np



from .seq import revcomp



K = 11

Match = namedtuple("Match", "q_start q_end r_start strand")





def threshold(k=K):

    return 2 * k - 1





def pack_kmers(seq, k=K):

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

    if k > 16:

        raise ValueError("k must be <= 16 for uint32 packing")

    chunk = chunk or k * 1_000_000

    chunk -= chunk % k

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

        rk, rok = rk[::k], rok[::k]

        lo = np.searchsorted(q_sort, rk, "left")

        hi = np.searchsorted(q_sort, rk, "right")

        for gi in np.flatnonzero((hi > lo) & rok):

            r = start + int(gi) * k

            for qi in q_pos[lo[gi]:hi[gi]]:

                found.add(_extend(query, ref, int(qi), r, k))

    return found





def scan_both_strands(query, ref, k=K, chunk=None):

    n = query.size

    out = [Match(a, b, r, "+") for a, b, r in find_matches(query, ref, k, chunk)]

    for a, b, r in find_matches(revcomp(query), ref, k, chunk):

        out.append(Match(n - b, n - a, r, "-"))

    return sorted(out)





def coverage_profile(intervals, n):

    cov = np.zeros(n, np.int32)

    for a, b in intervals:

        np.maximum(cov[a:b], b - a, out=cov[a:b])

    return cov





def window_longest(intervals, n, L):

    nw = n - L + 1

    ws = np.arange(nw)

    best = np.zeros(nw, np.int32)

    for a, b in intervals:

        lo, hi = max(0, a - L + 1), min(nw, b)

        if hi <= lo:

            continue

        w = ws[lo:hi]

        ov = np.minimum(b, w + L) - np.maximum(a, w)

        np.maximum(best[lo:hi], ov, out=best[lo:hi])

    return best
