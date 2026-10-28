import numpy as np

import pytest



from dsguard.brute import brute_matches

from dsguard.engine import find_matches, scan_both_strands, threshold, window_longest

from dsguard.seq import encode, revcomp





def rand(rng, n):

    return rng.integers(0, 4, n).astype(np.uint8)





@pytest.mark.parametrize("k", [3, 4, 5])

def test_exact_vs_bruteforce(k):

    rng = np.random.default_rng(k)

    for _ in range(300):

        q = rand(rng, int(rng.integers(20, 60)))

        r = rand(rng, int(rng.integers(30, 200)))

        if rng.random() < 0.5:

            L = int(rng.integers(2 * k - 1, 20))

            qs = int(rng.integers(0, q.size - L + 1))

            rs = int(rng.integers(0, r.size - L + 1))

            r[rs:rs + L] = q[qs:qs + L]

        thr = threshold(k)

        got = find_matches(q, r, k=k, chunk=k * 7)

        assert {m for m in got if m[1] - m[0] >= thr} == brute_matches(q, r, thr)

        assert got <= brute_matches(q, r, k)





def test_k11_threshold_21():

    rng = np.random.default_rng(1)

    q = rand(rng, 500)

    for L in (20, 21, 22, 40):

        r = rand(rng, 5000)

        r[1234:1234 + L] = q[100:100 + L]

        found = [m for m in find_matches(q, r, k=11) if m[1] - m[0] >= 21]

        assert bool(found) == (L >= 21)





def test_ambiguity_never_matches():

    q = encode("ACGTACGTACGTACGTACGTACGTAC")

    r = encode("ACGTACGTACGNACGTACGTACGTAC")

    assert all(m[1] - m[0] < 21 for m in find_matches(q, r, k=11))





def test_antisense_mapped_to_sense_coords():

    rng = np.random.default_rng(2)

    q = rand(rng, 300)

    r = rand(rng, 2000)

    r[500:530] = revcomp(q[100:130])

    ms = [m for m in scan_both_strands(q, r) if m.q_end - m.q_start >= 21]

    assert ms and ms[0].strand == "-" and (ms[0].q_start, ms[0].q_end) == (100, 130)





def test_window_longest_clips_to_window():



    best = window_longest([(100, 130)], 1000, 400)

    assert best[0] == 30 and best[100] == 30 and best[110] == 20 and best[129] == 1 and best[130] == 0
