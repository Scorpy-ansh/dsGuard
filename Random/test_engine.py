import numpy as np
import pytest

from dsguard.brute import brute_matches
from dsguard.engine import find_matches, scan_both_strands, threshold, window_longest
from dsguard.seq import encode, revcomp


def rand(rng, n):
    return rng.integers(0, 4, n).astype(np.uint8)


@pytest.mark.parametrize("k", [3, 4, 5])
def test_exact_vs_bruteforce(k):
    """Gate: every match >= 2k-1 found, and nothing reported that isn't real."""
    rng = np.random.default_rng(k)
    for _ in range(300):
        q = rand(rng, int(rng.integers(20, 60)))
        r = rand(rng, int(rng.integers(30, 200)))
        if rng.random() < 0.5:  # plant a long match so long runs are exercised
            L = int(rng.integers(2 * k - 1, 20))
            qs = int(rng.integers(0, q.size - L + 1))
            rs = int(rng.integers(0, r.size - L + 1))
            r[rs:rs + L] = q[qs:qs + L]
        thr = threshold(k)
        got = find_matches(q, r, k=k, chunk=k * 7)  # tiny chunks stress the boundaries
        assert {m for m in got if m[1] - m[0] >= thr} == brute_matches(q, r, thr)
        assert got <= brute_matches(q, r, k)  # soundness for short matches too


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
    # one 30-nt match at [100,130); a 400-window starting at 110 only holds 20 of it
    best = window_longest([(100, 130)], 1000, 400)
    assert best[0] == 30 and best[100] == 30 and best[110] == 20 and best[129] == 1 and best[130] == 0
import pytest

@pytest.mark.parametrize("k", [0, -1, 17, 1.5, True])
def test_invalid_k_rejected(k):
    q = encode("ACGTACGT")
    r = encode("ACGTACGT")

    with pytest.raises(ValueError, match="between 1 and 16"):
        find_matches(q, r, k=k)


@pytest.mark.parametrize("k", [0, -1, 17, 1.5, True])
def test_invalid_threshold_k_rejected(k):
    with pytest.raises(ValueError, match="between 1 and 16"):
        threshold(k)

def test_k11_randomized_against_bruteforce():
    """Independently verify production k=11 matches against brute force."""
    rng = np.random.default_rng(11011)
    k = 11
    min_len = 21

    for _ in range(100):
        q = rand(rng, int(rng.integers(50, 150)))
        r = rand(rng, int(rng.integers(100, 800)))

        # Plant a guaranteed match, including cases near chunk boundaries.
        length = int(rng.integers(min_len, 45))
        q_start = int(rng.integers(0, q.size - length + 1))
        r_start = int(rng.integers(0, r.size - length + 1))
        r[r_start:r_start + length] = q[q_start:q_start + length]

        # Compare all qualifying matches with the independent implementation.
        expected = brute_matches(q, r, min_len)
        got = {
            match for match in find_matches(q, r, k=k, chunk=77)
            if match[1] - match[0] >= min_len
        }

        assert got == expected

@pytest.mark.parametrize("k", [3, 5, 11])
def test_long_matches_are_invariant_across_chunk_sizes(k):
    """Qualifying matches must not depend on reference chunk size."""
    rng = np.random.default_rng(20261023 + k)
    query = rand(rng, 120)
    reference = rand(rng, 900)

    min_len = 2 * k - 1
    q_start = 37
    match_len = max(min_len, 30)
    reference_start = 5 * k + 2

    reference[
        reference_start:reference_start + match_len
    ] = query[q_start:q_start + match_len]

    expected = brute_matches(query, reference, min_len)

    for chunk in (k, 2 * k, 7 * k, 31 * k):
        got = {
            match
            for match in find_matches(
                query, reference, k=k, chunk=chunk
            )
            if match[1] - match[0] >= min_len
        }

        assert got == expected, (
            f"Results differed for k={k}, chunk={chunk}"
        )


@pytest.mark.parametrize(
    "interval, expected",
    [
        ((10, 20), [10, 9, 8, 7, 6]),
        ((15, 25), [5, 6, 7, 8, 9]),
        ((20, 30), [0, 1, 2, 3, 4]),
    ],
)
def test_window_longest_boundary_cases(interval, expected):
    """Check match clipping at window boundaries."""
    best = window_longest([interval], n=30, L=10)
    assert best[10:15].tolist() == expected


def test_window_longest_match_at_sequence_end():
    """A match ending at the sequence boundary is handled correctly."""
    best = window_longest([(20, 30)], n=30, L=10)

    assert best[20] == 10
    assert best[-1] == 10


def test_window_longest_rejects_invalid_window_length():
    """A window longer than the sequence produces no windows."""
    best = window_longest([(2, 5)], n=10, L=11)

    assert best.size == 0
@pytest.mark.parametrize("window_length", [0, -1, -10])
def test_window_longest_rejects_nonpositive_window_length(window_length):
    with pytest.raises(ValueError, match="positive"):
        window_longest([(2, 5)], n=10, L=window_length)

@pytest.mark.parametrize("sequence_length", [-1, -10, 1.5, True])
def test_window_longest_rejects_invalid_sequence_length(sequence_length):
    with pytest.raises(ValueError, match="n must be a nonnegative integer"):
        window_longest([], n=sequence_length, L=5)


@pytest.mark.parametrize(
    "interval",
    [
        (-1, 5),   # starts before the sequence
        (2, 11),   # ends after the sequence
        (5, 2),    # end before start
        (0, -1),   # negative end
    ],
)
def test_window_longest_rejects_invalid_intervals(interval):
    with pytest.raises(ValueError, match="interval"):
        window_longest([interval], n=10, L=5)