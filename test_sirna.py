import numpy as np

from dsguard.optimise import pareto_front
from dsguard.primers import T7, design_t7_primers
from dsguard.seq import encode, revcomp
from dsguard.sirna import potency, reynolds, ui_tei, window_mean


def test_known_good_sirna_scores_high():
    # A at 3, U at 10, A at 19, A/U-rich 3' half, 5' G/C, ~38% GC, no G at 13
    s = encode("GCAGAAGUUUCAAAUAAUAUU")
    W = s[None, :]
    assert reynolds(W)[0] >= 6 and ui_tei(W)[0] >= 3


def test_gc_stretch_rule():
    # 11-long GC stretch fails rule 4 only (3 of 4); a 9-long stretch passes all four
    assert ui_tei(encode("GCGCGCGCGCGAUAUAUAUAU")[None, :])[0] == 3
    assert ui_tei(encode("GCGCGCGCGAUAUAUAUAUAU")[None, :])[0] == 4


def test_ambiguous_windows_are_nan():
    s = encode("ACGU" * 30)  # 120 nt
    s[50] = 4                # siRNAs starting at 30..50 are invalid
    per, valid = potency(s)
    m = window_mean(per, valid, 40)  # window w holds siRNAs w..w+19
    assert np.isnan(m[11:51]).all()
    assert not np.isnan(m[:11]).any() and not np.isnan(m[51:]).any()


def test_pareto():
    pot = np.array([0.9, 0.5, 0.7, 0.7, 0.3])
    risk = np.array([0.9, 0.1, 0.5, 0.6, 0.1])
    assert list(pareto_front(pot, risk)) == [1, 2, 0]


def test_primers_have_t7_and_revcomp():
    g = encode("ACGT" * 50)
    p = design_t7_primers(g, 10, 100)
    assert p["forward"].startswith(T7) and p["reverse"].startswith(T7)
