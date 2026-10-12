"""On-target potency from published siRNA design rules (sense-strand 21-mers).

CAVEAT (must stay in the docs): these rules were derived from synthetic siRNA in
mammalian cells. Whether they predict insect Dicer-2 processing of long dsRNA is
unproven. Treat potency as a weak prior and test it against your lab data.
"""
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view


def sirna_matrix(seq):
    """(n-20, 21) view of every 21-mer on the sense strand."""
    return sliding_window_view(seq, 21)


def reynolds(W):
    """Reynolds 2004, 7 of the 8 criteria (0-7). 'No internal repeats' is NOT
    implemented yet -- see TODO. Positions are 1-based in the paper, 0-based here."""
    g = W[:, :19]
    gc = ((g == 1) | (g == 2)).mean(axis=1)
    au = ((W[:, 14:19] == 0) | (W[:, 14:19] == 3)).sum(axis=1)
    s = ((gc >= 0.30) & (gc <= 0.52)).astype(int)
    s += au >= 3                      # >=3 A/U at positions 15-19
    s += W[:, 18] == 0                # A at 19
    s += W[:, 2] == 0                 # A at 3
    s += W[:, 9] == 3                 # U at 10
    s += ~((W[:, 18] == 1) | (W[:, 18] == 2))  # not G/C at 19
    s += W[:, 12] != 2                # not G at 13
    return s  # TODO: internal-repeat criterion


def ui_tei(W):
    """Ui-Tei 2004, four rules (0-4). Antisense 5' end pairs with sense position 19."""
    au = (W == 0) | (W == 3)
    gc = (~au) & (W < 4)
    r1 = au[:, 18]                    # A/U at antisense 5' end
    r2 = gc[:, 0]                     # G/C at sense 5' end
    r3 = au[:, 12:19].sum(axis=1) >= 4  # >=4 A/U in antisense 5' 7 bp
    run = np.zeros(W.shape[0], int)
    best = np.zeros(W.shape[0], int)
    for j in range(19):
        run = np.where(gc[:, j], run + 1, 0)
        best = np.maximum(best, run)
    r4 = best <= 9                    # no GC stretch longer than 9
    return r1.astype(int) + r2 + r3 + r4


def potency(seq):
    """Composite potency in [0, 1] per sense-strand 21-mer, plus validity mask."""
    W = sirna_matrix(seq)
    valid = (W < 4).all(axis=1)
    score = (reynolds(W) / 7.0 + ui_tei(W) / 4.0) / 2.0
    return np.where(valid, score, 0.0), valid


def window_mean(per_sirna, valid, L):
    """Mean potency of the L-20 siRNAs in every window of length L (nan if the
    window contains an ambiguous base)."""
    m = L - 20
    nw = per_sirna.size - m + 1
    c = np.concatenate(([0.0], np.cumsum(per_sirna)))
    bad = np.concatenate(([0], np.cumsum(~valid)))
    mean = (c[m:m + nw] - c[:nw]) / m
    nbad = bad[m:m + nw] - bad[:nw]
    return np.where(nbad == 0, mean, np.nan)
