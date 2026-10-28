import numpy as np

from numpy.lib.stride_tricks import sliding_window_view





def sirna_matrix(seq):

    return sliding_window_view(seq, 21)





def reynolds(W):

    g = W[:, :19]

    gc = ((g == 1) | (g == 2)).mean(axis=1)

    au = ((W[:, 14:19] == 0) | (W[:, 14:19] == 3)).sum(axis=1)

    s = ((gc >= 0.30) & (gc <= 0.52)).astype(int)

    s += au >= 3

    s += W[:, 18] == 0

    s += W[:, 2] == 0

    s += W[:, 9] == 3

    s += ~((W[:, 18] == 1) | (W[:, 18] == 2))

    s += W[:, 12] != 2

    return s





def ui_tei(W):

    au = (W == 0) | (W == 3)

    gc = (~au) & (W < 4)

    r1 = au[:, 18]

    r2 = gc[:, 0]

    r3 = au[:, 12:19].sum(axis=1) >= 4

    run = np.zeros(W.shape[0], int)

    best = np.zeros(W.shape[0], int)

    for j in range(19):

        run = np.where(gc[:, j], run + 1, 0)

        best = np.maximum(best, run)

    r4 = best <= 9

    return r1.astype(int) + r2 + r3 + r4





def potency(seq):

    W = sirna_matrix(seq)

    valid = (W < 4).all(axis=1)

    score = (reynolds(W) / 7.0 + ui_tei(W) / 4.0) / 2.0

    return np.where(valid, score, 0.0), valid





def window_mean(per_sirna, valid, L):

    m = L - 20

    nw = per_sirna.size - m + 1

    c = np.concatenate(([0.0], np.cumsum(per_sirna)))

    bad = np.concatenate(([0], np.cumsum(~valid)))

    mean = (c[m:m + nw] - c[:nw]) / m

    nbad = bad[m:m + nw] - bad[:nw]

    return np.where(nbad == 0, mean, np.nan)
