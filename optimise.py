import numpy as np


def pareto_front(potency, risk):
    """Indices not dominated by any other (maximise potency, minimise risk),
    sorted by increasing risk. nan potency is ignored."""
    idx = np.flatnonzero(~np.isnan(potency))
    order = idx[np.lexsort((-potency[idx], risk[idx]))]
    keep, best = [], -np.inf
    for i in order:
        if potency[i] > best:
            keep.append(int(i))
            best = potency[i]
    return np.array(keep, dtype=int)


def recommend(potency, hazard, candidates, max_hazard=0.05):
    """Highest-potency candidate with hazard <= max_hazard, else the lowest hazard."""
    c = np.asarray(candidates)
    ok = c[hazard[c] <= max_hazard]
    if ok.size:
        return int(ok[np.argmax(potency[ok])])
    return int(c[np.argmin(hazard[c])])
