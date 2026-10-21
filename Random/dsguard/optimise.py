
import numpy as np


def _as_1d_numeric(values, name):
    """Convert an input to a finite, one-dimensional numeric array."""
    try:
        arr = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a one-dimensional numeric array") from exc

    if arr.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")

    return arr


def pareto_front(potency, risk):
    """Return non-dominated indices: maximise potency, minimise risk.

    Results are sorted by increasing risk. For identical risk, higher
    potency is considered first. Inputs must be finite 1-D arrays of equal
    length.
    """
    pot = _as_1d_numeric(potency, "potency")
    hazard = _as_1d_numeric(risk, "risk")

    if pot.shape != hazard.shape:
        raise ValueError("potency and risk must have the same length")

    order = np.lexsort((-pot, hazard))
    keep = []
    best_potency = -np.inf

    for i in order:
        if pot[i] > best_potency:
            keep.append(int(i))
            best_potency = pot[i]

    return np.asarray(keep, dtype=int)


def recommend(potency, hazard, candidates, max_hazard=0.05):
    """Choose highest-potency eligible candidate, otherwise lowest risk.

    Candidate indices must be valid, non-empty, unique integers.
    """
    pot = _as_1d_numeric(potency, "potency")
    risk = _as_1d_numeric(hazard, "hazard")

    if pot.shape != risk.shape:
        raise ValueError("potency and hazard must have the same length")

    try:
        raw_candidates = np.asarray(candidates)
    except (TypeError, ValueError) as exc:
        raise ValueError("candidates must be a one-dimensional index array") from exc

    if raw_candidates.ndim != 1:
        raise ValueError("candidates must be one-dimensional")
    if raw_candidates.size == 0:
        raise ValueError("candidates must not be empty")
    if not np.issubdtype(raw_candidates.dtype, np.integer):
        raise ValueError("candidate indices must be integers")

    c = raw_candidates.astype(np.int64, copy=False)

    if np.any(c < 0) or np.any(c >= pot.size):
        raise ValueError("candidate index out of range")
    if np.unique(c).size != c.size:
        raise ValueError("candidate indices must be unique")

    try:
        limit = float(max_hazard)
    except (TypeError, ValueError) as exc:
        raise ValueError("max_hazard must be finite") from exc

    if not np.isfinite(limit):
        raise ValueError("max_hazard must be finite")

    eligible = c[risk[c] <= limit]

    if eligible.size:
        return int(eligible[np.argmax(pot[eligible])])

    return int(c[np.argmin(risk[c])])
