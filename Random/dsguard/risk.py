
"""Graded non-target risk heuristic.

IMPORTANT:
L50 and kappa are placeholder parameters, not fitted biological values.
This score is not a validated probability of biological harm.
"""

import numpy as np

L50_DEFAULT = 21.0
KAPPA_DEFAULT = 2.0


def graded_hazard(longest, l50=L50_DEFAULT, kappa=KAPPA_DEFAULT):
    """Return heuristic hazard scores from longest exact-match lengths."""

    if not np.isfinite(l50):
        raise ValueError("l50 must be finite")

    if not np.isfinite(kappa) or kappa <= 0:
        raise ValueError("kappa must be finite and greater than zero")

    values = np.asarray(longest, dtype=float)

    if not np.all(np.isfinite(values)):
        raise ValueError("longest must contain only finite values")

    x = (values - l50) / kappa
    result = np.empty_like(x, dtype=float)

    positive = x >= 0
    result[positive] = 1.0 / (1.0 + np.exp(-x[positive]))

    exp_x = np.exp(x[~positive])
    result[~positive] = exp_x / (1.0 + exp_x)

    return float(result) if result.ndim == 0 else result
