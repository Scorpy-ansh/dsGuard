"""Graded non-target risk.

The binary 21-nt screen is always reported unchanged (it is what a dossier is
assessed against). The graded hazard runs alongside it:

    hazard(m) = 1 / (1 + exp(-(m - L50) / kappa))

PLACEHOLDER DEFAULTS: L50 = 21, kappa = 2 are NOT fitted values. Replace them
with values fitted to published non-target outcome data (and report the
confidence interval) before claiming a result.
"""
import numpy as np

L50_DEFAULT = 21.0
KAPPA_DEFAULT = 2.0


def graded_hazard(longest, l50=L50_DEFAULT, kappa=KAPPA_DEFAULT):
    x = (np.asarray(longest, dtype=float) - l50) / kappa
    return 1.0 / (1.0 + np.exp(-x))
