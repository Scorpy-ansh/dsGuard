import numpy as np



L50_DEFAULT = 21.0

KAPPA_DEFAULT = 2.0





def graded_hazard(longest, l50=L50_DEFAULT, kappa=KAPPA_DEFAULT):

    x = (np.asarray(longest, dtype=float) - l50) / kappa

    return 1.0 / (1.0 + np.exp(-x))
