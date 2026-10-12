import numpy as np

from dsguard.batch import concat_genes, split_matches
from dsguard.engine import scan_both_strands
from dsguard.seq import revcomp


def test_batch_equals_one_by_one():
    """Scanning 4 genes at once must give exactly the per-gene results."""
    rng = np.random.default_rng(5)
    genes = [rng.integers(0, 4, n).astype(np.uint8) for n in (300, 120, 500, 80)]
    ref = rng.integers(0, 4, 20000).astype(np.uint8)
    ref[1000:1030] = genes[0][50:80]
    ref[5000:5025] = revcomp(genes[2][200:225])
    ref[9000:9022] = genes[1][10:32]
    ref[12000:12040] = genes[3][20:60]
    Q, starts = concat_genes(genes)
    batch = split_matches(scan_both_strands(Q, ref, k=11), starts)
    for g, arr in enumerate(genes):
        single = {m for m in scan_both_strands(arr, ref, k=11) if m.q_end - m.q_start >= 21}
        got = {m for m in batch[g] if m.q_end - m.q_start >= 21}
        assert got == single and (single or g == 99)
