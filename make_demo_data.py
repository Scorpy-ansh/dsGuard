import os

import sys



import numpy as np



sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from dsguard.seq import decode, revcomp



rng = np.random.default_rng(0)

out = sys.argv[1] if len(sys.argv) > 1 else "demo"

os.makedirs(out, exist_ok=True)





def rand(n):

    return rng.integers(0, 4, n).astype(np.uint8)





def write(path, records):

    with open(path, "w") as fh:

        for name, arr in records:

            s = decode(arr)

            fh.write(f">{name}\n")

            for i in range(0, len(s), 70):

                fh.write(s[i:i + 70] + "\n")





genes = {"gene_A": rand(2000), "gene_B": rand(1800), "gene_C": rand(2000),

         "gene_D": rand(300), "gene_E": rand(2200)}

write(f"{out}/genes.fa", list(genes.items()))

write(f"{out}/target.fa", [("gene_A", genes["gene_A"])])





def organism(plants=(), mb=3_000_000, ntx=300):

    size = mb // ntx

    recs = [[f"tx{i}", rand(size)] for i in range(ntx)]

    for j, plant in enumerate(plants):

        recs[7 + j][1][1000:1000 + plant.size] = plant

    return [(n, s) for n, s in recs]





A, C = genes["gene_A"], genes["gene_C"]

write(f"{out}/crop.fa", organism([C[o:o + 25] for o in range(100, 1950, 300)]))

write(f"{out}/predator.fa", organism([A[600:630]]))

write(f"{out}/pollinator.fa", organism([A[1200:1222]]))

write(f"{out}/human.fa", organism([revcomp(A[300:325])]))

print("demo data written to", out)
