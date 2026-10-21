"""T7 primers for a chosen window. Gene-specific 3' part + T7 promoter 5' tail."""
from .seq import decode, revcomp

T7 = "TAATACGACTCACTATAGGG"


def _gc(s):
    return 100.0 * sum(c in "GC" for c in s) / len(s)


def _tm(s):  # basic Tm for ~20-mers; screening only, check with Primer3 later
    return 64.9 + 41.0 * (sum(c in "GC" for c in s) - 16.4) / len(s)


def design_t7_primers(gene, start, length, plen=20):
    g = decode(gene[start:start + length])
    fs, rs = g[:plen], decode(revcomp(gene[start + length - plen:start + length]))
    return {
        "amplicon_bp": length,
        "forward": T7 + fs, "forward_gene_specific_gc": round(_gc(fs), 1), "forward_tm": round(_tm(fs), 1),
        "reverse": T7 + rs, "reverse_gene_specific_gc": round(_gc(rs), 1), "reverse_tm": round(_tm(rs), 1),
        "note": "Basic GC/Tm only. Check dimers/hairpins and amplicon uniqueness before ordering.",
    }
