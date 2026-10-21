"""dsguard command line. The target FASTA may hold one gene or many (batch mode).

  dsguard --target genes.fa --panel crop=rice.fa --panel human=human.fa --out results/

Outputs:
  summary.csv            one row per gene (ranked by status, then hazard)
  manifest.json          versions, parameters, input hashes, timings
  <gene>/report.md       readable detail for the gene
  <gene>/windows.csv     every candidate window (1-based inclusive coordinates)
  <gene>/pareto.csv      non-dominated windows (potency vs hazard), binary-screen-clean only
  <gene>/hits.csv        every match >= guaranteed threshold, with transcript IDs
  <gene>/profiles.csv    per-position longest covering match per organism (for plots)
  <gene>/primers.json    T7 primers for the recommended window
"""
import argparse
import csv
import datetime
import json
import os
import platform
import re
import sys
import time

import numpy as np

from . import __version__
from .batch import concat_genes, split_matches
from .engine import K, coverage_profile, scan_both_strands, threshold, window_longest
from .optimise import pareto_front, recommend
from .primers import design_t7_primers
from .risk import KAPPA_DEFAULT, L50_DEFAULT, graded_hazard
from .seq import Reference, decode, encode, read_fasta, sha256_file
from .sirna import potency, window_mean

LIMITS = ("dsGuard ranks hazard by perfect sequence similarity only. It does not model "
          "uptake, nuclease stability, RISC loading or dose, and cannot detect miRNA-like "
          "partial-complementarity effects. A low score is necessary, not sufficient, for "
          "non-target safety.")


def _slug(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s)[:60] or "gene"


def _write_csv(path, header, rows, comment=None):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        if comment:
            w.writerow([comment])
        w.writerow(header)
        w.writerows(rows)


def analyse_gene(name, q, orgs, gene_matches, a, thr, outdir, ref_names):
    """gene_matches: {organism: [(Match, transcript_id, transcript_offset), ...]}"""
    n, L = q.size, a.window
    gc = 100.0 * float(np.isin(q, (1, 2)).sum()) / max(1, int((q < 4).sum()))
    hits = [[o, tid, off + 1, m.q_start + 1, m.q_end, m.q_end - m.q_start, m.strand]
            for o in orgs for (m, tid, off) in gene_matches[o] if m.q_end - m.q_start >= thr]
    row = {"gene": name, "length_bp": n, "gc_pct": round(gc, 1)}
    for o in orgs:
        row[f"hits_{o}"] = sum(1 for h in hits if h[0] == o)
    row["total_hits"] = len(hits)

    gdir = os.path.join(outdir, _slug(name))
    os.makedirs(gdir, exist_ok=True)
    _write_csv(os.path.join(gdir, "hits.csv"),
               ["organism", "transcript", "transcript_pos_1based", "gene_start", "gene_end", "length", "strand"], hits)

    if n < L:
        row.update(status=f"too_short(<{L}bp)", windows=0, clean_windows=0)
        _report(gdir, name, row, None, hits, None, a, thr)
        return row

    nw = n - L + 1
    longest = np.array([window_longest([(m.q_start, m.q_end) for (m, _, _) in gene_matches[o]], n, L) for o in orgs])
    cover = [coverage_profile([(m.q_start, m.q_end) for (m, _, _) in gene_matches[o]], n) for o in orgs]
    max_match = longest.max(axis=0)
    hazard = graded_hazard(max_match, a.l50, a.kappa)
    flagged = max_match >= thr
    per_sirna, valid = potency(q)
    pot = window_mean(per_sirna, valid, L)

    _write_csv(os.path.join(gdir, "windows.csv"),
               ["start", "end", "potency"] + [f"longest_{o}" for o in orgs] + ["max_match", "hazard", f"binary_flag_ge{thr}"],
               [[i + 1, i + L, "" if np.isnan(pot[i]) else round(float(pot[i]), 4)] + [int(x) for x in longest[:, i]]
                + [int(max_match[i]), round(float(hazard[i]), 6), int(flagged[i])] for i in range(nw)],
               comment=f"# 1-based inclusive; max_match exact when >= {thr}, otherwise a lower bound")
    _write_csv(os.path.join(gdir, "profiles.csv"), ["position_1based"] + list(orgs),
               [[i + 1] + [int(c[i]) for c in cover] for i in range(n)])

    clean = np.flatnonzero(~flagged & ~np.isnan(pot))
    pool = clean if clean.size else np.flatnonzero(~np.isnan(pot))
    row.update(windows=nw, clean_windows=int(clean.size))
    if pool.size == 0:
        row.update(status="no_valid_window(ambiguous bases)")
        _report(gdir, name, row, None, hits, None, a, thr)
        return row
    front = pool[pareto_front(pot[pool], hazard[pool])]
    _write_csv(os.path.join(gdir, "pareto.csv"), ["start", "end", "potency", "max_match", "hazard"],
               [[i + 1, i + L, round(float(pot[i]), 4), int(max_match[i]), round(float(hazard[i]), 6)] for i in front])
    best = recommend(pot, hazard, front, a.max_hazard)
    primers = design_t7_primers(q, best, L)
    primers.update(window_start_1based=best + 1, window_end_1based=best + L)
    with open(os.path.join(gdir, "primers.json"), "w") as fh:
        json.dump(primers, fh, indent=2)
    threshold_met = bool(hazard[best] <= a.max_hazard)
    recommendation_status = (
        "within_threshold" if threshold_met else "fallback_threshold_exceeded"
    )
    rec = {"start": best + 1, "end": best + L, "potency": round(float(pot[best]), 4),
           "max_match": int(max_match[best]), "hazard": round(float(hazard[best]), 6),
           "recommendation_status": recommendation_status,
           "max_hazard": float(a.max_hazard)}
    row.update(status="ok" if clean.size else f"no_window_free_of_{thr}nt_match",
               recommendation_status=recommendation_status,
               rec_start=rec["start"], rec_end=rec["end"], rec_potency=rec["potency"],
               rec_max_match=rec["max_match"], rec_hazard=rec["hazard"])
    _report(gdir, name, row, rec, hits, primers, a, thr)
    return row


def _report(gdir, name, row, rec, hits, primers, a, thr):
    L = ["# dsGuard report: %s" % name, "",
         f"- Length: {row['length_bp']} bp, GC {row['gc_pct']}%",
         f"- Status: **{row['status']}**",
         f"- Candidate windows ({a.window} bp): {row.get('windows', 0)}; "
         f"free of any >= {thr} nt match to the panel: {row.get('clean_windows', 0)}", ""]
    if rec:
        threshold_met = rec["recommendation_status"] == "within_threshold"
        recommendation_label = (
            "WITHIN CONFIGURED THRESHOLD"
            if threshold_met
            else "FALLBACK — CONFIGURED THRESHOLD EXCEEDED"
        )
        L += ["## Recommended window",
              f"- Recommendation status: **{recommendation_label}**",
              f"- Position {rec['start']}-{rec['end']} (1-based), potency {rec['potency']}",
              f"- Longest panel match inside it: {rec['max_match']} nt (exact if >= {thr}; otherwise a lower bound)",
              f"- Graded hazard: {rec['hazard']}",
              f"- Configured maximum hazard: {rec['max_hazard']}",
              "- Hazard is a heuristic score, not a validated biological probability.",
              f"- T7 forward: `{primers['forward']}`", f"- T7 reverse: `{primers['reverse']}`", ""]
    L += [f"## Non-target matches >= {thr} nt across the whole gene"]
    if hits:
        L += ["", "| organism | transcript | transcript pos | gene pos | length | strand |", "|---|---|---|---|---|---|"]
        L += [f"| {h[0]} | {h[1]} | {h[2]} | {h[3]}-{h[4]} | {h[5]} | {h[6]} |" for h in hits]
    else:
        L += ["", "None found (exact for >= %d nt)." % thr]
    L += ["", "## Method notes",
          f"- Risk model L50={a.l50}, kappa={a.kappa}: " + ("PLACEHOLDER defaults, not fitted to data." if (a.l50, a.kappa) == (L50_DEFAULT, KAPPA_DEFAULT) else "user-supplied."),
          "- Potency (Reynolds + Ui-Tei) was derived from mammalian siRNA data and is not validated for insects.",
          "", "**Limits:** " + LIMITS]
    with open(os.path.join(gdir, "report.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


def main(argv=None):
    p = argparse.ArgumentParser(prog="dsguard", description="dsRNA window design against an organism panel")
    p.add_argument("--target", required=True, help="FASTA of target gene CDS(s); one or many records")
    p.add_argument("--panel", action="append", default=[], metavar="NAME=FASTA",
                   help="non-target transcriptome; repeat for each organism")
    p.add_argument("--window", type=int, default=400, help="dsRNA length in bp (default 400)")
    p.add_argument("--k", type=int, default=K, help="seed length; exact for matches >= 2k-1 (default 11 -> 21 nt)")
    p.add_argument("--l50", type=float, default=L50_DEFAULT, help="PLACEHOLDER risk midpoint (not fitted)")
    p.add_argument("--kappa", type=float, default=KAPPA_DEFAULT, help="PLACEHOLDER risk steepness")
    p.add_argument("--max-hazard", type=float, default=0.05, help="recommend best-potency window under this hazard")
    p.add_argument("--out", default="dsguard_out")
    a = p.parse_args(argv)
    if not a.panel:
        p.error("give at least one --panel NAME=FASTA")
    os.makedirs(a.out, exist_ok=True)
    thr = threshold(a.k)

    names, arrays, seen = [], [], {}
    for nm, raw in read_fasta(a.target):
        seen[nm] = seen.get(nm, 0) + 1
        names.append(nm if seen[nm] == 1 else f"{nm}_{seen[nm]}")
        arrays.append(encode(raw))
    if not arrays:
        p.error("no records in target FASTA")
    Q, starts = concat_genes(arrays)
    print(f"{len(arrays)} gene(s); window {a.window} bp; exact for matches >= {thr} nt (k={a.k})")

    manifest = {"dsguard_version": __version__, "python": sys.version.split()[0], "numpy": np.__version__,
                "platform": platform.platform(),
                "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                "params": {"window": a.window, "k": a.k, "threshold": thr, "l50": a.l50, "kappa": a.kappa,
                           "max_hazard": a.max_hazard, "l50_kappa_fitted": False},
                "target": {"path": a.target, "genes": len(arrays), "sha256": sha256_file(a.target)},
                "panel": {}, "limits": LIMITS}

    orgs = []
    per_gene = [{} for _ in arrays]  # per gene: organism -> [(Match, transcript, offset)]
    for spec in a.panel:
        o, path = spec.split("=", 1)
        t0 = time.time()
        ref = Reference.from_fasta(path)
        by_gene = split_matches(scan_both_strands(Q, ref.seq, a.k), starts)
        dt = time.time() - t0
        orgs.append(o)
        for gi, ms in enumerate(by_gene):
            per_gene[gi][o] = [(m, *ref.locate(m.r_start)) for m in ms]
        manifest["panel"][o] = {"path": path, "sha256": sha256_file(path), "Mb": round(ref.seq.size / 1e6, 2),
                                "seconds": round(dt, 2)}
        print(f"  scanned {o}: {ref.seq.size / 1e6:.1f} Mb in {dt:.1f}s (all genes at once)")

    rows = [analyse_gene(names[i], arrays[i], orgs, per_gene[i], a, thr, a.out, None) for i in range(len(arrays))]
    order = {"ok": 0}
    rows.sort(key=lambda r: (order.get(r["status"], 1), r.get("rec_hazard", 9), -r.get("rec_potency", 0)))
    cols = ["gene", "status", "recommendation_status", "length_bp", "gc_pct",
            "windows", "clean_windows", "rec_start", "rec_end", "rec_potency",
            "rec_max_match", "rec_hazard", "total_hits"] + [f"hits_{o}" for o in orgs]
    _write_csv(os.path.join(a.out, "summary.csv"), cols, [[r.get(c, "") for c in cols] for r in rows])
    manifest["genes"] = names
    with open(os.path.join(a.out, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"\n{'gene':<22}{'status':<28}{'rec window':<14}{'potency':<9}{'maxmatch':<9}{'hits':<5}")
    for r in rows:
        win = f"{r['rec_start']}-{r['rec_end']}" if "rec_start" in r else "-"
        print(f"{r['gene'][:21]:<22}{r['status'][:27]:<28}{win:<14}{r.get('rec_potency', '-')!s:<9}"
              f"{r.get('rec_max_match', '-')!s:<9}{r['total_hits']:<5}")
    print(f"\nResults in {a.out}/ (summary.csv, <gene>/report.md)\nLIMITS: {LIMITS}")
    return 0
