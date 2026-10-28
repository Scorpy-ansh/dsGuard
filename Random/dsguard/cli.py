import argparse

import csv

import datetime

import json

import math

import os

import platform

import re

import sys

import time

import uuid



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







CANDIDATE_LIMITS = ("The reported candidate sequence is a real window copied from the "

                    "target CDS: selecting it is not evidence of biological efficacy, "

                    "uptake, stability or non-target safety, and no experimental "

                    "validation is implied by this tool.")



SCORING = {

    "hazard": {

        "formula": "hazard(m) = 1 / (1 + exp(-(m - L50) / kappa))",

        "input": "m = longest perfect (gapless) match, in nt, between the candidate "

                 "window and any reference record on either strand; m is clipped to the window.",

        "range": "[0, 1]",

        "meaning": "Monotone heuristic that rises towards 1 as the longest perfect match "

                   "grows. It is NOT a calibrated probability of biological harm.",

        "calibrated": False,

        "caveat": "L50/kappa are placeholder constants, not values fitted to published "

                  "non-target outcome data. No calibration or validation evidence is "

                  "bundled with this build.",

    },

    "potency": {

        "formula": "potency = mean over the sense-strand 21-mers in the window of "

                   "0.5 * (reynolds/7 + ui_tei/4)",

        "range": "[0, 1]",

        "meaning": "Composite of published mammalian siRNA design rules (Reynolds 2004, "

                   "7 of 8 criteria; Ui-Tei 2004, 4 rules). Used as a weak ordering prior "

                   "only.",

        "calibrated": False,

        "caveat": "Derived from synthetic siRNA in mammalian cells; not validated for "

                  "insect Dicer-2 processing of long dsRNA.",

    },

    "matching_method": {

        "name": "strided seed-and-extend exact match",

        "description": "A k-mer seed is placed on every k-th reference position; seeds "

                       "are extended in both directions while bases are identical. Only "

                       "perfect, gapless matches are reported, so identity is 100% and "

                       "mismatch/gap counts are 0 for every reported row. No gapped or "

                       "approximate alignment is performed.",

        "completeness": "Every match of length >= 2k-1 is guaranteed to be found. Shorter "

                        "matches are best-effort and the enumeration may be incomplete; "

                        "for those rows the reported length is a lower bound.",

        "reverse_complement": "Both strands of the candidate are screened; antisense hits "

                              "are mapped back onto sense coordinates.",

        "ambiguous_bases": "Any base outside A/C/G/T/U (including N) is treated as a "

                           "separate symbol that can never match, so matches cannot span "

                           "ambiguous positions.",

    },

}









CATEGORY_BANDS = ((0.5, "Critical"), (0.25, "High"), (0.05, "Moderate"), (0.0, "Low"))









PANEL_ACCESSION_SAMPLE = 25





def hazard_equivalent_match(max_hazard, l50=L50_DEFAULT, kappa=KAPPA_DEFAULT):

    p = float(max_hazard)

    if not 0.0 < p < 1.0:

        return None

    return float(l50) + float(kappa) * math.log(p / (1.0 - p))





def risk_category(hazard):

    h = float(hazard)

    for floor, name in CATEGORY_BANDS:

        if h >= floor:

            return name

    return "Low"





def threshold_status(hazard, max_hazard):

    if hazard is None or not math.isfinite(float(hazard)):

        return "not_assessable"

    return "meets_threshold" if float(hazard) <= float(max_hazard) else "exceeds_threshold"





def panel_provenance(name, path, ref, seconds):

    sample = [str(n) for n in list(ref.names)[:PANEL_ACCESSION_SAMPLE]]

    return {

        "name": name,

        "path": str(path),

        "sha256": sha256_file(path),

        "Mb": round(ref.seq.size / 1e6, 2),

        "records": len(ref.names),

        "seconds": round(seconds, 2),

        "records_without_annotation": ref.records_without_annotation,

        "record_accessions": sample,

        "record_accessions_omitted": max(0, len(ref.names) - len(sample)),

        "annotation_example": (ref.descriptions[0] if ref.descriptions else ""),

    }





def _slug(s):

    return re.sub(r"[^A-Za-z0-9._-]+", "_", s)[:60] or "gene"





def _write_csv(path, header, rows, comment=None):

    with open(path, "w", newline="", encoding="utf-8") as fh:

        w = csv.writer(fh)

        if comment:

            w.writerow([comment])

        w.writerow(header)

        w.writerows(rows)





def _json_default(value):

    if isinstance(value, np.integer):

        return int(value)

    if isinstance(value, np.floating):

        return float(value)

    if isinstance(value, np.bool_):

        return bool(value)

    if isinstance(value, np.ndarray):

        return value.tolist()

    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")





def _write_json(path, payload):

    with open(path, "w", encoding="utf-8") as fh:

        json.dump(payload, fh, indent=2, ensure_ascii=False, sort_keys=False,

                  default=_json_default)

        fh.write("\n")













OFFTARGET_COLUMNS = [

    "rank",

    "target_id",

    "target_window",

    "matched_organism",

    "transcript_or_reference_id",

    "reference_panel",

    "transcript_pos_1based",

    "target_start",

    "target_end",

    "match_length_bp",

    "identity_pct",

    "mismatches",

    "gaps",

    "strand",

    "potency",

    "hazard",

    "risk_category",

    "threshold_status",

    "hazard_is_exact",

    "evidence",

    "in_recommended_window",

    "governs_recommendation",

]





def rank_offtargets(rows):

    seen = set()

    unique = []

    for row in rows:

        key = tuple(row.get(column) for column in OFFTARGET_COLUMNS if column != "rank")

        if key in seen:

            continue

        seen.add(key)

        unique.append(row)

    unique.sort(

        key=lambda row: (

            -float(row["hazard"]),

            -int(row["match_length_bp"]),

            str(row["reference_panel"]),

            str(row["matched_organism"]),

            str(row["transcript_or_reference_id"]),

            int(row["target_start"]),

            str(row["strand"]),

        )

    )

    for position, row in enumerate(unique, start=1):

        row["rank"] = position

    return unique





def _window_potency(per_window, a, b, n, L):

    nw = per_window.size

    if nw <= 0 or b <= a:

        return None

    lo = max(0, b - L)

    hi = min(nw - 1, a)

    if hi < lo:

        return None

    block = per_window[lo:hi + 1]

    if block.size == 0 or np.isnan(block).all():

        return None

    return round(float(np.nanmax(block)), 4)





def build_offtarget_rows(gene, orgs, gene_matches, thr, max_hazard, pot, best, window_length, n,

                        min_report_match, rec_max_match=None):

    rows = []

    total = 0

    below_floor = 0

    for organism in orgs:

        for match, transcript_id, offset in gene_matches[organism]:

            total += 1

            length = int(match.q_end - match.q_start)

            inside_recommended = (best is not None

                                  and best <= match.q_start

                                  and match.q_end <= best + window_length)

            governs_recommendation = bool(

                inside_recommended and rec_max_match is not None and length == rec_max_match)

            if length < min_report_match and not governs_recommendation:

                below_floor += 1

                continue

            hazard = float(graded_hazard(length))

            exact = length >= thr

            potency_value = _window_potency(pot, int(match.q_start), int(match.q_end),

                                            n, window_length)

            if exact:

                evidence = ("confirmed exact match (enumeration is complete for matches "

                            f">= {thr} nt)")

            elif governs_recommendation:

                evidence = ("best-effort exact match; it sets the hazard of the recommended "

                            "window, but the full enumeration is only guaranteed for matches "

                            f">= {thr} nt, so its length is a lower bound")

            else:

                evidence = ("best-effort exact match (shorter than the guaranteed "

                            f"{thr} nt detection threshold; length is a lower bound)")

            rows.append({

                "target_id": gene,

                "target_window": f"{int(match.q_start) + 1}-{int(match.q_end)}",

                "matched_organism": organism,

                "transcript_or_reference_id": transcript_id,

                "reference_panel": organism,

                "transcript_pos_1based": int(offset) + 1,

                "target_start": int(match.q_start) + 1,

                "target_end": int(match.q_end),

                "match_length_bp": length,

                "identity_pct": 100.0,

                "mismatches": 0,

                "gaps": 0,

                "strand": match.strand,

                "potency": potency_value,

                "hazard": round(hazard, 6),

                "risk_category": risk_category(hazard),

                "threshold_status": threshold_status(hazard, max_hazard),

                "hazard_is_exact": bool(exact),

                "evidence": evidence,

                "in_recommended_window": int(inside_recommended),

                "governs_recommendation": governs_recommendation,

            })

    return rank_offtargets(rows), below_floor, total













def preventive_measures(gene_results, panels_searched, panels_skipped, max_hazard,

                        equivalent_match, invalid_bases=0, thr=None,

                        input_type=None, input_label=""):

    findings = []



    def add(finding, why, action, priority=None, evidence=""):

        entry = {"finding": finding, "why": why, "action": action,

                 "evidence": evidence}

        if priority is not None:

            entry["priority"] = priority

        findings.append(entry)





    skipped = []

    for panel in panels_skipped or []:

        if isinstance(panel, dict):

            skipped.append({"name": panel.get("name", "?"),

                            "detail": panel.get("detail") or "panel was not available"})

        else:

            skipped.append({"name": str(panel), "detail": "panel was not available"})



    if thr is not None and equivalent_match is not None and equivalent_match < thr:

        add(

            f"The configured maximum hazard of {max_hazard} corresponds to a perfect "

            f"match of only ~{equivalent_match:.2f} nt, which is shorter than the "

            f"{thr} nt guaranteed-detection threshold.",

            "At 15-17 nt, a perfect match against a genome-scale reference is expected by "

            "chance (roughly 2*n*R/4**m occurrences for a query of n nt against R nt of "

            "reference). A strict interpretation of this threshold will therefore often "

            "report 'exceeds threshold' even for a target with no real homology. This is a "

            "property of the placeholder L50/kappa and the chosen maximum hazard, not "

            "evidence of biological off-target risk.",

            "Decide the maximum hazard deliberately. If a 0.05 cap is required, fit "

            "L50/kappa to real non-target outcome data first; otherwise choose a "

            "maximum hazard whose equivalent match length is at or above the "

            "guaranteed-detection threshold, and record the choice. Do NOT simply "

            "raise the cap to make results look acceptable.",

            "High",

            f"max_hazard={max_hazard}, equivalent perfect match ~{equivalent_match:.2f} nt, "

            f"guaranteed detection threshold={thr} nt (2k-1)")



    if not panels_searched:

        add(

            "No reference panel was searched.",

            "Off-target similarity cannot be assessed at all without a reference set.",

            "Install at least one curated non-target reference FASTA and re-run.",

            "Critical",

            "reference_panels.searched is empty for this run")

    else:

        coverage = "; ".join(

            f"{p.get('name')} ({p.get('records')} records, {p.get('Mb')} Mb, "

            f"sha256 {str(p.get('sha256') or '')[:12]}\u2026)"

            for p in panels_searched

        )

        add(

            f"Reference coverage is limited to the {len(panels_searched)} searched "

            f"panel(s): {coverage}.",

            "Only these FASTA files were compared against the target. Any organism, "

            "transcript, isoform or variant that is not present in them was never "

            "searched, so no match being reported is not evidence that no match exists "

            "in the wider environment.",

            "Confirm the searched panels cover every organism relevant to this "

            "application and add the missing FASTA files (recording their accessions "

            "and versions) before relying on a negative result.",

            "Info",

            coverage)



    for panel in skipped:

        add(

            f"Reference panel not installed: {panel['name']}.",

            f"The candidate was never screened against this organism, so any match to it "

            f"is invisible in these results. Panel status: {panel['detail']}.",

            f"Install the FASTA for '{panel['name']}' (see the README) and re-run. Until "

            "then, treat coverage for that organism as an unresolved risk.",

            "High",

            f"reference_panels.skipped: name={panel['name']}, detail={panel['detail']}")



    under_annotated = [

        p for p in panels_searched if p.get("records_without_annotation")

    ]

    for panel in under_annotated:

        add(

            f"Reference record(s) without annotation in panel "

            f"'{panel.get('name')}': {panel.get('records_without_annotation')} of "

            f"{panel.get('records')} record(s) have only an identifier in their FASTA "

            f"header.",

            "A match to such a record can be located and scored, but the record cannot "

            "be confirmed from the results alone to be the intended organism or "

            "transcript, so the match cannot be attributed with confidence to a "

            "sequence whose biology is known.",

            "Check the affected FASTA headers against the source database and add the "

            "accession and annotation, or re-download the panel with complete headers, "

            "before relying on matches to those records.",

            "Medium",

            f"panel '{panel.get('name')}': records={panel.get('records')}, "

            f"records_without_annotation={panel.get('records_without_annotation')}, "

            f"path={panel.get('path')}")



    if invalid_bases:

        add(

            f"{invalid_bases} ambiguous/invalid base(s) in the target sequence.",

            "Ambiguous positions can never match, so any match spanning them is "

            "missed, and windows covering them are excluded from candidacy.",

            "Confirm the input sequence is the intended CDS and replace ambiguous "

            "bases before relying on the result.",

            "Medium",

            f"input.ambiguous_bases={invalid_bases} (target sequence)")



    if str(input_type or "").lower() == "ncbi":

        add(

            f"The target sequence was retrieved from a sequence database at run time "

            f"({input_type}: {input_label}).",

            "Database records can be updated, re-annotated or withdrawn, and the "

            "screening is only valid for the exact sequence that was fetched. A failed "

            "or ambiguous lookup aborts the request before any screening, so no results "

            "are produced for it and none are fabricated.",

            "Record the accession and version shown here with the results, and re-run "

            "the screening if the database record changes.",

            "Info",

            f"input.type={input_type}, input.label={input_label}")



    for gene in gene_results:

        row = gene["row"]

        name = gene["row"]["gene"]

        if row.get("status", "").startswith("too_short"):

            add(

                f"{name}: sequence shorter than one candidate window.",

                "No {0} bp window can be placed, so no candidate was evaluated."

                        .format(gene.get("window", "")),

                "Provide the full-length CDS or reduce the window length and re-run.",

                "Medium",

                f"{name}: length_bp={row.get('length_bp')}, window={gene.get('window')}")

            continue



        confirmed = [m for m in gene["ranked_offtargets"] if m["hazard_is_exact"]]

        if confirmed:

            worst = confirmed[0]

            add(

                f"{name}: confirmed exact match of {worst['match_length_bp']} nt to "

                f"{worst['matched_organism']} reference "

                f"{worst['transcript_or_reference_id']} "

                f"(target {worst['target_start']}-{worst['target_end']}, "

                f"strand {worst['strand']}).",

                "Matches at or above the guaranteed detection threshold are the "

                "criterion a non-target dossier is assessed against. This is an "

                "observed perfect sequence match; the hazard score derived from it is a "

                "heuristic and does not by itself demonstrate a biological effect.",

                "Avoid candidate windows containing this match; confirm the "

                "matched reference record and its annotation before discarding "

                "or accepting the target.",

                "Critical",

                f"{name}: match_length_bp={worst['match_length_bp']} "

                f"(>= guaranteed threshold), reference={worst['transcript_or_reference_id']}, "

                f"target={worst['target_start']}-{worst['target_end']}, "

                f"strand={worst['strand']}, hazard={worst['hazard']}")



        if row.get("recommendation_status") in {

            "fallback_threshold_exceeded",

            "within_threshold",

        }:

            add(

                f"{name}: selected candidate hazard is "

                f"{'above' if row.get('recommendation_status') == 'fallback_threshold_exceeded' else 'at or below'} "

                f"the configured reference value of {max_hazard}.",

                "This is an informational comparison only; candidate selection is "

                "based on the documented potency/hazard ordering and not a biological "

                "safety verdict.",

                "Review the match evidence and reference-panel coverage before "

                "making any biological decision.",

                "Medium",

                f"{name}: rec window {row.get('rec_start')}-{row.get('rec_end')}, "

                f"rec_hazard={row.get('rec_hazard')}, max_hazard={max_hazard}, "

                f"rec_max_match={row.get('rec_max_match')} nt")



    categories = {

        "before_accepting_a_candidate": [

            "Investigate the highest-ranked predicted matches first; they are ordered by the "

            "documented hazard score.",

            "Confirm the matched organism and the exact transcript/accession identifier in the "

            "reference FASTA header.",

            "Verify the reference sequence and its annotation actually correspond to that "

            "organism and transcript.",

            "Review match length, identity and the fact that only perfect, gapless matches are "

            "reported by this build.",

            "Check whether every organism you care about is present in the installed reference "

            "panels; missing panels are listed above and are unresolved risks.",

        ],

        "reducing_potential_off_target_risk": [

            "Re-assess any candidate with a concerning predicted match using an independent, "

            "validated method (e.g. a gapped aligner or a curated off-target database).",

            "Compare alternative candidate windows using the same reference set and the same "

            "scoring parameters so the comparison is fair.",

            "Treat the configured hazard comparison as an informational sequence-screening "

            "signal, not as validated biological safety or efficacy.",

            "Where a real-world decision depends on the result, obtain specialist review or "

            "independent experimental validation before use.",

            "Add the reference organisms/transcripts that are missing from the panels to "

            "improve coverage.",

        ],

        "validation_and_quality_assurance": [

            "Confirm the reference datasets are current and traceable (record accession and "

            "version for each panel).",

            "Re-run the analysis after any change to the target sequence, parameters, or "

            "reference panels; results are only valid for the run recorded here.",

            "Where available, compare against positive and negative controls.",

            "Record the software version, scoring configuration, reference checksums and the "

            "analysis timestamp (all stored in this run's manifest.json / results.json).",

            "Keep the three levels distinct: an observed perfect sequence match is a fact about "

            "the sequences compared; the hazard score is a documented heuristic derived from "

            "that match length against placeholder L50/kappa; a biological off-target effect "

            "would require experimental validation, which this tool does not perform.",

        ],

    }

    return {"findings": findings, "categories": categories}













def candidate_from_window(gene, seq, start, requested_length):

    if isinstance(requested_length, bool) or not isinstance(

            requested_length, (int, np.integer)) or requested_length <= 0:

        raise ValueError("requested length must be a positive integer number of bp")

    if isinstance(start, bool) or not isinstance(start, (int, np.integer)) or start < 0:

        raise ValueError("window start must be a non-negative integer index")



    requested = int(requested_length)

    first = int(start)

    n = int(np.asarray(seq).size)



    if requested > n:

        raise ValueError(

            f"target '{gene}' is {n} nt, shorter than the requested {requested} bp "

            "window; the engine does not pad or truncate the target, so no candidate "

            "can be produced at this length"

        )

    if first + requested > n:

        raise ValueError(

            f"window [{first + 1}, {first + requested}] (1-based) runs past the end of "

            f"the {n} nt target '{gene}'"

        )



    window = np.asarray(seq)[first:first + requested]

    ambiguous = int((window > 3).sum())

    if ambiguous:

        raise ValueError(

            f"window [{first + 1}, {first + requested}] of '{gene}' contains {ambiguous} "

            "ambiguous base(s); such windows are excluded from candidacy and the engine "

            "does not substitute bases, so no candidate can be reported"

        )



    sequence = decode(window)

    return {

        "candidate_id": f"{_slug(gene)}:{first + 1}-{first + len(sequence)}",

        "gene": gene,

        "requested_length_bp": requested,

        "actual_length_bp": len(sequence),

        "sequence": sequence,

        "start_1based": first + 1,

        "end_1based": first + len(sequence),

        "strand": "sense",

        "ambiguous_bases": 0,

    }













def analyse_gene(name, q, orgs, gene_matches, a, thr, outdir, analysis_meta):

    n, L = q.size, a.window

    gc = 100.0 * float(np.isin(q, (1, 2)).sum()) / max(1, int((q < 4).sum()))



    all_matches = [(m, tid, off, o) for o in orgs for (m, tid, off) in gene_matches[o]]

    hits = [[o, tid, off + 1, m.q_start + 1, m.q_end, m.q_end - m.q_start, m.strand]

            for (m, tid, off, o) in all_matches if m.q_end - m.q_start >= thr]

    row = {"gene": name, "length_bp": n, "gc_pct": round(gc, 1)}

    for o in orgs:

        row[f"hits_{o}"] = sum(1 for h in hits if h[0] == o)

    row["total_hits"] = len(hits)

    row["max_hazard"] = float(a.max_hazard)

    row["hazard_threshold_equivalent_match_bp"] = (

        None if a.hazard_equivalent_match is None

        else round(float(a.hazard_equivalent_match), 2)

    )



    gdir = os.path.join(outdir, _slug(name))

    os.makedirs(gdir, exist_ok=True)

    _write_csv(os.path.join(gdir, "hits.csv"),

               ["organism", "transcript", "transcript_pos_1based", "gene_start",

                "gene_end", "length", "strand"], hits)



    result = {

        "row": row,

        "ranked_offtargets": [],

        "pareto": [],

        "recommendation": None,

        "primers": None,

        "candidate": None,

        "candidate_error": None,

        "window": L,

        "total_matches_found": len(all_matches),

        "matches_below_report_floor": 0,

        "screening_status": "screening_incomplete",

        "selection_policy": "pareto_highest_potency_then_lowest_hazard",

    }



    if n < L:

        row.update(status=f"too_short(<{L}bp)",

                   screening_status="screened_no_candidate",

                   binary_status="not_evaluated",

                   windows=0, clean_windows=0, offtargets_reported=0,

                   matches_below_report_floor=0)





        result["candidate_error"] = (

            f"target '{name}' is {n} nt, shorter than the requested {L} bp window. The "

            "engine places a full-length window against the target and never pads, "

            "truncates or repeats bases, so no candidate sequence can be produced at "

            "this length. Provide a longer target or select a shorter window."

        )

        result["screening_status"] = "screened_no_candidate"

        _write_csv(os.path.join(gdir, "offtarget_ranked.csv"), OFFTARGET_COLUMNS, [])

        _report(gdir, name, row, None, hits, None, a, thr, [], None)

        return result



    nw = n - L + 1

    longest = np.array([

        window_longest([(m.q_start, m.q_end) for (m, _, _) in gene_matches[o]], n, L)

        for o in orgs

    ])

    cover = [coverage_profile([(m.q_start, m.q_end) for (m, _, _) in gene_matches[o]], n)

             for o in orgs]

    max_match = longest.max(axis=0)

    hazard = graded_hazard(max_match, a.l50, a.kappa)

    flagged = max_match >= thr

    per_sirna, valid = potency(q)

    pot = window_mean(per_sirna, valid, L)



    _write_csv(os.path.join(gdir, "windows.csv"),

               ["start", "end", "potency"] + [f"longest_{o}" for o in orgs]

               + ["max_match", "hazard", f"binary_flag_ge{thr}"],

               [[i + 1, i + L, "" if np.isnan(pot[i]) else round(float(pot[i]), 4)]

                + [int(x) for x in longest[:, i]]

                + [int(max_match[i]), round(float(hazard[i]), 6), int(flagged[i])]

                for i in range(nw)],

               comment=(f"# 1-based inclusive; max_match is exact when >= {thr}, "

                        f"otherwise a lower bound; hazard is a heuristic, not a probability"))

    _write_csv(os.path.join(gdir, "profiles.csv"), ["position_1based"] + list(orgs),

               [[i + 1] + [int(c[i]) for c in cover] for i in range(n)])



    clean = np.flatnonzero(~flagged & ~np.isnan(pot))

    pool = clean if clean.size else np.flatnonzero(~np.isnan(pot))

    row.update(windows=nw, clean_windows=int(clean.size),

               binary_status=("all_windows_clean" if clean.size == nw

                              else ("no_window_clean" if clean.size == 0 else "partially_clean")))

    if pool.size == 0:

        row.update(status="no_valid_window(ambiguous bases)",

                   screening_status="screened_no_candidate",

                   recommendation_status="not_assessable", offtargets_reported=0,

                   matches_below_report_floor=0)



        result["candidate_error"] = (

            f"none of the {nw} candidate window(s) of {L} bp is free of an ambiguous "

            f"base, so no window in '{name}' is eligible. The engine excludes such "

            "windows from candidacy instead of substituting bases, so no candidate "

            "sequence can be produced."

        )

        result["screening_status"] = "screened_no_candidate"

        _write_csv(os.path.join(gdir, "offtarget_ranked.csv"), OFFTARGET_COLUMNS, [])

        _report(gdir, name, row, None, hits, None, a, thr, [], None)

        return result



    front = pool[pareto_front(pot[pool], hazard[pool])]

    pareto_rows = [[i + 1, i + L, round(float(pot[i]), 4), int(max_match[i]),

                    round(float(hazard[i]), 6),

                    "meets_threshold" if float(hazard[i]) <= a.max_hazard else "exceeds_threshold",

                    int(max_match[i] >= thr)]

                   for i in front]

    _write_csv(os.path.join(gdir, "pareto.csv"),

               ["start", "end", "potency", "max_match", "hazard", "threshold_status",

                "hazard_is_exact"], pareto_rows)



    best = int(recommend(pot, hazard, front, a.max_hazard))









    candidate = candidate_from_window(name, q, best, L)

    primers = design_t7_primers(q, best, L)

    primers.update(window_start_1based=best + 1, window_end_1based=best + L)

    _write_json(os.path.join(gdir, "primers.json"), primers)

    threshold_met = bool(hazard[best] <= a.max_hazard)

    recommendation_status = ("within_threshold" if threshold_met

                             else "fallback_threshold_exceeded")

    rec = {"start": best + 1, "end": best + L, "potency": round(float(pot[best]), 4),

           "max_match": int(max_match[best]),

           "max_match_is_exact": bool(int(max_match[best]) >= thr),

           "hazard": round(float(hazard[best]), 6),

           "hazard_is_calibrated": False,

           "recommendation_status": recommendation_status,

           "threshold_status": threshold_status(float(hazard[best]), a.max_hazard),

           "risk_category": risk_category(float(hazard[best])),

           "max_hazard": float(a.max_hazard),

           "hazard_threshold_equivalent_match_bp": row["hazard_threshold_equivalent_match_bp"],





           "is_fallback": False,

           "screening_status": "screened_candidate_available",

           "selection_policy": "pareto_highest_potency_then_lowest_hazard",

           "candidate": candidate,

           "t7_forward": primers["forward"],

           "t7_reverse": primers["reverse"]}

    row.update(status="ok" if clean.size else f"no_window_free_of_{thr}nt_match",

               screening_status="screened_candidate_available",

               recommendation_status=recommendation_status,

               rec_start=rec["start"], rec_end=rec["end"], rec_potency=rec["potency"],

               rec_max_match=rec["max_match"], rec_hazard=rec["hazard"])



    ranked, below_floor, total_matches = build_offtarget_rows(

        name, orgs, gene_matches, thr, a.max_hazard, pot, best, L, n,

        a.min_report_match, rec_max_match=int(max_match[best]))

    _write_csv(os.path.join(gdir, "offtarget_ranked.csv"), OFFTARGET_COLUMNS,

               [[r[c] for c in OFFTARGET_COLUMNS] for r in ranked],

               comment=(f"# Ranked highest off-target hazard first. Only perfect gapless matches "

                        f">= {a.min_report_match} nt are listed, so identity_pct is 100 and "

                        f"mismatches/gaps are 0. Matches below {a.min_report_match} nt are "

                        f"ordinary random k-mer collisions against a large reference and are "

                        f"counted in matches_below_report_floor, not listed here; the longest "

                        f"match inside the recommended window is always listed. hazard is a "

                        f"heuristic score, not a probability. potency is the best on-target "

                        f"potency among candidate windows that fully contain the match "

                        f"(blank when no window can contain it)."))

    row["offtargets_reported"] = len(ranked)

    row["matches_below_report_floor"] = below_floor



    result.update({

        "ranked_offtargets": ranked,

        "matches_below_report_floor": below_floor,

        "total_matches_found": total_matches,

        "min_report_match": int(a.min_report_match),

        "pareto": [dict(zip(["start", "end", "potency", "max_match", "hazard",

                             "threshold_status", "hazard_is_exact"], values))

                   for values in pareto_rows],

        "recommendation": rec,

        "primers": primers,

        "candidate": candidate,

        "screening_status": "screened_candidate_available",

        "selection_policy": "pareto_highest_potency_then_lowest_hazard",

    })

    _report(gdir, name, row, rec, hits, primers, a, thr, ranked, analysis_meta)

    return result





def _report(gdir, name, row, rec, hits, primers, a, thr, ranked, analysis_meta):

    lines = ["# dsGuard report: %s" % name, ""]

    if analysis_meta:

        lines += [

            f"- Project: dsGuard {analysis_meta.get('dsguard_version', '')}",

            f"- Analysis ID: {analysis_meta.get('analysis_id', 'n/a')}",

            f"- Generated (UTC): {analysis_meta.get('timestamp_utc', 'n/a')}",

            f"- Input type: {analysis_meta.get('input_type', 'n/a')}",

            f"- Input identifier: {analysis_meta.get('input_label', 'n/a')}",

            "",

        ]

    lines += [

        f"- Length: {row['length_bp']} bp, GC {row['gc_pct']}%",

        f"- Screening status: **{row.get('screening_status', 'screening_incomplete')}**",

        f"- Candidate windows ({a.window} bp): {row.get('windows', 0)}; "

        f"free of any >= {thr} nt perfect match to the panel: {row.get('clean_windows', 0)}",

        f"- Windows: {a.window}; k-mer: {a.k} (exact for matches >= {thr} nt); "

        f"L50: {a.l50}; kappa: {a.kappa}; maximum hazard: {a.max_hazard}",

        "",

    ]

    if analysis_meta and analysis_meta.get("panels"):

        lines += ["## Reference panels", ""]

        for panel in analysis_meta["panels"]:

            lines.append(f"- {panel['name']}: {panel['status']} ({panel.get('detail', '')})")

        lines.append("")

    if rec:

        if rec["max_match_is_exact"]:

            match_note = "exact; enumeration is complete at or above this length"

        else:

            match_note = (f"a lower bound; the true longest match may be longer, but any "

                          f"match of {thr} nt or more would have been found")

        lines += [

            "## Recommended window",

            f"- Screening status: **{rec.get('screening_status', 'screened_candidate_available')}**",

            f"- Position {rec['start']}-{rec['end']} (1-based), potency {rec['potency']}",

            f"- Longest panel match inside it: {rec['max_match']} nt ({match_note})",

            f"- Graded hazard (heuristic): {rec['hazard']} "

            f"[{rec['risk_category']}, {rec['threshold_status']}]",

            f"- Configured maximum hazard reference: {rec['max_hazard']} "

            f"({rec['threshold_status']}; informational only)",

        ]

        if rec["hazard_threshold_equivalent_match_bp"] is not None:

            lines.append(

                f"- Equivalent-match reference: {rec['hazard_threshold_equivalent_match_bp']} nt "

                "under the configured hazard model; this does not gate candidate selection.")

        lines += [

            "- Hazard is a heuristic score, not a validated biological probability.",

            f"- T7 forward: `{primers['forward']}`",

            f"- T7 reverse: `{primers['reverse']}`", "",

        ]

        candidate = rec.get("candidate")

        if candidate:

            lines += [

                "## Candidate dsRNA sequence",

                f"- Candidate ID: `{candidate['candidate_id']}`",

                f"- Requested length: {candidate['requested_length_bp']} bp; "

                f"actual length: {candidate['actual_length_bp']} bp "

                f"(measured on the sequence below)",

                f"- Position {candidate['start_1based']}-{candidate['end_1based']} "

                f"(1-based, {candidate['strand']} strand of the target)",

                "- Target gene sequence and candidate: the sequence below is a window "

                "copied from the target CDS; the full target sequence is not replaced by "

                "it.",

                "",

                "```",

                candidate["sequence"],

                "```",

                "",

                "- This is an on-target candidate sequence only. Selecting or emitting it "

                "is not evidence of biological efficacy, uptake, stability or non-target "

                "safety, and no experimental validation is implied.",

                "",

            ]

    lines += [f"## Non-target matches >= {thr} nt across the whole gene"]

    if hits:

        lines += ["", "| organism | transcript | transcript pos | gene pos | length | strand |",

                  "|---|---|---|---|---|---|"]

        lines += [f"| {h[0]} | {h[1]} | {h[2]} | {h[3]}-{h[4]} | {h[5]} | {h[6]} |"

                  for h in hits[:100]]

        if len(hits) > 100:

            lines.append(f"| … | {len(hits) - 100} further matches | see hits.csv | | | |")

    else:

        lines += ["", f"None found (enumeration is complete for >= {thr} nt)."]



    lines += ["", "## Ranked off-target results (highest hazard first)"]

    if ranked:

        lines += ["", "| rank | organism | transcript/reference | target | length | identity | "

                      "strand | hazard | category | threshold |",

                  "|---|---|---|---|---|---|---|---|---|---|"]

        lines += [f"| {r['rank']} | {r['matched_organism']} | {r['transcript_or_reference_id']} | "

                  f"{r['target_start']}-{r['target_end']} | {r['match_length_bp']} | "

                  f"{r['identity_pct']}% | {r['strand']} | {r['hazard']} | {r['risk_category']} | "

                  f"{r['threshold_status']} |" for r in ranked[:50]]

        if len(ranked) > 50:

            lines.append(f"| … | {len(ranked) - 50} further rows | see offtarget_ranked.csv | | | | | | | |")

    else:

        lines += ["", "No off-target matches at or above the reporting floor were found against "

                      "the searched reference panel(s)."]

    floor = row.get("matches_below_report_floor")

    if floor:

        lines += ["",

                  f"_{floor} additional perfect match(es) shorter than {a.min_report_match} nt were "

                  f"found and are not listed: at this length a perfect match against a large "

                  f"reference is a random k-mer collision, and it cannot change whether a window "

                  f"meets the default hazard threshold. They are counted in "

                  f"`matches_below_report_floor`._"]



    lines += ["", "## Method notes and score definitions", "",

              f"- Hazard: {SCORING['hazard']['formula']}. {SCORING['hazard']['meaning']}",

              f"  {SCORING['hazard']['caveat']}",

              f"- Potency: {SCORING['potency']['formula']}. {SCORING['potency']['caveat']}",

              f"- Matching: {SCORING['matching_method']['description']}",

              f"  {SCORING['matching_method']['completeness']}",

              "",

              "**Limitations:** " + LIMITS]

    with open(os.path.join(gdir, "report.md"), "w", encoding="utf-8") as fh:

        fh.write("\n".join(lines) + "\n")





def main(argv=None):

    p = argparse.ArgumentParser(prog="dsguard",

                               description="dsRNA window design against an organism panel")

    p.add_argument("--target", required=True, help="FASTA of target gene CDS(s); one or many records")

    p.add_argument("--panel", action="append", default=[], metavar="NAME=FASTA",

                   help="non-target transcriptome; repeat for each organism")

    p.add_argument("--skipped-panel", action="append", default=[], metavar="NAME=REASON",

                   help="a configured panel that was NOT searched (missing/unreadable FASTA); "

                        "recorded in results.json and surfaced as an unresolved risk")

    p.add_argument("--window", type=int, default=400, help="dsRNA length in bp (default 400)")

    p.add_argument("--k", type=int, default=K,

                   help="seed length; exact for matches >= 2k-1 (default 11 -> 21 nt)")

    p.add_argument("--l50", type=float, default=L50_DEFAULT, help="PLACEHOLDER risk midpoint (not fitted)")

    p.add_argument("--kappa", type=float, default=KAPPA_DEFAULT, help="PLACEHOLDER risk steepness")

    p.add_argument("--max-hazard", type=float, default=0.05,

                   help="recommend the best-potency window at or under this hazard")

    p.add_argument("--min-report-match", type=int, default=None, metavar="NT",

                   help="only list off-target matches of at least this many nt "

                        "(default: the guaranteed detection threshold 2k-1 = %(default)s"

                        " overridden below by --k). Below that length a perfect match against a "

                        "genome-scale reference is expected by chance, so such matches are "

                        "counted (matches_below_report_floor) but not listed. The match that "

                        "sets the recommended window's hazard is always listed regardless.")

    p.add_argument("--out", default="dsguard_out")

    p.add_argument("--input-type", default="fasta_file",

                   help="how the target was supplied (fasta_file, pasted_sequence, ncbi)")

    p.add_argument("--input-label", default="",

                   help="human-readable identifier for the target input")

    a = p.parse_args(argv)

    if not a.panel:

        p.error("give at least one --panel NAME=FASTA")

    if a.window < 1:

        p.error("--window must be a positive integer")

    if not 1 <= a.k <= 16:

        p.error("--k must be an integer between 1 and 16")

    if a.kappa <= 0:

        p.error("--kappa must be greater than zero")

    if not 0.0 <= a.max_hazard <= 1.0:

        p.error("--max-hazard must be between 0 and 1")

    if a.min_report_match is not None and a.min_report_match < 1:

        p.error("--min-report-match must be a positive integer")

    os.makedirs(a.out, exist_ok=True)

    thr = threshold(a.k)

    if a.min_report_match is None:

        a.min_report_match = thr

    equivalent_match = hazard_equivalent_match(a.max_hazard, a.l50, a.kappa)

    a.hazard_equivalent_match = equivalent_match



    names, arrays, seen = [], [], {}

    invalid_bases = 0

    for nm, raw in read_fasta(a.target):

        seen[nm] = seen.get(nm, 0) + 1

        names.append(nm if seen[nm] == 1 else f"{nm}_{seen[nm]}")

        arr = encode(raw)

        invalid_bases += int((arr > 3).sum())

        arrays.append(arr)

    if not arrays:

        p.error("no records in target FASTA")

    Q, starts = concat_genes(arrays)

    print(f"{len(arrays)} gene(s); window {a.window} bp; exact for matches >= {thr} nt (k={a.k})")

    if invalid_bases:

        print(f"  WARNING: {invalid_bases} ambiguous/invalid base(s) in the target sequence")



    analysis_id = str(uuid.uuid4())

    timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()

    manifest = {

        "dsguard_version": __version__, "python": sys.version.split()[0],

        "numpy": np.__version__, "platform": platform.platform(),

        "analysis_id": analysis_id, "timestamp_utc": timestamp,

        "input_type": a.input_type, "input_label": a.input_label or a.target,

        "params": {"window": a.window, "k": a.k, "threshold": thr, "l50": a.l50,

                   "kappa": a.kappa, "max_hazard": a.max_hazard,

                   "min_report_match": a.min_report_match,

                   "hazard_threshold_equivalent_match_bp": (

                       None if equivalent_match is None else round(equivalent_match, 2)),

                   "l50_kappa_fitted": False},

        "target": {"path": a.target, "genes": len(arrays), "sha256": sha256_file(a.target),

                   "ambiguous_bases": invalid_bases},

        "panel": {}, "limits": LIMITS, "scoring": SCORING,

    }



    orgs = []

    per_gene = [{} for _ in arrays]

    for spec in a.panel:

        o, path = spec.split("=", 1)

        t0 = time.time()

        ref = Reference.from_fasta(path)

        by_gene = split_matches(scan_both_strands(Q, ref.seq, a.k), starts)

        dt = time.time() - t0

        orgs.append(o)

        for gi, ms in enumerate(by_gene):

            per_gene[gi][o] = [(m, *ref.locate(m.r_start)) for m in ms]

        manifest["panel"][o] = panel_provenance(o, path, ref, dt)

        print(f"  scanned {o}: {ref.seq.size / 1e6:.1f} Mb, {len(ref.names)} records "

              f"in {dt:.1f}s (all genes at once)")

        if ref.records_without_annotation:

            print(f"    WARNING: {ref.records_without_annotation} of {len(ref.names)} "

                  f"record(s) in '{o}' have no annotation in their FASTA header")



    skipped = []

    for spec in a.skipped_panel:

        nm, _, reason = spec.partition("=")

        skipped.append({"name": nm, "status": "not_searched",

                        "detail": reason or "configured panel was not available"})



    analysis_meta = {

        "analysis_id": analysis_id, "timestamp_utc": timestamp,

        "input_type": a.input_type, "input_label": a.input_label or a.target,

        "dsguard_version": __version__,

        "panels": [{"name": o, "status": "searched",

                    "detail": f"{manifest['panel'][o]['Mb']} Mb, "

                              f"{manifest['panel'][o]['records']} records"}

                   for o in orgs] + skipped,

    }



    gene_results = []

    for i in range(len(arrays)):

        gene_results.append(

            analyse_gene(names[i], arrays[i], orgs, per_gene[i], a, thr, a.out, analysis_meta)

        )

    rows = [g["row"] for g in gene_results]

    order = {"ok": 0}

    rows.sort(key=lambda r: (order.get(r["status"], 1), r.get("rec_hazard", 9),

                             -r.get("rec_potency", 0)))

    cols = ["gene", "status", "recommendation_status", "binary_status", "length_bp", "gc_pct",

            "windows", "clean_windows", "rec_start", "rec_end", "rec_potency", "rec_max_match",

            "rec_hazard", "max_hazard", "hazard_threshold_equivalent_match_bp", "total_hits",

            "offtargets_reported", "matches_below_report_floor"] \
        + [f"hits_{o}" for o in orgs]

    _write_csv(os.path.join(a.out, "summary.csv"), cols,

               [[r.get(c, "") if r.get(c) is not None else "" for c in cols] for r in rows])



    preventive = preventive_measures(gene_results,

                                     [{"name": o, **manifest["panel"][o]} for o in orgs],

                                     skipped, a.max_hazard,

                                     equivalent_match, invalid_bases, thr,

                                     input_type=a.input_type,

                                     input_label=a.input_label or a.target)

    results_payload = {

        "project": "dsGuard",

        "analysis_id": analysis_id,

        "created_utc": timestamp,

        "software": {"dsguard_version": __version__, "python": sys.version.split()[0],

                     "numpy": np.__version__, "platform": platform.platform()},

        "input": {"path": a.target, "type": a.input_type,

                  "label": a.input_label or a.target, "sha256": sha256_file(a.target),

                  "genes": len(arrays), "ambiguous_bases": invalid_bases,

                  "records": [{"name": n, "length_bp": int(arr.size)}

                              for n, arr in zip(names, arrays)]},

        "parameters": manifest["params"],

        "reporting": {"min_report_match_bp": a.min_report_match,

                      "guaranteed_detection_threshold_bp": thr,

                      "note": "Only perfect gapless matches of at least min_report_match_bp are "

                              "listed in ranked_offtargets. By default that floor is the "

                              "guaranteed detection threshold 2k-1, at which a perfect match "

                              "against a genome-scale reference is not expected by chance and "

                              "the enumeration is provably complete. Shorter best-effort matches "

                              "are counted in matches_below_report_floor and are not listed, "

                              "except the match that sets the recommended window's hazard, which "

                              "is always listed so the recommendation is traceable. Lower "

                              "--min-report-match to list more best-effort matches."},

        "scoring": SCORING,

        "reference_panels": {"searched": [{"name": o, **manifest["panel"][o]} for o in orgs],

                             "skipped": skipped},

        "summary": rows,

        "selection_policy": "pareto_highest_potency_then_lowest_hazard",

        "genes": [{"gene": g["row"]["gene"], "row": g["row"],

                   "length_bp": g["row"]["length_bp"],

                   "target_sequence": decode(arrays[i]),

                   "window": g["window"],

                   "screening_status": g["screening_status"],

                   "selection_policy": g["selection_policy"],

                   "total_matches_found": g["total_matches_found"],

                   "matches_below_report_floor": g["matches_below_report_floor"],

                   "candidate": g["candidate"],

                   "candidate_error": g["candidate_error"],

                   "recommendation": g["recommendation"],

                   "primers": g["primers"],

                   "pareto": g["pareto"],

                   "ranked_offtargets": g["ranked_offtargets"]} for g in gene_results],

        "preventive_measures": preventive,

        "limitations": [LIMITS, CANDIDATE_LIMITS],

    }

    _write_json(os.path.join(a.out, "results.json"), results_payload)

    manifest["genes"] = names

    _write_json(os.path.join(a.out, "manifest.json"), manifest)



    print(f"\n{'gene':<22}{'status':<28}{'rec window':<14}{'potency':<9}{'maxmatch':<9}{'hits':<5}")

    for r in rows:

        win = f"{r['rec_start']}-{r['rec_end']}" if "rec_start" in r else "-"

        print(f"{r['gene'][:21]:<22}{r['status'][:27]:<28}{win:<14}"

              f"{r.get('rec_potency', '-')!s:<9}{r.get('rec_max_match', '-')!s:<9}"

              f"{r['total_hits']:<5}")

    print(f"\nResults in {a.out}/ (summary.csv, results.json, <gene>/report.md)\nLIMITS: {LIMITS}")

    return 0
