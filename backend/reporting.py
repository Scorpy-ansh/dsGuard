

from __future__ import annotations



import csv

import io

import json

import re

from datetime import datetime, timezone

from xml.sax.saxutils import escape as xml_escape



from reportlab.lib import colors

from reportlab.lib.enums import TA_LEFT

from reportlab.lib.pagesizes import A4, landscape

from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet

from reportlab.lib.units import mm

from reportlab.platypus import (

    KeepTogether,

    PageBreak,

    Paragraph,

    SimpleDocTemplate,

    Spacer,

    Table,

    TableStyle,

)



PROJECT = "dsGuard"









def _as_int(value):

    value = int(value or 0)

    return value





RANKED_FIELDS = [

    ("rank", "Rank", lambda r: r.get("rank")),

    ("target_id", "Target ID", lambda r: r.get("target_id")),

    ("target_window", "Target window", lambda r: r.get("target_window")),

    ("matched_organism", "Matched organism", lambda r: r.get("matched_organism")),

    ("transcript_or_reference_id", "Transcript / reference ID",

     lambda r: r.get("transcript_or_reference_id")),

    ("reference_panel", "Panel", lambda r: r.get("reference_panel")),

    ("match_length_bp", "Match (bp)", lambda r: r.get("match_length_bp")),

    ("identity_pct", "Identity (%)", lambda r: r.get("identity_pct")),

    ("mismatches", "Mismatch / gap",

     lambda r: f"{_as_int(r.get('mismatches'))}/{_as_int(r.get('gaps'))}"),

    ("strand", "Strand", lambda r: r.get("strand")),

    ("target_start", "Target coords",

     lambda r: f"{r.get('target_start')}-{r.get('target_end')}"),

    ("transcript_pos_1based", "Ref pos", lambda r: r.get("transcript_pos_1based")),

    ("potency", "Potency",

     lambda r: "" if r.get("potency") is None else r.get("potency")),

    ("hazard", "Risk score (hazard)", lambda r: r.get("hazard")),

    ("risk_category", "Risk category", lambda r: r.get("risk_category")),

    ("threshold_status", "Threshold status", lambda r: r.get("threshold_status")),

    ("evidence", "Evidence", lambda r: r.get("evidence")),

]



RANKED_HEADERS = [header for _, header, _ in RANKED_FIELDS]





def _specimen_rank(rows):

    seen = {}

    counter = 0

    for row in rows:

        key = (row.get("matched_organism"), row.get("transcript_or_reference_id"),

               row.get("strand"))

        if key not in seen:

            counter += 1

            seen[key] = counter

        yield seen[key]





def flatten_ranked_offtargets(results):

    base_ranked = []

    for gene in results.get("genes", []):

        gene_name = gene.get("gene") or gene.get("row", {}).get("gene", "")

        for row in gene.get("ranked_offtargets", []):

            base_ranked.append({"gene": gene_name, **row})





    base_ranked.sort(

        key=lambda r: (

            -float(r.get("hazard") or 0.0),

            -(int(r.get("match_length_bp") or 0)),

            str(r.get("matched_organism") or ""),

            str(r.get("transcript_or_reference_id") or ""),

            int(r.get("target_start") or 0),

        )

    )



    specimen = list(_specimen_rank(base_ranked))

    out = []

    for index, (row, spec_rank) in enumerate(zip(base_ranked, specimen), start=1):





        out.append({**row, "rank": index, "species_rank": spec_rank})

    return out





def build_offtarget_csv(results) -> str:

    rows = flatten_ranked_offtargets(results)



    output = io.StringIO(newline="")

    writer = csv.writer(output, lineterminator="\r\n")

    writer.writerow(["# dsGuard ranked off-target results"])

    writer.writerow([f"# analysis_id,{results.get('analysis_id', '')}"])

    writer.writerow([f"# created_utc,{results.get('created_utc', '')}"])

    writer.writerow([f"# input_label,{results.get('input', {}).get('label', '')}"])

    threshold = results.get("parameters", {}).get("max_hazard", "")

    writer.writerow([f"# configured_max_hazard,{threshold}"])

    writer.writerow([

        "# Ranking is by the documented off-target risk score (hazard), highest "

        "first; the eligibility decision imposed by the configured maximum hazard "

        "is reported separately in the 'threshold_status' column."

    ])

    writer.writerow([

        "# Only perfect, gapless matches >= min_report_match are listed, so "

        "identity is 100% and mismatch/gap counts are 0. hazard is a heuristic "

        "score, not a calibrated biological probability."

    ])

    writer.writerow(["gene"] + RANKED_HEADERS)



    for row in rows:

        writer.writerow(

            [row.get("gene", "")]

            + [formatter(row) if formatter(row) is not None else ""

               for _, _, formatter in RANKED_FIELDS]

        )



    return output.getvalue()





def build_results_json(results) -> str:

    payload = dict(results)

    payload["export"] = {

        "generated_utc": datetime.now(timezone.utc).isoformat(),

        "format": "dsguard-results-v1",

        "note": "Generated from the current analysis result object.",

    }

    return json.dumps(payload, indent=2, ensure_ascii=False, default=str) + "\n"

























FASTA_LINE_WIDTH = 60





def _fasta_field(value, limit=180) -> str:

    text = "" if value is None else str(value).replace(">", "")

    return " ".join(text.split())[:limit]





def clean_candidate_sequence(value) -> str:

    text = re.sub(r"\s+", "", "" if value is None else str(value))



    if not text:

        raise ValueError("the recorded candidate sequence is empty")



    invalid = re.search(r"[^A-Za-z]", text)

    if invalid:

        raise ValueError(

            "the recorded candidate sequence contains a non-nucleotide "

            f"character ({invalid.group(0)!r})"

        )



    return text





def candidate_gene_rows(results):

    rows = []



    for gene in results.get("genes", []) or []:

        candidate = gene.get("candidate")

        sequence = candidate.get("sequence") if isinstance(candidate, dict) else None



        if isinstance(sequence, str) and sequence.strip():

            rows.append((gene, candidate, clean_candidate_sequence(sequence)))

        else:

            rows.append((gene, None, ""))



    return rows





def candidate_records(results):

    return [

        (gene, candidate, sequence)

        for gene, candidate, sequence in candidate_gene_rows(results)

        if candidate is not None

    ]





def candidate_organism(results, gene) -> str:

    for value in (

        gene.get("organism"),

        results.get("organism"),

        (results.get("input") or {}).get("organism"),

    ):

        if isinstance(value, str) and value.strip():

            return value.strip()



    return ""





def build_candidate_fasta(results) -> str:

    lines = []



    for gene, candidate, sequence in candidate_records(results):

        start = candidate.get("start_1based")

        end = candidate.get("end_1based")

        position = f"{start}-{end}" if start is not None and end is not None else ""



        fields = [

            _fasta_field(

                gene.get("gene") or candidate.get("gene") or "candidate",

                limit=80,

            )

        ]

        for key, value in (

            ("candidate_id", candidate.get("candidate_id")),

            ("organism", candidate_organism(results, gene)),

            ("requested_length_bp", candidate.get("requested_length_bp")),





            ("length_bp", len(sequence)),

            ("target_position", position),

            ("strand", candidate.get("strand")),

        ):

            field = _fasta_field(value)

            if field:

                fields.append(f"{key}={field}")



        lines.append(">" + "|".join(fields))

        lines.extend(

            sequence[index:index + FASTA_LINE_WIDTH]

            for index in range(0, len(sequence), FASTA_LINE_WIDTH)

        )



    if not lines:

        return ""



    return "\n".join(lines) + "\n"





def wrapped_sequence(sequence, width=FASTA_LINE_WIDTH) -> str:

    return "<br/>".join(

        xml_escape(sequence[index:index + width])

        for index in range(0, len(sequence), width)

    )













def _styles():

    styles = getSampleStyleSheet()

    styles.add(ParagraphStyle(

        name="DsH1", parent=styles["Heading1"], fontSize=19, leading=23,

        spaceAfter=4, textColor=colors.HexColor("#111111"),

    ))

    styles.add(ParagraphStyle(

        name="DsH2", parent=styles["Heading2"], fontSize=13.5, leading=17,

        spaceBefore=13, spaceAfter=5, textColor=colors.HexColor("#1a1a1a"),

    ))

    styles.add(ParagraphStyle(

        name="DsH3", parent=styles["Heading3"], fontSize=11.5, leading=14,

        spaceBefore=9, spaceAfter=3, textColor=colors.HexColor("#222222"),

    ))

    styles.add(ParagraphStyle(

        name="DsBody", parent=styles["BodyText"], fontSize=9, leading=12.5,

        alignment=TA_LEFT, spaceAfter=3,

    ))

    styles.add(ParagraphStyle(

        name="DsCell", parent=styles["BodyText"], fontSize=7, leading=8.6,

        spaceAfter=0,

    ))

    styles.add(ParagraphStyle(

        name="DsSmall", parent=styles["BodyText"], fontSize=8, leading=11,

        textColor=colors.HexColor("#555555"), spaceAfter=3,

    ))





    styles.add(ParagraphStyle(

        name="DsSeq", parent=styles["BodyText"], fontName="Courier", fontSize=7,

        leading=8.6, spaceBefore=2, spaceAfter=6,

    ))

    return styles





def _kv_table(pairs, styles, key_width=52 * mm):

    body = [[Paragraph(f"<b>{k}</b>", styles["DsCell"]),

             Paragraph(str(v), styles["DsCell"])] for k, v in pairs]

    table = Table(body, colWidths=[key_width, None], hAlign="LEFT")

    table.setStyle(TableStyle([

        ("VALIGN", (0, 0), (-1, -1), "TOP"),

        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),

        ("TOPPADDING", (0, 0), (-1, -1), 2.5),

        ("LEFTPADDING", (0, 0), (-1, -1), 0),

        ("LINEBELOW", (0, 0), (-1, -2), 0.2, colors.HexColor("#dddddd")),

    ]))

    return table





def _data_table(header, rows, widths, styles, max_rows=None):

    body = [[Paragraph(f"<b>{h}</b>", styles["DsCell"]) for h in header]]

    body += [[Paragraph("" if c is None else str(c), styles["DsCell"]) for c in row]

             for row in (rows if max_rows is None else rows[:max_rows])]

    table = Table(body, colWidths=widths, repeatRows=1, hAlign="LEFT")

    table.setStyle(TableStyle([

        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#f0f0f4")),

        ("GRID", (0, 0), (-1, -1), 0.2, colors.HexColor("#cccccc")),

        ("VALIGN", (0, 0), (-1, -1), "TOP"),

        ("TOPPADDING", (0, 0), (-1, -1), 2),

        ("BOTTOMPADDING", (0, 0), (-1, -1), 2),

        ("LEFTPADDING", (0, 0), (-1, -1), 3),

        ("RIGHTPADDING", (0, 0), (-1, -1), 3),

    ]))

    return table





def _risk_colour(category):

    return {

        "Critical": colors.HexColor("#fde2e1"),

        "High": colors.HexColor("#fdeacd"),

        "Moderate": colors.HexColor("#fdf6d3"),

        "Low": colors.HexColor("#e8f5e9"),

    }.get(str(category), colors.white)





def build_pdf_report(results) -> bytes:

    styles = _styles()

    story = []



    params = results.get("parameters", {}) or {}

    input_info = results.get("input", {}) or {}

    reporting = results.get("reporting", {}) or {}

    panels = results.get("reference_panels", {}) or {}

    searched = panels.get("searched", []) or []

    skipped = panels.get("skipped", []) or []

    summary = results.get("summary", []) or []

    genes = results.get("genes", []) or []

    preventive = results.get("preventive_measures", {}) or {}

    scoring = results.get("scoring", {}) or {}

    max_hazard = params.get("max_hazard")



    story.append(Paragraph(f"{PROJECT} — Off-Target Risk Report", styles["DsH1"]))

    story.append(Paragraph(

        "Sequence-similarity screening of a candidate dsRNA target against a "

        "non-target reference panel. Research prototype, not a validated "

        "biological safety predictor.", styles["DsSmall"]))

    story.append(Spacer(1, 4))



    created = results.get("created_utc") or datetime.now(timezone.utc).isoformat()

    story.append(_kv_table([

        ("Project", PROJECT),

        ("Analysis ID", results.get("analysis_id", "n/a")),

        ("Generated (UTC)", created),

        ("Input identifier", input_info.get("label", "n/a")),

        ("Input type", input_info.get("type", "n/a")),

        ("Input records", input_info.get("genes", len(input_info.get("records", [])))),

        ("Ambiguous bases", input_info.get("ambiguous_bases", 0)),

        ("Reference SHA-256", str(input_info.get("sha256", ""))[:32]),

    ], styles))



    story.append(Paragraph("Analysis parameters", styles["DsH2"]))

    story.append(_kv_table([

        ("Window length", f"{params.get('window')} bp"),

        ("Seed length k", f"{params.get('k')} (exact for matches >= {params.get('threshold')} nt)"),

        ("L50 (placeholder)", params.get("l50")),

        ("Kappa (placeholder)", params.get("kappa")),

        ("Configured maximum hazard", max_hazard),

        ("Equivalent match length",

         f"{params.get('hazard_threshold_equivalent_match_bp')} nt "

         "(no perfect match longer than this can satisfy the threshold)"),

        ("Reporting floor", f"{reporting.get('min_report_match_bp')} nt"),

    ], styles))

    story.append(Paragraph(

        "L50 and kappa are placeholder constants and are NOT fitted to published "

        "non-target outcome data. The hazard value is therefore a monotone "

        "heuristic, not a validated probability.", styles["DsSmall"]))



    story.append(Paragraph("Reference panels", styles["DsH2"]))

    panel_rows = []

    for panel in searched:

        digest = str(panel.get("sha256", ""))

        panel_rows.append([

            panel.get("name", ""), "searched",

            f"{panel.get('Mb', '?')} Mb, {panel.get('records', '?')} records"

            + (f", sha256 {digest[:12]}…" if digest else ", sha256 not recorded"),

            str(panel.get("path", "")),

        ])

    for panel in skipped:

        panel_rows.append([

            panel.get("name", ""), "NOT SEARCHED",

            panel.get("detail", "configured panel was not available"), "-",

        ])

    if panel_rows:

        story.append(_data_table(

            ["Panel", "Status", "Detail", "Path"],

            panel_rows, [30 * mm, 22 * mm, 60 * mm, None], styles))





        for panel in searched:

            annotation = str(panel.get("annotation_example", "")).strip()

            if annotation:

                story.append(Paragraph(

                    f"{panel.get('name', '')} header annotation: {annotation}",

                    styles["DsSmall"]))

        if skipped:

            story.append(Paragraph(

                "Panels listed as NOT SEARCHED were never compared against the "

                "target. A missing panel is an unresolved gap in coverage and must "

                "not be read as absence of off-target risk.", styles["DsSmall"]))

    else:

        story.append(Paragraph(

            "No reference panel was searched. Off-target similarity cannot be "

            "assessed from this run.", styles["DsBody"]))



    story.append(Paragraph("Per-gene screening status", styles["DsH2"]))

    gene_rows = []

    for row in summary:

        verdict = row.get("screening_status")

        if not verdict:

            status = str(row.get("recommendation_status", ""))

            verdict = {

                "within_threshold": "screened_candidate_available",

                "fallback_threshold_exceeded": "screened_candidate_available",

            }.get(status, "legacy_or_unavailable")

        gene_rows.append([

            row.get("gene", ""),

            row.get("status", ""),

            verdict,

            f"{row.get('rec_start', '')}-{row.get('rec_end', '')}",

            row.get("rec_potency", ""),

            row.get("rec_max_match", ""),

            row.get("rec_hazard", ""),

            row.get("total_hits", ""),

        ])

    if gene_rows:

        story.append(_data_table(

            ["Gene", "Engine status", "Screening status", "Window", "Potency",

             "Max match (nt)", "Risk score", "Hits"],

            gene_rows, [38 * mm, 34 * mm, 55 * mm, 22 * mm, 15 * mm, 18 * mm, 17 * mm, 12 * mm],

            styles))

    else:

        story.append(Paragraph("No gene results were recorded.", styles["DsBody"]))



    story.append(Paragraph("Highest-risk findings", styles["DsH2"]))

    worst = []

    for gene in genes:

        rows = gene.get("ranked_offtargets", [])

        if rows:

            worst.append((gene.get("gene", ""), rows[0]))

    if worst:

        story.append(_data_table(

            ["Gene", "Matched organism", "Transcript / reference", "Match (bp)",

             "Risk score", "Category", "Threshold"],

            [[g, r.get("matched_organism", ""), r.get("transcript_or_reference_id", ""),

              r.get("match_length_bp", ""), r.get("hazard", ""), r.get("risk_category", ""),

              r.get("threshold_status", "")] for g, r in worst],

            [34 * mm, 34 * mm, 45 * mm, 20 * mm, 20 * mm, 22 * mm, 30 * mm], styles))

    else:

        story.append(Paragraph(

            "No off-target matches at or above the reporting floor were recorded "

            "for any gene against the searched panel(s).", styles["DsBody"]))



    story.append(PageBreak())



    ranked = flatten_ranked_offtargets(results)

    story.append(Paragraph("Ranked off-target matches", styles["DsH2"]))

    story.append(Paragraph(

        f"Ranked by the documented off-target risk score (hazard), highest first. "

        f"{len(ranked)} row(s) at or above the reporting floor of "

        f"{reporting.get('min_report_match_bp')} nt. The comparison with the "

        f"configured maximum hazard reference ({max_hazard}) is reported "

        f"separately in the Threshold status column and does not affect "

        f"candidate eligibility or ordering.", styles["DsSmall"]))

    if ranked:

        display = [[

            r.get("rank"), r.get("gene", ""), r.get("target_window", ""),

            r.get("matched_organism", ""), r.get("transcript_or_reference_id", ""),

            r.get("match_length_bp", ""), r.get("identity_pct", ""),

            f"{r.get('mismatches', 0)}/{r.get('gaps', 0)}", r.get("strand", ""),

            r.get("potency", ""), r.get("hazard", ""), r.get("risk_category", ""),

            r.get("threshold_status", ""),

        ] for r in ranked]

        table = _data_table(

            ["Rk", "Target ID", "Target window", "Matched organism",

             "Transcript / reference", "bp", "Ident %", "Mism/gap", "Strand",

             "Potency", "Hazard", "Category", "Threshold"],

            display,

            [9 * mm, 26 * mm, 22 * mm, 30 * mm, 32 * mm, 10 * mm, 13 * mm,

             14 * mm, 11 * mm, 13 * mm, 15 * mm, 17 * mm, 30 * mm],

            styles, max_rows=250)



        style_cmds = []

        for i, row in enumerate(display[:250], start=1):

            style_cmds.append(("BACKGROUND", (11, i), (11, i),

                               _risk_colour(row[11])))

            if str(row[12]) == "exceeds_threshold":

                style_cmds.append(("TEXTCOLOR", (12, i), (12, i),

                                   colors.HexColor("#a30000")))

        table.setStyle(TableStyle(style_cmds))

        story.append(table)

        if len(ranked) > 250:

            story.append(Paragraph(

                f"{len(ranked) - 250} further row(s) are in the CSV/JSON export "

                f"for this analysis.", styles["DsSmall"]))

    else:

        story.append(Paragraph(

            "No off-target matches were reported for the searched panel(s).",

            styles["DsBody"]))



    story.append(Paragraph("dsRNA candidate sequences", styles["DsH2"]))

    story.append(Paragraph(

        "The window this run selected inside each target gene, exactly as it was "

        "recorded in the stored analysis. These are candidate dsRNA sequences, "

        "NOT the submitted target gene sequences, and none of them was rebuilt "

        "from coordinates during export. Selecting a sequence is not evidence of "

        "biological efficacy, uptake, stability or non-target safety.",

        styles["DsSmall"]))

    candidate_rows = candidate_gene_rows(results)

    if candidate_rows:

        for gene, candidate, sequence in candidate_rows:

            name = xml_escape(str(gene.get("gene") or "unknown gene"))



            if candidate is None:





                error = gene.get("candidate_error")

                reason = f" Engine reason: {error}" if error else ""

                story.append(Paragraph(

                    f"<b>{name}</b> — no candidate sequence was generated for "

                    f"this gene.{xml_escape(reason)}",

                    styles["DsBody"]))

                continue



            start = candidate.get("start_1based")

            end = candidate.get("end_1based")

            target_length = gene.get("length_bp")

            target_note = (

                f", not the {target_length} bp target gene"

                if target_length else ", not the target gene itself"

            )



            story.append(Paragraph(name, styles["DsH3"]))

            story.append(_kv_table([

                ("Candidate ID", xml_escape(str(candidate.get("candidate_id", "n/a")))),

                ("Requested length", f"{candidate.get('requested_length_bp', 'n/a')} bp"),

                ("Actual length", f"{candidate.get('actual_length_bp', 'n/a')} bp "

                                  "(measured on the sequence below)"),

                ("Target window",

                 f"{start}-{end}" if start is not None and end is not None else "n/a"),

                ("Strand", xml_escape(str(candidate.get("strand", "n/a")))),

                ("Ambiguous bases", candidate.get("ambiguous_bases", "n/a")),

            ], styles))

            story.append(Paragraph(

                f"Candidate sequence, {len(sequence)} bp (the dsRNA window"

                f"{target_note}):", styles["DsSmall"]))

            story.append(Paragraph(wrapped_sequence(sequence), styles["DsSeq"]))

    else:

        story.append(Paragraph(

            "No gene in this analysis recorded a candidate sequence.",

            styles["DsBody"]))



    if genes and not any("candidate" in gene for gene in genes):

        story.append(Paragraph(

            "This analysis predates candidate-sequence output, so no candidate "

            "sequence could be reported. Re-run the screening to generate one.",

            styles["DsSmall"]))



    story.append(Paragraph("Interpretation of the scores", styles["DsH2"]))

    hazard = scoring.get("hazard", {})

    potency = scoring.get("potency", {})

    matching = scoring.get("matching_method", {})

    story.append(_kv_table([

        ("Risk score (hazard)", hazard.get("formula", "n/a")),

        ("Risk score input", hazard.get("input", "n/a")),

        ("Risk score range", hazard.get("range", "n/a")),

        ("Calibrated probability?", "No — heuristic. " + hazard.get("caveat", "")),

        ("Potency", potency.get("formula", "n/a")),

        ("Potency meaning", potency.get("meaning", "n/a")),

        ("Potency caveat", potency.get("caveat", "")),

        ("Matching method", matching.get("description", "n/a")),

        ("Enumeration completeness", matching.get("completeness", "n/a")),

    ], styles))

    story.append(Paragraph(

        f"Identity is reported as 100% and mismatch/gap counts as 0 because this "

        f"build reports only perfect, gapless matches. Partial and "

        f"reverse-complement matches: reverse-complement matches are reported on "

        f"the '-' strand; positions containing ambiguous bases never match and "

        f"shorten any match spanning them. A high heuristic score is evidence of a "

        f"perfect sequence match of the reported length, not proof of biological "

        f"harm.", styles["DsSmall"]))



    story.append(Paragraph("Preventive measures & risk mitigation", styles["DsH2"]))

    findings = preventive.get("findings", []) or []

    if findings:

        story.append(_data_table(

            ["Finding / risk", "Why it matters", "Supporting evidence",

             "Recommended preventive action", "Priority"],

            [[f.get("finding", ""), f.get("why", ""), f.get("evidence", ""),

              f.get("action", ""), f.get("priority", "")] for f in findings],

            [42 * mm, 46 * mm, 52 * mm, 58 * mm, 16 * mm], styles))

    else:

        story.append(Paragraph(

            "No specific findings were generated for this run.", styles["DsBody"]))



    categories = preventive.get("categories", {}) or {}

    for key, title in (

        ("before_accepting_a_candidate", "Before accepting a candidate"),

        ("reducing_potential_off_target_risk", "Reducing potential off-target risk"),

        ("validation_and_quality_assurance", "Validation and quality assurance"),

    ):

        items = categories.get(key) or []

        if not items:

            continue

        story.append(Paragraph(title, styles["DsH3"]))

        for item in items:

            story.append(Paragraph(f"• {item}", styles["DsBody"]))



    story.append(Paragraph("Scientific limitations and missing data", styles["DsH2"]))

    for limitation in results.get("limitations", []) or []:

        story.append(Paragraph(str(limitation), styles["DsBody"]))

    if skipped:

        story.append(Paragraph(

            "Missing data: reference panel(s) "

            + ", ".join(p.get("name", "?") for p in skipped)

            + " were not installed and were not searched.", styles["DsBody"]))

    story.append(Paragraph(

        "This report contains computational predictions only. No experimental "

        "validation was performed and no validated model performance is claimed.",

        styles["DsSmall"]))



    buffer = io.BytesIO()

    doc = SimpleDocTemplate(

        buffer,

        pagesize=landscape(A4),

        leftMargin=12 * mm,

        rightMargin=12 * mm,

        topMargin=12 * mm,

        bottomMargin=12 * mm,

        title=f"{PROJECT} off-target risk report",

        author=PROJECT,

        subject=str(results.get("analysis_id", "")),

    )

    doc.build(story)

    return buffer.getvalue()
