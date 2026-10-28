

from dsguard.cli import preventive_measures





def gene_row(name="G1", **overrides):

    row = {

        "gene": name,

        "status": "ok",

        "recommendation_status": "within_threshold",

        "rec_start": 1,

        "rec_end": 30,

        "rec_hazard": 0.01,

        "rec_max_match": 12,

        "length_bp": 40,

    }

    row.update(overrides)

    return row





def gene_result(row=None, offtargets=None, window=30):

    return {

        "row": row or gene_row(),

        "ranked_offtargets": offtargets or [],

        "window": window,

    }





def match(**overrides):

    row = {

        "hazard_is_exact": True,

        "match_length_bp": 25,

        "matched_organism": "rice",

        "transcript_or_reference_id": "NC_008394.4",

        "target_start": 100,

        "target_end": 124,

        "strand": "+",

        "hazard": 0.880797,

    }

    row.update(overrides)

    return row





def panel(name="rice", **overrides):

    record = {

        "name": name,

        "path": f"reference_panels/{name}.fasta",

        "sha256": "abcdef0123456789",

        "Mb": 382.78,

        "records": 12,

    }

    record.update(overrides)

    return record





def find(findings, text):

    return [f for f in findings if text.lower() in f["finding"].lower()]





def build(**kwargs):

    return preventive_measures(

        kwargs.get("gene_results", []),

        kwargs.get("panels_searched", []),

        kwargs.get("panels_skipped", []),

        kwargs.get("max_hazard", 0.05),

        kwargs.get("equivalent_match"),

        kwargs.get("invalid_bases", 0),

        kwargs.get("thr"),

        input_type=kwargs.get("input_type"),

        input_label=kwargs.get("input_label", ""),

    )["findings"]













def test_no_panel_searched_is_critical():

    findings = build()



    hits = find(findings, "No reference panel was searched")

    assert len(hits) == 1

    assert hits[0]["priority"] == "Critical"

    assert "reference_panels.searched is empty" in hits[0]["evidence"]





def test_searched_panel_states_limited_coverage_with_provenance():

    findings = build(panels_searched=[panel()])



    hits = find(findings, "Reference coverage is limited to the 1 searched panel")

    assert len(hits) == 1

    assert hits[0]["priority"] == "Info"



    assert "rice (12 records, 382.78 Mb" in hits[0]["evidence"]

    assert "sha256 abcdef012345" in hits[0]["evidence"]



    assert "not evidence that no match exists" in hits[0]["why"]





def test_skipped_panel_is_reported_as_high_priority_gap():

    findings = build(

        panels_searched=[panel()],

        panels_skipped=[{"name": "honeybee", "detail": "FASTA not installed"}],

    )



    hits = find(findings, "Reference panel not installed: honeybee")

    assert len(hits) == 1

    assert hits[0]["priority"] == "High"

    assert "detail=FASTA not installed" in hits[0]["evidence"]





def test_skipped_panel_accepts_a_plain_name():

    findings = build(panels_skipped=["honeybee"])



    hits = find(findings, "Reference panel not installed: honeybee")

    assert len(hits) == 1

    assert "name=honeybee" in hits[0]["evidence"]





def test_missing_annotation_is_reported_with_panel_metadata():

    findings = build(

        panels_searched=[panel(records=12, records_without_annotation=3)],

    )



    hits = find(findings, "without annotation")

    assert len(hits) == 1

    assert hits[0]["priority"] == "Medium"

    evidence = hits[0]["evidence"]

    assert "records=12" in evidence

    assert "records_without_annotation=3" in evidence

    assert "path=reference_panels/rice.fasta" in evidence





def test_fully_annotated_panel_produces_no_annotation_finding():

    findings = build(panels_searched=[panel(records_without_annotation=0)])



    assert find(findings, "without annotation") == []













def test_threshold_below_detection_floor_is_flagged_high():

    findings = build(equivalent_match=15.11, thr=21)



    hits = find(findings, "shorter than the 21 nt guaranteed-detection threshold")

    assert len(hits) == 1

    assert hits[0]["priority"] == "High"

    assert "equivalent perfect match ~15.11 nt" in hits[0]["evidence"]



    assert "not evidence of biological off-target risk" in hits[0]["why"]





def test_threshold_above_detection_floor_is_not_flagged():

    findings = build(equivalent_match=25.0, thr=21)



    assert find(findings, "guaranteed-detection threshold") == []





def test_confirmed_exact_match_is_critical_and_keeps_levels_distinct():

    findings = build(

        gene_results=[gene_result(offtargets=[match()])],

        panels_searched=[panel()],

    )



    hits = find(findings, "confirmed exact match")

    assert len(hits) == 1

    assert hits[0]["priority"] == "Critical"

    assert "NC_008394.4" in hits[0]["finding"]

    assert "target 100-124" in hits[0]["finding"]



    assert "observed perfect sequence match" in hits[0]["why"]

    assert "heuristic" in hits[0]["why"]

    assert "does not by itself demonstrate a biological effect" in hits[0]["why"]





def test_sub_floor_match_only_is_not_a_confirmed_finding():

    findings = build(

        gene_results=[

            gene_result(

                offtargets=[match(hazard_is_exact=False, match_length_bp=17,

                                   hazard=0.119203)]

            )

        ],

    )



    assert find(findings, "confirmed exact match") == []





def test_hazard_above_reference_is_informational_and_cites_the_threshold():

    findings = build(

        panels_searched=[panel()],

        equivalent_match=15.11,

        thr=21,

        gene_results=[

            gene_result(

                row=gene_row(

                    recommendation_status="fallback_threshold_exceeded",

                    rec_start=1374,

                    rec_end=1773,

                    rec_hazard=0.119203,

                    rec_max_match=17,

                )

            )

        ],

    )



    hits = find(findings, "selected candidate hazard is above")

    assert len(hits) == 1

    assert hits[0]["priority"] == "Medium"

    assert "rec_hazard=0.119203, max_hazard=0.05" in hits[0]["evidence"]

    assert "rec_max_match=17 nt" in hits[0]["evidence"]

    assert "informational" in hits[0]["why"]





def test_hazard_at_or_below_reference_is_informational():

    findings = build(

        gene_results=[

            gene_result(row=gene_row(rec_hazard=0.01, rec_start=5, rec_end=34))

        ],

    )



    hits = find(findings, "selected candidate hazard is at or below")

    assert len(hits) == 1

    assert hits[0]["priority"] == "Medium"

    assert "rec_hazard=0.01, max_hazard=0.05" in hits[0]["evidence"]

    assert "informational" in hits[0]["why"]





def test_too_short_gene_reports_no_candidate_and_no_recommendation():

    findings = build(

        gene_results=[

            gene_result(row=gene_row(status="too_short(<30bp)", length_bp=20))

        ],

    )



    assert find(findings, "shorter than one candidate window")

    assert find(findings, "selected candidate hazard") == []

    assert find(findings, "satisfied the configured maximum hazard") == []













def test_ambiguous_bases_reported_with_count():

    findings = build(invalid_bases=4)



    hits = find(findings, "ambiguous/invalid base(s) in the target sequence")

    assert len(hits) == 1

    assert hits[0]["priority"] == "Medium"

    assert "input.ambiguous_bases=4" in hits[0]["evidence"]





def test_database_lookup_finding_only_for_ncbi_input():

    ncbi = find(

        build(input_type="ncbi", input_label="NM_000546.6"),

        "retrieved from a sequence database",

    )

    assert len(ncbi) == 1

    assert ncbi[0]["priority"] == "Info"

    assert "input.label=NM_000546.6" in ncbi[0]["evidence"]



    for input_type in ("fasta_file", "pasted_sequence"):

        assert find(build(input_type=input_type), "retrieved from a sequence database") == []













def test_every_finding_has_finding_why_action_and_evidence():

    findings = build(

        gene_results=[

            gene_result(offtargets=[match()],

                        row=gene_row(recommendation_status="fallback_threshold_exceeded"))

        ],

        panels_searched=[panel(records_without_annotation=2)],

        panels_skipped=[{"name": "honeybee", "detail": "FASTA not installed"}],

        max_hazard=0.05,

        equivalent_match=15.11,

        thr=21,

        invalid_bases=2,

        input_type="ncbi",

        input_label="NM_000546.6",

    )



    assert findings

    for entry in findings:

        for key in ("finding", "why", "action", "evidence"):

            assert entry.get(key), f"{key} missing from {entry}"





def test_qa_category_keeps_observed_match_heuristic_and_validation_apart():

    preventive = preventive_measures([], [panel()], [], 0.05, None, 0, 21)



    bullets = " ".join(preventive["categories"]["validation_and_quality_assurance"])

    assert "observed perfect sequence match" in bullets

    assert "heuristic" in bullets

    assert "experimental validation, which this tool does not perform" in bullets





def test_missing_panels_are_listed_as_unresolved_risk_in_categories():

    preventive = preventive_measures([], [panel()], [], 0.05, None, 0, 21)



    before = " ".join(preventive["categories"]["before_accepting_a_candidate"])

    assert "missing panels are listed above and are unresolved risks" in before
