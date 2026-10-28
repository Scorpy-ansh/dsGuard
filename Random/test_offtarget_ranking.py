

import numpy as np

from types import SimpleNamespace



from dsguard.cli import build_offtarget_rows, rank_offtargets





def row(**overrides):

    base = {

        "target_id": "G1",

        "target_window": "1-25",

        "matched_organism": "rice",

        "transcript_or_reference_id": "rice.1",

        "reference_panel": "rice",

        "transcript_pos_1based": 5,

        "target_start": 1,

        "target_end": 25,

        "match_length_bp": 25,

        "identity_pct": 100.0,

        "mismatches": 0,

        "gaps": 0,

        "strand": "+",

        "potency": 0.5,

        "hazard": 0.880797,

        "risk_category": "Critical",

        "threshold_status": "exceeds_threshold",

        "hazard_is_exact": True,

        "evidence": "confirmed exact match",

        "in_recommended_window": 1,

        "governs_recommendation": False,

    }

    base.update(overrides)

    return base





def seg(start, end, strand="+"):

    return SimpleNamespace(q_start=start, q_end=end, strand=strand)













def test_rank_orders_by_hazard_descending():

    ranked = rank_offtargets([

        row(hazard=0.1, match_length_bp=17),

        row(hazard=0.99, match_length_bp=25),

        row(hazard=0.5, match_length_bp=21),

    ])



    assert [r["hazard"] for r in ranked] == [0.99, 0.5, 0.1]

    assert [r["rank"] for r in ranked] == [1, 2, 3]





def test_rank_tie_break_is_deterministic():

    ranked = rank_offtargets([

        row(reference_panel="rice", matched_organism="rice",

            transcript_or_reference_id="NC_008394.4", target_start=900, hazard=0.5),

        row(reference_panel="rice", matched_organism="rice",

            transcript_or_reference_id="NC_008394.4", target_start=100, hazard=0.5),

        row(reference_panel="crop", matched_organism="maize",

            transcript_or_reference_id="zm.1", target_start=10, hazard=0.5),

        row(reference_panel="rice", matched_organism="rice",

            transcript_or_reference_id="NC_008400.2", target_start=400, hazard=0.5),

    ])





    assert [r["reference_panel"] for r in ranked] == ["crop", "rice", "rice", "rice"]

    assert [r["transcript_or_reference_id"] for r in ranked] == [

        "zm.1", "NC_008394.4", "NC_008394.4", "NC_008400.2",

    ]

    assert [r["target_start"] for r in ranked] == [10, 100, 900, 400]





def test_rank_collapses_exact_duplicates():

    ranked = rank_offtargets([row(), row()])



    assert len(ranked) == 1













def build(matches, **overrides):

    n, window = overrides.get("n", 60), overrides.get("window", 30)

    return build_offtarget_rows(

        overrides.get("gene", "G1"),

        overrides.get("orgs", ["rice"]),

        overrides.get("gene_matches", {"rice": matches}),

        overrides.get("thr", 21),

        overrides.get("max_hazard", 0.05),

        np.full(n - window + 1, 0.5),

        overrides.get("best", 0),

        window,

        n,

        overrides.get("min_report_match", 21),

        overrides.get("rec_max_match"),

    )





def test_matches_above_floor_are_listed_and_below_floor_counted():

    rows, below_floor, total = build(

        [(seg(0, 25), "rice.1", 4),

         (seg(40, 55), "rice.2", 90)],

    )



    assert total == 2

    assert below_floor == 1

    assert len(rows) == 1

    assert rows[0]["match_length_bp"] == 25

    assert rows[0]["hazard_is_exact"] is True





def test_sub_floor_match_that_governs_the_recommendation_is_still_listed():

    rows, below_floor, total = build(

        [(seg(5, 22), "rice.1", 4)],

        rec_max_match=17,

    )



    assert total == 1

    assert below_floor == 0

    assert len(rows) == 1

    assert rows[0]["match_length_bp"] == 17

    assert rows[0]["hazard_is_exact"] is False

    assert rows[0]["governs_recommendation"] is True

    assert "lower bound" in rows[0]["evidence"]



    assert rows[0]["match_length_bp"] < 21





def test_hits_count_uses_the_threshold_while_ranked_rows_add_sub_floor():

    matches = [(seg(0, 30), "rice.1", 0),

               (seg(40, 61), "rice.2", 80),

               (seg(5, 22), "rice.3", 200)]

    rows, below_floor, total = build(matches, rec_max_match=17)



    counted_hits = [r for r in rows if r["match_length_bp"] >= 21]



    assert total == 3

    assert len(counted_hits) == 2

    assert len(rows) == 3

    assert below_floor == 0





def test_rows_are_ranked_highest_hazard_first():

    rows, _, _ = build([(seg(0, 25), "rice.1", 4), (seg(0, 30), "rice.2", 50)])



    hazards = [r["hazard"] for r in rows]

    assert hazards == sorted(hazards, reverse=True)

    assert [r["rank"] for r in rows] == [1, 2]
