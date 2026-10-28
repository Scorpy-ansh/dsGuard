import numpy as np



from dsguard.optimise import pareto_front

from dsguard.primers import T7, design_t7_primers

from dsguard.seq import encode, revcomp

from dsguard.sirna import potency, reynolds, ui_tei, window_mean





def test_known_good_sirna_scores_high():



    s = encode("GCAGAAGUUUCAAAUAAUAUU")

    W = s[None, :]

    assert reynolds(W)[0] >= 6 and ui_tei(W)[0] >= 3





def test_gc_stretch_rule():



    assert ui_tei(encode("GCGCGCGCGCGAUAUAUAUAU")[None, :])[0] == 3

    assert ui_tei(encode("GCGCGCGCGAUAUAUAUAUAU")[None, :])[0] == 4





def test_ambiguous_windows_are_nan():

    s = encode("ACGU" * 30)

    s[50] = 4

    per, valid = potency(s)

    m = window_mean(per, valid, 40)

    assert np.isnan(m[11:51]).all()

    assert not np.isnan(m[:11]).any() and not np.isnan(m[51:]).any()





def test_pareto():

    pot = np.array([0.9, 0.5, 0.7, 0.7, 0.3])

    risk = np.array([0.9, 0.1, 0.5, 0.6, 0.1])

    assert list(pareto_front(pot, risk)) == [1, 2, 0]





def test_primers_have_t7_and_revcomp():

    g = encode("ACGT" * 50)

    p = design_t7_primers(g, 10, 100)

    assert p["forward"].startswith(T7) and p["reverse"].startswith(T7)



import pytest



from dsguard.optimise import recommend





def test_recommend_highest_potency_under_limit():

    potency_values = np.array([0.9, 0.5, 0.7])

    hazards = np.array([0.9, 0.02, 0.04])



    assert recommend(potency_values, hazards, [0, 1, 2], 0.05) == 2





def test_recommend_falls_back_to_lowest_hazard():

    potency_values = np.array([0.9, 0.5, 0.7])

    hazards = np.array([0.9, 0.8, 0.7])



    assert recommend(potency_values, hazards, [0, 1, 2], 0.05) == 2





@pytest.mark.parametrize(

    "potency_values,hazards,candidates",

    [

        ([0.1], [], [0]),

        ([0.1, 0.2], [0.1], [0]),

        ([0.1], [0.2], []),

        ([0.1], [0.2], [-1]),

        ([0.1], [0.2], [1]),

        ([0.1], [0.2], [0, 0]),

        ([np.nan], [0.2], [0]),

        ([0.1], [np.inf], [0]),

    ],

)

def test_recommend_rejects_invalid_inputs(potency_values, hazards, candidates):

    with pytest.raises(ValueError):

        recommend(potency_values, hazards, candidates)





def test_pareto_rejects_mismatched_lengths():

    from dsguard.optimise import pareto_front



    with pytest.raises(ValueError, match="same length"):

        pareto_front([0.2, 0.3], [0.1])





def test_pareto_rejects_non_finite_values():

    from dsguard.optimise import pareto_front



    with pytest.raises(ValueError, match="finite"):

        pareto_front([0.2, np.nan], [0.1, 0.2])
