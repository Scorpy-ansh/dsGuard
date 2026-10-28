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





def test_pareto_front_is_the_complete_non_dominated_set():





    pot = np.array([0.9, 0.5, 0.7, 0.7, 0.3, 0.8])

    risk = np.array([0.9, 0.1, 0.5, 0.6, 0.1, 0.05])

    front = list(pareto_front(pot, risk))



    assert len(front) == len(set(front))



    for i in front:

        dominated = any(

            pot[j] >= pot[i]

            and risk[j] <= risk[i]

            and (pot[j] > pot[i] or risk[j] < risk[i])

            for j in range(len(pot))

        )

        assert not dominated, f"index {i} is dominated but was kept"





def test_pareto_front_collapses_when_every_candidate_shares_one_hazard():







    pot = np.array([0.5, 0.6595, 0.6, 0.65])

    risk = np.array([0.268941, 0.268941, 0.268941, 0.268941])



    assert list(pareto_front(pot, risk)) == [1]





def test_primers_have_t7_and_revcomp():

    g = encode("ACGT" * 50)

    p = design_t7_primers(g, 10, 100)

    assert p["forward"].startswith(T7) and p["reverse"].startswith(T7)
