

import json



import numpy as np

import pytest



from dsguard.cli import candidate_from_window, main

from dsguard.optimise import recommend

from dsguard.seq import encode













def write_fasta(path, records):

    with open(path, "w", encoding="utf-8") as fh:

        for header, seq in records:

            fh.write(f">{header}\n")

            for i in range(0, len(seq), 60):

                fh.write(seq[i:i + 60] + "\n")





def random_seq(seed, n):

    rng = np.random.default_rng(seed)

    return "".join("ACGT"[b] for b in rng.integers(0, 4, n))





def run_engine(tmp_path, target, window, panel=None, header="G1 candidate gene"):

    target_fa = tmp_path / "target.fa"

    panel_fa = tmp_path / "panel.fa"

    write_fasta(target_fa, [(header, target)])

    write_fasta(panel_fa, [("NC_TEST.1 reference test record", panel or "ACGT" * 50)])



    out = tmp_path / "out"

    code = main([

        "--target", str(target_fa),

        "--panel", f"rice={panel_fa}",

        "--window", str(window),

        "--out", str(out),

    ])



    assert code == 0

    return json.loads((out / "results.json").read_text(encoding="utf-8"))













@pytest.mark.parametrize("length", [21, 22, 50, 100, 201, 300, 400, 501])

def test_candidate_is_exactly_the_requested_window(length):

    target = random_seq(7, 800)

    seq = encode(target)



    candidate = candidate_from_window("G1", seq, 13, length)



    assert candidate["sequence"] == target[13:13 + length]

    assert candidate["requested_length_bp"] == length

    assert candidate["actual_length_bp"] == length

    assert candidate["start_1based"] == 14

    assert candidate["end_1based"] == 13 + length

    assert candidate["strand"] == "sense"





def test_actual_length_is_measured_from_the_sequence_not_hard_coded():

    target = random_seq(11, 900)

    seq = encode(target)



    lengths = [21, 60, 250, 501]

    actual = [

        candidate_from_window("G1", seq, 0, length)["actual_length_bp"]

        for length in lengths

    ]







    assert actual == lengths

    for length, value in zip(lengths, actual):

        candidate = candidate_from_window("G1", seq, 0, length)

        assert value == len(candidate["sequence"])





def test_candidate_id_is_deterministic_and_locates_the_window():

    seq = encode(random_seq(3, 500))



    first = candidate_from_window("G1", seq, 100, 150)

    again = candidate_from_window("G1", seq, 100, 150)



    assert first["candidate_id"] == again["candidate_id"]

    assert first["candidate_id"] == "G1:101-250"

    assert candidate_from_window("G1", seq, 200, 150)["candidate_id"] != first["candidate_id"]





def test_max_hazard_does_not_gate_candidate_selection():

    potency = np.array([0.9, 0.8, 0.7])

    hazard = np.array([0.9, 0.01, 0.2])

    candidates = np.array([0, 1, 2])



    assert recommend(potency, hazard, candidates, max_hazard=0.05) == 0

    assert recommend(potency, hazard, candidates, max_hazard=0.95) == 0





def test_recommend_uses_hazard_then_coordinate_for_equal_potency():

    potency = np.array([0.7, 0.7, 0.7])

    hazard = np.array([0.2, 0.1, 0.1])

    candidates = np.array([0, 1, 2])



    assert recommend(potency, hazard, candidates, max_hazard=0.05) == 1













def test_short_target_is_rejected_and_never_padded():

    seq = encode(random_seq(5, 40))



    with pytest.raises(ValueError) as excinfo:

        candidate_from_window("G1", seq, 0, 400)



    message = str(excinfo.value)

    assert "40 nt" in message and "400 bp" in message

    assert "pad" in message and "truncate" in message





def test_target_exactly_the_requested_length_is_accepted():

    seq = encode(random_seq(5, 100))



    candidate = candidate_from_window("G1", seq, 0, 100)



    assert candidate["actual_length_bp"] == 100

    assert candidate["start_1based"] == 1

    assert candidate["end_1based"] == 100





def test_window_running_past_the_target_end_is_rejected():

    seq = encode(random_seq(5, 300))



    with pytest.raises(ValueError, match="runs past the end"):

        candidate_from_window("G1", seq, 200, 150)





def test_last_valid_window_start_is_accepted():

    seq = encode(random_seq(5, 300))



    candidate = candidate_from_window("G1", seq, 150, 150)



    assert candidate["start_1based"] == 151

    assert candidate["end_1based"] == 300













def test_window_containing_an_ambiguous_base_is_rejected():

    seq = encode("ACGTACGTAC" * 6)

    seq[10] = 4



    with pytest.raises(ValueError, match="ambiguous"):

        candidate_from_window("G1", seq, 0, 21)





def test_valid_window_next_to_an_ambiguous_base_is_still_accepted():

    target = "ACGTACGTAC" * 6

    seq = encode(target)

    seq[10] = 4



    candidate = candidate_from_window("G1", seq, 11, 21)



    assert candidate["ambiguous_bases"] == 0

    assert candidate["actual_length_bp"] == 21

    assert candidate["sequence"] == target[11:32]





@pytest.mark.parametrize("length", [0, -3, 1.5, "50", None, True])

def test_non_positive_or_non_integer_length_is_rejected(length):

    seq = encode(random_seq(9, 200))



    with pytest.raises(ValueError):

        candidate_from_window("G1", seq, 0, length)





@pytest.mark.parametrize("start", [-1, 1.5, "0", None, True])

def test_invalid_window_start_is_rejected(start):

    seq = encode(random_seq(9, 200))



    with pytest.raises(ValueError):

        candidate_from_window("G1", seq, start, 50)













def test_candidate_is_a_copy_of_the_target_not_a_view():

    target = "ACGTACGTAC" * 6

    seq = encode(target)



    candidate = candidate_from_window("G1", seq, 0, 21)

    before = candidate["sequence"]



    seq[0] = 4



    assert isinstance(candidate["sequence"], str)

    assert candidate["sequence"] == before == target[:21]













def test_results_carry_the_selected_candidate_and_keep_the_target_separate(tmp_path, capsys):

    target = random_seq(20261027, 600)



    results = run_engine(tmp_path, target, window=60)

    gene = results["genes"][0]

    candidate = gene["candidate"]



    capsys.readouterr()



    assert candidate is not None

    assert gene["candidate_error"] is None







    assert gene["window"] == 60

    assert candidate["requested_length_bp"] == 60

    assert candidate["actual_length_bp"] == 60 == len(candidate["sequence"])





    assert candidate["sequence"] in target

    assert candidate["sequence"] == target[candidate["start_1based"] - 1:candidate["end_1based"]]





    assert gene["length_bp"] == 600

    assert gene["target_sequence"] == target

    assert candidate["actual_length_bp"] != gene["length_bp"]

    assert candidate["sequence"] != target





    assert candidate["start_1based"] == gene["recommendation"]["start"]

    assert candidate["end_1based"] == gene["recommendation"]["end"]

    assert gene["recommendation"]["candidate"] == candidate





@pytest.mark.parametrize("window", [21, 50, 100, 250, 501])

def test_requested_length_flows_through_the_cli_to_the_candidate(tmp_path, window, capsys):

    target = random_seq(4, 700)



    results = run_engine(tmp_path, target, window=window)

    candidate = results["genes"][0]["candidate"]



    capsys.readouterr()



    assert candidate["requested_length_bp"] == window

    assert candidate["actual_length_bp"] == window

    assert len(candidate["sequence"]) == window





def test_target_shorter_than_the_window_reports_no_candidate(tmp_path, capsys):

    target = random_seq(31, 120)



    results = run_engine(tmp_path, target, window=400)

    gene = results["genes"][0]



    capsys.readouterr()



    assert gene["candidate"] is None

    assert gene["recommendation"] is None

    assert gene["candidate_error"] is not None

    assert "120 nt" in gene["candidate_error"]

    assert "400 bp" in gene["candidate_error"]

    assert "never pads" in gene["candidate_error"]





def test_target_with_no_eligible_window_reports_no_candidate(tmp_path, capsys):

    target = "N" * 200



    results = run_engine(tmp_path, target, window=50)

    gene = results["genes"][0]



    capsys.readouterr()



    assert gene["candidate"] is None

    assert gene["recommendation"] is None

    assert gene["candidate_error"] is not None

    assert "ambiguous" in gene["candidate_error"]





def test_run_states_that_the_candidate_is_not_experimental_validation(tmp_path, capsys):

    results = run_engine(tmp_path, random_seq(2, 300), window=100)



    capsys.readouterr()



    limitations = " ".join(results["limitations"])

    assert "not evidence of biological efficacy" in limitations

    assert "no experimental validation is implied" in limitations
