

import copy

import json

from pathlib import Path



import numpy as np



from dsguard.cli import main





def write_fasta(path, records):

    with open(path, "w", encoding="utf-8") as fh:

        for header, seq in records:

            fh.write(f">{header}\n")

            for i in range(0, len(seq), 60):

                fh.write(seq[i:i + 60] + "\n")





def random_seq(rng, n):

    return "".join("ACGT"[b] for b in rng.integers(0, 4, n))





def build_inputs(tmp_path):

    rng = np.random.default_rng(20261024)





    target = random_seq(rng, 600)







    panel = list(random_seq(rng, 900))

    panel[100:130] = target[50:80]

    panel[300:317] = target[150:167]

    panel = "".join(panel)



    target_fa = tmp_path / "target.fa"

    panel_fa = tmp_path / "panel.fa"

    write_fasta(target_fa, [("G1 candidate gene", target)])

    write_fasta(panel_fa, [

        ("NC_TEST.1 Oryza sativa Japonica Group test chromosome 1", panel),

    ])

    return target_fa, panel_fa





def run_engine(tmp_path, out_name, target_fa, panel_fa, capsys):

    out = tmp_path / out_name

    code = main([

        "--target", str(target_fa),

        "--panel", f"rice={panel_fa}",

        "--window", "400",

        "--k", "11",

        "--max-hazard", "0.05",

        "--out", str(out),

        "--input-type", "fasta_file",

        "--input-label", "G1 candidate gene",

    ])

    capsys.readouterr()

    assert code == 0

    return out





def load(out, name):

    return json.loads((Path(out) / name).read_text(encoding="utf-8"))





def strip_volatile(panel):

    panel = copy.deepcopy(panel)

    panel.pop("seconds", None)

    return panel





def test_repeated_runs_are_identical(tmp_path, capsys):

    target_fa, panel_fa = build_inputs(tmp_path)



    r1 = run_engine(tmp_path, "run1", target_fa, panel_fa, capsys)

    r2 = run_engine(tmp_path, "run2", target_fa, panel_fa, capsys)





    assert (r1 / "summary.csv").read_bytes() == (r2 / "summary.csv").read_bytes()

    assert (r1 / "G1" / "offtarget_ranked.csv").read_bytes() == \

        (r2 / "G1" / "offtarget_ranked.csv").read_bytes()



    res1, res2 = load(r1, "results.json"), load(r2, "results.json")





    assert res1["analysis_id"] != res2["analysis_id"]



    assert res1["summary"] == res2["summary"]

    assert res1["software"] == res2["software"]

    assert res1["parameters"] == res2["parameters"]

    assert res1["preventive_measures"] == res2["preventive_measures"]

    for g1, g2 in zip(res1["genes"], res2["genes"]):

        assert g1["recommendation"] == g2["recommendation"]

        assert g1["ranked_offtargets"] == g2["ranked_offtargets"]

        assert g1["pareto"] == g2["pareto"]





    p1 = strip_volatile(res1["reference_panels"]["searched"][0])

    p2 = strip_volatile(res2["reference_panels"]["searched"][0])

    assert p1 == p2

    assert len(p1["sha256"]) == 64 and all(c in "0123456789abcdef" for c in p1["sha256"])

    assert p1["record_accessions"] == ["NC_TEST.1"]

    assert "Oryza sativa" in p1["annotation_example"]





    assert res1["genes"][0]["ranked_offtargets"], "expected at least one ranked match"
