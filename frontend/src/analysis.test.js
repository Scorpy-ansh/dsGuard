



import test from "node:test";

import assert from "node:assert/strict";



import {

  DEFAULT_WINDOW,

  WINDOW_MAX,

  WINDOW_MIN,

  WINDOW_PRESETS,

  candidateSequence,

  candidateView,

  catalogGenesForOrganism,

  csvField,

  formatHazard,

  geneDisplayName,

  geneSummaryCsv,

  getHazardLimit,

  getRecommendationLabel,

  isThresholdExceeded,

  screeningStatus,

  screeningStatusLabel,

  jobDisplayName,

  paretoPlotPoints,

  presetForLength,

  targetGeneMetadata,

  validateWindowLength,

} from "./analysis.js";



const results = {

  analysis: { parameters: { max_hazard: 0.05 } },

  manifest: { params: { max_hazard: 0.9 } },

};



test("getHazardLimit prefers the analysis object over the manifest", () => {

  assert.equal(getHazardLimit(results), 0.05);

});



test("getHazardLimit falls back to the manifest for older jobs", () => {

  const older = { manifest: { params: { max_hazard: 0.05 } } };



  assert.equal(getHazardLimit(older), 0.05);

});



test("getHazardLimit skips empty values instead of reading them as zero", () => {

  const partial = {

    analysis: { parameters: { max_hazard: "" } },

    manifest: { params: { max_hazard: 0.05 } },

  };



  assert.equal(getHazardLimit(partial), 0.05);

  assert.equal(getHazardLimit({}), null);

});



test("isThresholdExceeded trusts an explicit fallback status", () => {

  const row = { recommendation_status: "fallback_threshold_exceeded" };



  assert.equal(isThresholdExceeded(row, results), true);

});



test("isThresholdExceeded catches a score above the configured threshold", () => {

  const row = { recommendation_status: "within_threshold", rec_hazard: 0.119203 };



  assert.equal(isThresholdExceeded(row, results), true);

});



test("isThresholdExceeded is false when the score is within the threshold", () => {

  const row = { recommendation_status: "within_threshold", rec_hazard: 0.01 };



  assert.equal(isThresholdExceeded(row, results), false);

});



test("isThresholdExceeded never treats a missing threshold as passing silently", () => {

  const row = { recommendation_status: "within_threshold", rec_hazard: 0.9 };



  assert.equal(isThresholdExceeded(row, {}), false);

  assert.equal(isThresholdExceeded(row, { analysis: { parameters: {} } }), false);

});



test("formatHazard renders six decimals and a dash for missing values", () => {

  assert.equal(formatHazard(0.119203), "0.119203");

  assert.equal(formatHazard(null), "—");

  assert.equal(formatHazard(undefined), "—");

  assert.equal(formatHazard("abc"), "—");

});



test("getRecommendationLabel reports neutral status for a legacy fallback", () => {

  const row = { recommendation_status: "fallback_threshold_exceeded", rec_hazard: 0.12 };



  assert.equal(getRecommendationLabel(row, results), "Screened — candidate available");

});



test("getRecommendationLabel reports neutral status for a legacy passing result", () => {

  const row = { recommendation_status: "within_threshold", rec_hazard: 0.01 };



  assert.equal(getRecommendationLabel(row, results), "Screened — candidate available");

});



test("getRecommendationLabel passes unknown statuses through and defaults honestly", () => {

  assert.equal(getRecommendationLabel({ recommendation_status: "queued" }, results), "Legacy result — status unavailable");

  assert.equal(getRecommendationLabel({}, results), "Legacy result — status unavailable");

});



test("screening status distinguishes new result states", () => {

  assert.equal(

    screeningStatus({ screening_status: "screened_candidate_available" }, {}),

    "screened_candidate_available"

  );

  assert.equal(

    screeningStatus({ screening_status: "screened_no_candidate" }, {}),

    "screened_no_candidate"

  );

  assert.equal(

    screeningStatus({ screening_status: "screening_incomplete" }, {}),

    "screening_incomplete"

  );

  assert.equal(

    screeningStatusLabel({ screening_status: "screened_no_candidate" }, {}),

    "Screened — no candidate"

  );

});



test("target catalog filters genes by target organism", () => {

  const catalog = [

    { scientific_name: "Chilo suppressalis", genes: [{ gene_name: "CsEcR" }] },

    { scientific_name: "Sogatella furcifera", genes: [{ gene_name: "SfEcR" }] },

  ];



  assert.deepEqual(catalogGenesForOrganism(catalog, "Chilo suppressalis"), [

    { gene_name: "CsEcR" },

  ]);

  assert.deepEqual(catalogGenesForOrganism(catalog, "Nilaparvata lugens"), []);

});



test("target metadata uses persisted accession and organism without conflating candidate data", () => {

  const metadata = targetGeneMetadata(

    {

      gene: "NlRan",

      accession: "KT313028.1",

      organism: "Nilaparvata lugens",

      length_bp: 645,

      target_sequence: "ACGT",

      candidate: { sequence: "TGCA" },

    },

    { input: { type: "ncbi" } }

  );



  assert.equal(metadata.accession, "KT313028.1");

  assert.equal(metadata.organism, "Nilaparvata lugens");

  assert.equal(metadata.sequence, "ACGT");

  assert.notEqual(metadata.sequence, "TGCA");

});



test("target metadata reports unavailable values instead of inferring them", () => {

  const metadata = targetGeneMetadata(

    { gene: "NlRan", target_sequence: "ACGT" },

    { input: { type: "pasted_sequence" } }

  );



  assert.equal(metadata.accession, "");

  assert.equal(metadata.organism, "");

});







test("paretoPlotPoints returns exactly one marker per candidate", () => {

  const points = [

    { start: 1, end: 400, potency: 0.5, hazard: 0.05, threshold_status: "meets_threshold" },

    { start: 101, end: 500, potency: 0.6, hazard: 0.1, threshold_status: "exceeds_threshold" },

    { start: 201, end: 600, potency: 0.7, hazard: 0.2, threshold_status: "exceeds_threshold" },

    { start: 301, end: 700, potency: 0.55, hazard: 0.02, threshold_status: "meets_threshold" },

  ];



  const plotted = paretoPlotPoints(points);



  assert.equal(plotted.length, points.length);

  assert.deepEqual(

    plotted.map((point) => point.start),

    points.map((point) => point.start)

  );

});



test("paretoPlotPoints keeps distinct candidates at distinct coordinates", () => {

  const plotted = paretoPlotPoints([

    { start: 1, end: 400, potency: 0.5, hazard: 0.05, threshold_status: "meets_threshold" },

    { start: 101, end: 500, potency: 0.6, hazard: 0.1, threshold_status: "exceeds_threshold" },

  ]);



  assert.notEqual(plotted[0].x, plotted[1].x);

  assert.notEqual(plotted[0].y, plotted[1].y);



  assert.ok(plotted[1].y < plotted[0].y);

});



test("paretoPlotPoints centres a single-point front instead of clipping it", () => {

  const single = paretoPlotPoints([

    { start: 1039, end: 1438, potency: 0.6595, hazard: 0.268941, threshold_status: "exceeds_threshold" },

  ]);



  assert.equal(single.length, 1);



  assert.ok(single[0].x > 36 && single[0].x < 520 - 36);

  assert.ok(single[0].y > 36 && single[0].y < 220 - 36);

});



test("paretoPlotPoints sketches nothing for an empty front", () => {

  assert.deepEqual(paretoPlotPoints([]), []);

  assert.deepEqual(paretoPlotPoints(undefined), []);

});



test("paretoPlotPoints ignores records without actual numeric hazard and potency", () => {

  const plotted = paretoPlotPoints([

    { start: 1, end: 400, potency: 0.5, hazard: 0.05 },

    { start: 101, end: 500, potency: "not-recorded", hazard: 0.1 },

  ]);



  assert.equal(plotted.length, 1);

  assert.equal(plotted[0].start, 1);

  assert.equal(plotted[0].hazard, 0.05);

  assert.equal(plotted[0].potency, 0.5);

});









test("geneDisplayName shows the alias for a pipe-delimited project header", () => {

  assert.equal(

    geneDisplayName("NlRan|accession=KT313028.1|cds_length=645"),

    "NlRan"

  );

  assert.equal(geneDisplayName("Nlsid-1|accession=MH184552.1"), "Nlsid-1");

});



test("geneDisplayName leaves a plain identifier untouched", () => {

  assert.equal(geneDisplayName("KT313028.1"), "KT313028.1");

  assert.equal(geneDisplayName("NlRan"), "NlRan");

  assert.equal(geneDisplayName("NlRan some annotation"), "NlRan some annotation");

});



test("geneDisplayName never renders an empty cell", () => {

  assert.equal(geneDisplayName(""), "Unknown gene");

  assert.equal(geneDisplayName("   "), "Unknown gene");

  assert.equal(geneDisplayName(null), "Unknown gene");

  assert.equal(geneDisplayName(undefined), "Unknown gene");

});



test("jobDisplayName cleans a pasted FASTA header stored as the target label", () => {

  assert.equal(

    jobDisplayName({ target_filename: "NlRan|accession=KT313028.1|cds_length=645" }),

    "NlRan"

  );

});



test("jobDisplayName keeps file names and NCBI labels readable", () => {

  assert.equal(

    jobDisplayName({ target_filename: "nilaparvata_lugens_13_cds.fasta" }),

    "nilaparvata_lugens_13_cds.fasta"

  );

  assert.equal(jobDisplayName({ target_filename: "my genes.fasta" }), "my genes.fasta");

  assert.equal(

    jobDisplayName({ target_filename: "NCBI:KT313028.1 (NlRan)" }),

    "NCBI:KT313028.1 (NlRan)"

  );

});



test("jobDisplayName defaults honestly when a job has no label yet", () => {

  assert.equal(jobDisplayName({}), "Screening job");

  assert.equal(jobDisplayName(null), "Screening job");

});











test("the five preset labels are exactly the specified text", () => {

  assert.deepEqual(

    WINDOW_PRESETS.map((preset) => preset.label),

    [

      "50\u2013101 bp",

      "100\u2013201 bp",

      "200\u2013301 bp",

      "300\u2013401 bp",

      "400\u2013501 bp",

    ]

  );

});



test("every preset sends an in-range, unique lower bound of its own band", () => {

  const sent = WINDOW_PRESETS.map((preset) => preset.value);



  for (const preset of WINDOW_PRESETS) {

    assert.equal(preset.value, preset.min, `${preset.label} must send its lower bound`);

    assert.ok(preset.min >= WINDOW_MIN, `${preset.label} starts below the supported minimum`);

    assert.ok(preset.max <= WINDOW_MAX, `${preset.label} ends above the supported maximum`);

    assert.ok(preset.min <= preset.value && preset.value <= preset.max);

  }

  assert.equal(new Set(sent).size, sent.length);

});



test("WINDOW_MIN is the engine's own floor and WINDOW_MAX is the preset range top", () => {



  assert.equal(WINDOW_MIN, 21);

  assert.equal(WINDOW_MAX, Math.max(...WINDOW_PRESETS.map((preset) => preset.max)));

});



test("the default length is still 400 and belongs to the last preset", () => {

  assert.equal(DEFAULT_WINDOW, 400);

  assert.equal(

    WINDOW_PRESETS.find((preset) => preset.value === DEFAULT_WINDOW).label,

    "400\u2013501 bp"

  );

});



test("presetForLength resolves overlapping boundaries by the first band", () => {



  assert.equal(presetForLength(100).label, "50\u2013101 bp");

  assert.equal(presetForLength(101).label, "50\u2013101 bp");

  assert.equal(presetForLength(102).label, "100\u2013201 bp");

  assert.equal(presetForLength(250).label, "200\u2013301 bp");

  assert.equal(presetForLength(401).label, "300\u2013401 bp");

  assert.equal(presetForLength(501).label, "400\u2013501 bp");



  assert.equal(presetForLength(30), null);

  assert.equal(presetForLength(999), null);

  assert.equal(presetForLength(""), null);

});



test("validateWindowLength accepts valid lengths, including shared boundaries", () => {

  for (const value of [21, 50, 100, 101, 200, 201, 300, 301, 400, 401, 501]) {

    assert.deepEqual(validateWindowLength(value), { ok: true, value });

  }

  assert.deepEqual(validateWindowLength("250"), { ok: true, value: 250 });

  assert.deepEqual(validateWindowLength(" 400 "), { ok: true, value: 400 });

});



test("validateWindowLength rejects invalid lengths instead of substituting 400", () => {

  for (const bad of ["", "   ", "abc", "40 0", "40.5", "1e3", "-5", "0", "20", "502", "5000"]) {

    const result = validateWindowLength(bad);

    assert.equal(result.ok, false, `${JSON.stringify(bad)} should be rejected`);

    assert.ok(result.error.length > 0);

    assert.equal(result.value, undefined);

  }

});



test("a rejected length reports the supported range rather than 400", () => {

  assert.match(validateWindowLength(502).error, /between 21 and 501 bp/);

  assert.match(validateWindowLength(20).error, /between 21 and 501 bp/);

  assert.doesNotMatch(validateWindowLength("abc").error, /400/);

});













const targetSequence =

  "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACAC";



function geneWithCandidate(overrides = {}) {

  const sequence = overrides.sequence ?? targetSequence.slice(11, 32);

  return {

    gene: "NlRan|cds|accession=XM_TEST1",

    length_bp: targetSequence.length,

    candidate: {

      candidate_id: "NlRan:12-32",

      gene: "NlRan|cds|accession=XM_TEST1",

      requested_length_bp: 21,

      actual_length_bp: sequence.length,

      sequence,

      start_1based: 12,

      end_1based: 32,

      strand: "sense",

      ambiguous_bases: 0,

    },

    candidate_error: null,

    ...overrides,

  };

}



test("candidateSequence returns the recorded candidate, never the target", () => {

  const gene = geneWithCandidate();

  const sequence = candidateSequence(gene);



  assert.equal(sequence, gene.candidate.sequence);

  assert.notEqual(sequence, targetSequence);

  assert.equal(sequence.length, gene.candidate.actual_length_bp);

});



test("candidateSequence is empty when a run recorded no candidate", () => {

  assert.equal(candidateSequence({ candidate: null }), "");

  assert.equal(candidateSequence({ candidate: { sequence: "" } }), "");

  assert.equal(candidateSequence({}), "");

  assert.equal(candidateSequence(undefined), "");

});



test("candidateView exposes the exact fields the API returned", () => {

  const gene = geneWithCandidate();

  const view = candidateView(gene);



  assert.equal(view.available, true);

  assert.equal(view.candidateId, "NlRan:12-32");

  assert.equal(view.sequence, gene.candidate.sequence);

  assert.equal(view.requestedLengthBp, 21);

  assert.equal(view.actualLengthBp, gene.candidate.actual_length_bp);

  assert.equal(view.start1Based, 12);

  assert.equal(view.end1Based, 32);

  assert.equal(view.strand, "sense");

  assert.equal(view.ambiguousBases, 0);

  assert.equal(view.message, "");

});



test("candidateView keeps the candidate and target lengths distinct", () => {

  const gene = geneWithCandidate();

  const view = candidateView(gene);



  assert.notEqual(view.actualLengthBp, gene.length_bp);

  assert.equal(view.actualLengthBp, view.sequence.length);

});



test("candidateView surfaces the backend candidate_error for a null candidate", () => {

  const view = candidateView({

    gene: "shortgene",

    length_bp: 15,

    candidate: null,

    candidate_error: "target is 15 bp, shorter than the requested 21 bp window",

  });



  assert.equal(view.available, false);

  assert.equal(view.sequence, "");

  assert.match(view.message, /shorter than the requested 21 bp/);

});



test("candidateView distinguishes a missing candidate from an older result", () => {

  const recorded = candidateView({ gene: "g", candidate: null });

  const older = candidateView({ gene: "g" });



  assert.equal(recorded.available, false);

  assert.equal(older.available, false);

  assert.notEqual(recorded.message, older.message);

  assert.match(older.message, /predates candidate output/);

});



test("candidateView leaves missing length fields null instead of zero", () => {

  const view = candidateView({

    gene: "g",

    candidate: { candidate_id: "g:1-21", sequence: "ACGTACGTACACGTACGTACAC" },

  });



  assert.equal(view.available, true);

  assert.equal(view.requestedLengthBp, null);

  assert.equal(view.actualLengthBp, view.sequence.length);

  assert.equal(view.start1Based, null);

  assert.equal(view.end1Based, null);

  assert.equal(view.ambiguousBases, null);

});



test("candidateView never invents a sequence for a malformed candidate", () => {

  for (const bad of [null, undefined, 42, {}, { sequence: 123 }, "ACGT"]) {

    const view = candidateView({ gene: "g", candidate: bad });

    assert.equal(view.available, false, `${JSON.stringify(bad)} must not be shown`);

    assert.equal(view.sequence, "");

  }

});



test("targetGeneMetadata reports the target separately from the candidate", () => {

  const gene = geneWithCandidate();

  const metadata = targetGeneMetadata(gene, {

    input: { label: "demo.fasta", type: "fasta_file" },

  });



  assert.equal(metadata.displayName, "NlRan");

  assert.equal(metadata.accession, "XM_TEST1");

  assert.equal(metadata.lengthBp, targetSequence.length);

  assert.equal(metadata.inputLabel, "demo.fasta");

  assert.equal(metadata.inputType, "fasta_file");

  assert.notEqual(metadata.lengthBp, gene.candidate.actual_length_bp);

});



test("targetGeneMetadata leaves absent target fields empty rather than guessed", () => {

  const metadata = targetGeneMetadata({ gene: "plain_gene" }, {});



  assert.equal(metadata.gene, "plain_gene");

  assert.equal(metadata.displayName, "plain_gene");

  assert.equal(metadata.accession, "");

  assert.equal(metadata.lengthBp, null);

  assert.equal(metadata.inputLabel, "");

  assert.equal(metadata.inputType, "");

});













function parseCsv(text) {

  const rows = [];

  let row = [];

  let field = "";

  let quoted = false;



  for (let index = 0; index < text.length; index += 1) {

    const char = text[index];



    if (quoted) {

      if (char === '"' && text[index + 1] === '"') {

        field += '"';

        index += 1;

      } else if (char === '"') {

        quoted = false;

      } else {

        field += char;

      }

    } else if (char === '"') {

      quoted = true;

    } else if (char === ",") {

      row.push(field);

      field = "";

    } else if (char === "\r") {
      continue;

    } else if (char === "\n") {

      row.push(field);

      rows.push(row);

      row = [];

      field = "";

    } else {

      field += char;

    }

  }



  if (field !== "" || row.length) {

    row.push(field);

    rows.push(row);

  }



  return rows;

}



function exportResults(genes, summary) {

  return {

    job_id: "job-1",

    analysis: { parameters: { max_hazard: 0.05 }, genes },

    summary,

  };

}



const summaryRow = {

  gene: "NlRan|cds|accession=XM_TEST1",

  status: "ok",

  recommendation_status: "fallback_threshold_exceeded",

  rec_hazard: 0.377541,

  length_bp: 63,

  rec_max_match: 20,

  total_hits: 1,

};



test("csvField quotes a value so its content cannot split a row", () => {

  assert.equal(csvField("a,b"), '"a,b"');

  assert.equal(csvField('say "hi"'), '"say ""hi"""');

  assert.equal(csvField(null), '""');

  assert.equal(csvField(21), '"21"');

});



test("geneSummaryCsv exports the candidate, never the target sequence", () => {

  const gene = geneWithCandidate();

  const rows = parseCsv(

    geneSummaryCsv(exportResults([gene], [{ ...summaryRow, gene: gene.gene }]))

  );



  assert.equal(rows.length, 2, "one header row plus one row per gene");

  assert.equal(rows[0].length, rows[1].length);



  const header = rows[0];

  const value = (name) => rows[1][header.indexOf(name)];



  assert.equal(value("Gene"), gene.gene);

  assert.equal(value("Candidate ID"), "NlRan:12-32");

  assert.equal(value("Candidate Requested Length (bp)"), "21");

  assert.equal(value("Candidate Actual Length (bp)"), "21");

  assert.equal(value("Candidate Sequence"), gene.candidate.sequence);

  assert.equal(value("Candidate Sequence").length, gene.candidate.actual_length_bp);

  assert.notEqual(value("Candidate Sequence"), targetSequence);

  assert.equal(value("Candidate Note"), "");



  assert.equal(value("Recommended Hazard Score"), "0.377541");

  assert.equal(value("Configured Maximum Hazard"), "0.05");

});



test("geneSummaryCsv keeps a null candidate empty and carries its reason", () => {

  const gene = {

    gene: "shortgene",

    length_bp: 15,

    candidate: null,

    candidate_error: "target is 15 bp, shorter than the requested 30 bp window",

  };

  const rows = parseCsv(

    geneSummaryCsv(exportResults([gene], [{ ...summaryRow, gene: gene.gene }]))

  );

  const header = rows[0];

  const value = (name) => rows[1][header.indexOf(name)];



  assert.equal(value("Candidate Sequence"), "");

  assert.equal(value("Candidate Actual Length (bp)"), "");

  assert.match(value("Candidate Note"), /shorter than the requested 30 bp/);

});



test("geneSummaryCsv marks a run that predates candidate output without crashing", () => {

  const gene = { gene: "legacy", length_bp: 40 };

  const rows = parseCsv(

    geneSummaryCsv(exportResults([gene], [{ ...summaryRow, gene: gene.gene }]))

  );

  const header = rows[0];

  const value = (name) => rows[1][header.indexOf(name)];



  assert.equal(value("Candidate Sequence"), "");

  assert.match(value("Candidate Note"), /predates candidate output/);

});



test("geneSummaryCsv escapes values containing commas, quotes and breaks", () => {

  const gene = geneWithCandidate();

  gene.gene = 'NlRan,accession="XM_TEST1"';

  gene.candidate.sequence = 'ACG"T,ACGT\nACGT';



  const csv = geneSummaryCsv(

    exportResults([gene], [{ ...summaryRow, gene: gene.gene }])

  );

  const rows = parseCsv(csv);

  const header = rows[0];



  assert.equal(rows.length, 2, "a quoted value must not add a CSV row");

  assert.equal(rows[1].length, header.length, "every row must match the header");

  assert.equal(rows[1][header.indexOf("Gene")], 'NlRan,accession="XM_TEST1"');

  assert.equal(

    rows[1][header.indexOf("Candidate Sequence")],

    'ACG"T,ACGT\nACGT',

    "the sequence must stay one field and round-trip unchanged"

  );

});



test("geneSummaryCsv returns an empty string when the run has no summary", () => {

  assert.equal(geneSummaryCsv({ analysis: { genes: [] } }), "");

  assert.equal(geneSummaryCsv(null), "");

});
