import assert from "node:assert/strict";

import test from "node:test";

import { readFileSync } from "node:fs";

import { dirname, join } from "node:path";

import { fileURLToPath } from "node:url";



const srcDir = dirname(fileURLToPath(import.meta.url));

const app = readFileSync(join(srcDir, "App.jsx"), "utf8");

const css = readFileSync(join(srcDir, "App.css"), "utf8");



const requiredHeadings = [

  "dsRNA Candidate Sequence",

  "Target Information and Original Target Sequence",

  "Screening Summary and Hazard Evidence",

  "Target Metadata and Run Parameters",

  "Detailed Off-Target Hits",

  "Preventive Measures and Risk Mitigation",

  "Gene report viewer",

  "Downloads and Exports",

];



test("Results page names the required sections and orders the candidate first", () => {

  for (const heading of requiredHeadings) {

    assert.ok(app.includes(heading), `missing Results heading: ${heading}`);

  }



  const summaryOrder = css.indexOf(".results-panel > .results-summary { order: 1; }");

  const candidateOrder = css.lastIndexOf(".results-panel > .results-candidate { order: 3; }");

  const targetOrder = css.indexOf(".results-panel > .results-target { order: 4; }");



  assert.ok(summaryOrder >= 0, "summary section must be near the top");

  assert.ok(candidateOrder > summaryOrder, "candidate section must follow the summary");

  assert.ok(targetOrder > candidateOrder, "target section must follow candidate section");

});



test("Hazard comparison notice is neutral and explicitly heuristic", () => {

  assert.ok(app.includes("Informational heuristic hazard comparison"));

  assert.ok(app.includes("not a pass/fail verdict"));

  assert.ok(app.includes("not proof of biological safety"));

  assert.ok(css.includes("background: var(--info-soft)"));

  assert.ok(css.includes("border: 1px solid var(--info-border)"));

  assert.ok(css.includes("color: var(--info)"));

});
