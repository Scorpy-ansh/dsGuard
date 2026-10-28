













import test from "node:test";

import assert from "node:assert/strict";



import {

  ORGANISM_NAMES,

  findOrganism,

  organismDisplay,

  organismLabel,

  organismScientificName,

} from "./organisms.js";



test("verified common names are displayed for known mappings", () => {

  const expected = [

    ["Apis mellifera", "Honey bee"],

    ["Aedes aegypti", "Yellow fever mosquito"],

    ["Oryza sativa", "Rice"],

    ["Nilaparvata lugens", "Brown planthopper"],

    ["Cyrtorhinus lividipennis", "Mirid bug"],

  ];



  for (const [scientific, common] of expected) {

    assert.equal(organismLabel(scientific), common, scientific);

    assert.equal(organismLabel(common), common, `${common} (by common name)`);

  }

});



test("engine panel keys resolve to the same verified names", () => {



  assert.equal(organismLabel("honeybee"), "Honey bee");

  assert.equal(organismLabel("rice"), "Rice");

  assert.equal(organismLabel("cyrtorhinus_lividipennis"), "Mirid bug");

  assert.equal(organismScientificName("honeybee"), "Apis mellifera");

  assert.equal(organismScientificName("rice"), "Oryza sativa");

});



test("composite reference-panel labels are resolved through their binomial", () => {

  assert.equal(organismLabel("Honeybee (Apis mellifera)"), "Honey bee");

  assert.equal(organismLabel("Rice (Oryza sativa)"), "Rice");

  assert.equal(organismScientificName("Rice (Oryza sativa)"), "Oryza sativa");

});



test("scientific names and identifiers are preserved, never rewritten", () => {



  for (const record of ORGANISM_NAMES) {

    assert.equal(organismScientificName(record.scientific), record.scientific);

    for (const key of record.panelKeys) {

      assert.equal(organismScientificName(key), record.scientific);

      assert.equal(findOrganism(key).common, record.common);

    }

  }

});



test("an unknown organism falls back to its scientific name, never a guess", () => {

  const unknown = "Danio rerio";



  assert.equal(findOrganism(unknown), null);

  assert.equal(organismScientificName(unknown), "");

  assert.equal(organismLabel(unknown), unknown, "the given name must survive unchanged");

  assert.deepEqual(organismDisplay(unknown), { primary: unknown, secondary: "" });

});



test("missing metadata falls back to the caller's placeholder", () => {

  for (const empty of ["", "   ", null, undefined, 42, {}]) {

    assert.equal(findOrganism(empty), null);

    assert.equal(organismLabel(empty), "—");

    assert.equal(organismLabel(empty, "Not reported"), "Not reported");

    assert.deepEqual(organismDisplay(empty), { primary: "—", secondary: "" });

  }

});



test("a panel that names no species is not given one", () => {





  assert.equal(findOrganism("spider"), null);

  assert.equal(organismScientificName("spider"), "");

  assert.equal(organismLabel("spider"), "Spider (reference species panel)");

  assert.deepEqual(organismDisplay("spider"), {

    primary: "Spider (reference species panel)",

    secondary: "",

  });

});



test("a mapped organism never repeats its scientific name as secondary", () => {



  assert.deepEqual(organismDisplay("ApiS MELlifera"), {

    primary: "Honey bee",

    secondary: "Apis mellifera",

  });

  assert.deepEqual(organismDisplay("rice"), {

    primary: "Rice",

    secondary: "Oryza sativa",

  });

  assert.deepEqual(organismDisplay(""), { primary: "—", secondary: "" });

});



test("lookups are case- and whitespace-insensitive", () => {

  assert.equal(organismLabel("  HONEYBEE  "), "Honey bee");

  assert.equal(organismLabel("apis mellifera"), "Honey bee");

  assert.equal(organismLabel("Nilaparvata Lugens"), "Brown planthopper");

});
