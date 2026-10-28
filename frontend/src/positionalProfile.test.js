import assert from "node:assert/strict";

import test from "node:test";



import {

  contiguousSegments,

  numericRange,

  positionalProfileData,

} from "./positionalProfile.js";



test("positionalProfileData maps windows and per-position profiles without changing coordinates", () => {

  const data = positionalProfileData({

    target_length: "4306",

    windows: [

      { start: "1", end: "400", potency: "0.5096", max_match: "18" },

      { start: "2", end: "401", potency: "0.5104", max_match: "19" },

    ],

    profiles: [

      { position_1based: "1", rice: "14" },

      { position_1based: "2", rice: "17" },

    ],

  });



  assert.equal(data.available, true);

  assert.deepEqual(data.potencyPoints, [

    { position: 1, value: 0.5096 },

    { position: 2, value: 0.5104 },

  ]);

  assert.deepEqual(data.profilePoints("rice"), [

    { position: 1, value: 14 },

    { position: 2, value: 17 },

  ]);

  assert.deepEqual(numericRange(data.profilePoints("rice")), { min: 14, max: 17 });

});



test("invalid and missing values remain gaps instead of becoming zero", () => {

  const data = positionalProfileData({

    target_length: 10,

    windows: [

      { start: 1, end: 4, potency: "0.5" },

      { start: 2, end: 5, potency: "not recorded" },

      { start: 3, end: 6, potency: "" },

    ],

    profiles: [

      { position_1based: 1, rice: "18" },

      { position_1based: 2, rice: "invalid" },

      { position_1based: 3, rice: "" },

    ],

  });



  assert.deepEqual(data.potencyPoints.map((point) => point.value), [0.5, null, null]);

  assert.deepEqual(data.profilePoints("rice").map((point) => point.value), [18, null, null]);

  assert.deepEqual(contiguousSegments(data.profilePoints("rice")), [[

    { position: 1, value: 18 },

  ]]);

});



test("positionalProfileData reports unavailable data when no valid rows exist", () => {

  assert.deepEqual(

    positionalProfileData({ target_length: 10, windows: [], profiles: [] }),

    { available: false, reason: "This run did not record valid positional data." }

  );

  assert.equal(positionalProfileData({ windows: [{ start: 1 }] }).available, false);

});
