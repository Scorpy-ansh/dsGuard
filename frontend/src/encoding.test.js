







import test from "node:test";

import assert from "node:assert/strict";

import { readFileSync, readdirSync } from "node:fs";

import { dirname, join } from "node:path";

import { fileURLToPath } from "node:url";



const srcDir = dirname(fileURLToPath(import.meta.url));







const sourceFiles = readdirSync(srcDir).filter(

  (name) => /\.(jsx?|css)$/.test(name) && !name.endsWith(".test.js")

);





const MOJIBAKE_SENTINELS = ["â€", "â†", "â‰", "Â·", "Ã¢", "Ã©", "Ã¨", "Ã¶", "Ã¼", "Ã¤"];



test("no dashboard source file contains mojibake byte sequences", () => {

  for (const name of sourceFiles) {

    const text = readFileSync(join(srcDir, name), "utf8");

    for (const sentinel of MOJIBAKE_SENTINELS) {

      assert.ok(

        !text.includes(sentinel),

        `${name} still contains the mojibake sequence ${JSON.stringify(sentinel)}`

      );

    }

  }

});



test("App.jsx has no UTF-8 byte-order mark", () => {

  const text = readFileSync(join(srcDir, "App.jsx"), "utf8");



  assert.notEqual(text.charCodeAt(0), 0xfeff);

  assert.ok(text.startsWith("import"), "App.jsx should start with its first import");

});



test("App.jsx renders the threshold symbols and chart units as real Unicode", () => {

  const text = readFileSync(join(srcDir, "App.jsx"), "utf8");





  assert.ok(text.includes("off-target hazard \u2192"), "hazard axis arrow missing");

  assert.ok(text.includes("potency \u2192"), "potency axis arrow missing");

  assert.ok(text.includes("\u2265"), "the >= threshold symbol is missing or mojibaked");

  assert.ok(text.includes("\u2014"), "the em dash is missing or mojibaked");

});
