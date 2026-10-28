





import test from "node:test";

import assert from "node:assert/strict";



import { copyFailureMessage, copyText } from "./clipboard.js";



function fakeClipboard({ reject = false } = {}) {

  const writes = [];

  return {

    writes,

    api: {

      writeText: async (value) => {

        if (reject) throw new Error("denied by user");

        writes.push(value);

      },

    },

  };

}



test("copyText writes exactly the candidate sequence it is given", async () => {

  const { writes, api } = fakeClipboard();

  const candidate = "ACGTACGTACACGTACGTACAC";

  const target = `${candidate}TTTTTTTTTTTTTTTTTTTT`;



  const result = await copyText(candidate, { clipboard: api });



  assert.deepEqual(result, { ok: true, length: candidate.length });

  assert.deepEqual(writes, [candidate], "the target sequence must never be copied");

  assert.notEqual(writes[0], target);

});



test("copyText reports an empty value and writes nothing", async () => {

  const { writes, api } = fakeClipboard();



  assert.deepEqual(await copyText("", { clipboard: api }), {

    ok: false,

    reason: "empty",

  });

  assert.deepEqual(await copyText(null, { clipboard: api }), {

    ok: false,

    reason: "empty",

  });

  assert.deepEqual(await copyText(undefined, { clipboard: api }), {

    ok: false,

    reason: "empty",

  });

  assert.deepEqual(writes, []);

});



test("copyText reports an unsupported clipboard instead of throwing", async () => {

  assert.deepEqual(await copyText("ACGT", { clipboard: null }), {

    ok: false,

    reason: "unsupported",

  });

  assert.deepEqual(await copyText("ACGT", { clipboard: {} }), {

    ok: false,

    reason: "unsupported",

  });

});



test("copyText surfaces a denied clipboard permission", async () => {

  const { api } = fakeClipboard({ reject: true });



  const result = await copyText("ACGT", { clipboard: api });



  assert.equal(result.ok, false);

  assert.equal(result.reason, "denied");

  assert.ok(result.error instanceof Error);

});



test("copyFailureMessage explains each failure mode distinctly", () => {

  const messages = ["empty", "unsupported", "denied", "anything-else"].map(

    copyFailureMessage

  );



  assert.equal(new Set(messages).size, messages.length);

  for (const message of messages) {

    assert.ok(message.trim().length > 0);





    assert.match(message, /copy|clipboard/i);

  }

});
