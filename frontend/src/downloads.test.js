import assert from "node:assert/strict";

import test from "node:test";



import { downloadFromBackend } from "./downloads.js";



function fakeDom() {

  const clicked = [];

  const link = {

    href: "",

    download: "",

    click: () => clicked.push({ href: link.href, download: link.download }),

    remove: () => {},

  };

  const doc = {

    body: { appendChild: () => {} },

    createElement: (tag) => {

      assert.equal(tag, "a");

      return link;

    },

  };

  const urlApi = {

    createObjectURL: (blob) => `blob:${blob}`,

    revokeObjectURL: () => {},

  };

  return { doc, urlApi, clicked };

}



test("downloadFromBackend saves the file under the given name", async () => {

  const { doc, urlApi, clicked } = fakeDom();

  const errors = [];

  const result = await downloadFromBackend(

    "http://127.0.0.1:8000/api/jobs/j1/export/report.pdf",

    "dsguard-report-j1.pdf",

    (message) => errors.push(message),

    {

      document: doc,

      URL: urlApi,

      fetch: async () => ({

        ok: true,

        blob: async () => "pdf-bytes",

      }),

    }

  );



  assert.deepEqual(result, { ok: true });

  assert.deepEqual(errors, []);

  assert.deepEqual(clicked, [

    { href: "blob:pdf-bytes", download: "dsguard-report-j1.pdf" },

  ]);

});



test("downloadFromBackend surfaces the backend message and saves nothing", async () => {

  const { doc, urlApi, clicked } = fakeDom();

  const errors = [];

  const result = await downloadFromBackend(

    "http://127.0.0.1:8000/api/jobs/legacy/export/report.pdf",

    "dsguard-report-legacy.pdf",

    (message) => errors.push(message),

    {

      document: doc,

      URL: urlApi,

      fetch: async () => ({

        ok: false,

        status: 404,

        json: async () => ({

          detail: "Re-run the screening to generate results.json.",

        }),

      }),

    }

  );



  assert.equal(result.ok, false);

  assert.equal(result.status, 404);

  assert.deepEqual(errors, ["Re-run the screening to generate results.json."]);

  assert.deepEqual(clicked, [], "an error must not be saved as a .pdf/.csv file");

});



test("downloadFromBackend falls back to a status message for a non-JSON error", async () => {

  const { doc, urlApi, clicked } = fakeDom();

  const errors = [];

  await downloadFromBackend("u", "f.csv", (message) => errors.push(message), {

    document: doc,

    URL: urlApi,

    fetch: async () => ({

      ok: false,

      status: 500,

      json: async () => {

        throw new Error("not json");

      },

    }),

  });



  assert.deepEqual(errors, ["Download failed (HTTP 500)."]);

  assert.deepEqual(clicked, []);

});



test("downloadFromBackend reports an unreachable backend", async () => {

  const { doc, urlApi, clicked } = fakeDom();

  const errors = [];

  const result = await downloadFromBackend(

    "u",

    "f.json",

    (message) => errors.push(message),

    {

      document: doc,

      URL: urlApi,

      fetch: async () => {

        throw new Error("fetch failed");

      },

    }

  );



  assert.equal(result.reason, "network");

  assert.deepEqual(errors, ["Could not reach the backend: fetch failed"]);

  assert.deepEqual(clicked, []);

});
