import assert from "node:assert/strict";

import test from "node:test";

import { getApiBaseUrl } from "./apiConfig.js";

test("uses the configured frontend API base URL", () => {
  assert.equal(
    getApiBaseUrl({ VITE_API_BASE_URL: "https://api.example.com/" }),
    "https://api.example.com",
  );
});

test("falls back to the local API during development", () => {
  assert.equal(getApiBaseUrl({}), "http://127.0.0.1:8000");
});
