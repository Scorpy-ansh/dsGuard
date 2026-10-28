











export async function downloadFromBackend(url, filename, onError, deps = {}) {

  const fetchImpl = deps.fetch || globalThis.fetch;

  const doc = deps.document || globalThis.document;

  const urlApi = deps.URL || globalThis.URL;

  let response;



  try {

    response = await fetchImpl(url);

  } catch (networkError) {

    onError?.(`Could not reach the backend: ${networkError.message}`);

    return { ok: false, reason: "network" };

  }



  if (!response.ok) {

    let detail = `Download failed (HTTP ${response.status}).`;

    try {

      const payload = await response.json();

      if (payload?.detail) detail = payload.detail;

    } catch {



    }

    onError?.(detail);

    return { ok: false, reason: "http", status: response.status, detail };

  }



  const blob = await response.blob();

  const objectUrl = urlApi.createObjectURL(blob);

  const link = doc.createElement("a");

  link.href = objectUrl;

  link.download = filename;

  doc.body.appendChild(link);

  link.click();

  link.remove();

  urlApi.revokeObjectURL(objectUrl);

  return { ok: true };

}
