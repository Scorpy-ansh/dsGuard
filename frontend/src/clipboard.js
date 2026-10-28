













export async function copyText(value, deps = {}) {

  const text = value === null || value === undefined ? "" : String(value);

  if (!text) return { ok: false, reason: "empty" };



  const clipboard = deps.clipboard || globalThis.navigator?.clipboard;

  if (!clipboard || typeof clipboard.writeText !== "function") {

    return { ok: false, reason: "unsupported" };

  }



  try {

    await clipboard.writeText(text);

    return { ok: true, length: text.length };

  } catch (error) {

    return { ok: false, reason: "denied", error };

  }

}







export function copyFailureMessage(reason) {

  switch (reason) {

    case "empty":

      return "There is no candidate sequence to copy.";

    case "unsupported":

      return "This browser blocked clipboard access. Select the sequence and copy it manually.";

    case "denied":

      return "Clipboard permission was denied. Select the sequence and copy it manually.";

    default:

      return "The candidate sequence could not be copied. Select it and copy manually.";

  }

}
