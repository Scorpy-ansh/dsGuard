













export function getHazardLimit(results) {

  const sources = [

    results?.analysis?.parameters?.max_hazard,

    results?.manifest?.params?.max_hazard,

  ];



  for (const value of sources) {

    if (value === null || value === undefined || value === "") continue;

    const number = Number(value);

    if (Number.isFinite(number)) return number;

  }



  return null;

}



export function isThresholdExceeded(row, results) {

  const status = String(row?.recommendation_status || "").toLowerCase();

  if (status === "fallback_threshold_exceeded") return true;



  const hazard = Number(row?.rec_hazard);

  const limit = getHazardLimit(results);



  return limit !== null && Number.isFinite(hazard) && hazard > limit;

}



export function screeningStatus(row, analysis) {

  const explicit = String(

    row?.screening_status || analysis?.screening_status || ""

  ).trim();

  if (

    explicit === "screened_candidate_available" ||

    explicit === "screened_no_candidate" ||

    explicit === "screening_incomplete"

  ) {

    return explicit;

  }

  if (row?.candidate || row?.recommendation_status === "within_threshold" ||

      row?.recommendation_status === "fallback_threshold_exceeded") {

    return "screened_candidate_available";

  }

  if (row?.status && String(row.status).startsWith("too_short")) {

    return "screened_no_candidate";

  }

  return "legacy_or_unavailable";

}



export function screeningStatusLabel(row, analysis) {

  return {

    screened_candidate_available: "Screened — candidate available",

    screened_no_candidate: "Screened — no candidate",

    screening_incomplete: "Screening incomplete",

    legacy_or_unavailable: "Legacy result — status unavailable",

  }[screeningStatus(row, analysis)];

}





















export function geneDisplayName(header) {

  if (header === null || header === undefined) return "Unknown gene";

  const alias = String(header).split("|")[0].trim();

  return alias || "Unknown gene";

}









export function jobDisplayName(job) {

  const raw = job?.target_filename || job?.target_label || "";

  return String(raw).trim() ? geneDisplayName(raw) : "Screening job";

}



export function formatHazard(value) {





  if (value === null || value === undefined || value === "") return "—";

  const number = Number(value);

  return Number.isFinite(number) ? number.toFixed(6) : "—";

}



export function getRecommendationLabel(row, results) {

  return screeningStatusLabel(row, results?.analysis);

}

































export const WINDOW_MIN = 21;

export const WINDOW_MAX = 501;

export const DEFAULT_WINDOW = 400;















export const WINDOW_PRESETS = [

  { label: "50\u2013101 bp", min: 50, max: 101, value: 50 },

  { label: "100\u2013201 bp", min: 100, max: 201, value: 100 },

  { label: "200\u2013301 bp", min: 200, max: 301, value: 200 },

  { label: "300\u2013401 bp", min: 300, max: 401, value: 300 },

  { label: "400\u2013501 bp", min: 400, max: 501, value: 400 },

];







export function presetForLength(length) {

  const value = Number(length);

  if (!Number.isFinite(value)) return null;

  return WINDOW_PRESETS.find((preset) => value >= preset.min && value <= preset.max) || null;

}







export function validateWindowLength(raw) {

  const text = String(raw ?? "").trim();

  if (!text) return { ok: false, error: "Enter a sequence length in bp." };

  if (!/^\d+$/.test(text)) {

    return {

      ok: false,

      error: `Sequence length must be a whole number of bp (got "${text}").`,

    };

  }



  const value = Number(text);

  if (value < WINDOW_MIN || value > WINDOW_MAX) {

    return {

      ok: false,

      error: `Sequence length must be between ${WINDOW_MIN} and ${WINDOW_MAX} bp (got ${value}).`,

    };

  }

  return { ok: true, value };

}



export function catalogGenesForOrganism(catalog, scientificName) {

  const organism = (catalog || []).find(

    (entry) => entry?.scientific_name === scientificName

  );

  return Array.isArray(organism?.genes) ? organism.genes : [];

}











export function paretoPlotPoints(

  points,

  { width = 520, height = 220, pad = 36 } = {}

) {

  const rows = (points || [])

    .map((point) => ({

      hazard: Number(point?.hazard),

      potency: Number(point?.potency),

      start: point?.start,

      end: point?.end,

      threshold_status: point?.threshold_status,

    }))

    .filter((row) => Number.isFinite(row.hazard) && Number.isFinite(row.potency));



  if (!rows.length) return [];



  const hazards = rows.map((row) => row.hazard);

  const potencies = rows.map((row) => row.potency);

  const minHazard = Math.min(...hazards);

  const maxHazard = Math.max(...hazards);

  const minPotency = Math.min(...potencies);

  const maxPotency = Math.max(...potencies);









  const scale = (value, min, max, low, high) =>

    max > min

      ? low + ((value - min) / (max - min)) * (high - low)

      : (low + high) / 2;



  return rows.map((row) => ({

    ...row,

    x: scale(row.hazard, minHazard, maxHazard, pad, width - pad),

    y: scale(row.potency, minPotency, maxPotency, height - pad, pad),

  }));

}



























export function candidateSequence(gene) {

  const sequence = gene?.candidate?.sequence;

  return typeof sequence === "string" && sequence.length > 0 ? sequence : "";

}













export function candidateView(gene) {

  const candidate = gene?.candidate || null;

  const sequence = candidateSequence(gene);



  if (sequence) {

    const requested = Number(candidate.requested_length_bp);



    return {

      available: true,

      candidateId: candidate.candidate_id || "",

      sequence,





      requestedLengthBp: Number.isFinite(requested) ? requested : null,

      actualLengthBp: sequence.length,

      start1Based: candidate.start_1based ?? null,

      end1Based: candidate.end_1based ?? null,

      strand: candidate.strand || "",

      ambiguousBases: candidate.ambiguous_bases ?? null,

      message: "",

    };

  }



  const recorded =

    Boolean(gene) && Object.prototype.hasOwnProperty.call(gene, "candidate");

  const error =

    typeof gene?.candidate_error === "string" ? gene.candidate_error.trim() : "";



  return {

    available: false,

    candidateId: "",

    sequence: "",

    requestedLengthBp: null,

    actualLengthBp: null,

    start1Based: null,

    end1Based: null,

    strand: "",

    ambiguousBases: null,

    message:

      error ||

      (recorded

        ? "The engine did not produce a candidate sequence for this gene in this run."

        : "This run does not record a candidate sequence (it predates candidate output). Re-run the screening to generate one."),

  };

}











export function targetGeneMetadata(gene, analysis) {

  const header = String(gene?.gene ?? "");

  const row = gene?.row || {};

  const accessionMatch = header.match(/(?:^|\|)accession=([^|]+)/i);

  const length = Number(gene?.length_bp);

  const sourceMetadata = gene?.source_metadata || gene?.metadata || {};

  const accession =

    gene?.accession ||

    gene?.accession_version ||

    row.accession ||

    row.accession_version ||

    sourceMetadata.accession ||

    sourceMetadata.accession_version ||

    (accessionMatch ? accessionMatch[1] : "");

  const organism =

    gene?.organism ||

    row.organism ||

    sourceMetadata.organism ||

    analysis?.input?.organism ||

    analysis?.target_metadata?.organism ||

    "";



  return {

    gene: header,

    displayName: header ? geneDisplayName(header) : "Unknown gene",

    accession: String(accession || "").trim(),

    lengthBp: Number.isFinite(length) ? length : null,

    inputLabel: analysis?.input?.label || "",

    inputType: analysis?.input?.type || "",

    organism: String(organism || "").trim(),

    sequence: typeof gene?.target_sequence === "string" ? gene.target_sequence : "",

  };

}







export function csvField(value) {

  return `"${String(value ?? "").replace(/"/g, '""')}"`;

}



















export function geneSummaryCsv(results) {

  const rows = results?.summary || [];

  if (!rows.length) return "";



  const maxHazard = getHazardLimit(results) ?? "";

  const genes = results?.analysis?.genes || results?.genes || [];

  const geneEntries = new Map(

    genes.map((gene) => [gene.gene || gene.row?.gene || "", gene])

  );

  const candidateOf = (row) => geneEntries.get(row.gene)?.candidate || {};



  const columns = [







    ["Gene", (row) => row.gene],

    ["Status", (row) => row.status],

    ["Recommendation Status", (row) => row.recommendation_status],

    ["Recommended Hazard Score", (row) => row.rec_hazard],

    ["Configured Maximum Hazard", () => maxHazard],

    ["Length (bp)", (row) => row.length_bp],

    ["GC (%)", (row) => row.gc_pct],

    ["Maximum Match (nt)", (row) => row.rec_max_match],

    ["Exact-Match Hits", (row) => row.total_hits],

    ["Candidate ID", (row) => candidateOf(row).candidate_id],

    ["Candidate Requested Length (bp)", (row) => candidateOf(row).requested_length_bp],

    ["Candidate Actual Length (bp)", (row) => candidateOf(row).actual_length_bp],

    ["Candidate Sequence", (row) => candidateOf(row).sequence],

    ["Candidate Note", (row) => {

      const entry = geneEntries.get(row.gene);





      return entry ? candidateView(entry).message : "";

    }],

  ];



  return [

    columns.map(([label]) => csvField(label)).join(","),

    ...rows.map((row) =>

      columns.map(([, value]) => csvField(value(row))).join(",")

    ),

  ].join("\r\n");

}
