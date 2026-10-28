













































export const ORGANISM_NAMES = [

  {

    scientific: "Oryza sativa",

    common: "Rice",

    panelKeys: ["rice"],

    source: "NCBI Taxonomy 4530 (GenBank common name: Asian cultivated rice; also 'rice')",

  },

  {

    scientific: "Apis mellifera",

    common: "Honey bee",

    panelKeys: ["honeybee"],

    source: "NCBI Taxonomy 7460 (GenBank common name: honey bee)",

  },

  {

    scientific: "Nilaparvata lugens",

    common: "Brown planthopper",

    panelKeys: ["nilaparvata_lugens"],

    source: "NCBI Taxonomy 108931; EFSA Journal 2023;21(6):7999",

  },

  {

    scientific: "Cyrtorhinus lividipennis",

    common: "Mirid bug",

    panelKeys: ["cyrtorhinus_lividipennis"],

    source: "ICAR-CRRI Technology Bulletin 273 ('Common name: Mirid bug')",

  },

  {

    scientific: "Aedes aegypti",

    common: "Yellow fever mosquito",

    panelKeys: ["aedes_aegypti"],

    source: "NCBI Taxonomy 7159; ITIS TSN 126240; ECDC vector factsheet",

  },

];









const PANEL_LABELS = {

  spider: "Spider (reference species panel)",

};





const byPanelKey = new Map();

const byScientificName = new Map();



for (const record of ORGANISM_NAMES) {

  for (const key of record.panelKeys || []) {

    byPanelKey.set(key.toLowerCase(), record);

  }

  byScientificName.set(record.scientific.toLowerCase(), record);

}



function normalise(value) {

  return typeof value === "string" ? value.trim() : "";

}







function binomialIn(label) {

  const match = normalise(label).match(/\(([^()]+)\)/);

  return match ? match[1].trim() : "";

}







export function findOrganism(value) {

  const raw = normalise(value);

  if (!raw) return null;



  const lowered = raw.toLowerCase();

  if (byPanelKey.has(lowered)) return byPanelKey.get(lowered);

  if (byScientificName.has(lowered)) return byScientificName.get(lowered);



  const binomial = binomialIn(raw);

  if (binomial && byScientificName.has(binomial.toLowerCase())) {

    return byScientificName.get(binomial.toLowerCase());

  }



  return null;

}









export function organismLabel(value, fallback = "—") {

  const raw = normalise(value);

  if (!raw) return fallback;



  const record = findOrganism(raw);

  if (record) return record.common || record.scientific;



  return PANEL_LABELS[raw.toLowerCase()] || raw;

}









export function organismScientificName(value) {

  const record = findOrganism(value);

  return record ? record.scientific : "";

}









export function organismDisplay(value, fallback = "—") {

  const raw = normalise(value);

  const record = findOrganism(raw);



  return {

    primary: organismLabel(raw, fallback),

    secondary: record && record.common ? record.scientific : "",

  };

}
