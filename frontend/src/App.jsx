import { useEffect, useRef, useState } from "react";

import {

  Activity,

  ArrowUpRight,

  BarChart3,

  CheckCircle2,

  CircleHelp,

  Clock3,

  Copy,

  Dna,

  Download,

  FileSearch,

  FlaskConical,

  LayoutDashboard,

  Settings2,

  ShieldCheck,

  Upload,

  RefreshCw,

  LoaderCircle,

  AlertCircle,

  FileText,

} from "lucide-react";

import "./App.css";

import {

  DEFAULT_WINDOW,

  WINDOW_MAX,

  WINDOW_MIN,

  WINDOW_PRESETS,

  candidateView,

  catalogGenesForOrganism,

  formatHazard,

  geneDisplayName,

  geneSummaryCsv,

  getHazardLimit,

  getRecommendationLabel,

  isThresholdExceeded,

  screeningStatusLabel,

  jobDisplayName,

  paretoPlotPoints,

  presetForLength,

  targetGeneMetadata,

  validateWindowLength,

} from "./analysis";

import { copyFailureMessage, copyText } from "./clipboard";

import { downloadFromBackend } from "./downloads";
import { getApiBaseUrl } from "./apiConfig";
import { organismDisplay, organismLabel } from "./organisms";

import {

  contiguousSegments,

  numericRange,

  positionalProfileData,

} from "./positionalProfile";









const API = getApiBaseUrl();

const API_HOST = API.replace(/^https?:\/\//, "");



function downloadResultsCSV(results) {

  const csv = geneSummaryCsv(results);

  if (!csv) return;



  const blob = new Blob(["\uFEFF", csv], {

    type: "text/csv;charset=utf-8;",

  });



  const url = URL.createObjectURL(blob);

  const link = document.createElement("a");

  link.href = url;

  link.download = `dsguard-gene-summary-${results.job_id || "screening"}.csv`;

  document.body.appendChild(link);

  link.click();

  link.remove();

  URL.revokeObjectURL(url);

}







function downloadDetailedHits(jobId, onError) {

  if (!jobId) return;



  downloadFromBackend(

    `${API}/api/jobs/${encodeURIComponent(jobId)}/hits/download`,

    `dsguard-hits-${jobId}.csv`,

    onError

  );

}







function downloadExport(jobId, kind, filename, onError) {

  if (!jobId) return;



  downloadFromBackend(

    `${API}/api/jobs/${encodeURIComponent(jobId)}/export/${kind}`,

    filename,

    onError

  );

}







function downloadCandidateFasta(jobId, onError) {

  if (!jobId) return;



  downloadFromBackend(

    `${API}/api/jobs/${encodeURIComponent(jobId)}/export/candidates.fasta`,

    `dsguard-candidates-${jobId}.fasta`,

    onError

  );

}







function flattenRankedOfftargets(analysis) {

  const rows = [];



  (analysis?.genes || []).forEach((gene) => {

    const geneName = gene.gene || gene.row?.gene || "";

    (gene.ranked_offtargets || []).forEach((row) => {

      rows.push({ gene: geneName, ...row });

    });

  });



  rows.sort(

    (a, b) =>

      (Number(b.hazard) || 0) - (Number(a.hazard) || 0) ||

      (Number(b.match_length_bp) || 0) - (Number(a.match_length_bp) || 0) ||

      String(a.reference_panel || "").localeCompare(

        String(b.reference_panel || "")

      ) ||

      String(a.matched_organism || "").localeCompare(

        String(b.matched_organism || "")

      ) ||

      String(a.transcript_or_reference_id || "").localeCompare(

        String(b.transcript_or_reference_id || "")

      ) ||

      (Number(a.target_start) || 0) - (Number(b.target_start) || 0)

  );



  return rows.map((row, index) => ({ ...row, job_rank: index + 1 }));

}



function ParetoScatter({ points }) {





  const plotted = paretoPlotPoints(points);



  if (!plotted.length) return null;



  const width = 520;

  const height = 220;

  const pad = 36;



  return (

    <svg

      className="pareto-plot"

      viewBox={`0 0 ${width} ${height}`}

      role="img"

      aria-label="Pareto front: potency against off-target hazard"

    >

      <line x1={pad} y1={height - pad} x2={width - pad} y2={height - pad} />

      <line x1={pad} y1={pad} x2={pad} y2={height - pad} />

      <text x={width / 2} y={height - 6} textAnchor="middle">

        off-target hazard →

      </text>

      <text x={12} y={height / 2} transform={`rotate(-90 12 ${height / 2})`}>

        potency →

      </text>

      {plotted.map((point, index) => (

        <circle

          key={`${point.start}-${index}`}

          cx={point.x}

          cy={point.y}

          r={4}

          className="pareto-point"

        >

          <title>

            {`${point.start}-${point.end}: potency ${point.potency}, hazard ${point.hazard} (${point.threshold_status})`}

          </title>

        </circle>

      ))}

    </svg>

  );

}



function PositionalProfile({ profile, candidate, mode, species, layers, onSpeciesChange, onModeChange, onToggleLayer }) {

  const data = positionalProfileData(profile);

  if (!data.available) {

    return (

      <div className="notice">

        <AlertCircle size={17} />

        <p>{data.reason || profile?.reason || "Positional data is unavailable for this run."}</p>

      </div>

    );

  }



  const panelSpecies = data.species;

  const selectedSpecies = species !== "all" && panelSpecies.includes(species) ? species : null;

  const riskPoints = data.profilePoints(selectedSpecies || panelSpecies[0]).map((point, index) => {

    if (selectedSpecies) return point;

    const values = panelSpecies

      .map((name) => data.profilePoints(name)[index]?.value)

      .filter((value) => value !== null && value !== undefined);

    return { ...point, value: values.length ? Math.max(...values) : null };

  });

  const potencyRange = numericRange(data.potencyPoints);

  const rawRiskRange = numericRange(riskPoints);

  const riskRange = rawRiskRange && {

    min: Math.min(0, rawRiskRange.min),

    max: Math.max(21, rawRiskRange.max),

  };

  const selectedStart = Number(candidate?.start1Based);

  const selectedEnd = Number(candidate?.end1Based);

  const xCoordinate = (position) => Math.max(0, Math.min(100, ((position - 1) / Math.max(1, data.targetLength - 1)) * 100));

  const x = (position) => `${xCoordinate(position)}%`;

  const y = (value, range) => {

    if (value === null || !range) return null;

    if (range.max === range.min) return "50%";

    return `${100 - ((value - range.min) / (range.max - range.min)) * 100}%`;

  };

  const line = (points, range, className) => contiguousSegments(points).map((segment, index) => (

    <polyline

      key={index}

      className={className}

      points={segment.map((point) => `${xCoordinate(point.position) * 5.2},${y(point.value, range).replace("%", "")}`).join(" ")}

      vectorEffect="non-scaling-stroke"

    />

  ));

  const profileAxis = `${data.targetLength} nt`;



  return (

    <div className="profile-content">

      <div className="profile-controls">

        <div className="mode-buttons" role="group" aria-label="Profile display mode">

          {["dual", "combined"].map((value) => (

            <button type="button" key={value} className={`mode-button ${mode === value ? "active" : ""}`} onClick={() => onModeChange(value)}>

              {value === "dual" ? "Dual Tracks" : "Combined"}

            </button>

          ))}

        </div>

        <label className="profile-select">

          Species

          <select value={species} onChange={(event) => onSpeciesChange(event.target.value)}>

            <option value="all">All Panel Species (Worst-Case Peak)</option>

            {panelSpecies.map((name) => <option key={name} value={name}>{panelOptionLabel(name)}</option>)}

          </select>

        </label>

      </div>

      <div className={`profile-chart ${mode}`}>

        <div className="profile-axis">

          <span>1</span><span>{Math.round(data.targetLength / 2)}</span><span>{profileAxis}</span>

        </div>

        <div className="profile-track">

          <div className="profile-track-label">On-target efficacy</div>

          <svg className="profile-lines" viewBox="0 0 520 100" preserveAspectRatio="none" aria-label="On-target efficacy by window start">

            {layers.potency && line(data.potencyPoints, potencyRange, "potency-line")}

          </svg>

          {!potencyRange && <p className="profile-unavailable">On-target efficacy is unavailable.</p>}

          <div className="profile-points" aria-hidden="true">

            {layers.potency && data.potencyPoints.map((point, index) => point.value === null ? null : (

              <span key={index} className="profile-point potency-point" style={{ left: x(point.position), bottom: y(point.value, potencyRange) }} title={`Window starting ${point.position}: potency ${point.value}`} />

            ))}

          </div>

        </div>

        {mode === "dual" && <div className="profile-track">

          <div className="profile-track-label">Non-target continuous match risk</div>

          <svg className="profile-lines" viewBox="0 0 520 100" preserveAspectRatio="none" aria-label="Non-target match risk by target position">

            {layers.risk && line(riskPoints, riskRange, "risk-line")}

          </svg>

          {!riskRange && <p className="profile-unavailable">Per-position match risk is unavailable.</p>}

          <div className="profile-threshold" style={{ bottom: y(21, riskRange) || "0%" }}><b>21 nt</b></div>

        </div>}

        {mode === "combined" && <div className="profile-track combined-track">

          <div className="profile-track-label">Combined: potency (blue) · match risk (red)</div>

          <svg className="profile-lines" viewBox="0 0 520 100" preserveAspectRatio="none" aria-label="Combined positional profile">

            {layers.risk && line(riskPoints, riskRange, "risk-line")}

            {layers.potency && line(data.potencyPoints, potencyRange, "potency-line")}

          </svg>

          {!riskRange && !potencyRange && <p className="profile-unavailable">No valid positional series are available.</p>}

          <div className="profile-threshold" style={{ bottom: y(21, riskRange) || "0%" }}><b>21 nt</b></div>

        </div>}

        {Number.isFinite(selectedStart) && Number.isFinite(selectedEnd) && (

          <div className="profile-window" style={{ left: x(selectedStart), width: `${Math.max(0.5, ((selectedEnd - selectedStart + 1) / Math.max(1, data.targetLength)) * 100)}%` }} title={`Selected construct: ${selectedStart}-${selectedEnd}`} />

        )}



      </div>

      <div className="profile-legend">

        {[["potency", "On-target potency", "potency-key"], ["accessibility", "mRNA accessibility", "accessibility-key"], ["risk", "Non-target match risk", "risk-key"]].map(([key, label, className]) => (

          <button type="button" className={`legend-toggle ${layers[key] ? "enabled" : ""}`} key={key} onClick={() => onToggleLayer(key)}>

            <span className={`legend-swatch ${className}`} />{label}

          </button>

        ))}

        <span className="legend-note">Selected construct</span>

      </div>

      <p className="panel-note">Values are the engine&apos;s recorded positional outputs. Accessibility is unavailable unless the engine provides it; no substitute curve is shown.</p>

    </div>

  );

}



function SequenceViewer({ target, candidate, hits }) {

  if (!target?.sequence) {

    return <div className="notice"><AlertCircle size={17} /><p>The complete target sequence is unavailable for this historical run.</p></div>;

  }

  const sequence = target.sequence.toUpperCase();

  const start = Number(candidate?.start1Based);

  const end = Number(candidate?.end1Based);

  const matchRanges = (hits || [])

    .filter((hit) => hit.gene_start && hit.gene_end)

    .map((hit) => ({ start: Number(hit.gene_start), end: Number(hit.gene_end) }));

  const lineLength = 80;

  const lines = [];

  for (let offset = 0; offset < sequence.length; offset += lineLength) {

    const line = sequence.slice(offset, offset + lineLength);

    lines.push(

      <div className="sequence-line" key={offset}>

        <span className="sequence-coordinate">{String(offset + 1).padStart(7, " ")}</span>

        <span>{[...line].map((base, index) => {

          const position = offset + index + 1;

          const inCandidate = position >= start && position <= end;

          const inHit = matchRanges.some((range) => position >= range.start && position <= range.end);

          return <span key={position} className={inCandidate ? "candidate-base" : inHit ? "hit-base" : ""} title={`Position ${position}`}>{base}</span>;

        })}</span>

      </div>

    );

  }

  return (

    <div className="sequence-viewer">

      <div className="sequence-viewer-meta">{sequence.length} nt · selected construct {Number.isFinite(start) ? `${start}-${end}` : "Unavailable"}</div>

      <div className="sequence-viewer-scroll">{lines}</div>

      <div className="sequence-legend"><span><i className="candidate-key" /> Selected construct</span><span><i className="hit-key" /> Off-target match</span></div>

    </div>

  );

}



function getPanelLabel(job) {

  const rawNames = job?.panel_names;



  const names = Array.isArray(rawNames)

    ? rawNames

    : typeof rawNames === "string"

      ? [rawNames]

      : [];



  const validNames = names.filter(

    (name) =>

      typeof name === "string" &&

      name.trim() !== "" &&

      name.trim().toLowerCase() !== "string"

  );



  if (validNames.length > 0) {





    return validNames.map((name) => organismLabel(name)).join(", ");

  }



  const manifestPanel = job?.manifest?.panel;



  if (

    manifestPanel &&

    typeof manifestPanel === "object" &&

    Object.keys(manifestPanel).length > 0

  ) {

    return Object.keys(manifestPanel).map((name) => organismLabel(name)).join(", ");

  }



  return job?.panel_name || job?.panel_filename || "—";

}











function OrganismCell({ value, fallback = "—" }) {

  const { primary, secondary } = organismDisplay(value, fallback);



  return (

    <>

      <span>{primary}</span>

      {secondary ? (

        <>

          <br />

          <span className="detail-note">{secondary}</span>

        </>

      ) : null}

    </>

  );

}







function panelOptionLabel(value) {

  const { primary, secondary } = organismDisplay(value, value);

  return secondary ? `${primary} (${secondary})` : primary;

}



function App() {

  const [activeTab, setActiveTab] = useState("Overview");

  const [health, setHealth] = useState(null);

  const [jobs, setJobs] = useState([]);

  const [loading, setLoading] = useState(true);

  const [error, setError] = useState("");



  const [selectedJob, setSelectedJob] = useState(null);

  const [results, setResults] = useState(null);



  const [hits, setHits] = useState([]);

  const [hitsLoading, setHitsLoading] = useState(false);

  const [profile, setProfile] = useState(null);

  const [profileMode, setProfileMode] = useState("dual");

  const [profileSpecies, setProfileSpecies] = useState("all");

  const [profileLayers, setProfileLayers] = useState({

    potency: true,

    accessibility: true,

    risk: true,

  });

  const [hitsPage, setHitsPage] = useState(1);



  const [geneReports, setGeneReports] = useState([]);

  const [selectedReportGene, setSelectedReportGene] = useState("");

  const [geneReport, setGeneReport] = useState("");

  const [reportLoading, setReportLoading] = useState(false);



  const [submitting, setSubmitting] = useState(false);

  const [notice, setNotice] = useState("");



  const [targetMode, setTargetMode] = useState("file");

  const [targetFile, setTargetFile] = useState(null);

  const [targetQuery, setTargetQuery] = useState("");

  const [targetSequence, setTargetSequence] = useState("");

  const [targetOrganism, setTargetOrganism] = useState("Nilaparvata lugens");

  const [targetCatalog, setTargetCatalog] = useState([]);

  const [panelKey, setPanelKey] = useState("rice");

  const [referencePanels, setReferencePanels] = useState([]);

  const [paretoGene, setParetoGene] = useState("");

  const [candidateGene, setCandidateGene] = useState("");

  const [targetCopyNotice, setTargetCopyNotice] = useState("");









  const [copyNotice, setCopyNotice] = useState({ key: "", text: "" });

  const copyNoticeTimerRef = useRef(0);









  const resultsRequestRef = useRef(0);

  const reportRequestRef = useRef(0);



  const [params, setParams] = useState({

    window: DEFAULT_WINDOW,

    k: 11,

    l50: 21,

    kappa: 2,

    max_hazard: 0.05,

  });







  const [windowMode, setWindowMode] = useState(String(DEFAULT_WINDOW));







  const windowBand = presetForLength(params.window);

  const targetGenes = catalogGenesForOrganism(targetCatalog, targetOrganism);

  const selectedCatalogGene = targetGenes.find(

    (gene) => gene.gene_name === targetQuery || gene.accession === targetQuery

  );



  async function loadDashboard() {

    try {

      const [healthResponse, jobsResponse] = await Promise.all([

        fetch(`${API}/api/health`),

        fetch(`${API}/api/jobs?limit=50`),

      ]);



      if (!healthResponse.ok || !jobsResponse.ok) {

        throw new Error("The backend returned an error.");

      }



      setHealth(await healthResponse.json());

      setJobs(await jobsResponse.json());

      setError("");

    } catch {

      setError(

        `Cannot reach the backend. Check that FastAPI is running at ${API_HOST}.`

      );

    } finally {

      setLoading(false);

    }

  }



  useEffect(() => {





    const initial = setTimeout(loadDashboard, 0);

    fetch(`${API}/api/panels`)

      .then((response) => response.ok ? response.json() : [])

      .then((data) => setReferencePanels(Array.isArray(data) ? data : []))

      .catch(() => setReferencePanels([]));

    fetch(`${API}/api/targets`)

      .then((response) => response.ok ? response.json() : { organisms: [] })

      .then((data) => setTargetCatalog(Array.isArray(data.organisms) ? data.organisms : []))

      .catch(() => setTargetCatalog([]));

    const timer = setInterval(loadDashboard, 5000);

    return () => {

      clearTimeout(initial);

      clearInterval(timer);

    };

  }, []);





  useEffect(() => () => clearTimeout(copyNoticeTimerRef.current), []);



  async function loadGeneReport(jobId, geneName) {

    if (!jobId || !geneName) return;



    const requestId = (reportRequestRef.current += 1);



    setSelectedReportGene(geneName);

    setGeneReport("");

    setReportLoading(true);

    setNotice("");



    try {

      const response = await fetch(

        `${API}/api/jobs/${encodeURIComponent(jobId)}/reports/${encodeURIComponent(geneName)}`

      );

      const data = await response.json();







      if (requestId !== reportRequestRef.current) return;



      if (!response.ok) {

        throw new Error(data.detail || "Could not load gene report.");

      }









      setGeneReport(

        String(data.report || "").replaceAll("\u00e2\u20ac\u201d", "\u2014")

      );

    } catch (err) {

      if (requestId === reportRequestRef.current) {

        setNotice(err.message);

      }

    } finally {

      if (requestId === reportRequestRef.current) {

        setReportLoading(false);

      }

    }

  }



  async function deleteJob(job) {

    const confirmed = window.confirm(

      `Delete screening "${jobDisplayName(job)}"?\n\nThis will permanently remove its history record, uploaded files, and generated results.`

    );



    if (!confirmed) return;



    setError("");

    setNotice("");



    try {

      const response = await fetch(

        `${API}/api/jobs/${encodeURIComponent(job.job_id)}`,

        { method: "DELETE" }

      );



      const data = await response.json().catch(() => ({}));



      if (!response.ok) {

        throw new Error(data.detail || "Could not delete this screening.");

      }









      if (selectedJob === job.job_id) {

        resultsRequestRef.current += 1;

        setSelectedJob(null);

        setResults(null);

        setHits([]);

        setProfile(null);

        setHitsPage(1);

        setGeneReports([]);

        setSelectedReportGene("");

        setGeneReport("");

      }



      await loadDashboard();

      setNotice("Screening history entry deleted successfully.");

    } catch (err) {

      setError(err.message || "Could not delete this screening.");

    }

  }





  async function openJob(jobId) {





    const requestId = (resultsRequestRef.current += 1);



    setSelectedJob(jobId);

    setResults(null);

    setHits([]);

    setGeneReports([]);

    setSelectedReportGene("");

    setGeneReport("");

    setHitsLoading(true);

    setActiveTab("Results");

    setNotice("");



    try {

      const response = await fetch(

        `${API}/api/jobs/${encodeURIComponent(jobId)}/results`

      );

      const data = await response.json();



      if (requestId !== resultsRequestRef.current) return;



      if (!response.ok) {

        throw new Error(data.detail || "Could not load results.");

      }



      setResults(data);



      if (data.status === "completed") {

        const [hitsResult, reportsResult] = await Promise.allSettled([

          fetch(`${API}/api/jobs/${encodeURIComponent(jobId)}/hits`),

          fetch(`${API}/api/jobs/${encodeURIComponent(jobId)}/reports`),

        ]);



        if (requestId !== resultsRequestRef.current) return;



        if (hitsResult.status === "fulfilled") {

          const hitsResponse = hitsResult.value;

          const hitsData = await hitsResponse.json();



          if (hitsResponse.ok) {

            setHits(Array.isArray(hitsData.hits) ? hitsData.hits : []);

          } else {

            setNotice(

              hitsData.detail || "Could not load detailed hit records."

            );

          }

        } else {

          setNotice("Could not connect to the detailed-hits endpoint.");

        }



        if (reportsResult.status === "fulfilled") {

          const reportsResponse = reportsResult.value;

          const reportsData = await reportsResponse.json();



          if (reportsResponse.ok) {

            setGeneReports(

              Array.isArray(reportsData.reports) ? reportsData.reports : []

            );

          } else {

            setNotice(

              reportsData.detail || "Could not load the gene report list."

            );

          }

        } else {

          setNotice("Could not connect to the gene reports endpoint.");

        }

      }

    } catch (err) {

      if (requestId === resultsRequestRef.current) {

        setNotice(err.message);

      }

    } finally {





      if (requestId === resultsRequestRef.current) {

        setHitsLoading(false);

      }

    }

  }



  async function submitScreening(event) {

    event.preventDefault();

    setNotice("");

    if (targetMode === "file" && !targetFile) {

      setNotice("Choose a target FASTA file.");

      return;

    }

    if (targetMode === "ncbi" && !targetQuery.trim()) {

      setNotice("Enter an NCBI accession ID or gene name.");

      return;

    }

    if (targetMode === "ncbi" && selectedCatalogGene &&

        selectedCatalogGene.accession_status !== "verified") {

      setNotice(

        `${selectedCatalogGene.gene_name} has no verified accession. Supply a FASTA sequence or resolve the target before screening.`

      );

      return;

    }

    if (targetMode === "sequence" && !targetSequence.trim()) {

      setNotice("Paste a FASTA or nucleotide sequence.");

      return;

    }

    const selectedPanel = referencePanels.find((panel) => panel.key === panelKey);

    if (selectedPanel && !selectedPanel.available) {

      setNotice(`The ${panelOptionLabel(selectedPanel.label)} reference panel is not installed in the backend yet.`);

      return;

    }





    const windowCheck = validateWindowLength(params.window);

    if (!windowCheck.ok) {

      setNotice(windowCheck.error);

      return;

    }



    const form = new FormData();

    form.append("target_mode", targetMode);

    if (targetMode === "file" && targetFile) form.append("target_file", targetFile);

    if (targetMode === "ncbi") {

      form.append("target_query", targetQuery.trim());

      form.append("target_organism", targetOrganism.trim());

    }

    if (targetMode === "sequence") form.append("target_sequence", targetSequence);

    form.append("panel_key", panelKey);

    Object.entries(params).forEach(([key, value]) => form.append(key, String(value)));

    setSubmitting(true);

    try {

      const response = await fetch(`${API}/api/screen`, { method: "POST", body: form });

      const data = await response.json();

      if (!response.ok) throw new Error(data.detail || "Screening submission failed.");

      setNotice(`Screening submitted. Job ID: ${data.job_id}`);

      setTargetFile(null);

      setTargetQuery("");

      setTargetSequence("");

      event.target.reset();

      setActiveTab("Screening history");

      await loadDashboard();

    } catch (err) {

      setNotice(err.message);

    } finally {

      setSubmitting(false);

    }

  }





  const navItems = [

    { label: "Overview", icon: LayoutDashboard },

    { label: "New screening", icon: FlaskConical },

    { label: "Screening history", icon: Clock3 },

    { label: "Results", icon: BarChart3 },

    { label: "Detailed hits", icon: FileSearch },

  ];



  const latestCompleted = jobs.find((job) => job.status === "completed");

  const latestSummary = latestCompleted?.summary || [];

  const reportedHits = latestSummary.reduce(

    (total, row) => total + (Number(row.total_hits) || 0),

    0

  );



  const analysis = results?.analysis || null;

  const rankedOfftargets = analysis ? flattenRankedOfftargets(analysis) : [];

  const searchedPanels = analysis?.reference_panels?.searched || [];

  const skippedPanels = analysis?.reference_panels?.skipped || [];

  const preventive = analysis?.preventive_measures || null;

  const runThreshold = getHazardLimit(results);

  const paretoGenes = (analysis?.genes || []).filter(

    (gene) => (gene.pareto || []).length > 0

  );

  const activeParetoGene =

    paretoGenes.find((gene) => gene.gene === paretoGene) ||

    paretoGenes[0] ||

    null;













  const candidateGenes = analysis?.genes || [];

  const activeCandidateGene =

    candidateGenes.find((gene) => gene.gene === candidateGene) ||

    candidateGenes[0] ||

    null;

  const activeCandidate = activeCandidateGene

    ? candidateView(activeCandidateGene)

    : null;

  const activeTarget = activeCandidateGene

    ? targetGeneMetadata(activeCandidateGene, analysis)

    : null;

  const activeGeneHits = hits.filter((hit) => geneDisplayName(hit.gene) === geneDisplayName(activeCandidateGene?.gene));

  const hitsPageSize = 25;

  const hitsPageCount = Math.max(1, Math.ceil(hits.length / hitsPageSize));

  const visibleHits = hits.slice((hitsPage - 1) * hitsPageSize, hitsPage * hitsPageSize);







  const copyNoticeKey = `${selectedJob ?? ""}|${activeCandidateGene?.gene ?? ""}`;

  const copyNoticeText =

    activeCandidateGene && copyNotice.key === copyNoticeKey ? copyNotice.text : "";



  useEffect(() => {

    if (!selectedJob || !activeCandidateGene?.gene || activeTab !== "Results") return;

    let cancelled = false;

    fetch(`${API}/api/jobs/${encodeURIComponent(selectedJob)}/profiles/${encodeURIComponent(activeCandidateGene.gene)}`)

      .then((response) => response.json().then((data) => ({ response, data })))

      .then(({ response, data }) => {

        if (cancelled) return;

        setProfile(response.ok ? data : { available: false, reason: data.detail || "Positional data unavailable." });

      })

      .catch(() => {

        if (!cancelled) setProfile({ available: false, reason: "Positional data could not be loaded." });

      });

    return () => { cancelled = true; };

  }, [selectedJob, activeCandidateGene?.gene, activeTab]);











  async function copyCandidateSequence(gene) {

    const view = candidateView(gene);

    const key = `${selectedJob ?? ""}|${gene?.gene ?? ""}`;



    if (!view.available) {

      setCopyNotice({ key, text: copyFailureMessage("empty") });

      return;

    }



    const outcome = await copyText(view.sequence);

    setCopyNotice({

      key,

      text: outcome.ok

        ? `Copied ${view.sequence.length} bp candidate sequence.`

        : copyFailureMessage(outcome.reason),

    });



    clearTimeout(copyNoticeTimerRef.current);

    copyNoticeTimerRef.current = setTimeout(

      () => setCopyNotice({ key: "", text: "" }),

      5000

    );

  }



  async function copyTargetSequence(gene) {

    const sequence = typeof gene?.target_sequence === "string" ? gene.target_sequence : "";

    if (!sequence) {

      setTargetCopyNotice("The original target sequence is not stored for this historical run.");

      return;

    }

    const outcome = await copyText(sequence);

    setTargetCopyNotice(

      outcome.ok

        ? `Copied ${sequence.length} bp original target sequence.`

        : copyFailureMessage(outcome.reason)

    );

  }



  function goTo(tab) {

    setActiveTab(tab);

    setNotice("");



    if (tab !== "Results" && tab !== "Detailed hits") {

      setSelectedJob(null);

      setResults(null);

      setHits([]);

      setGeneReports([]);

      setSelectedReportGene("");

      setGeneReport("");

    }

  }



  return (

    <div className="app-shell">

      <aside className="sidebar">

        <div className="brand">

          <div className="brand-mark">

            <Dna size={23} />

          </div>

          <div>

            <h1>dsGuard</h1>

            <span>SEQUENCE SAFETY PLATFORM</span>

          </div>

        </div>



        <div className="workspace-label">WORKSPACE</div>



        <nav className="nav-list">

          {navItems.map(({ label, icon: Icon }) => (

            <button

              key={label}

              className={`nav-item ${activeTab === label ? "active" : ""}`}

              onClick={() => goTo(label)}

            >

              <Icon size={18} />

              <span>{label}</span>

              {label === "Screening history" && (

                <span className="nav-count">{jobs.length}</span>

              )}

            </button>

          ))}

        </nav>



        <div className="sidebar-bottom">

          <div className="connection-indicator">

            <span

              className={`status-dot ${health?.status !== "ok" ? "offline" : ""}`}

            />

            <div>

              <strong>Local engine</strong>

              <small>

                {health?.engine_exists

                  ? "Executable found"

                  : "Checking connection…"}

              </small>

            </div>

          </div>

          <div className="sidebar-footer">

            dsGuard v0.3.0 · Local development

          </div>

        </div>

      </aside>



      <main className="main-area">

        <header className="topbar">

          <div className="breadcrumb">

            Workspace <span>/</span> {activeTab}

          </div>

          <div className="topbar-right">

            <span className="environment-badge">

              <span className="status-dot" /> LOCAL

            </span>

            <div className="avatar">DS</div>

          </div>

        </header>



        <div className="page-content">

          {error && (

            <div className="notice error-notice">

              <AlertCircle size={18} />

              <p>{error}</p>

              <button className="text-button" onClick={loadDashboard}>

                Retry

              </button>

            </div>

          )}



          {notice && (

            <div className="notice">

              <CircleHelp size={18} />

              <p>{notice}</p>

              <button

                className="text-button"

                onClick={() => setNotice("")}

              >

                Dismiss

              </button>

            </div>

          )}



          {activeTab === "Overview" && (

            <>

              <section className="welcome-row">

                <div>

                  <div className="eyebrow">

                    SEQUENCE INTELLIGENCE WORKSPACE

                  </div>

                  <h2>Screen with confidence.</h2>

                  <p>

                    Design candidates. Screen non-target similarity. Review

                    the evidence.

                  </p>

                </div>

                <button

                  className="primary-button"

                  onClick={() => goTo("New screening")}

                >

                  <FlaskConical size={17} /> New screening{" "}

                  <ArrowUpRight size={16} />

                </button>

              </section>



              <section className="metric-grid">

                <article className="metric-card">

                  <div className="metric-icon violet">

                    <Dna size={20} />

                  </div>

                  <div className="metric-label">TARGET GENES PROCESSED</div>

                  <div className="metric-value">

                    {latestSummary.length || "—"}

                  </div>

                  <div className="metric-foot">

                    <CheckCircle2 size={14} /> Latest completed job

                  </div>

                </article>



                <article className="metric-card">

                  <div className="metric-icon blue">

                    <FileSearch size={20} />

                  </div>

                  <div className="metric-label">REFERENCE PANEL</div>

                  <div className="metric-value">

                    {latestCompleted?.manifest?.panel?.[

                      getPanelLabel(latestCompleted)

                    ]?.Mb ?? "—"}{" "}

                    <small>MB</small>

                  </div>

                  <div className="metric-foot">

                    Latest completed screening

                  </div>

                </article>



                <article className="metric-card">

                  <div className="metric-icon amber">

                    <Activity size={20} />

                  </div>

                  <div className="metric-label">REPORTED EXACT-MATCH HITS</div>

                  <div className="metric-value">

                    {latestSummary.length ? reportedHits : "—"}

                  </div>

                  <div className="metric-foot">Latest completed job</div>

                </article>



                <article className="metric-card">

                  <div className="metric-icon green">

                    <ShieldCheck size={20} />

                  </div>

                  <div className="metric-label">ENGINE STATUS</div>

                  <div className="metric-value status-value">

                    {loading

                      ? "Checking"

                      : health?.status === "ok" && health?.engine_exists

                        ? "Ready"

                        : "Offline"}

                  </div>

                  <div className="metric-foot">

                    {health?.mongodb || "Checking database"}

                  </div>

                </article>

              </section>



              <section className="content-grid">

                <article className="panel recent-panel">

                  <div className="panel-heading">

                    <div>

                      <h3>Recent screening jobs</h3>

                      <p>Live data from MongoDB</p>

                    </div>

                    <button

                      className="text-button"

                      onClick={() => goTo("Screening history")}

                    >

                      View history <ArrowUpRight size={15} />

                    </button>

                  </div>



                  <div className="table-wrap">

                    <table>

                      <thead>

                        <tr>

                          <th>SCREENING JOB</th>

                          <th>TARGET</th>

                          <th>STATUS</th>

                          <th>RESULTS</th>

                        </tr>

                      </thead>

                      <tbody>

                        {jobs.slice(0, 5).map((job) => (

                          <tr key={job.job_id}>

                            <td>

                              <button

                                className="job-link"

                                onClick={() => openJob(job.job_id)}

                              >

                                {jobDisplayName(job)}

                              </button>

                              <div className="job-subtitle">

                                #{job.job_id.slice(0, 8)}

                              </div>

                            </td>

                            <td>{getPanelLabel(job)}</td>

                            <td>

                              <span className={`job-status ${job.status}`}>

                                <span className="status-dot" />

                                {job.status}

                              </span>

                            </td>

                            <td>

                              <button

                                className="text-button"

                                onClick={() => openJob(job.job_id)}

                              >

                                View <ArrowUpRight size={13} />

                              </button>

                            </td>

                          </tr>

                        ))}

                        {!jobs.length && (

                          <tr>

                            <td colSpan="4">

                              {loading ? "Loading jobs…" : "No jobs found yet."}

                            </td>

                          </tr>

                        )}

                      </tbody>

                    </table>

                  </div>

                </article>



                <article className="panel system-panel">

                  <div className="panel-heading">

                    <div>

                      <h3>System health</h3>

                      <p>Live service status</p>

                    </div>

                    <Settings2 size={19} />

                  </div>



                  {[

                    ["FastAPI backend", !error],

                    ["MongoDB Atlas", health?.mongodb === "connected"],

                    ["dsGuard CLI", health?.engine_exists],

                  ].map(([label, ok]) => (

                    <div className="health-row" key={label}>

                      <span>

                        <span

                          className={`status-dot ${ok ? "" : "offline"}`}

                        />

                        {label}

                      </span>

                      <strong>{ok ? "Connected" : "Unavailable"}</strong>

                    </div>

                  ))}



                  <div className="notice">

                    <CircleHelp size={17} />

                    <p>

                      Sequence similarity is a screening signal, not proof of

                      biological safety.

                    </p>

                  </div>

                </article>

              </section>



              <section className="bottom-banner">

                <div className="banner-icon">

                  <Upload size={21} />

                </div>

                <div>

                  <h3>Ready to screen a new target?</h3>

                  <p>

                    Upload a target FASTA and a non-target reference panel to

                    begin.

                  </p>

                </div>

                <button

                  className="secondary-button"

                  onClick={() => goTo("New screening")}

                >

                  Open screening form <ArrowUpRight size={15} />

                </button>

              </section>

            </>

          )}



          {activeTab === "New screening" && (

            <section className="panel form-panel">

              <div className="panel-heading">

                <div>

                  <div className="eyebrow">NEW ANALYSIS</div>

                  <h3>Configure a screening job</h3>

                  <p>

                    Upload FASTA files to run the existing dsGuard engine.

                  </p>

                </div>

                <FlaskConical size={22} />

              </div>



              <form onSubmit={submitScreening} className="screen-form">

                <div className="input-mode-switch" role="group" aria-label="Target input method">

                  <button type="button" className={targetMode === "file" ? "mode-button active" : "mode-button"} onClick={() => setTargetMode("file")}>Upload FASTA</button>

                  <button type="button" className={targetMode === "ncbi" ? "mode-button active" : "mode-button"} onClick={() => setTargetMode("ncbi")}>Gene / NCBI search</button>

                  <button type="button" className={targetMode === "sequence" ? "mode-button active" : "mode-button"} onClick={() => setTargetMode("sequence")}>Paste sequence</button>

                </div>

                {targetMode === "file" && (

                  <label>

                    Target gene FASTA <span>Required</span>

                    <input type="file" accept=".fa,.fasta,.fna,.fas,.txt" required onChange={(e) => setTargetFile(e.target.files?.[0] || null)} />

                    <small>Upload the candidate gene or dsRNA-region nucleotide sequence.</small>

                  </label>

                )}

                {targetMode === "ncbi" && (

                  <>

                    <label>

                      Target organism

                      <select

                        value={targetOrganism}

                        onChange={(event) => {

                          setTargetOrganism(event.target.value);

                          setTargetQuery("");

                        }}

                      >

                        {(targetCatalog.length ? targetCatalog : [

                          { scientific_name: "Nilaparvata lugens", genes: [] },

                          { scientific_name: "Chilo suppressalis", genes: [] },

                          { scientific_name: "Sogatella furcifera", genes: [] },

                        ]).map((organism) => (

                          <option key={organism.scientific_name} value={organism.scientific_name}>

                            {organism.scientific_name}

                          </option>

                        ))}

                      </select>

                    </label>

                    <label>

                      Catalog target gene

                      <select

                        value={selectedCatalogGene?.gene_name || ""}

                        onChange={(event) => {

                          const gene = targetGenes.find(

                            (entry) => entry.gene_name === event.target.value

                          );

                          if (!gene) {

                            setTargetQuery("");

                            return;

                          }

                          setTargetQuery(gene.accession_status === "verified"

                            ? gene.accession

                            : gene.gene_name);

                        }}

                      >

                        <option value="">Manual accession or gene name</option>

                        {targetGenes.map((gene) => (

                          <option key={`${gene.source_workbook}-${gene.source_row}-${gene.gene_name}`} value={gene.gene_name}>

                            {gene.gene_name} — {gene.accession_status === "verified"

                              ? gene.accession

                              : gene.accession_status}

                          </option>

                        ))}

                      </select>

                      <small>

                        Catalog entries select target metadata only; they are not non-target reference panels.

                      </small>

                    </label>

                    <label>

                      NCBI accession ID or gene name <span>Required</span>

                      <input value={targetQuery} onChange={(e) => setTargetQuery(e.target.value)} placeholder="e.g. NlRan or an accession ID" required />

                      <small>An accession, or a project gene alias such as NlRan, resolves to the exact record curated for it; other gene names are searched in NCBI Nucleotide.</small>

                    </label>

                  </>

                )}

                {targetMode === "sequence" && (

                  <label>

                    Paste target FASTA or nucleotide sequence <span>Required</span>

                    <textarea className="sequence-input" value={targetSequence} onChange={(e) => setTargetSequence(e.target.value)} placeholder={">NlRan candidate gene\nATG..."} rows={7} required />

                    <small>Accepts FASTA header + sequence, or raw nucleotide bases.</small>

                  </label>

                )}

                <label>

                  Non-target organism reference panel <span>Required</span>

                  <select value={panelKey} onChange={(e) => setPanelKey(e.target.value)} required>

                    {referencePanels.length ? referencePanels.map((panel) => (

                      <option key={panel.key} value={panel.key} disabled={!panel.available}>

                        {panelOptionLabel(panel.label)}{panel.available ? "" : " — reference FASTA not installed"}

                      </option>

                    )) : <>

                      <option value="rice">Rice (Oryza sativa)</option>

                      <option value="cyrtorhinus_lividipennis">Mirid bug (Cyrtorhinus lividipennis)</option>

                      <option value="honeybee">Honey bee (Apis mellifera)</option>

                      <option value="spider">Spider (reference species panel)</option>

                    </>}

                  </select>

                  <small>Each organism must have a curated FASTA panel installed on the backend.</small>

                </label>



                <div className="parameter-grid">

                  <label>

                    Sequence length (dsRNA window)

                    <select

                      value={windowMode}

                      onChange={(e) => {

                        const choice = e.target.value;

                        setWindowMode(choice);

                        if (choice !== "custom") {

                          setParams({ ...params, window: choice });

                        }

                      }}

                      required

                    >

                      {WINDOW_PRESETS.map((preset) => (

                        <option key={preset.label} value={String(preset.value)}>

                          {preset.label}

                        </option>

                      ))}

                      <option value="custom">Custom length{"\u2026"}</option>

                    </select>

                    {windowMode === "custom" && (

                      <input

                        type="number"

                        min={WINDOW_MIN}

                        max={WINDOW_MAX}

                        step="1"

                        value={params.window}

                        onChange={(e) =>

                          setParams({ ...params, window: e.target.value })

                        }

                        placeholder={`${WINDOW_MIN}-${WINDOW_MAX}`}

                        required

                      />

                    )}

                    <small>

                      {windowMode === "custom"

                        ? `A whole number between ${WINDOW_MIN} and ${WINDOW_MAX} bp.` +

                          (windowBand ? ` ${params.window} bp is inside the ${windowBand.label} band.` : "")

                        : `Sends ${params.window} bp. The presets overlap at their boundaries (100/101, 200/201, 300/301, 400/401); each sends the lower end of its band.`}

                    </small>

                  </label>



                  <label>

                    k-mer length

                    <input

                      type="number"

                      min="2"

                      max="31"

                      value={params.k}

                      onChange={(e) =>

                        setParams({ ...params, k: e.target.value })

                      }

                      required

                    />

                  </label>



                  <label>

                    L50

                    <input

                      type="number"

                      min="0.1"

                      step="any"

                      value={params.l50}

                      onChange={(e) =>

                        setParams({ ...params, l50: e.target.value })

                      }

                      required

                    />

                  </label>



                  <label>

                    Kappa

                    <input

                      type="number"

                      min="0.1"

                      step="any"

                      value={params.kappa}

                      onChange={(e) =>

                        setParams({ ...params, kappa: e.target.value })

                      }

                      required

                    />

                  </label>



                  <label>

                    Maximum hazard

                    <input

                      type="number"

                      min="0"

                      max="1"

                      step="any"

                      value={params.max_hazard}

                      onChange={(e) =>

                        setParams({ ...params, max_hazard: e.target.value })

                      }

                      required

                    />

                  </label>

                </div>



                <div className="form-actions">

                  <button

                    className="primary-button"

                    disabled={submitting}

                    type="submit"

                  >

                    {submitting ? (

                      <LoaderCircle className="spin" size={17} />

                    ) : (

                      <FlaskConical size={17} />

                    )}

                    {submitting ? "Submitting…" : "Start screening"}

                  </button>

                  <button

                    className="secondary-button"

                    type="button"

                    onClick={() => goTo("Overview")}

                  >

                    Cancel

                  </button>

                </div>

              </form>

            </section>

          )}



          {activeTab === "Screening history" && (

            <section className="panel">

              <div className="panel-heading">

                <div>

                  <h3>Screening history</h3>

                  <p>{jobs.length} jobs retrieved from MongoDB</p>

                </div>

                <div className="form-actions">

                  <button className="text-button" onClick={loadDashboard}>

                    <RefreshCw size={15} /> Refresh

                  </button>

                  <button

                    className="delete-button"

                    type="button"

                    disabled={!jobs.length}

                    onClick={async () => {

                      const confirmed = window.confirm(

                        "Clear all screening history?\n\nThis permanently removes all application-scoped screening records, uploads, and generated reports. Shared reference panels, target catalogs, and source datasets are not affected."

                      );

                      if (!confirmed) return;

                      setError("");

                      setNotice("");

                      try {

                        const response = await fetch(`${API}/api/jobs`, { method: "DELETE" });

                        const data = await response.json().catch(() => ({}));

                        if (!response.ok) throw new Error(data.detail || "Could not clear screening history.");

                        if (selectedJob) {

                          resultsRequestRef.current += 1;

                          setSelectedJob(null);

                          setResults(null);

                          setHits([]);

                        }

                        await loadDashboard();

                        setNotice(`Cleared ${data.deleted_count ?? 0} screening history record(s).`);

                      } catch (err) {

                        setError(err.message || "Could not clear screening history.");

                      }

                    }}

                  >

                    Clear History

                  </button>

                </div>

              </div>



              <div className="table-wrap">

                <table>

                  <thead>

                    <tr>

                      <th>JOB ID</th>

                      <th>TARGET FILE</th>

                      <th>PANEL</th>

                      <th>STATUS</th>

                      <th>CREATED</th>

                      <th></th>

                    </tr>

                  </thead>

                  <tbody>

                    {jobs.map((job) => (

                      <tr key={job.job_id}>

                        <td>{job.job_id.slice(0, 8)}</td>

                        <td>

                          {job.target_filename ? jobDisplayName(job) : "—"}

                        </td>

                        <td>{getPanelLabel(job)}</td>

                        <td>

                          <span className={`job-status ${job.status}`}>

                            <span className="status-dot" />

                            {job.status}

                          </span>

                        </td>

                        <td>

                          {job.created_at

                            ? new Date(job.created_at).toLocaleString()

                            : "—"}

                        </td>

                        <td>



<div className="history-actions">

  <button

    className="text-button"

    onClick={() => openJob(job.job_id)}

  >

    Open

  </button>



  <button

    className="delete-button"

    onClick={() => deleteJob(job)}

    disabled={job.status === "queued" || job.status === "running"}

    title={

      job.status === "queued" || job.status === "running"

        ? "Wait until this screening finishes before deleting it."

        : "Delete this screening"

    }

  >

    Delete

  </button>

</div>



                        </td>

                      </tr>

                    ))}

                    {!jobs.length && (

                      <tr>

                        <td colSpan="6">No jobs found.</td>

                      </tr>

                    )}

                  </tbody>

                </table>

              </div>

            </section>

          )}



          {activeTab === "Results" && (

            <section className="panel results-panel">

              <div className="panel-heading">

                <div>

                  <h3>Screening results</h3>

                  <p>

                    {selectedJob ||

                      "Select a job from history to inspect results."}

                  </p>

                </div>



                {results?.summary?.length > 0 && (

                  <div className="form-actions">

                    <button

                      className="secondary-button"

                      type="button"

                      title="One row per target gene (screening status, recommendation, hazard, exact-match hit count). This is a per-gene summary, not the complete ranked off-target results."

                      onClick={() =>

                        downloadResultsCSV({

                          ...results,

                          job_id: selectedJob,

                        })

                      }

                    >

                      <FileText size={16} />

                      Download gene-summary CSV (per-gene)

                    </button>



                    <button

                      className="primary-button"

                      type="button"

                      onClick={() => downloadDetailedHits(selectedJob, setNotice)}

                    >

                      <FileSearch size={16} />

                      Download Detailed Hits (≥{analysis?.parameters?.threshold ?? 21}nt)

                    </button>

                    <button

                      className="secondary-button"

                      type="button"

                      onClick={() =>

                        downloadExport(

                          selectedJob,

                          "report.pdf",

                          `dsguard-full-report-${selectedJob}.pdf`,

                          setNotice

                        )

                      }

                    >

                      <FileText size={16} />

                      Download Full PDF Report

                    </button>

                  </div>

                )}

              </div>



              {!results && (

                <p>

                  {selectedJob

                    ? "Loading results…"

                    : "Choose a screening job from the Overview or Screening history page."}

                </p>

              )}



              {results && (

                <>

                  <div className={`job-status results-status ${results.status}`}>

                    <span className="status-dot" />

                    {results.status}

                  </div>



                  {results.error && (

                    <div className="notice error-notice">

                      <AlertCircle size={17} />

                      <p>{results.error}</p>

                    </div>

                  )}



                  {results.status === "completed" && !analysis && (

                    <div className="notice error-notice">

                      <AlertCircle size={17} />

                      <p>

                        This screening was recorded before dsGuard stored structured

                        results (<code>results.json</code>), so the ranked off-target

                        matches, run parameters, reference coverage and the

                        downloadable PDF/ranked-CSV/JSON reports are unavailable for

                        it. The gene summary and detailed hits shown here are this

                        job&apos;s own recorded output. Re-run the screening from its

                        original target and reference panels to regenerate the full

                        structured results.

                      </p>

                    </div>

                  )}



                  {Array.isArray(results.summary) &&

                    results.summary.some((row) => isThresholdExceeded(row, results)) && (

                      <div className="threshold-warning results-threshold-warning" role="alert">

                        <AlertCircle size={20} />

                        <div>

                          <strong>Informational heuristic hazard comparison</strong>

                          <p>

                          One or more selected candidates exceed the configured maximum hazard reference

                            {getHazardLimit(results) !== null

                              ? ` of ${formatHazard(getHazardLimit(results))}`

                            : ""}. This heuristic comparison does not determine candidate eligibility,

                          is not a pass/fail verdict, and is not proof of biological safety. Review

                          each gene's match evidence, panel coverage, and heuristic limitations.

                          </p>

                        </div>

                      </div>

                    )}



                  {results.summary && (

                    <section className="results-summary">

                      <div className="panel-heading">

                        <div>

                          <h3>Screening Summary and Hazard Evidence</h3>

                          <p>Screening status is separate from the informational hazard comparison.</p>

                        </div>

                        <ShieldCheck size={20} />

                      </div>

                      <div className="table-wrap">

                      <table>

                        <thead>

                          <tr>

                            <th>GENE</th>

                            <th>SCREENING STATUS</th>

                            <th>HAZARD REFERENCE</th>

                            <th>HAZARD SCORE</th>

                            <th>MAX HAZARD</th>

                            <th>LENGTH</th>

                            <th>GC %</th>

                            <th>MAX MATCH</th>

                            <th title={`Confirmed exact matches >= ${analysis?.parameters?.threshold ?? 21} nt (the guaranteed-detection threshold). This is narrower than the reported off-target matches listed further down.`}>

                              HITS (&ge;{analysis?.parameters?.threshold ?? 21} nt)

                            </th>

                          </tr>

                        </thead>

                        <tbody>

                          {results.summary.map((row, index) => (

                            <tr key={`${row.gene}-${index}`}>

                              <td>{geneDisplayName(row.gene)}</td>

                              <td>{screeningStatusLabel(row, analysis)}</td>

                              <td>

                                <span className="recommendation-badge neutral">

                                  {getRecommendationLabel(row, results)}

                                </span>

                              </td>

                              <td>{formatHazard(row.rec_hazard)}</td>

                              <td>{getHazardLimit(results) === null ? "—" : formatHazard(getHazardLimit(results))}</td>

                              <td>{row.length_bp}</td>

                              <td>{row.gc_pct}</td>

                              <td>{row.rec_max_match}</td>

                              <td>{row.total_hits}</td>

                            </tr>

                          ))}

                        </tbody>

                      </table>

                      </div>

                    </section>

                  )}



                  {analysis && (

                    <section className="panel results-section profile-section" style={{ marginTop: "20px" }}>

                      <div className="panel-heading">

                        <div>

                          <h3>Positional Efficacy &amp; Non-Target Match Risk Profile</h3>

                          <p>Engine-calculated positional values for the selected target gene and reference panel.</p>

                        </div>

                        <Activity size={20} />

                      </div>

                      {profile === null ? <p>Loading positional data...</p> : (

                        <PositionalProfile

                          profile={profile}

                          candidate={activeCandidate}

                          mode={profileMode}

                          species={profileSpecies}

                          layers={profileLayers}

                          onSpeciesChange={setProfileSpecies}

                          onModeChange={setProfileMode}

                          onToggleLayer={(key) => setProfileLayers((current) => ({ ...current, [key]: !current[key] }))}

                        />

                      )}

                    </section>

                  )}



                  {analysis && (

                    <section className="panel results-section results-metadata" style={{ marginTop: "20px" }}>

                      <div className="panel-heading">

                        <div>

                          <h3>Target Metadata and Run Parameters</h3>

                          <p>

                            {analysis.analysis_id

                              ? `Analysis ${String(analysis.analysis_id).slice(0, 8)} · ${analysis.created_utc || ""}`

                              : "Current analysis"}

                          </p>

                        </div>

                        <Settings2 size={20} />

                      </div>



                      <div className="coverage-grid">

                        <div>

                          <div className="eyebrow">PARAMETERS</div>

                          <ul className="detail-list">

                            <li>

                              Window:{" "}

                              <strong>{analysis.parameters?.window} bp</strong>

                            </li>

                            <li>

                              Seed k:{" "}

                              <strong>

                                {analysis.parameters?.k} (exact at ≥{" "}

                                {analysis.parameters?.threshold} nt)

                              </strong>

                            </li>

                            <li>

                              L50 / kappa:{" "}

                              <strong>

                                {analysis.parameters?.l50} /{" "}

                                {analysis.parameters?.kappa}

                              </strong>{" "}

                              (placeholder, not fitted)

                            </li>

                            <li>

                              Maximum hazard used:{" "}

                              <strong>{analysis.parameters?.max_hazard}</strong>

                              {analysis.parameters

                                ?.hazard_threshold_equivalent_match_bp != null && (

                                <>

                                  {" "}

                                  (≈{" "}

                                  {

                                    analysis.parameters

                                      .hazard_threshold_equivalent_match_bp

                                  }{" "}

                                  nt equivalent match)

                                </>

                              )}

                            </li>

                            <li>

                              Reporting floor:{" "}

                              <strong>

                                {analysis.reporting?.min_report_match_bp} nt

                              </strong>

                            </li>

                            <li>

                              Input: <strong>{analysis.input?.label}</strong> (

                              {analysis.input?.type})

                            </li>

                          </ul>

                        </div>



                        <div>

                          <div className="eyebrow">REFERENCE PANELS SEARCHED</div>

                          {searchedPanels.length ? (

                            <ul className="detail-list">

                              {searchedPanels.map((panel) => (

                                <li key={panel.name}>

                                  <strong>{panelOptionLabel(panel.name)}</strong> — {panel.Mb} Mb,{" "}

                                  {panel.records} records

                                  {panel.annotation_example ? (

                                    <>

                                      <br />

                                      <span className="detail-note">

                                        {panel.annotation_example}

                                      </span>

                                    </>

                                  ) : null}

                                  {panel.sha256 ? (

                                    <>

                                      <br />

                                      <span className="detail-note">

                                        sha256 {String(panel.sha256).slice(0, 16)}…

                                      </span>

                                    </>

                                  ) : null}

                                </li>

                              ))}

                            </ul>

                          ) : (

                            <p>No panel was searched for this run.</p>

                          )}



                          {skippedPanels.length > 0 && (

                            <>

                              <div

                                className="eyebrow"

                                style={{ marginTop: "12px" }}

                              >

                                PANELS NOT SEARCHED (UNRESOLVED RISK)

                              </div>

                              <ul className="detail-list skipped">

                                {skippedPanels.map((panel) => (

                                  <li key={panel.name}>

                                    <strong>{panelOptionLabel(panel.name)}</strong> — {panel.detail}

                                  </li>

                                ))}

                              </ul>

                            </>

                          )}

                        </div>

                      </div>

                    </section>

                  )}



                  {analysis && (

                    <section className="panel results-section results-hits" style={{ marginTop: "20px" }}>

                      <div className="panel-heading">

                        <div>

                          <h3>Detailed Off-Target Hits</h3>

                          <p>

                            {rankedOfftargets.length

                              ? `${rankedOfftargets.length} reported match(es), highest predicted off-target risk first`

                              : "No matches at or above the reporting floor"}

                          </p>

                          {rankedOfftargets.length > 0 && (

                            <p className="panel-note">

                              These are the <em>reported off-target matches</em>: perfect, gapless

                              matches at or above the reporting floor

                              ({analysis?.reporting?.min_report_match_bp ?? analysis?.parameters?.threshold ?? 21}{" "}

                              nt), plus any shorter match inside the recommended window that sets

                              its hazard score (listed so the recommendation stays traceable).

                              The <em>HITS</em> column in the gene summary is a different, narrower

                              count: only confirmed exact matches at or above the

                              guaranteed-detection threshold

                              (&ge;{analysis?.parameters?.threshold ?? 21} nt), whose enumeration

                              is provably complete. The two counts are therefore expected to

                              differ and are not meant to be equal.

                            </p>

                          )}

                        </div>

                        <ShieldCheck size={20} />

                      </div>



                      {rankedOfftargets.length > 0 ? (

                        <div className="table-wrap table-wrap--wide">

                          <table>

                            <thead>

                              <tr>

                                <th>RANK</th>

                                <th>TARGET ID</th>

                                <th>TARGET WINDOW</th>

                                <th>MATCHED ORGANISM</th>

                                <th>TRANSCRIPT / REFERENCE ID</th>

                                <th>PANEL</th>

                                <th>MATCH (BP)</th>

                                <th>IDENTITY</th>

                                <th>MISMATCH / GAP</th>

                                <th>STRAND</th>

                                <th>TARGET COORDS</th>

                                <th>REF POS</th>

                                <th>POTENCY</th>

                                <th>RISK SCORE (HAZARD)</th>

                                <th>CATEGORY</th>

                                <th>THRESHOLD</th>

                                <th>EVIDENCE</th>

                              </tr>

                            </thead>

                            <tbody>

                              {rankedOfftargets.slice(0, 5).map((row) => (

                                <tr key={`${row.job_rank}-${row.gene}-${row.transcript_pos_1based}`}>

                                  <td>{row.job_rank}</td>

                                  <td>{geneDisplayName(row.gene)}</td>

                                  <td>{row.target_window}</td>

                                  <td className="wrap-cell">

                                    <OrganismCell value={row.matched_organism} />

                                  </td>

                                  <td className="wrap-cell">{row.transcript_or_reference_id}</td>

                                  <td>{organismLabel(row.reference_panel)}</td>

                                  <td>{row.match_length_bp}</td>

                                  <td>

                                    {row.identity_pct != null

                                      ? `${row.identity_pct}%`

                                      : "—"}

                                  </td>

                                  <td>

                                    {row.mismatches}/{row.gaps}

                                  </td>

                                  <td>{row.strand}</td>

                                  <td>

                                    {row.target_start}-{row.target_end}

                                  </td>

                                  <td>{row.transcript_pos_1based}</td>

                                  <td>{row.potency ?? "—"}</td>

                                  <td>{formatHazard(row.hazard)}</td>

                                  <td>

                                    <span

                                      className={`risk-chip ${String(

                                        row.risk_category || ""

                                      ).toLowerCase()}`}

                                    >

                                      {row.risk_category}

                                    </span>

                                  </td>

                                  <td>

                                    {row.threshold_status}

                                  </td>

                                  <td className="evidence-cell">

                                    {row.evidence}

                                  </td>

                                </tr>

                              ))}

                            </tbody>

                          </table>

                        </div>

                      ) : (

                        <p>

                          No off-target matches at or above the reporting floor

                          were found against the searched reference panel(s).

                          Absence of a match here is not proof of non-target

                          safety, and a missing panel is not zero risk.

                        </p>

                      )}



                      {rankedOfftargets.length > 5 && (

                        <button

                          className="secondary-button show-all-hits"

                          type="button"

                          onClick={() => setActiveTab("Detailed hits")}

                        >

                          Show all {rankedOfftargets.length} hits

                        </button>

                      )}



                      {analysis.genes?.some(

                        (gene) => (gene.matches_below_report_floor || 0) > 0

                      ) && (

                        <div className="notice">

                          <CircleHelp size={17} />

                          <p>

                            Perfect matches shorter than the reporting floor were

                            also found and are counted per gene

                            (matches_below_report_floor) but not listed, because

                            at that length a perfect match against a large

                            reference is expected by chance.

                          </p>

                        </div>

                      )}

                    </section>

                  )}



                  {paretoGenes.length > 0 && (

                    <section className="panel results-section results-pareto" style={{ marginTop: "20px" }}>

                      <div className="panel-heading">

                        <div>

                          <h3>Pareto front (potency vs hazard)</h3>

                          <p>

                            Non-dominated candidate windows: maximise potency,

                            minimise off-target hazard.

                          </p>

                        </div>

                        <BarChart3 size={20} />

                      </div>



                      <label className="pareto-select">

                        Target gene

                        <select

                          value={activeParetoGene?.gene || ""}

                          onChange={(event) => setParetoGene(event.target.value)}

                        >

                          {paretoGenes.map((gene) => (

                            <option key={gene.gene} value={gene.gene}>

                              {geneDisplayName(gene.gene)}

                            </option>

                          ))}

                        </select>

                      </label>



                      {activeParetoGene && (

                        <>

                          <ParetoScatter points={activeParetoGene.pareto} />

                          <p className="detail-note">

                            {activeParetoGene.pareto.length} non-dominated

                            {activeParetoGene.pareto.length === 1

                              ? " window"

                              : " windows"}{" "}

                            plotted for {geneDisplayName(activeParetoGene.gene)} —

                            potency on the y-axis, off-target hazard on the

                            x-axis. A single marker means this gene's front has

                            one member, not a missing chart.

                          </p>

                          <div className="table-wrap">

                            <table>

                              <thead>

                                <tr>

                                  <th>START</th>

                                  <th>END</th>

                                  <th>POTENCY</th>

                                  <th>MAX MATCH (NT)</th>

                                  <th>RISK SCORE (HAZARD)</th>

                                  <th>THRESHOLD</th>

                                  <th>MATCH IS EXACT</th>

                                </tr>

                              </thead>

                              <tbody>

                                {activeParetoGene.pareto.map((point) => (

                                  <tr key={`${point.start}-${point.end}`}>

                                    <td>{point.start}</td>

                                    <td>{point.end}</td>

                                    <td>{point.potency}</td>

                                    <td>{point.max_match}</td>

                                    <td>{formatHazard(point.hazard)}</td>

                                    <td>

                                      {point.threshold_status}

                                    </td>

                                    <td>{point.hazard_is_exact ? "yes" : "no"}</td>

                                  </tr>

                                ))}

                              </tbody>

                            </table>

                          </div>

                        </>

                      )}

                    </section>

                  )}



                  {candidateGenes.length > 0 && activeCandidateGene && (

                    <section className="panel results-section results-candidate" style={{ marginTop: "20px" }}>

                      <div className="panel-heading">

                        <div>

                          <h3>dsRNA Candidate Sequence</h3>

                          <p>

                            The window this run selected, copied out of the target

                            gene. Selecting a sequence is not evidence of biological

                            efficacy or safety.

                          </p>

                        </div>

                        <Dna size={20} />

                      </div>



                      <label className="pareto-select">

                        Target gene

                        <select

                          value={activeCandidateGene.gene || ""}

                          onChange={(event) => setCandidateGene(event.target.value)}

                        >

                          {candidateGenes.map((gene) => (

                            <option key={gene.gene} value={gene.gene}>

                              {geneDisplayName(gene.gene)}

                            </option>

                          ))}

                        </select>

                      </label>



                      {activeCandidate.available ? (

                        <>

                          <div className="coverage-grid">

                            <div>

                              <div className="eyebrow">CANDIDATE (dsRNA WINDOW)</div>

                              <ul className="detail-list">

                                <li>

                                  Candidate ID:{" "}

                                  <strong>{activeCandidate.candidateId || "—"}</strong>

                                </li>

                                <li>

                                  Requested length:{" "}

                                  <strong>

                                    {activeCandidate.requestedLengthBp ?? "—"} bp

                                  </strong>

                                </li>

                                <li>

                                  Actual length:{" "}

                                  <strong>

                                    {activeCandidate.actualLengthBp ?? "—"} bp

                                  </strong>{" "}

                                  (measured on the sequence below)

                                </li>

                                <li>

                                  Position:{" "}

                                  <strong>

                                    {activeCandidate.start1Based ?? "—"}–

                                    {activeCandidate.end1Based ?? "—"}

                                  </strong>{" "}

                                  (1-based

                                  {activeCandidate.strand

                                    ? `, ${activeCandidate.strand} strand`

                                    : ""}

                                  )

                                </li>

                                <li>

                                  Ambiguous bases:{" "}

                                  <strong>{activeCandidate.ambiguousBases ?? "—"}</strong>

                                </li>

                              </ul>

                            </div>



                            <div>

                              <div className="eyebrow">

                                TARGET GENE (SOURCE, KEPT SEPARATE)

                              </div>

                              <ul className="detail-list">

                                <li>

                                  Gene: <strong>{activeTarget.displayName}</strong>

                                </li>

                                <li>

                                  Target sequence length:{" "}

                                  <strong>{activeTarget.lengthBp ?? "—"} bp</strong> —

                                  the submitted target, not the dsRNA window

                                </li>

                                {activeTarget.accession ? (

                                  <li>

                                    Accession:{" "}

                                    <strong>{activeTarget.accession}</strong>

                                  </li>

                                ) : null}

                                {activeTarget.inputLabel ? (

                                  <li>

                                    Input: <strong>{activeTarget.inputLabel}</strong>

                                    {activeTarget.inputType

                                      ? ` (${activeTarget.inputType})`

                                      : ""}

                                  </li>

                                ) : null}

                                <li>

                                  <span className="detail-note">

                                    {activeTarget.gene}

                                  </span>

                                </li>

                              </ul>

                            </div>

                          </div>



                          <div className="sequence-block">

                            <div className="sequence-head">

                              <span className="eyebrow">

                                ACTUAL CANDIDATE SEQUENCE ·{" "}

                                {activeCandidate.sequence.length} bp

                              </span>

                              <div className="sequence-actions">

                                <button

                                  type="button"

                                  className="secondary-button"

                                  onClick={() =>

                                    copyCandidateSequence(activeCandidateGene)

                                  }

                                >

                                  <Copy size={15} />

                                  Copy Sequence

                                </button>

                                <button

                                  type="button"

                                  className="secondary-button"

                                  title="FASTA of the candidate dsRNA sequence(s) recorded by this run — the target gene sequence is never exported."

                                  onClick={() =>

                                    downloadCandidateFasta(selectedJob, setNotice)

                                  }

                                >

                                  <Download size={15} />

                                  Download FASTA

                                </button>

                                {copyNoticeText ? (

                                  <span className="copy-status">{copyNoticeText}</span>

                                ) : null}

                              </div>

                            </div>

                            <div className="sequence-text">

                              {activeCandidate.sequence}

                            </div>

                          </div>

                        </>

                      ) : (

                        <div className="notice">

                          <AlertCircle size={17} />

                          <p>{activeCandidate.message}</p>

                        </div>

                      )}

                    </section>

                  )}



                  {candidateGenes.length > 0 && activeCandidateGene && (

                    <section className="panel results-section results-target" style={{ marginTop: "20px" }}>

                      <div className="panel-heading">

                        <div>

                          <h3>Target Information and Original Target Sequence</h3>

                          <p>The submitted target is shown separately from the generated dsRNA candidate.</p>

                        </div>

                        <FileSearch size={20} />

                      </div>

                      <div className="coverage-grid">

                        <ul className="detail-list">

                          <li>Gene: <strong>{activeTarget.displayName}</strong></li>

                          <li>Organism: <strong>{activeTarget.organism || "Unavailable"}</strong></li>

                          <li>Accession/version: <strong>{activeTarget.accession || "Unavailable"}</strong></li>

                          <li>Input source: <strong>{activeTarget.inputType || "Unavailable"}</strong></li>

                          <li>Original target length: <strong>{activeTarget.lengthBp ?? "Unavailable"}{activeTarget.lengthBp != null ? " bp" : ""}</strong></li>

                        </ul>

                        <div className="sequence-block">

                          <div className="sequence-head">

                            <span className="eyebrow">TARGET SEQUENCE (ORIGINAL)</span>

                            <button

                              type="button"

                              className="secondary-button"

                              onClick={() => copyTargetSequence(activeCandidateGene)}

                            >

                              <Copy size={15} /> Copy Original Sequence

                            </button>

                          </div>

                          {activeTarget.sequence ? (

                            <div className="sequence-text target-sequence-text">{activeTarget.sequence}</div>

                          ) : (

                            <p className="detail-note">The original target sequence is unavailable because this historical result predates target-sequence persistence.</p>

                          )}

                          {targetCopyNotice ? <p className="copy-status">{targetCopyNotice}</p> : null}

                        </div>

                      </div>

                      <div className="interactive-viewer-block">

                        <div className="panel-heading">

                          <div>

                            <h3>dsGuard Interactive Sequence Result Viewer</h3>

                            <p>Target coordinates are 1-based and map directly to the stored target sequence.</p>

                          </div>

                          <Dna size={20} />

                        </div>

                        <SequenceViewer target={activeTarget} candidate={activeCandidate} hits={activeGeneHits} />

                      </div>

                    </section>

                  )}



                  {preventive && (

                    <section className="panel results-section results-preventive" style={{ marginTop: "20px" }}>

                      <div className="panel-heading">

                        <div>

                          <h3>Preventive Measures and Risk Mitigation</h3>

                          <p>

                            Derived from the findings of this analysis, not

                            generic advice.

                          </p>

                        </div>

                        <ShieldCheck size={20} />

                      </div>



                      {preventive.findings?.length > 0 ? (

                        <div className="table-wrap table-wrap--text">

                          <table>

                            <thead>

                              <tr>

                                <th>PRIORITY</th>

                                <th>FINDING / RISK</th>

                                <th>WHY IT MATTERS</th>

                                <th>SUPPORTING EVIDENCE</th>

                                <th>RECOMMENDED ACTION</th>

                              </tr>

                            </thead>

                            <tbody>

                              {preventive.findings.map((finding, index) => (

                                <tr key={`${finding.finding}-${index}`}>

                                  <td>

                                    <span

                                      className={`priority-chip ${String(

                                        finding.priority || ""

                                      ).toLowerCase()}`}

                                    >

                                      {finding.priority}

                                    </span>

                                  </td>

                                  <td className="wrap-cell">{finding.finding}</td>

                                  <td className="wrap-cell">{finding.why}</td>

                                  <td className="wrap-cell">{finding.evidence}</td>

                                  <td className="wrap-cell">{finding.action}</td>

                                </tr>

                              ))}

                            </tbody>

                          </table>

                        </div>

                      ) : (

                        <p>No specific findings were generated for this run.</p>

                      )}



                      {[

                        [

                          "before_accepting_a_candidate",

                          "Before accepting a candidate",

                        ],

                        [

                          "reducing_potential_off_target_risk",

                          "Reducing potential off-target risk",

                        ],

                        [

                          "validation_and_quality_assurance",

                          "Validation and quality assurance",

                        ],

                      ].map(([key, title]) => {

                        const items = preventive.categories?.[key] || [];

                        if (!items.length) return null;

                        return (

                          <div key={key} className="measure-block">

                            <div className="eyebrow">{title.toUpperCase()}</div>

                            <ul className="detail-list">

                              {items.map((item, index) => (

                                <li key={`${key}-${index}`}>{item}</li>

                              ))}

                            </ul>

                          </div>

                        );

                      })}

                    </section>

                  )}



                  {analysis && (

                    <section className="panel results-section results-exports" style={{ marginTop: "20px" }}>

                      <div className="panel-heading">

                        <div>

                          <h3>Downloads and Exports</h3>

                          <p>

                            Generated by the backend from this run's analysis

                            object. Threshold used:{" "}

                            {runThreshold === null || runThreshold === undefined

                              ? "—"

                              : runThreshold}

                            . The PDF, ranked off-target CSV and JSON below are

                            the complete results; the gene-summary CSV is a

                            per-gene overview only.

                          </p>

                        </div>

                        <FileText size={20} />

                      </div>



                      <div className="form-actions">

                        <button

                          className="primary-button"

                          type="button"

                          onClick={() =>

                            downloadExport(

                              selectedJob,

                              "report.pdf",

                              `dsguard-report-${selectedJob}.pdf`,

                              setNotice

                            )

                          }

                        >

                          <FileText size={16} />

                          Download PDF report

                        </button>

                        <button

                          className="secondary-button"

                          type="button"

                          title="Complete ranked off-target results: rank, target gene/window/coordinates, matched organism, reference accession/record ID, match length, identity/mismatch/gap, risk score and category, and threshold status."

                          onClick={() =>

                            downloadExport(

                              selectedJob,

                              "results.csv",

                              `dsguard-offtargets-${selectedJob}.csv`,

                              setNotice

                            )

                          }

                        >

                          <FileSearch size={16} />

                          Download ranked off-target CSV (complete)

                        </button>

                        <button

                          className="secondary-button"

                          type="button"

                          onClick={() =>

                            downloadExport(

                              selectedJob,

                              "results.json",

                              `dsguard-results-${selectedJob}.json`,

                              setNotice

                            )

                          }

                        >

                          <FileText size={16} />

                          Download JSON

                        </button>

                      </div>

                    </section>

                  )}



                  {activeTab === "__legacy-detailed-hits-disabled__" && <section

                    className="panel results-section results-confirmed-hits"

                    style={{ marginTop: "20px" }}

                  >

                    <div className="panel-heading">

                      <div>

                        <h3>Detailed screening hits</h3>

                        <p>

                          {hitsLoading

                            ? "Loading detailed records..."

                            : `${hits.length} confirmed exact-match hit record(s) at or above the ${analysis?.parameters?.threshold ?? 21} nt guaranteed-detection threshold`}

                        </p>

                        <p className="panel-note">

                          These are the confirmed exact matches counted by the gene summary&apos;s

                          HITS column. <em>Reference record</em> is the identifier (accession) of

                          the matched record as written in the panel FASTA header; it can be a

                          genomic accession (e.g. NC_…) or a transcript (e.g. XM_/NM_), depending

                          on how the panel was built.

                        </p>

                      </div>

                      <FileSearch size={20} />

                    </div>



                    {hitsLoading ? (

                      <p>Loading hit records...</p>

                    ) : hits.length > 0 ? (

                      <div className="table-wrap">

                        <table>

                          <thead>

                            <tr>

                              <th>GENE</th>

                              <th>ORGANISM</th>

                              <th>REFERENCE RECORD</th>

                              <th>REFERENCE POSITION</th>

                              <th>TARGET START</th>

                              <th>TARGET END</th>

                              <th>MATCH LENGTH</th>

                              <th>STRAND</th>

                            </tr>

                          </thead>

                          <tbody>

                            {hits.map((hit, index) => (

                              <tr

                                key={`${hit.gene}-${hit.transcript}-${index}`}

                              >

                                <td>{geneDisplayName(hit.gene)}</td>

                                <td className="wrap-cell">

                                  <OrganismCell value={hit.organism} />

                                </td>

                                <td>{hit.transcript}</td>

                                <td>{hit.transcript_pos_1based}</td>

                                <td>{hit.gene_start}</td>

                                <td>{hit.gene_end}</td>

                                <td>{hit.length} nt</td>

                                <td>{hit.strand}</td>

                              </tr>

                            ))}

                          </tbody>

                        </table>

                      </div>

                    ) : (

                      <p>

                        No detailed hits are available for this screening job.

                      </p>

                    )}

                  </section>}



                  <section

                    className="panel results-section results-report"

                    style={{ marginTop: "20px" }}

                  >

                    <div className="panel-heading">

                      <div>

                        <h3>Gene report viewer</h3>

                        <p>

                          Inspect the report generated for an individual

                          target gene.

                        </p>

                      </div>

                      <FileText size={20} />

                    </div>



                    <label className="field-block">

                      Select target gene

                      <select

                        className="select-block"

                        value={selectedReportGene}

                        onChange={(event) =>

                          loadGeneReport(selectedJob, event.target.value)

                        }

                        disabled={!geneReports.length || reportLoading}

                      >

                        <option value="">Choose a gene report...</option>

                        {geneReports.map((report) => (

                          <option

                            key={report.gene}

                            value={report.gene}

                          >

                            {report.display_name}

                          </option>

                        ))}

                      </select>

                    </label>



                    {reportLoading ? (

                      <p>Loading gene report...</p>

                    ) : geneReport ? (

                      <pre className="report-view">{geneReport}</pre>

                    ) : (

                      <p>

                        {geneReports.length

                          ? "Choose a gene above to inspect its report."

                          : "No gene reports are available for this job."}

                      </p>

                    )}

                  </section>



                </>

              )}

            </section>

          )}



          {activeTab === "Detailed hits" && (

            <section className="panel results-confirmed-hits">

              <div className="panel-heading">

                <div>

                  <h3>Detailed Off-Target Hits</h3>

                  <p>

                    {hitsLoading

                      ? "Loading detailed records..."

                      : `${hits.length} reported exact-match hit record(s). Showing ${hits.length ? (hitsPage - 1) * hitsPageSize + 1 : 0}–${Math.min(hitsPage * hitsPageSize, hits.length)} of ${hits.length}.`}

                  </p>

                  <p className="panel-note">

                    The reporting floor is {analysis?.reporting?.min_report_match_bp ?? analysis?.parameters?.threshold ?? 21} nt;

                    the guaranteed-detection threshold is {analysis?.reporting?.guaranteed_detection_threshold_bp ?? analysis?.parameters?.threshold ?? 21} nt.

                    These counts are not necessarily the same as all matches found.

                  </p>

                </div>

                <FileSearch size={20} />

              </div>

              {!selectedJob || !results ? (

                <p>Select a completed screening job from Results before opening this view.</p>

              ) : hitsLoading ? (

                <p>Loading hit records...</p>

              ) : hits.length > 0 ? (

                <>

                  <div className="table-wrap table-wrap--wide">

                    <table>

                      <thead><tr>

                        <th>GENE</th><th>ORGANISM</th><th>REFERENCE RECORD</th>

                        <th>REFERENCE POSITION</th><th>TARGET START</th><th>TARGET END</th>

                        <th>MATCH LENGTH</th><th>STRAND</th>

                      </tr></thead>

                      <tbody>{visibleHits.map((hit, index) => (

                        <tr key={`${hit.gene}-${hit.transcript}-${index}`}>

                          <td>{geneDisplayName(hit.gene)}</td>

                          <td className="wrap-cell"><OrganismCell value={hit.organism} /></td>

                          <td>{hit.transcript}</td><td>{hit.transcript_pos_1based}</td>

                          <td>{hit.gene_start}</td><td>{hit.gene_end}</td>

                          <td>{hit.length} nt</td><td>{hit.strand}</td>

                        </tr>

                      ))}</tbody>

                    </table>

                  </div>

                  <div className="pagination">

                    <button className="secondary-button" type="button" disabled={hitsPage <= 1} onClick={() => setHitsPage((page) => page - 1)}>Previous</button>

                    <span>Page {hitsPage} of {hitsPageCount}</span>

                    <button className="secondary-button" type="button" disabled={hitsPage >= hitsPageCount} onClick={() => setHitsPage((page) => page + 1)}>Next</button>

                  </div>

                </>

              ) : <p>No detailed hits are available for this screening job.</p>}

            </section>

          )}

        </div>

      </main>

    </div>

  );

}



export default App;
