
import { useEffect, useState } from "react";
import {
  Activity,
  ArrowUpRight,
  BarChart3,
  CheckCircle2,
  CircleHelp,
  Clock3,
  Dna,
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

const API = "http://127.0.0.1:8000";

function getGeneName(header) {
  if (!header) return "Unknown gene";
  return String(header).split("|")[0].trim();
}

function downloadResultsCSV(results) {
  if (!results?.summary?.length) return;

  const columns = [
    ["Gene", "gene"],
    ["Status", "status"],
    ["Length (bp)", "length_bp"],
    ["GC (%)", "gc_pct"],
    ["Maximum Match (nt)", "rec_max_match"],
    ["Exact-Match Hits", "total_hits"],
  ];

  const escapeCSV = (value) => {
    const text = String(value ?? "");
    return `"${text.replace(/"/g, '""')}"`;
  };

  const csv = [
    columns.map(([label]) => escapeCSV(label)).join(","),
    ...results.summary.map((row) =>
      columns
        .map(([, key]) =>
          escapeCSV(
            key === "gene" ? getGeneName(row[key]) : row[key]
          )
        )
        .join(",")
    ),
  ].join("\r\n");

  const blob = new Blob(["\uFEFF", csv], {
    type: "text/csv;charset=utf-8;",
  });

  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = `dsguard-results-${results.job_id || "screening"}.csv`;
  document.body.appendChild(link);
  link.click();
  link.remove();
  URL.revokeObjectURL(url);
}

function downloadDetailedHits(jobId) {
  if (!jobId) return;

  const link = document.createElement("a");
  link.href = `${API}/api/jobs/${encodeURIComponent(jobId)}/hits/download`;
  link.download = `dsguard-hits-${jobId}.csv`;
  document.body.appendChild(link);
  link.click();
  link.remove();
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
    return validNames.join(", ");
  }

  const manifestPanel = job?.manifest?.panel;

  if (
    manifestPanel &&
    typeof manifestPanel === "object" &&
    Object.keys(manifestPanel).length > 0
  ) {
    return Object.keys(manifestPanel).join(", ");
  }

  return job?.panel_name || job?.panel_filename || "—";
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

  const [geneReports, setGeneReports] = useState([]);
  const [selectedReportGene, setSelectedReportGene] = useState("");
  const [geneReport, setGeneReport] = useState("");
  const [reportLoading, setReportLoading] = useState(false);

  const [submitting, setSubmitting] = useState(false);
  const [notice, setNotice] = useState("");

  const [targetFile, setTargetFile] = useState(null);
  const [panelFile, setPanelFile] = useState(null);
  const [panelName, setPanelName] = useState("rice");

  const [params, setParams] = useState({
    window: 400,
    k: 11,
    l50: 21,
    kappa: 2,
    max_hazard: 0.05,
  });

  async function loadDashboard() {
    setError("");

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
    } catch {
      setError(
        "Cannot reach the backend. Check that FastAPI is running at 127.0.0.1:8000."
      );
    } finally {
      setLoading(false);
    }
  }

  useEffect(() => {
    loadDashboard();
    const timer = setInterval(loadDashboard, 5000);
    return () => clearInterval(timer);
  }, []);

  async function loadGeneReports(jobId) {
    try {
      const response = await fetch(
        `${API}/api/jobs/${encodeURIComponent(jobId)}/reports`
      );
      const data = await response.json();

      if (!response.ok) {
        throw new Error(data.detail || "Could not load gene reports.");
      }

      setGeneReports(Array.isArray(data.reports) ? data.reports : []);
    } catch (err) {
      setGeneReports([]);
      setNotice(err.message);
    }
  }

  async function loadGeneReport(jobId, geneName) {
    if (!jobId || !geneName) return;

    setSelectedReportGene(geneName);
    setGeneReport("");
    setReportLoading(true);
    setNotice("");

    try {
      const response = await fetch(
        `${API}/api/jobs/${encodeURIComponent(jobId)}/reports/${encodeURIComponent(geneName)}`
      );
      const data = await response.json();

      if (!response.ok) {
        throw new Error(data.detail || "Could not load gene report.");
      }

      setGeneReport(data.report || "");
    } catch (err) {
      setNotice(err.message);
    } finally {
      setReportLoading(false);
    }
  }

  async function deleteJob(job) {
    const confirmed = window.confirm(
      `Delete screening "${job.target_filename || job.job_id}"?\n\nThis will permanently remove its history record, uploaded files, and generated results.`
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

      // Clear the results view if this was the selected job.
      if (selectedJob === job.job_id) {
        setSelectedJob(null);
        setResults(null);
        setHits([]);
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

      if (!response.ok) {
        throw new Error(data.detail || "Could not load results.");
      }

      setResults(data);

      if (data.status === "completed") {
        const [hitsResult, reportsResult] = await Promise.allSettled([
          fetch(`${API}/api/jobs/${encodeURIComponent(jobId)}/hits`),
          fetch(`${API}/api/jobs/${encodeURIComponent(jobId)}/reports`),
        ]);

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
      setNotice(err.message);
    } finally {
      setHitsLoading(false);
    }
  }

  async function submitScreening(event) {
    event.preventDefault();
    setNotice("");

    if (!targetFile || !panelFile) {
      setNotice("Please select both the target FASTA and panel FASTA files.");
      return;
    }

    const form = new FormData();
    form.append("target_file", targetFile);
    form.append("panel_file", panelFile);
    form.append("panel_name", panelName);

    Object.entries(params).forEach(([key, value]) => {
      form.append(key, String(value));
    });

    setSubmitting(true);

    try {
      const response = await fetch(`${API}/api/screen`, {
        method: "POST",
        body: form,
      });

      const data = await response.json();

      if (!response.ok) {
        throw new Error(data.detail || "Screening submission failed.");
      }

      setNotice(`Screening submitted. Job ID: ${data.job_id}`);
      setTargetFile(null);
      setPanelFile(null);
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
  ];

  const latestCompleted = jobs.find((job) => job.status === "completed");
  const latestSummary = latestCompleted?.summary || [];
  const reportedHits = latestSummary.reduce(
    (total, row) => total + (Number(row.total_hits) || 0),
    0
  );

  function goTo(tab) {
    setActiveTab(tab);
    setNotice("");

    if (tab !== "Results") {
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
                                {job.target_filename || "Screening job"}
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
                <label>
                  Target FASTA <span>Required</span>
                  <input
                    type="file"
                    accept=".fa,.fasta,.fna,.fas,.txt"
                    required
                    onChange={(e) =>
                      setTargetFile(e.target.files?.[0] || null)
                    }
                  />
                </label>

                <label>
                  Non-target panel FASTA <span>Required</span>
                  <input
                    type="file"
                    accept=".fa,.fasta,.fna,.fas,.txt"
                    required
                    onChange={(e) =>
                      setPanelFile(e.target.files?.[0] || null)
                    }
                  />
                </label>

                <label>
                  Panel name
                  <input
                    value={panelName}
                    onChange={(e) => setPanelName(e.target.value)}
                    required
                    maxLength={40}
                  />
                </label>

                <div className="parameter-grid">
                  <label>
                    Window (bp)
                    <input
                      type="number"
                      min="1"
                      value={params.window}
                      onChange={(e) =>
                        setParams({ ...params, window: e.target.value })
                      }
                      required
                    />
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
                <button className="text-button" onClick={loadDashboard}>
                  <RefreshCw size={15} /> Refresh
                </button>
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
                        <td>{job.target_filename || "—"}</td>
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
            <section className="panel">
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
                      onClick={() =>
                        downloadResultsCSV({
                          ...results,
                          job_id: selectedJob,
                        })
                      }
                    >
                      <FileText size={16} />
                      Export Summary CSV
                    </button>

                    <button
                      className="primary-button"
                      type="button"
                      onClick={() => downloadDetailedHits(selectedJob)}
                    >
                      <FileSearch size={16} />
                      Download Detailed Hits
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
                  <div className={`job-status ${results.status}`}>
                    <span className="status-dot" />
                    {results.status}
                  </div>

                  {results.error && (
                    <div className="notice error-notice">
                      <AlertCircle size={17} />
                      <p>{results.error}</p>
                    </div>
                  )}

                  {results.summary && (
                    <div className="table-wrap">
                      <table>
                        <thead>
                          <tr>
                            <th>GENE</th>
                            <th>STATUS</th>
                            <th>LENGTH</th>
                            <th>GC %</th>
                            <th>MAX MATCH</th>
                            <th>HITS</th>
                          </tr>
                        </thead>
                        <tbody>
                          {results.summary.map((row, index) => (
                            <tr key={`${row.gene}-${index}`}>
                              <td>{getGeneName(row.gene)}</td>
                              <td>{row.status}</td>
                              <td>{row.length_bp}</td>
                              <td>{row.gc_pct}</td>
                              <td>{row.rec_max_match}</td>
                              <td>{row.total_hits}</td>
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </div>
                  )}

                  <section
                    className="panel"
                    style={{ marginTop: "20px" }}
                  >
                    <div className="panel-heading">
                      <div>
                        <h3>Detailed screening hits</h3>
                        <p>
                          {hitsLoading
                            ? "Loading detailed records..."
                            : `${hits.length} exact-match hit records`}
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
                              <th>REFERENCE TRANSCRIPT</th>
                              <th>TRANSCRIPT POSITION</th>
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
                                <td>{hit.gene}</td>
                                <td>{hit.organism}</td>
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
                  </section>

                  <section
                    className="panel"
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

                    <label
                      style={{
                        display: "block",
                        marginBottom: "16px",
                        fontWeight: 600,
                      }}
                    >
                      Select target gene
                      <select
                        value={selectedReportGene}
                        onChange={(event) =>
                          loadGeneReport(selectedJob, event.target.value)
                        }
                        disabled={!geneReports.length || reportLoading}
                        style={{
                          display: "block",
                          width: "100%",
                          marginTop: "8px",
                          padding: "12px",
                          borderRadius: "8px",
                          border: "1px solid var(--border-color, #d1d5db)",
                          background: "var(--panel-bg, #ffffff)",
                          color: "inherit",
                        }}
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
                      <pre
                        style={{
                          whiteSpace: "pre-wrap",
                          overflowWrap: "anywhere",
                          fontFamily: "inherit",
                          fontSize: "14px",
                          lineHeight: 1.7,
                          padding: "16px",
                          borderRadius: "10px",
                          background: "var(--report-bg, #f7f8fa)",
                          color: "inherit",
                          maxHeight: "600px",
                          overflowY: "auto",
                        }}
                      >
                        {geneReport}
                      </pre>
                    ) : (
                      <p>
                        {geneReports.length
                          ? "Choose a gene above to inspect its report."
                          : "No gene reports are available for this job."}
                      </p>
                    )}
                  </section>

                  {results.manifest?.limits && (
                    <div className="notice">
                      <CircleHelp size={17} />
                      <p>{results.manifest.limits}</p>
                    </div>
                  )}
                </>
              )}
            </section>
          )}
        </div>
      </main>
    </div>
  );
}

export default App;
