/*
 * app.js - Virtual Knockdown Explorer (static site, no server).
 *
 * Data files (made by preprocess.py):
 *   data/index.json              gene names + per-condition perturbation stats
 *   data/<Condition>/<L>.json    measured targets for perturbed genes starting with letter L
 *   data/model_check.json        model check across many genes (scripts/model_check.js)
 */
(function () {
  "use strict";

  const CONDS = ["Rest", "Stim8hr", "Stim48hr"];
  const COND_LABEL = { Rest: "Rest", Stim8hr: "Stim 8hr", Stim48hr: "Stim 48hr" };
  const REASONS = ["a neighboring gene was also knocked down", "possible off-target guide binding", "the gene's own expression is too low to measure knockdown reliably"];
  // Demo genes (shown only if present in the data)
  const EXAMPLES = [
    { g: "GATA3", why: "T helper 2 master regulator: hundreds of targets in every condition" },
    { g: "CD3E", why: "T cell receptor part: 2 targets at rest, ~5,000 at Stim 8hr" },
    { g: "NFAT5", why: "The reverse: ~1,100 targets at rest, none at Stim 8hr" },
    { g: "IRF4", why: "Key transcription factor after activation" },
    { g: "ZAP70", why: "T cell receptor signaling kinase: huge effects only after activation" },
    { g: "STAT5B", why: "Relays the IL-2 signal" },
    { g: "IL2RA", why: "IL-2 receptor (CD25): strong knockdown, few targets" },
    { g: "FOXP3", why: "Reduced to ~5% of normal, yet only 2-3 genes change" },
  ];

  const $ = (id) => document.getElementById(id);
  const state = { gene: null, cond: "Stim8hr", mode: "measured", topN: 25, depth: 2, cutoff: 0.1, cond2: "Rest", up: false };
  const DB = { genes: [], byName: new Map(), pert: {}, removed: {}, shards: {}, check: null };
  let renderToken = 0;
  const chainCache = new Map();
  let cy1 = null, cy2 = null;
  let current = { rows: [], sortKey: null, sortDir: -1, columns: [] };

  // =====================================================================
  // Data loading
  // =====================================================================
  async function fetchJSON(url) {
    const r = await fetch(url);
    if (!r.ok) throw new Error(`${url}: HTTP ${r.status}`);
    return r.json();
  }

  async function loadIndex() {
    const ix = await fetchJSON("data/index.json");
    DB.genes = ix.genes;
    DB.params = ix.params;
    ix.genes.forEach((n, i) => DB.byName.set(n.toUpperCase(), i));
    const F = ix.pert_fields;
    const fi = Object.fromEntries(F.map((k, i) => [k, i]));
    for (const c of CONDS) {
      const m = new Map();
      for (const a of ix.pert[c] || []) {
        m.set(a[fi.gene], {
          nSig: a[fi.n_sig], nUp: a[fi.n_up], kd: a[fi.kd_lfc], kdz: a[fi.kd_z],
          kdSig: !!a[fi.kd_sig], nCells: a[fi.n_cells], shard: a[fi.shard],
        });
      }
      DB.pert[c] = m;
      DB.removed[c] = new Map((ix.removed[c] || []).map(([g, r]) => [g, r]));
      DB.shards[c] = {};
    }
    // autocomplete candidates: every gene that was perturbed in some condition
    const set = new Set();
    for (const c of CONDS) { for (const g of DB.pert[c].keys()) set.add(g); for (const g of DB.removed[c].keys()) set.add(g); }
    DB.perturbed = [...set].sort((a, b) => DB.genes[a].localeCompare(DB.genes[b]));
    fetchJSON("data/model_check.json").then((d) => { DB.check = d; if (state.mode === "chain") renderGlobalCheck(); }).catch(() => {});
  }

  function loadShard(cond, key) {
    if (!DB.shards[cond][key]) {
      DB.shards[cond][key] = fetchJSON(`data/${cond}/${key}.json`).then((d) => {
        // convert keys to numbers once
        const m = new Map();
        for (const k in d) m.set(+k, d[k]);
        DB.shards[cond][key] = m;
        return m;
      }, (err) => {
        delete DB.shards[cond][key]; // forget the failure so the next try downloads again
        throw err;
      });
    }
    return DB.shards[cond][key];
  }

  /** Make sure measured targets of these genes are loaded (for a condition). */
  async function ensure(cond, gis) {
    const keys = new Set();
    for (const g of gis) {
      const p = DB.pert[cond].get(g);
      if (p && p.shard && !(DB.shards[cond][p.shard] instanceof Map)) keys.add(p.shard);
    }
    await Promise.all([...keys].map((k) => loadShard(cond, k)));
  }

  function getTargets(cond, gi) {
    const p = DB.pert[cond].get(gi);
    if (!p || !p.shard) return p ? [] : null;
    const s = DB.shards[cond][p.shard];
    return s instanceof Map ? s.get(gi) || [] : null;
  }
  const getKd = (cond, gi) => { const p = DB.pert[cond].get(gi); return p ? p.kd : null; };
  const name = (gi) => DB.genes[gi];

  // =====================================================================
  // Helpers
  // =====================================================================
  const fmt = (x, d = 2) => (x == null || !isFinite(x) ? "–" : (x > 0 ? "+" : "") + x.toFixed(d));
  const fmtP = (p) => (p == null ? "–" : p === 0 ? "< 1e-300" : p < 1e-3 ? p.toExponential(1) : p.toPrecision(2));
  const pctLeft = (lfc) => Math.round(100 * Math.pow(2, lfc));
  /** Plain-language level: "23% of normal" for decreases, "3.5× normal" for big increases. */
  const level = (lfc) => {
    const f = Math.pow(2, lfc);
    if (f < 2) return `${Math.round(100 * f)}% of normal`;
    return `${f < 10 ? f.toFixed(1) : Math.round(f).toLocaleString()}× normal`;
  };
  const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const reduceMotion = () => window.matchMedia && matchMedia("(prefers-reduced-motion: reduce)").matches;
  const debounce = (fn, ms) => { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); }; };

  function dirTag(v) {
    if (!v) return '<span class="tag neutral">–</span>';
    return v < 0 ? '<span class="tag down">down</span>' : '<span class="tag up">up</span>';
  }

  // =====================================================================
  // URL hash (so a specific view can be shared / bookmarked)
  // =====================================================================
  let lastMajor = null;
  function writeHash() {
    if (state.gene == null) return;
    const h = new URLSearchParams({ g: name(state.gene), c: state.cond, m: state.mode, n: state.topN, d: state.depth, k: state.cutoff, c2: state.cond2, u: state.up ? 1 : 0 });
    // a new gene, condition or view gets its own history entry, so the browser's Back button
    // returns to it; slider tweaks just update the current entry
    const major = `${state.gene}|${state.cond}|${state.mode}|${state.cond2}`;
    const url = "#" + h.toString();
    if (lastMajor !== null && major !== lastMajor && url !== location.hash) history.pushState(null, "", url);
    else history.replaceState(null, "", url);
    lastMajor = major;
  }
  function readHash() {
    const h = new URLSearchParams(location.hash.slice(1));
    if (!h.get("g")) return null;
    // any setting missing from the link falls back to its default
    state.cond = CONDS.includes(h.get("c")) ? h.get("c") : "Stim8hr";
    state.cond2 = CONDS.includes(h.get("c2")) && h.get("c2") !== state.cond ? h.get("c2") : CONDS.find((x) => x !== state.cond);
    state.mode = ["measured", "chain", "compare"].includes(h.get("m")) ? h.get("m") : "measured";
    state.topN = Math.min(200, Math.max(5, Math.round((+h.get("n") || 25) / 5) * 5));
    state.depth = Math.min(3, Math.max(1, +h.get("d") || 2));
    state.cutoff = Math.min(0.5, Math.max(0.02, +h.get("k") || 0.1));
    state.up = h.get("u") === "1";
    const g = h.get("g") && DB.byName.get(h.get("g").toUpperCase());
    return g != null ? g : null;
  }

  // =====================================================================
  // Controls
  // =====================================================================
  function syncControls() {
    document.body.dataset.mode = state.mode;
    for (const b of $("condSeg").children) b.setAttribute("aria-pressed", b.dataset.cond === state.cond);
    for (const b of $("modeSeg").children) b.setAttribute("aria-pressed", b.dataset.mode === state.mode);
    $("topN").value = state.topN; $("topNVal").textContent = state.topN;
    $("depth").value = state.depth; $("depthVal").textContent = state.depth;
    $("cutoff").value = state.cutoff; $("cutoffVal").textContent = (+state.cutoff).toFixed(2);
    $("cond2").value = state.cond2;
    $("propUp").checked = state.up;
    for (const o of $("cond2").options) o.disabled = o.value === state.cond;
    if (state.gene != null && document.activeElement !== $("geneInput")) $("geneInput").value = name(state.gene);
  }

  function setupControls() {
    $("condSeg").addEventListener("click", (e) => {
      const c = e.target.closest("button")?.dataset.cond; if (!c) return;
      state.cond = c;
      if (state.cond2 === c) state.cond2 = CONDS.find((x) => x !== c);
      update();
    });
    $("modeSeg").addEventListener("click", (e) => {
      const m = e.target.closest("button")?.dataset.mode; if (!m) return;
      state.mode = m; update();
    });
    const rerender = debounce(update, 120);
    $("topN").addEventListener("input", (e) => { state.topN = +e.target.value; $("topNVal").textContent = state.topN; rerender(); });
    $("depth").addEventListener("input", (e) => { state.depth = +e.target.value; $("depthVal").textContent = state.depth; rerender(); });
    $("cutoff").addEventListener("input", (e) => { state.cutoff = +e.target.value; $("cutoffVal").textContent = state.cutoff.toFixed(2); rerender(); });
    $("cond2").addEventListener("change", (e) => { state.cond2 = e.target.value; update(); });
    $("propUp").addEventListener("change", (e) => { state.up = e.target.checked; update(); });
    $("fitBtn").addEventListener("click", () => { refit(cy1); refit(cy2); });
    $("themeBtn").addEventListener("click", () => {
      const dark = getComputedStyle(document.documentElement).colorScheme === "dark";
      document.documentElement.dataset.theme = dark ? "light" : "dark";
      try { localStorage.setItem("kd-theme", document.documentElement.dataset.theme); } catch (e) {}
      // recolor in place so the selected gene and zoom are kept
      for (const cy of [cy1, cy2]) if (cy) cy.style(cyStyle());
      if (state.mode === "chain" && lastCheck) renderCheck(lastCheck.chk, lastCheck.steps); else renderGlobalCheck();
    });
    window.addEventListener("resize", debounce(() => { for (const cy of [cy1, cy2]) if (cy) { cy.resize(); refit(cy); } }, 150));
    setupAutocomplete();
  }

  function setupAutocomplete() {
    const input = $("geneInput"), list = $("acList");
    let items = [], active = -1;
    const close = () => { list.hidden = true; input.setAttribute("aria-expanded", "false"); input.removeAttribute("aria-activedescendant"); active = -1; };
    const hint = (gi) => {
      const p = DB.pert[state.cond].get(gi);
      if (p) return `${p.nSig} targets`;
      if (DB.removed[state.cond].has(gi)) return "filtered here";
      return "other conditions";
    };
    const show = () => {
      const q = input.value.trim().toUpperCase();
      if (!q) { close(); return; }
      const pre = [], sub = [];
      for (const gi of DB.perturbed) {
        const n = DB.genes[gi].toUpperCase();
        if (n.startsWith(q)) pre.push(gi); else if (n.includes(q)) sub.push(gi);
        if (pre.length >= 12) break;
      }
      items = pre.concat(sub).slice(0, 12);
      if (!items.length) { list.innerHTML = '<li class="hint">No perturbed gene matches</li>'; list.hidden = false; return; }
      list.innerHTML = items.map((gi, i) => `<li role="option" id="ac-opt-${i}" data-i="${i}" aria-selected="${i === active}"><span>${esc(DB.genes[gi])}</span><span class="hint">${hint(gi)}</span></li>`).join("");
      list.hidden = false; input.setAttribute("aria-expanded", "true");
      if (active >= 0) input.setAttribute("aria-activedescendant", `ac-opt-${active}`); else input.removeAttribute("aria-activedescendant");
    };
    const pick = (gi) => { close(); input.value = DB.genes[gi]; input.blur(); selectGene(gi); };
    input.addEventListener("input", () => { active = -1; show(); });
    input.addEventListener("focus", () => input.select());
    input.addEventListener("keydown", (e) => {
      if (e.key === "ArrowDown" || e.key === "ArrowUp") {
        e.preventDefault(); if (list.hidden) show();
        if (!items.length) return;
        active = (active + (e.key === "ArrowDown" ? 1 : -1) + items.length) % items.length; show();
      } else if (e.key === "Enter") {
        e.preventDefault();
        if (active >= 0 && items[active] != null) pick(items[active]);
        else {
          const gi = DB.byName.get(input.value.trim().toUpperCase());
          if (gi != null) pick(gi); else if (items.length) pick(items[0]);
        }
      } else if (e.key === "Escape") { close(); }
    });
    list.addEventListener("mousedown", (e) => {
      const li = e.target.closest("li[data-i]"); if (!li) return;
      e.preventDefault(); pick(items[+li.dataset.i]);
    });
    input.addEventListener("blur", () => setTimeout(close, 100));
  }

  function renderExamples() {
    const box = $("exampleBtns");
    box.innerHTML = "";
    for (const ex of EXAMPLES) {
      const gi = DB.byName.get(ex.g);
      if (gi == null) continue;
      const any = CONDS.some((c) => DB.pert[c].has(gi));
      if (!any) continue;
      const b = document.createElement("button");
      b.className = "chip"; b.title = ex.why; b.textContent = ex.g;
      b.addEventListener("click", () => selectGene(gi));
      box.appendChild(b);
    }
  }

  function selectGene(gi) {
    state.gene = gi;
    // if the gene wasn't kept in the current condition, switch to the condition with most targets
    if (!DB.pert[state.cond].has(gi)) {
      const best = CONDS.filter((c) => DB.pert[c].has(gi)).sort((a, b) => DB.pert[b].get(gi).nSig - DB.pert[a].get(gi).nSig)[0];
      if (best) { state.cond = best; if (state.cond2 === best) state.cond2 = CONDS.find((x) => x !== best); }
    }
    update();
  }

  let lastAnnounce = "";
  function announce() {
    if (state.gene == null) return;
    const msg = `${name(state.gene)} knocked down in ${COND_LABEL[state.cond]}, ${{ measured: "measured effects", chain: "chain reaction", compare: "compare with " + COND_LABEL[state.cond2] }[state.mode]}`;
    if (msg !== lastAnnounce) { $("srStatus").textContent = msg; lastAnnounce = msg; }
  }

  // =====================================================================
  // Summary card
  // =====================================================================
  function renderSummary() {
    const el = $("summary"), gi = state.gene, c = state.cond;
    if (gi == null) { el.innerHTML = '<p class="muted">Pick a gene above to begin.</p>'; return; }
    const p = DB.pert[c].get(gi);
    const conds = CONDS.map((cc) => {
      const q = DB.pert[cc].get(gi);
      const label = q ? `${q.nSig}` : DB.removed[cc].has(gi) ? "filtered" : "not tested";
      return `<button class="chip" data-c="${cc}" ${cc === c ? 'style="border-color:var(--accent)"' : ""}>${COND_LABEL[cc]} <small>${label}</small></button>`;
    }).join("");
    let body;
    if (p) {
      const kdTxt = p.kd != null
        ? `<span class="v">${fmt(p.kd)} log₂FC</span><span class="k">reduced to ~${pctLeft(p.kd)}% of normal</span><div class="kdbar"><i style="width:${Math.max(0, Math.min(100, 100 - pctLeft(p.kd)))}%"></i></div>`
        : `<span class="v">z = ${fmt(p.kdz, 1)}</span><span class="k">gene not in measured panel; z-score shown</span>`;
      body = `
        <div class="stat"><span class="k">How much was it reduced?</span>${kdTxt}
          ${p.kdSig ? "" : '<span class="warn"><svg class="ico" viewBox="0 0 16 16" aria-hidden="true"><path fill="currentColor" d="M8 1.5 15 14H1L8 1.5zm-.8 4.5v4h1.6V6H7.2zm0 5.2v1.6h1.6v-1.6H7.2z"/></svg>knockdown not statistically significant</span>'}</div>
        <div class="stat"><span class="v">${p.nSig.toLocaleString()}</span><span class="k">significant targets (adj. p &lt; ${DB.params.padj_threshold})</span></div>
        <div class="stat"><span class="v"><span style="color:var(--down)">${(p.nSig - p.nUp).toLocaleString()}↓</span> · <span style="color:var(--up)">${p.nUp.toLocaleString()}↑</span></span><span class="k">down · up</span></div>
        <div class="stat"><span class="v">${p.nCells ? p.nCells.toLocaleString() : "–"}</span><span class="k">cells with this knockdown</span></div>`;
    } else if (DB.removed[c].has(gi)) {
      const r = DB.removed[c].get(gi);
      const why = REASONS.filter((_, i) => r & (1 << i)).join("; ");
      body = `<p class="warn">Removed in ${COND_LABEL[c]} by quality control: ${why}.</p>`;
    } else {
      body = `<p class="muted">${esc(name(gi))} was not tested in ${COND_LABEL[c]}.</p>`;
    }
    el.innerHTML = `<div><div class="gname">${esc(name(gi))}</div><div class="muted small">knocked down in ${COND_LABEL[c]}</div></div>${body}
      <div class="stat"><span class="k">Targets per condition</span><div class="chips">${conds}</div></div>`;
    el.querySelectorAll("button[data-c]").forEach((b) => b.addEventListener("click", () => {
      state.cond = b.dataset.c; if (state.cond2 === state.cond) state.cond2 = CONDS.find((x) => x !== state.cond); update();
    }));
  }

  // =====================================================================
  // Network drawing (Cytoscape)
  // =====================================================================
  // circle size: square-root scale relative to the strongest effect drawn (so differences stay visible)
  let sizeMax = 2;
  const nodeSize = (v) => 14 + 32 * Math.sqrt(Math.min(Math.abs(v) / sizeMax, 1));
  const setSizeMax = (vals) => { sizeMax = Math.max(0.25, ...vals.map(Math.abs)); };

  function cyStyle() {
    const down = css("--down"), up = css("--up"), src = css("--src"), text = css("--text"), edge = css("--edge"), only = css("--only"), bg = css("--bg");
    return [
      { selector: "node", style: {
        "background-color": (n) => (n.data("v") < 0 ? down : up),
        width: "data(size)", height: "data(size)", label: "data(label)",
        "font-size": 12, color: text, "text-valign": "bottom", "text-margin-y": 2,
        "text-outline-color": bg, "text-outline-width": 2, "min-zoomed-font-size": 7,
        "border-width": 0,
      } },
      { selector: "node[?only]", style: { "border-width": 4, "border-color": only } },
      { selector: "node[?src]", style: { "background-color": src, "font-size": 14, "font-weight": 700, "text-valign": "center", "text-margin-y": 0, color: text } },
      { selector: "edge", style: {
        width: "data(w)", "line-color": edge, "target-arrow-color": edge, "target-arrow-shape": "triangle",
        "arrow-scale": 0.7, "curve-style": "straight", opacity: 0.55,
      } },
      { selector: "edge[?pred]", style: { "line-style": "dashed", "line-dash-pattern": [5, 4] } },
      { selector: ".faded", style: { opacity: 0.12 } },
      { selector: "node:selected", style: { "border-width": 3, "border-color": text } },
      { selector: "edge.hl", style: { opacity: 1, "line-color": text, "target-arrow-color": text } },
    ];
  }

  function makeCy(container, elements) {
    const cy = cytoscape({
      container, elements, style: cyStyle(), layout: { name: "preset" },
      minZoom: 0.15, maxZoom: 4,
      textureOnViewport: true, hideEdgesOnViewport: elements.length > 600, pixelRatio: "auto",
      boxSelectionEnabled: false,
    });
    refit(cy); // fit, but don't blow up tiny networks
    cy.on("tap", "node", (e) => onNodeTap(cy, e.target));
    cy.on("mouseover", "node", (e) => showTip(e.target, e.originalEvent));
    cy.on("mousemove", "node", (e) => moveTip(e.originalEvent));
    cy.on("mouseout", "node", hideTip);
    cy.on("viewport", hideTip);
    cy.on("tap", (e) => { if (e.target === cy) { cy.elements().removeClass("faded hl"); } });
    return cy;
  }

  // hover tooltip: gene name + effect (labels get too small to read in big networks)
  const tipEl = document.createElement("div");
  tipEl.className = "node-tip"; tipEl.hidden = true; document.body.appendChild(tipEl);
  function showTip(node, ev) {
    const v = node.data("v");
    tipEl.innerHTML = node.data("src") ? `<b>${esc(node.data("label"))}</b> · knocked down`
      : `<b>${esc(node.data("label"))}</b> · ${v < 0 ? "down" : "up"} ${fmt(v)} log₂FC (~${level(v)})${node.data("only") ? "<br>affected in only this condition" : ""}`;
    tipEl.hidden = false; moveTip(ev);
  }
  function moveTip(ev) { if (!ev || tipEl.hidden) return; tipEl.style.left = ev.clientX + 14 + "px"; tipEl.style.top = ev.clientY + 14 + "px"; }
  function hideTip() { tipEl.hidden = true; }

  function onNodeTap(cy, node) {
    cy.elements().addClass("faded");
    const nb = node.closedNeighborhood();
    nb.removeClass("faded"); nb.edges().addClass("hl");
    cy.$("node[?src]").removeClass("faded");
    showDetails(+node.id(), cy === cy2 ? state.cond2 : state.cond);
  }

  /** Place groups of node ids on concentric rings (group 0 = center). */
  function ringPositions(groups, maxSize) {
    const pos = {};
    let radius = 0;
    groups.forEach((ids, gi) => {
      if (!ids.length) return;
      if (gi === 0 && ids.length === 1) { pos[ids[0]] = { x: 0, y: 0 }; return; }
      let i = 0;
      while (i < ids.length) {
        radius += Math.max(70, maxSize * 1.6);
        const cap = Math.max(6, Math.floor((2 * Math.PI * radius) / (maxSize * 1.5)));
        const chunk = ids.slice(i, i + cap);
        const off = (gi % 2) * 0.3;
        chunk.forEach((id, j) => {
          const a = off + (2 * Math.PI * j) / chunk.length - Math.PI / 2;
          pos[id] = { x: radius * Math.cos(a), y: radius * Math.sin(a) };
        });
        i += cap;
      }
      radius += maxSize * 0.6; // extra gap between different steps
    });
    return pos;
  }

  function showGraphMsg(msg) {
    const el = $("graphMsg");
    if (msg) { el.innerHTML = msg; el.hidden = false; } else el.hidden = true;
  }

  function refit(cy) { if (!cy) return; cy.fit(undefined, 20); if (cy.zoom() > 1.3) { cy.zoom(1.3); cy.center(); } }

  function destroyCys() {
    hideTip();
    if (cy1) { cy1.destroy(); cy1 = null; }
    if (cy2) { cy2.destroy(); cy2 = null; }
  }

  // ---------- Measured view ----------
  function buildMeasured(cond, gi, topN, positions, otherSet) {
    const targets = (getTargets(cond, gi) || []).slice(0, topN);
    if (!otherSet) setSizeMax(targets.map((t) => t[1]));
    const els = [{ data: { id: String(gi), label: name(gi), src: true, v: -1, size: 48 }, position: { x: 0, y: 0 } }];
    for (const [t, lfc] of targets) {
      const only = otherSet ? !otherSet.has(t) : false;
      els.push({ data: { id: String(t), label: name(t), v: lfc, size: nodeSize(lfc), only }, position: positions[t] });
      els.push({ data: { id: `e${gi}-${t}`, source: String(gi), target: String(t), w: 1 + 2.5 * Math.min(Math.abs(lfc) / 2, 1) } });
    }
    return { els, targets };
  }

  function measuredPositions(ids) {
    // strongest first; rings of increasing radius
    return ringPositions([[], ids], 40);
  }

  /** Show "Loading..." only if the download takes more than a moment. */
  async function withLoading(token, promise) {
    const t = setTimeout(() => { if (token === renderToken) showGraphMsg("<div>Loading data…</div>"); }, 150);
    try { return await promise; } finally { clearTimeout(t); if (token === renderToken && /Loading/.test($("graphMsg").textContent)) showGraphMsg(null); }
  }

  async function renderMeasured(token) {
    const gi = state.gene, c = state.cond;
    await withLoading(token, ensure(c, [gi]));
    if (token !== renderToken) return;
    const all = getTargets(c, gi);
    if (!all || !all.length) { showGraphMsg(noTargetsMsg(gi, c)); renderTable([], "measured"); return; }
    const ids = all.slice(0, state.topN).map((t) => t[0]);
    const pos = measuredPositions(ids);
    const { els } = buildMeasured(c, gi, state.topN, pos, null);
    cy1 = makeCy($("cy1"), els);
    renderTable(all.map(([t, lfc, p]) => ({ gi: t, lfc, padj: p, shown: ids.includes(t) })), "measured");
  }

  function noTargetsMsg(gi, c) {
    const p = DB.pert[c].get(gi);
    if (p) return `<div><b>${esc(name(gi))}</b> knockdown in ${COND_LABEL[c]} changed no other genes significantly.<br><span class="small">Try another condition, or a gene from the examples.</span></div>`;
    if (DB.removed[c].has(gi)) return `<div>${esc(name(gi))} was removed in ${COND_LABEL[c]} by quality control.</div>`;
    return `<div>${esc(name(gi))} was not tested in ${COND_LABEL[c]}.</div>`;
  }

  // ---------- Chain reaction view ----------
  const chainKey = (cond, gi, depth, cutoff, up) => `${cond}|${gi}|${depth}|${cutoff}|${up}`;
  async function runChain(cond, gi, depth, cutoff, up) {
    const key = chainKey(cond, gi, depth, cutoff, up);
    if (chainCache.has(key)) return chainCache.get(key);
    const res = await KDModel.predictChain({
      source: gi, depth, cutoff, propagateUp: up,
      getTargets: (g) => getTargets(cond, g), getKd: (g) => getKd(cond, g),
      ensure: (gs) => ensure(cond, gs),
    });
    if (chainCache.size > 30) chainCache.delete(chainCache.keys().next().value);
    chainCache.set(key, res);
    return res;
  }

  async function renderChain(token) {
    const gi = state.gene, c = state.cond;
    await withLoading(token, ensure(c, [gi]));
    if (token !== renderToken) return;
    if (!(getTargets(c, gi) || []).length) { showGraphMsg(noTargetsMsg(gi, c)); renderTable([], "chain"); renderCheck(null); return; }
    showGraphMsg("<div>Computing chain reaction… (first time in a condition may download a few MB)</div>");
    const res = await runChain(c, gi, state.depth, state.cutoff, state.up);
    if (token !== renderToken) return;
    showGraphMsg(null);
    const nodes = [...res.nodes.values()];

    // choose which nodes to draw: about half the budget for step 1, the rest for later steps
    const N = state.topN;
    const byStep = [[], [], [], []];
    for (const n of nodes) byStep[n.step].push(n);
    byStep[1].sort((a, b) => Math.abs(b.direct) - Math.abs(a.direct));
    for (let k = 2; k <= 3; k++) byStep[k].sort((a, b) => Math.abs(b.total) - Math.abs(a.total));
    const later = byStep[2].length + byStep[3].length;
    const n1 = later ? Math.max(Math.ceil(N / 2), N - later) : N;
    const chosen = new Set(byStep[1].slice(0, n1).map((n) => n.gi));
    let budget = N - chosen.size;
    const laterSorted = byStep[2].concat(byStep[3]).sort((a, b) => Math.abs(b.total) - Math.abs(a.total));
    for (const n of laterSorted) { if (budget <= 0) break; chosen.add(n.gi); budget--; }
    // make sure each drawn step 2+ gene has its strongest parent drawn (so paths are visible)
    for (const g of [...chosen]) {
      let n = res.nodes.get(g), guard = 0;
      while (n && n.step > 1 && guard++ < 3) {
        const parents = n.incoming.filter((e) => e.step === n.step && res.nodes.has(e.from)).sort((a, b) => Math.abs(b.w) - Math.abs(a.w));
        if (!parents.length) break;
        if (parents.some((e) => chosen.has(e.from))) break;
        chosen.add(parents[0].from);
        n = res.nodes.get(parents[0].from);
      }
    }

    const rings = [[String(gi)], [], [], []];
    const sorted = [...chosen].map((g) => res.nodes.get(g)).sort((a, b) => Math.abs(b.total) - Math.abs(a.total));
    for (const n of sorted) rings[n.step].push(String(n.gi));
    const pos = ringPositions(rings, 40);
    setSizeMax(sorted.map((n) => n.total));
    const els = [{ data: { id: String(gi), label: name(gi), src: true, v: -1, size: 48 }, position: pos[String(gi)] }];
    for (const n of sorted) els.push({ data: { id: String(n.gi), label: name(n.gi), v: n.total, size: nodeSize(n.total) }, position: pos[String(n.gi)] });
    for (const n of sorted) {
      for (const e of n.incoming) {
        if (e.from !== gi && !chosen.has(e.from)) continue;
        const pred = e.step > 1;
        els.push({ data: { id: `e${e.from}-${n.gi}-${e.step}`, source: String(e.from), target: String(n.gi), pred, w: 0.8 + 2.2 * Math.min(Math.abs(e.w) / 1.5, 1) } });
      }
    }
    cy1 = makeCy($("cy1"), els);
    renderTable(nodes.map((n) => ({ gi: n.gi, step: n.step, total: n.total, direct: n.direct || null, indirect: n.indirect || null, shown: chosen.has(n.gi) })), "chain");
    renderCheck(KDModel.modelCheck(res.nodes, getTargets(c, gi), gi), res.steps);
  }

  // ---------- Compare view ----------
  async function renderCompare(token) {
    const gi = state.gene, a = state.cond, b = state.cond2;
    await withLoading(token, Promise.all([ensure(a, [gi]), ensure(b, [gi])]));
    if (token !== renderToken) return;
    const ta = getTargets(a, gi) || [], tb = getTargets(b, gi) || [];
    const ma = new Map(ta.map((t) => [t[0], t])), mb = new Map(tb.map((t) => [t[0], t]));
    $("paneTitle1").textContent = `${COND_LABEL[a]} (${ta.length ? DB.pert[a].get(gi).nSig : 0} targets)`;
    $("paneTitle2").textContent = `${COND_LABEL[b]} (${tb.length ? DB.pert[b].get(gi).nSig : 0} targets)`;
    // shared layout: union of displayed genes, strongest first, same positions in both panes
    const showA = ta.slice(0, state.topN).map((t) => t[0]), showB = tb.slice(0, state.topN).map((t) => t[0]);
    const union = [...new Set(showA.concat(showB))].sort((x, y) => {
      const vx = Math.max(Math.abs(ma.get(x)?.[1] || 0), Math.abs(mb.get(x)?.[1] || 0));
      const vy = Math.max(Math.abs(ma.get(y)?.[1] || 0), Math.abs(mb.get(y)?.[1] || 0));
      return vy - vx;
    });
    const pos = measuredPositions(union);
    setSizeMax(ta.slice(0, state.topN).concat(tb.slice(0, state.topN)).map((t) => t[1]));
    const setA = new Set(ma.keys()), setB = new Set(mb.keys());
    if (ta.length) cy1 = makeCy($("cy1"), buildMeasured(a, gi, state.topN, pos, setB).els);
    if (tb.length) cy2 = makeCy($("cy2"), buildMeasured(b, gi, state.topN, pos, setA).els);
    if (!ta.length && !tb.length) showGraphMsg(noTargetsMsg(gi, a));
    else if (!ta.length || !tb.length) {
      const missing = !ta.length ? a : b;
      const box = !ta.length ? $("cy1") : $("cy2");
      box.innerHTML = `<div class="graph-msg" style="position:static;height:100%">${noTargetsMsg(gi, missing)}</div>`;
    }
    // table + stats
    const all = new Set([...ma.keys(), ...mb.keys()]);
    const rows = [...all].map((t) => ({
      gi: t, a: ma.get(t)?.[1] ?? null, b: mb.get(t)?.[1] ?? null,
      status: ma.has(t) && mb.has(t) ? "both" : ma.has(t) ? `only ${COND_LABEL[a]}` : `only ${COND_LABEL[b]}`,
      shown: showA.includes(t) || showB.includes(t),
    }));
    renderTable(rows, "compare");
    const shared = rows.filter((r) => r.status === "both");
    const r = KDModel.pearson(shared.map((x) => x.a), shared.map((x) => x.b));
    const sameDir = shared.filter((x) => Math.sign(x.a) === Math.sign(x.b)).length;
    const pa = DB.pert[a].get(gi), pb = DB.pert[b].get(gi);
    const capped = (pa && pa.nSig > ta.length) || (pb && pb.nSig > tb.length);
    $("compareStats").innerHTML = `
      <div class="cmp-row"><span>Affected in both</span><b>${shared.length}</b></div>
      <div class="cmp-row"><span>Only in ${COND_LABEL[a]}</span><b>${rows.length - shared.length - rows.filter((x) => x.status === `only ${COND_LABEL[b]}`).length}</b></div>
      <div class="cmp-row"><span>Only in ${COND_LABEL[b]}</span><b>${rows.filter((x) => x.status === `only ${COND_LABEL[b]}`).length}</b></div>
      <div class="cmp-row"><span>Overlap (shared ÷ all)</span><b>${rows.length ? Math.round((100 * shared.length) / rows.length) : 0}%</b></div>
      <div class="cmp-row"><span>Same direction (shared genes)</span><b>${shared.length ? Math.round((100 * sameDir) / shared.length) + "%" : "–"}</b></div>
      <div class="cmp-row"><span>Correlation of shared effects</span><b>${r == null ? "–" : r.toFixed(2)}</b></div>
      <p class="small muted">Green rings mark genes affected in only one condition. Same gene = same position in both panels.${capped ? " Note: only the top 200 targets per condition are stored, so some \"only\" genes may be weaker hits in the other condition." : ""}</p>`;
  }

  // =====================================================================
  // Details panel
  // =====================================================================
  function showDetails(gi, cond) {
    const el = $("details");
    const src = state.gene;
    if (gi === src) {
      el.innerHTML = `<h2>${esc(name(gi))}</h2><p class="small">This is the knocked-down gene (center). Click any other circle to see how it changed.</p>`;
      return;
    }
    const meas = (getTargets(cond, src) || []).find((t) => t[0] === gi);
    if (state.mode === "compare") return showCompareDetails(gi);
    let chainInfo = "";
    if (state.mode === "chain") {
      const res = chainCache.get(chainKey(state.cond, src, state.depth, state.cutoff, state.up));
      const n = res && res.nodes.get(gi);
      if (n) {
        const paths = n.incoming.slice().sort((a, b) => Math.abs(b.w) - Math.abs(a.w))
          .map((e) => `${esc(name(e.from))} (step ${e.step}, ${fmt(e.w)})`).join("<br>");
        chainInfo = `<dt>First reached</dt><dd>step ${n.step}</dd><dt>Predicted total</dt><dd>${fmt(n.total)}</dd><dt>…of which indirect</dt><dd>${fmt(n.indirect)}</dd><dt>Strongest paths from</dt><dd>${paths || "–"}</dd>`;
      }
    }
    const p = DB.pert[cond].get(gi);
    const anyCond = CONDS.find((c) => DB.pert[c].has(gi));
    el.innerHTML = `<h2>${esc(name(gi))}</h2>
      <dl>
        <dt>Measured change</dt><dd>${meas ? `${dirTag(meas[1])} ${fmt(meas[1])} log₂FC (~${level(meas[1])})` : "not significant"}</dd>
        ${meas ? `<dt>Adjusted p</dt><dd>${fmtP(meas[2])}</dd>` : ""}
        ${chainInfo}
        <dt>Own knockdown</dt><dd>${p ? `${p.nSig} targets in ${COND_LABEL[cond]}${p.kd != null ? `, reduced to ~${pctLeft(p.kd)}%` : ""}` : DB.removed[cond].has(gi) ? "filtered by QC here" : "not knocked down in this condition"}</dd>
      </dl>
      <button class="btn" id="kdThis" ${anyCond ? "" : "disabled"}>${anyCond ? `Knock ${esc(name(gi))} down instead` : "Never knocked down in this screen"}</button>`;
    if (anyCond) $("kdThis").addEventListener("click", () => { selectGene(gi); window.scrollTo({ top: 0, behavior: reduceMotion() ? "auto" : "smooth" }); });
  }

  function showCompareDetails(gi) {
    const src = state.gene, el = $("details");
    const row = (c) => {
      const m = (getTargets(c, src) || []).find((t) => t[0] === gi);
      return `<dt>${COND_LABEL[c]}</dt><dd>${m ? `${dirTag(m[1])} ${fmt(m[1])} log₂FC (~${level(m[1])}), adj. p ${fmtP(m[2])}` : "not significant"}</dd>`;
    };
    const anyCond = CONDS.find((c) => DB.pert[c].has(gi));
    el.innerHTML = `<h2>${esc(name(gi))}</h2><p class="small muted" style="margin:0">Change after knocking down ${esc(name(src))}:</p>
      <dl>${row(state.cond)}${row(state.cond2)}</dl>
      <button class="btn" id="kdThis" ${anyCond ? "" : "disabled"}>${anyCond ? `Knock ${esc(name(gi))} down instead` : "Never knocked down in this screen"}</button>`;
    if (anyCond) $("kdThis").addEventListener("click", () => { selectGene(gi); window.scrollTo({ top: 0, behavior: reduceMotion() ? "auto" : "smooth" }); });
  }

  // =====================================================================
  // Table
  // =====================================================================
  const TABLES = {
    measured: {
      title: "Measured effects",
      cols: [
        { k: "gene", t: "Gene", f: (r) => `<button type="button" class="genelink">${esc(name(r.gi))}</button>`, s: (r) => name(r.gi) },
        { k: "dir", t: "Direction", f: (r) => dirTag(r.lfc), s: (r) => Math.sign(r.lfc) },
        { k: "lfc", t: "Effect (log₂FC)", f: (r) => fmt(r.lfc), s: (r) => Math.abs(r.lfc) },
        { k: "pct", t: "Level vs normal", f: (r) => level(r.lfc), s: (r) => r.lfc },
        { k: "padj", t: "Adjusted p", f: (r) => fmtP(r.padj), s: (r) => -r.padj },
      ],
      sortKey: "lfc",
    },
    chain: {
      title: "Predicted chain reaction",
      cols: [
        { k: "gene", t: "Gene", f: (r) => `<button type="button" class="genelink">${esc(name(r.gi))}</button>`, s: (r) => name(r.gi) },
        { k: "step", t: "Step", f: (r) => r.step, s: (r) => -r.step },
        { k: "dir", t: "Direction", f: (r) => dirTag(r.total), s: (r) => Math.sign(r.total) },
        { k: "total", t: "Predicted total", f: (r) => fmt(r.total), s: (r) => Math.abs(r.total) },
        { k: "direct", t: "Measured (step 1)", f: (r) => fmt(r.direct), s: (r) => Math.abs(r.direct || 0) },
        { k: "indirect", t: "Indirect (step 2+)", f: (r) => fmt(r.indirect), s: (r) => Math.abs(r.indirect || 0) },
      ],
      sortKey: "total",
    },
    compare: {
      title: "Compare conditions",
      cols: [
        { k: "gene", t: "Gene", f: (r) => `<button type="button" class="genelink">${esc(name(r.gi))}</button>`, s: (r) => name(r.gi) },
        { k: "a", t: () => COND_LABEL[state.cond] + " (log₂FC)", f: (r) => fmt(r.a), s: (r) => Math.abs(r.a || 0) },
        { k: "b", t: () => COND_LABEL[state.cond2] + " (log₂FC)", f: (r) => fmt(r.b), s: (r) => Math.abs(r.b || 0) },
        { k: "status", t: "Affected in", f: (r) => (r.status === "both" ? '<span class="tag neutral">both</span>' : `<span class="tag only">${r.status}</span>`), s: (r) => r.status },
      ],
      sortKey: "a",
    },
  };

  function renderTable(rows, kind) {
    const T = TABLES[kind];
    current = { rows, kind, sortKey: T.sortKey, sortDir: -1 };
    $("tableTitle").textContent = T.title;
    drawTable();
  }

  function drawTable() {
    const T = TABLES[current.kind];
    const col = T.cols.find((c) => c.k === current.sortKey) || T.cols[0];
    const rows = current.rows.slice().sort((a, b) => {
      const x = col.s(a), y = col.s(b);
      return (x < y ? -1 : x > y ? 1 : 0) * current.sortDir;
    });
    const LIMIT = 500;
    $("tableCount").textContent = rows.length ? `${rows.length.toLocaleString()} genes${rows.length > LIMIT ? ` (showing first ${LIMIT})` : ""} · click a header to sort, a row to select` : "no genes";
    const thead = $("table").tHead;
    thead.innerHTML = "<tr>" + T.cols.map((c) => `<th data-k="${c.k}"${c.k === current.sortKey ? ` aria-sort="${current.sortDir < 0 ? "descending" : "ascending"}"` : ""}><button type="button" class="th-btn">${typeof c.t === "function" ? c.t() : c.t}${c.k === current.sortKey ? (current.sortDir < 0 ? " ▼" : " ▲") : ""}</button></th>`).join("") + "</tr>";
    thead.querySelectorAll("th").forEach((th) => th.addEventListener("click", () => {
      if (current.sortKey === th.dataset.k) current.sortDir *= -1; else { current.sortKey = th.dataset.k; current.sortDir = -1; }
      drawTable();
    }));
    const tb = $("table").tBodies[0];
    tb.innerHTML = rows.slice(0, LIMIT).map((r) => `<tr data-g="${r.gi}">${T.cols.map((c) => `<td>${c.f(r)}</td>`).join("")}</tr>`).join("");
  }

  $("table").tBodies[0].addEventListener("click", (e) => {
    const tr = e.target.closest("tr[data-g]"); if (!tr) return;
    const gi = +tr.dataset.g;
    let found = false;
    for (const cy of [cy1, cy2]) {
      const node = cy && cy.getElementById(String(gi));
      if (node && node.length) { cy.$(":selected").unselect(); node.select(); onNodeTap(cy, node); cy.animate({ center: { eles: node } }, { duration: reduceMotion() ? 0 : 300 }); found = true; }
    }
    if (!found) showDetails(gi, state.cond);
    $("details").scrollIntoView({ behavior: reduceMotion() ? "auto" : "smooth", block: "nearest" });
  });

  // =====================================================================
  // Model check panel (scatter of predicted indirect vs measured)
  // =====================================================================
  let lastCheck = null;
  const MIN_CHECK_GENES = 10;
  function renderCheck(chk, steps) {
    const box = $("scatter"), stats = $("checkStats");
    lastCheck = chk ? { chk, steps } : null;
    if (!chk) { box.innerHTML = ""; stats.innerHTML = ""; renderGlobalCheck(); return; }
    if (state.depth < 2) {
      box.innerHTML = '<p class="small muted">Set chain depth to 2 or 3: step 1 <em>is</em> the measurement, so only steps 2+ are a real prediction.</p>';
      stats.innerHTML = ""; renderGlobalCheck(); return;
    }
    if (steps < 2) {
      box.innerHTML = `<p class="small muted">The chain stops at step 1: none of the affected genes ${state.up ? "" : "that went <b>down</b> "}were themselves knocked down in this condition with a usable knockdown, so the model has nothing to pass on. Try a lower minimum effect${state.up ? "" : ", tick \"also pass effects through genes that went up\""}, or a gene with more targets.</p>`;
      stats.innerHTML = ""; renderGlobalCheck(); return;
    }
    const pts = chk.points;
    const W = 300, H = 230, m = { l: 34, r: 8, t: 8, b: 30 };
    const lim = Math.max(0.2, ...pts.map((d) => Math.max(Math.abs(d.pred), Math.abs(d.meas)))) * 1.05;
    const sx = (v) => m.l + ((v + lim) / (2 * lim)) * (W - m.l - m.r);
    const sy = (v) => H - m.b - ((v + lim) / (2 * lim)) * (H - m.t - m.b);
    const down = css("--down"), up = css("--up"), border = css("--border"), muted = css("--muted");
    const ticks = [-lim, -lim / 2, 0, lim / 2, lim].map((v) => +v.toFixed(1));
    let svg = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Scatter of predicted vs measured effects">`;
    svg += `<rect x="${m.l}" y="${m.t}" width="${W - m.l - m.r}" height="${H - m.t - m.b}" fill="none" stroke="${border}"/>`;
    svg += `<line x1="${sx(0)}" x2="${sx(0)}" y1="${m.t}" y2="${H - m.b}" stroke="${border}"/><line y1="${sy(0)}" y2="${sy(0)}" x1="${m.l}" x2="${W - m.r}" stroke="${border}"/>`;
    svg += `<line x1="${sx(-lim)}" y1="${sy(-lim)}" x2="${sx(lim)}" y2="${sy(lim)}" stroke="${muted}" stroke-dasharray="3 3" opacity=".6"/>`;
    for (const t of ticks) {
      svg += `<text x="${sx(t)}" y="${H - m.b + 12}" text-anchor="middle">${t}</text>`;
      svg += `<text x="${m.l - 4}" y="${sy(t) + 3}" text-anchor="end">${t}</text>`;
    }
    svg += `<text x="${(W + m.l) / 2}" y="${H - 3}" text-anchor="middle">predicted indirect effect (step 2+)</text>`;
    svg += `<text transform="translate(9 ${(H - m.b) / 2}) rotate(-90)" text-anchor="middle">measured effect</text>`;
    for (const d of pts) {
      svg += `<circle cx="${sx(d.pred).toFixed(1)}" cy="${sy(d.meas).toFixed(1)}" r="3.2" fill="${d.meas < 0 ? down : up}" fill-opacity=".7"><title>${esc(name(d.gi))}: predicted ${fmt(d.pred)}, measured ${fmt(d.meas)}</title></circle>`;
    }
    svg += "</svg>";
    box.innerHTML = pts.length ? svg : '<p class="small muted">No measured genes received an indirect prediction. Try a lower minimum effect or more depth.</p>';
    const few = chk.n < MIN_CHECK_GENES;
    const verdict = chk.r == null ? "" : few ? `too few genes (< ${MIN_CHECK_GENES}) to judge` : chk.r > 0.5 ? "good match" : chk.r > 0.2 ? "partial match" : chk.r > -0.1 ? "little or no match" : "opposite direction";
    stats.innerHTML = `
      <div><div class="big"${few ? ' style="opacity:.45"' : ""}>${chk.r == null ? "–" : chk.r.toFixed(2)}</div><div class="small muted">correlation (r)${verdict ? ` · ${verdict}` : ""}</div></div>
      <div><div class="big">${chk.signAgreement == null ? "–" : Math.round(chk.signAgreement * 100) + "%"}</div><div class="small muted">same direction</div></div>
      <div><div class="big">${chk.n}</div><div class="small muted">genes compared</div></div>`;
    renderGlobalCheck();
  }

  function renderGlobalCheck() {
    const el = $("globalCheck"), d = DB.check;
    if (!d) { el.innerHTML = ""; return; }
    const v = d[state.up ? "up_and_down" : "down_only"];
    const c = v && v.conditions[state.cond];
    if (!c) { el.innerHTML = ""; return; }
    // histogram of per-gene correlations
    const bins = c.hist, max = Math.max(...bins.map((b) => b.count));
    const W = 300, H = 90, P = 10, bw = (W - 2 * P) / bins.length;
    const acc = css("--accent"), muted = css("--muted");
    let svg = `<svg viewBox="0 0 ${W} ${H + 16}" role="img" aria-label="Histogram of model-check correlations">`;
    bins.forEach((b, i) => {
      const h = max ? (b.count / max) * H : 0;
      svg += `<rect x="${P + i * bw + 1}" y="${H - h}" width="${bw - 2}" height="${h}" fill="${acc}" opacity="${b.lo >= 0 ? 0.9 : 0.45}"><title>r ${b.lo.toFixed(1)} to ${b.hi.toFixed(1)}: ${b.count} genes</title></rect>`;
    });
    [-1, -0.5, 0, 0.5, 1].forEach((v) => { svg += `<text x="${P + ((v + 1) / 2) * (W - 2 * P)}" y="${H + 14}" text-anchor="middle" fill="${muted}" font-size="12">${v}</text>`; });
    svg += "</svg>";
    el.innerHTML = `<h3 style="font-size:.9rem">Across all genes (${COND_LABEL[state.cond]})</h3>
      <p class="small" style="margin:0 0 4px">Median correlation <b>${c.median_r.toFixed(2)}</b> over ${c.n_genes.toLocaleString()} knockdowns (depth ${d.settings.depth}, min effect ${d.settings.cutoff}, ${state.up ? "up + down paths" : "down paths only"}); ${Math.round(c.frac_positive * 100)}% positive.</p>${svg}
      <p class="small" style="margin:4px 0 0">Baseline: comparing each prediction with a <em>random other</em> knockdown gives <b>${c.baseline_median_r.toFixed(2)}</b>. The closer these two numbers are, the more the match comes from a generic response shared by many knockdowns rather than from the chosen gene.</p>`;
  }

  // =====================================================================
  // Main update
  // =====================================================================
  async function update() {
    const token = ++renderToken;
    document.body.dataset.busy = "1";
    syncControls();
    writeHash();
    renderSummary();
    announce();
    destroyCys();
    $("cy1").innerHTML = ""; $("cy2").innerHTML = "";
    showGraphMsg(null);
    $("details").innerHTML = '<h2>Gene details</h2><p class="muted">Click any circle in the network to see its details.</p>';
    if (state.gene == null) { showGraphMsg("Pick a gene to knock down."); return; }
    try {
      if (state.mode === "measured") await renderMeasured(token);
      else if (state.mode === "chain") await renderChain(token);
      else await renderCompare(token);
    } catch (err) {
      console.error(err);
      if (token === renderToken) showGraphMsg(`<div class="warn">Something went wrong loading data: ${esc(err.message)}</div>`);
    }
    if (token === renderToken) {
      // a fast view switch can leave a stale "Computing..." message behind; clear it
      if (state.mode !== "chain" && /Computing/.test($("graphMsg").textContent)) showGraphMsg(null);
      document.body.dataset.busy = "0";
    }
  }

  async function init() {
    try { const t = localStorage.getItem("kd-theme"); if (t) document.documentElement.dataset.theme = t; } catch (e) {}
    setupControls();
    syncControls();
    if (typeof cytoscape === "undefined") {
      $("summary").innerHTML = '<p class="warn">Could not load the network library (Cytoscape.js) from the CDN. Check your internet connection.</p>';
      return;
    }
    try { await loadIndex(); } catch (err) {
      $("summary").innerHTML = `<p class="warn">Could not load data/index.json (${esc(err.message)}). If you opened index.html directly from disk, run a local web server instead (see README).</p>`;
      return;
    }
    renderExamples();
    // react when someone edits or pastes a link with a different #hash
    const fromUrl = () => { const g = readHash(); if (g != null) state.gene = g; lastMajor = null; update(); };
    window.addEventListener("hashchange", fromUrl);
    window.addEventListener("popstate", fromUrl);
    const fromHash = readHash();
    state.gene = fromHash != null ? fromHash : DB.byName.get("GATA3") ?? DB.perturbed[0];
    if (fromHash == null && !DB.pert[state.cond].has(state.gene)) selectGene(state.gene);
    else update();
  }

  init();
})();
