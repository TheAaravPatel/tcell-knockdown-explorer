#!/usr/bin/env node
/*
 * model_check.js - run the chain-reaction model (js/model.js, the same code the
 * website uses) on every knockdown, and measure how well its INDIRECT (step 2+)
 * predictions match the measured effects.
 *
 *   node scripts/model_check.js [depth=2] [cutoff=0.1]
 *
 * Writes data/model_check.json (used by the website's "Across all genes" chart)
 * and prints a summary.
 *
 * Baseline: compare gene A's prediction with a RANDOM other knockdown's measured
 * effects. If the real match is not better than this, the model is only picking
 * up genes that respond to almost any knockdown, not anything specific to A.
 */
const fs = require("fs");
const path = require("path");
const KDModel = require("../js/model.js");

const DATA = path.join(__dirname, "..", "data");
const DEPTH = +(process.argv[2] || 2);
const CUTOFF = +(process.argv[3] || 0.1);
const MIN_TARGETS = 20;   // only test knockdowns with at least this many measured targets
const MIN_POINTS = 10;    // need at least this many genes in both sets to compute r

const ix = JSON.parse(fs.readFileSync(path.join(DATA, "index.json")));
const fi = Object.fromEntries(ix.pert_fields.map((k, i) => [k, i]));

// simple seeded random generator so results are reproducible
let seed = 42;
const rand = () => ((seed = (seed * 1664525 + 1013904223) >>> 0) / 2 ** 32);

const median = (a) => { const s = a.slice().sort((x, y) => x - y); const m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; };

(async () => {
  const out = { settings: { depth: DEPTH, cutoff: CUTOFF, min_targets: MIN_TARGETS, min_points: MIN_POINTS } };
  for (const variant of ["down_only", "up_and_down"]) {
  const propagateUp = variant === "up_and_down";
  out[variant] = { conditions: {} };
  console.log(`\n=== ${variant === "down_only" ? "Only genes that went DOWN pass effects on (default)" : "Genes that went UP or DOWN pass effects on"} ===`);
  for (const cond of ix.conditions) {
    // load every shard of this condition
    const targets = new Map(), kd = new Map(), nSig = new Map();
    for (const f of fs.readdirSync(path.join(DATA, cond))) {
      const d = JSON.parse(fs.readFileSync(path.join(DATA, cond, f)));
      for (const k in d) targets.set(+k, d[k]);
    }
    for (const a of ix.pert[cond]) { kd.set(a[fi.gene], a[fi.kd_lfc]); nSig.set(a[fi.gene], a[fi.n_sig]); }
    const sources = [...targets.keys()].filter((g) => targets.get(g).length >= MIN_TARGETS);

    const rows = [];
    for (const g of sources) {
      const res = await KDModel.predictChain({
        source: g, depth: DEPTH, cutoff: CUTOFF, propagateUp,
        getTargets: (x) => targets.get(x) || (kd.has(x) ? [] : null),
        getKd: (x) => (kd.has(x) ? kd.get(x) : null),
      });
      const chk = KDModel.modelCheck(res.nodes, targets.get(g), g);
      if (chk.n < MIN_POINTS || chk.r == null) continue;
      // baseline: same prediction vs a random other knockdown's measurements
      let other;
      do { other = sources[Math.floor(rand() * sources.length)]; } while (other === g);
      const base = KDModel.modelCheck(res.nodes, targets.get(other).filter((t) => t[0] !== g), g);
      rows.push({ gene: ix.genes[g], n_sig: nSig.get(g), n: chk.n, r: chk.r, sign: chk.signAgreement,
                  base_r: base.n >= MIN_POINTS ? base.r : null });
    }
    const rs = rows.map((r) => r.r);
    const bs = rows.map((r) => r.base_r).filter((x) => x != null);
    const hist = [];
    for (let i = 0; i < 20; i++) {
      const lo = -1 + i * 0.1, hi = lo + 0.1;
      hist.push({ lo, hi, count: rs.filter((r) => r >= lo && (r < hi || (i === 19 && r <= 1))).length });
    }
    const top = rows.filter((r) => r.n >= 30).sort((a, b) => b.r - a.r).slice(0, 15);
    const worst = rows.filter((r) => r.n >= 30).sort((a, b) => a.r - b.r).slice(0, 10);
    out[variant].conditions[cond] = {
      n_tested: sources.length, n_genes: rows.length,
      median_r: median(rs), mean_r: rs.reduce((a, b) => a + b, 0) / rs.length,
      frac_positive: rs.filter((r) => r > 0).length / rs.length,
      median_sign_agreement: median(rows.map((r) => r.sign)),
      baseline_median_r: bs.length ? median(bs) : null,
      hist,
      best: top.map((r) => ({ gene: r.gene, r: +r.r.toFixed(3), n: r.n, n_sig: r.n_sig })),
      worst: worst.map((r) => ({ gene: r.gene, r: +r.r.toFixed(3), n: r.n, n_sig: r.n_sig })),
    };
    const c = out[variant].conditions[cond];
    console.log(`${cond.padEnd(9)} tested=${c.n_genes}/${c.n_tested}  median r=${c.median_r.toFixed(3)}  ` +
      `positive=${(100 * c.frac_positive).toFixed(0)}%  sign agreement=${(100 * c.median_sign_agreement).toFixed(0)}%  ` +
      `random-pair baseline median r=${c.baseline_median_r?.toFixed(3)}`);
    console.log(`   best: ${c.best.slice(0, 8).map((b) => `${b.gene}(${b.r.toFixed(2)})`).join(" ")}`);
  }
  }
  fs.writeFileSync(path.join(DATA, "model_check.json"), JSON.stringify(out));
})();
