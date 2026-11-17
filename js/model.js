/*
 * model.js - the "predicted chain reaction" model and the model check.
 *
 * Shared by the website (window.KDModel) and Node scripts (require).
 *
 * Model (simple linear propagation):
 *   Step 1: knocking down gene A changes each measured target B by lfc(A->B).
 *   Step k+1: every gene B that changed at step k, and was itself knocked down
 *     in the experiment (same condition), passes on its own measured effects,
 *     scaled by how much it changed relative to its own knockdown:
 *         scale(B)   = change_k(B) / kd(B)     (both are log2 fold changes)
 *         change_{k+1}(C) += scale(B) * lfc(B->C)
 *   Contributions from different paths are added up. Totals are capped.
 *   Genes never perturbed in the screen cannot pass anything on.
 *   Propagation stops at the chosen depth, or for changes smaller than `cutoff`.
 *   By default only genes that went DOWN pass effects on (see propagateUp).
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.KDModel = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  const DEFAULTS = {
    kdMin: 0.3,      // B passes effects on only if its own knockdown was at least this strong (|log2FC|)
    scaleCap: 2,     // |scale(B)| capped: B can't act like more than 2x its own knockdown
    effectCap: 5,    // predicted total change per gene capped at +/-5 log2FC (32-fold)
    maxIncoming: 3,  // remember the 3 strongest incoming paths per gene (for drawing)
    // The screen only measured what happens when a gene goes DOWN. By default a
    // gene that went UP does not pass effects on (that would assume "up" is the
    // exact mirror image of a knockdown, which the model check shows is wrong).
    propagateUp: false,
  };

  const clamp = (x, c) => (x > c ? c : x < -c ? -c : x);

  /**
   * Run the chain reaction.
   * @param {object} p
   *   source     gene index that is knocked down
   *   depth      1..3
   *   cutoff     minimum |change| (log2FC) to keep propagating / show a gene
   *   getTargets gi -> array of [targetGi, lfc, padj] or null (must be loaded)
   *   getKd      gi -> on-target knockdown log2FC (negative) or null
   *   propagateUp also pass effects through genes that went up (default false)
   *   ensure     optional async (array of gi) -> loads data for those genes
   * @returns {Promise<{nodes: Map, steps: number}>}
   *   nodes: gi -> {step, total, direct, indirect, byStep[], incoming[]}
   *     step     first step at which the gene changed by >= cutoff
   *     direct   measured step-1 effect (or 0)
   *     indirect sum of step 2+ contributions (the model's real "prediction")
   *     incoming [{from, w, step}] strongest contributing paths
   */
  async function predictChain(p) {
    const o = Object.assign({}, DEFAULTS, p);
    const nodes = new Map();
    const node = (gi) => {
      let n = nodes.get(gi);
      if (!n) {
        n = { gi, step: 0, total: 0, direct: 0, indirect: 0, byStep: [0, 0, 0, 0], incoming: [] };
        nodes.set(gi, n);
      }
      return n;
    };
    const addIncoming = (n, from, w, step) => {
      const inc = n.incoming;
      if (inc.length < o.maxIncoming) inc.push({ from, w, step });
      else {
        let mi = 0;
        for (let i = 1; i < inc.length; i++) if (Math.abs(inc[i].w) < Math.abs(inc[mi].w)) mi = i;
        if (Math.abs(w) > Math.abs(inc[mi].w)) inc[mi] = { from, w, step };
      }
    };

    // Step 1 = measured effects of the source
    if (o.ensure) await o.ensure([o.source]);
    let frontier = new Map(); // gi -> change at the current step
    for (const [t, lfc] of o.getTargets(o.source) || []) {
      if (t === o.source) continue;
      frontier.set(t, lfc);
      const n = node(t);
      n.direct = lfc;
      n.byStep[1] = lfc;
      if (Math.abs(lfc) >= o.cutoff) { n.step = 1; addIncoming(n, o.source, lfc, 1); }
    }

    let steps = 1;
    for (let k = 2; k <= o.depth; k++) {
      // genes that can pass effects on: changed enough, perturbed, real knockdown
      const movers = [];
      for (const [b, ch] of frontier) {
        if (Math.abs(ch) < o.cutoff) continue;
        if (ch > 0 && !o.propagateUp) continue;
        const kd = o.getKd(b);
        if (kd == null || kd > -o.kdMin) continue;
        movers.push(b);
      }
      if (!movers.length) break;
      if (o.ensure) await o.ensure(movers);
      const next = new Map();
      for (const b of movers) {
        const targets = o.getTargets(b);
        if (!targets) continue;
        const scale = clamp(frontier.get(b) / o.getKd(b), o.scaleCap);
        for (const [c, lfc] of targets) {
          if (c === o.source) continue; // the source itself is held down by CRISPRi
          const w = scale * lfc;
          next.set(c, (next.get(c) || 0) + w);
          if (Math.abs(w) >= o.cutoff) addIncoming(node(c), b, w, k);
        }
      }
      for (const [c, ch] of next) {
        const n = node(c);
        n.byStep[k] += ch;
        n.indirect += ch;
        if (!n.step && Math.abs(ch) >= o.cutoff) n.step = k;
      }
      frontier = next;
      steps = k;
    }

    for (const n of nodes.values()) {
      n.total = clamp(n.direct + n.indirect, o.effectCap);
      n.indirect = clamp(n.indirect, o.effectCap);
    }
    // drop genes that never reached the cutoff
    for (const [gi, n] of nodes) if (!n.step) nodes.delete(gi);
    return { nodes, steps };
  }

  /** Pearson correlation of two equal-length arrays (null if < 3 points or no variance). */
  function pearson(x, y) {
    const n = x.length;
    if (n < 3) return null;
    let sx = 0, sy = 0;
    for (let i = 0; i < n; i++) { sx += x[i]; sy += y[i]; }
    const mx = sx / n, my = sy / n;
    let sxy = 0, sxx = 0, syy = 0;
    for (let i = 0; i < n; i++) {
      const dx = x[i] - mx, dy = y[i] - my;
      sxy += dx * dy; sxx += dx * dx; syy += dy * dy;
    }
    if (sxx === 0 || syy === 0) return null;
    return sxy / Math.sqrt(sxx * syy);
  }

  /**
   * Model check: compare the INDIRECT prediction (step 2+, i.e. what the model
   * adds beyond the measured direct effects) with the measured effect, for
   * genes that are both measured-significant and indirectly predicted.
   * (Comparing step 1 to the measurement would be circular - they are equal.)
   */
  function modelCheck(chainNodes, measuredTargets, source) {
    const meas = new Map();
    for (const [t, lfc] of measuredTargets || []) if (t !== source) meas.set(t, lfc);
    const pts = [];
    let predictedOnly = 0;
    for (const n of chainNodes.values()) {
      if (n.indirect === 0) continue;
      if (meas.has(n.gi)) pts.push({ gi: n.gi, pred: n.indirect, meas: meas.get(n.gi) });
      else predictedOnly++;
    }
    const r = pearson(pts.map((d) => d.pred), pts.map((d) => d.meas));
    let agree = 0;
    for (const d of pts) if (Math.sign(d.pred) === Math.sign(d.meas)) agree++;
    return {
      points: pts,
      r,
      n: pts.length,
      signAgreement: pts.length ? agree / pts.length : null,
      predictedOnly,
      measuredCount: meas.size,
    };
  }

  return { predictChain, modelCheck, pearson, DEFAULTS };
});
