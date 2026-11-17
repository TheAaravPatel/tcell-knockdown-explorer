#!/usr/bin/env python3
"""
preprocess.py - turn the 17 GB Perturb-seq DE file into small JSON shards
for the Virtual Knockdown Explorer website.

Source data
-----------
Zhu, Dann et al. (Cell, 2026) genome-scale CRISPRi Perturb-seq in primary
human CD4+ T cells. File: GWCD4i.DE_stats.h5ad in the public bucket
s3://genome-scale-tcell-perturb-seq/marson2025_data/ (anonymous access).

The file is an AnnData object:
  obs    (33,983 rows)  one perturbed gene in one culture condition
  var    (10,282 cols)  measured genes
  layers log_fc, zscore, p_value, adj_p_value, lfcSE, baseMean
         each a dense, uncompressed float64 matrix (~2.8 GB)

Because the layers are stored contiguously, we can stream blocks of rows
directly from S3 with h5py + s3fs without downloading the whole file.
Only `log_fc` and `adj_p_value` are read (~5.6 GB of network traffic).

Usage
-----
    pip install h5py s3fs numpy pandas
    python preprocess.py                 # stream from S3 (default)
    python preprocess.py --local FILE    # use a downloaded copy instead

Output (in ./data)
------------------
  index.json                 every gene name + per-condition stats per perturbation
  <Condition>/<Letter>.json  measured targets for perturbed genes whose name
                             starts with <Letter> (non-letters -> "0")
"""

import argparse
import json
import os
import shutil
import time
from concurrent.futures import ThreadPoolExecutor

import h5py
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Settings - change these to re-run with different thresholds
# ---------------------------------------------------------------------------
PADJ_THRESHOLD = 0.05      # keep targets with adjusted p-value below this
TOP_N_TARGETS = 200        # keep at most this many targets per perturbation (by |log_fc|)
CONDITIONS = ["Rest", "Stim8hr", "Stim48hr"]
# Quality flags in obs; a perturbation is dropped if ANY of these is True
QUALITY_FLAGS = ["neighboring_gene_KD", "distal_offtarget_flag", "low_target_gex"]

S3_PATH = "genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad"
OUT_DIR = "data"
BLOCK_ROWS = 256           # rows streamed per request (~21 MB per layer)
N_WORKERS = 6              # parallel S3 connections


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def open_h5(local_path=None):
    """Open the .h5ad as an h5py File, either locally or remotely from S3."""
    if local_path:
        return h5py.File(local_path, "r")
    import s3fs
    fs = s3fs.S3FileSystem(anon=True)
    # cache_type="none": every h5py read becomes one ranged HTTP GET
    fobj = fs.open(S3_PATH, "rb", cache_type="none", block_size=2**20)
    return h5py.File(fobj, "r")


def read_column(group, key):
    """Read one AnnData dataframe column (handles categorical encoding)."""
    d = group[key]
    if isinstance(d, h5py.Group):  # categorical: codes + categories
        cats = d["categories"][:].astype(str)
        return pd.Categorical.from_codes(d["codes"][:], cats).astype(str)
    v = d[:]
    return v.astype(str) if v.dtype.kind in "OS" else v


def read_dataframe(group):
    return pd.DataFrame({k: read_column(group, k) for k in group.keys()})


def shard_key(name):
    """Shard files are split by the first letter of the perturbed gene name."""
    c = name[:1].upper()
    return c if "A" <= c <= "Z" else "0"


def r3(x):
    """Round to 3 decimals to keep JSON small."""
    return float(round(float(x), 3))


def sig2(p):
    """Keep 2 significant digits of a p-value (e.g. 0.0034, 1.2e-20)."""
    return float(f"{float(p):.2g}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--local", help="path to a local copy of GWCD4i.DE_stats.h5ad")
    ap.add_argument("--limit", type=int, default=None, help="debug: only first N rows")
    args = ap.parse_args()
    t0 = time.time()

    f = open_h5(args.local)
    obs = read_dataframe(f["obs"])
    var = read_dataframe(f["var"])
    n_obs, n_var = f["layers/log_fc"].shape
    print(f"Opened file: {n_obs} perturbations x {n_var} measured genes")

    # --- 1. Quality filter ---------------------------------------------------
    flags = obs[QUALITY_FLAGS].astype(bool)
    obs["drop"] = flags.any(axis=1)
    # bitmask of reasons: 1=neighboring_gene_KD, 2=distal_offtarget, 4=low_target_gex
    obs["reason"] = sum(flags[c].astype(int) * (2 ** i) for i, c in enumerate(QUALITY_FLAGS))
    print("\nQuality filter (removed = any of " + ", ".join(QUALITY_FLAGS) + "):")
    for cond in CONDITIONS:
        sub = obs[obs.culture_condition == cond]
        print(f"  {cond:9s} total={len(sub):6d}  removed={int(sub['drop'].sum()):5d}  "
              f"kept={int((~sub['drop']).sum()):5d}")

    # --- 2. Unified gene list ------------------------------------------------
    # Genes are identified by Ensembl ID; the website uses integer indices into
    # a single alphabetically sorted list of names to keep files compact.
    id2name = dict(zip(var.gene_ids, var.gene_name))
    for gid, name in zip(obs.target_contrast, obs.target_contrast_gene_name):
        id2name.setdefault(gid, name)
    ids_sorted = sorted(id2name, key=lambda g: (id2name[g].upper(), g))
    gene_index = {gid: i for i, gid in enumerate(ids_sorted)}
    gene_names = [id2name[g] for g in ids_sorted]
    var_to_global = np.array([gene_index[g] for g in var.gene_ids])
    var_pos = {gid: j for j, gid in enumerate(var.gene_ids)}  # Ensembl ID -> column

    # --- 3. Stream rows and extract significant targets ---------------------
    rows_needed = np.where(~obs["drop"].values)[0]
    if args.limit:
        rows_needed = rows_needed[: args.limit]
    keep_set = set(rows_needed.tolist())
    blocks = sorted({r // BLOCK_ROWS for r in rows_needed})
    results = {}  # obs row -> dict

    def process_block(b, fh):
        lo, hi = b * BLOCK_ROWS, min((b + 1) * BLOCK_ROWS, n_obs)
        lfc_blk = fh["layers/log_fc"][lo:hi]
        padj_blk = fh["layers/adj_p_value"][lo:hi]
        out = {}
        for k in range(hi - lo):
            i = lo + k
            if i not in keep_set:
                continue
            lfc, padj = lfc_blk[k], padj_blk[k]
            self_col = var_pos.get(obs.target_contrast.iat[i])
            # on-target knockdown strength (log2FC of the gene in its own row)
            kd_lfc = None
            if self_col is not None and np.isfinite(lfc[self_col]):
                kd_lfc = r3(lfc[self_col])
            # significant targets (NaN p-values compare False -> not significant)
            with np.errstate(invalid="ignore"):
                sig = (padj < PADJ_THRESHOLD) & np.isfinite(lfc)
            if self_col is not None:
                sig[self_col] = False  # exclude effect on itself
            cols = np.where(sig)[0]
            n_sig = len(cols)
            n_up = int((lfc[cols] > 0).sum())
            # top-N by absolute effect size
            if n_sig > TOP_N_TARGETS:
                cols = cols[np.argsort(-np.abs(lfc[cols]))[:TOP_N_TARGETS]]
            else:
                cols = cols[np.argsort(-np.abs(lfc[cols]))]
            out[i] = dict(kd=kd_lfc, n_sig=n_sig, n_up=n_up,
                          t=[[int(var_to_global[c]), r3(lfc[c]), sig2(padj[c])] for c in cols])
        return out

    print(f"\nStreaming {len(blocks)} blocks of {BLOCK_ROWS} rows with {N_WORKERS} workers...")
    # one file handle per worker thread (h5py handles are not thread safe)
    import threading
    local = threading.local()

    def worker(b):
        if not hasattr(local, "fh"):
            local.fh = open_h5(args.local)
        return process_block(b, local.fh)

    done = 0
    with ThreadPoolExecutor(N_WORKERS) as pool:
        for res in pool.map(worker, blocks):
            results.update(res)
            done += 1
            if done % 10 == 0 or done == len(blocks):
                el = time.time() - t0
                print(f"  {done}/{len(blocks)} blocks  ({el:.0f}s elapsed)", flush=True)

    # --- 4. Write JSON shards -----------------------------------------------
    if os.path.isdir(OUT_DIR):
        shutil.rmtree(OUT_DIR)
    os.makedirs(OUT_DIR)
    index = {
        "source": "Zhu, Dann et al., Cell 2026 - GWCD4i.DE_stats.h5ad",
        "params": {"padj_threshold": PADJ_THRESHOLD, "top_n_targets": TOP_N_TARGETS,
                   "quality_flags_removed": QUALITY_FLAGS},
        "conditions": CONDITIONS,
        "genes": gene_names,
        # per condition: list of perturbations kept, as compact arrays:
        # [gene_idx, n_sig_targets, n_up, kd_log2fc|null, kd_zscore,
        #  ontarget_significant(0/1), n_cells, shard]
        "pert_fields": ["gene", "n_sig", "n_up", "kd_lfc", "kd_z", "kd_sig", "n_cells", "shard"],
        "pert": {},
        # per condition: perturbations removed by the quality filter
        # [gene_idx, reason_bitmask]  (1=neighboring_gene_KD, 2=distal_offtarget, 4=low_target_gex)
        "removed": {},
    }
    summary = []
    for cond in CONDITIONS:
        cond_dir = os.path.join(OUT_DIR, cond)
        os.makedirs(cond_dir)
        sub = obs[obs.culture_condition == cond]
        shards, perts, n_edges = {}, [], 0
        for i in sub.index:
            gi = gene_index[obs.target_contrast.iat[i]]
            if obs["drop"].iat[i]:
                index["removed"].setdefault(cond, []).append([gi, int(obs.reason.iat[i])])
                continue
            if i not in results:
                continue  # only possible with --limit
            r = results[i]
            key = shard_key(gene_names[gi])
            if r["t"]:
                shards.setdefault(key, {})[str(gi)] = r["t"]
                n_edges += len(r["t"])
            z = obs.ontarget_effect_size.iat[i]
            nc = obs.n_cells_target.iat[i]
            perts.append([gi, r["n_sig"], r["n_up"], r["kd"],
                          r3(z) if np.isfinite(z) else None,
                          int(bool(obs.ontarget_significant.iat[i])),
                          int(nc) if np.isfinite(nc) else None,
                          key if r["t"] else None])
        index["pert"][cond] = perts
        size = 0
        for key, content in shards.items():
            p = os.path.join(cond_dir, f"{key}.json")
            with open(p, "w") as fh:
                json.dump(content, fh, separators=(",", ":"))
            size += os.path.getsize(p)
        biggest = max((os.path.getsize(os.path.join(cond_dir, f"{k}.json")), k) for k in shards)
        summary.append((cond, len(perts), sum(1 for p in perts if p[1] > 0), n_edges,
                        len(shards), size, biggest))

    with open(os.path.join(OUT_DIR, "index.json"), "w") as fh:
        json.dump(index, fh, separators=(",", ":"))

    # --- 5. Summary ----------------------------------------------------------
    print("\n================ SUMMARY ================")
    print(f"Total genes in index: {len(gene_names)}  "
          f"(distinct perturbed genes kept: "
          f"{len({p[0] for c in CONDITIONS for p in index['pert'][c]})})")
    for cond, n_p, n_p_t, n_e, n_sh, size, (bsz, bk) in summary:
        print(f"{cond:9s} perturbations kept={n_p:5d}  with targets={n_p_t:5d}  "
              f"edges={n_e:7d}  shards={n_sh:2d}  total={size/1e6:6.1f} MB  "
              f"largest shard={bk}.json {bsz/1e6:.1f} MB")
    print(f"index.json: {os.path.getsize(os.path.join(OUT_DIR, 'index.json'))/1e6:.2f} MB")
    print(f"Done in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
