# Virtual Knockdown Explorer

An interactive, static website: pick a gene, "turn it down," and see what happens in human CD4+ T cells.

- **Measured effects**: which genes went up or down when that gene was knocked down in a real experiment.
- **Predicted chain reaction**: a simple model of how those changes might ripple further through the gene network.
- **Model check**: how well the chain-reaction prediction matches the real measurements, for one gene and across all genes.
- **Compare conditions**: the same knockdown in resting vs. activated T cells, side by side.

It is plain HTML/CSS/JavaScript with no server, so it can be hosted for free on GitHub Pages.

## The data

Zhu, Dann *et al.* "Genome-scale perturb-seq in primary human CD4+ T cells maps context-specific regulators of T cell programs and human immune traits." *Cell* (2026).
The authors used CRISPR interference (CRISPRi) to knock down every expressed gene in primary human CD4+ T cells, in three conditions: **Rest**, **Stim8hr**, and **Stim48hr** (8 or 48 hours after activation). They then measured how every other gene responded (Perturb-seq).

- Dataset page: <https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq>
- File used: `GWCD4i.DE_stats.h5ad` (17 GB), from the public bucket `s3://genome-scale-tcell-perturb-seq/marson2025_data/`
- Authors' code: <https://github.com/emdann/GWT_perturbseq_analysis_2025>

The file holds 33,983 perturbations (one gene in one condition) × 10,282 measured genes, with layers `log_fc`, `zscore`, `p_value`, `adj_p_value`, `lfcSE`, and `baseMean`.

## How the data was processed (`preprocess.py`)

The 17 GB file is **never downloaded**. Its layers are stored uncompressed and contiguous, so `preprocess.py` streams blocks of rows directly from S3 with `h5py` + `s3fs`. It reads only `log_fc` and `adj_p_value`, about 5.6 GB of network traffic, and takes about 5 minutes.

For each condition separately:

1. **Quality filter.** It drops perturbations flagged `neighboring_gene_KD` (a nearby gene was also knocked down), `distal_offtarget_flag` (possible off-target guide), or `low_target_gex` (target expressed too low to judge the knockdown).

   | Condition | Total | Removed | Kept | Kept with ≥1 significant target |
   |---|---|---|---|---|
   | Rest | 11,287 | 3,647 | 7,640 | 4,760 |
   | Stim8hr | 11,415 | 3,677 | 7,738 | 4,979 |
   | Stim48hr | 11,281 | 3,342 | 7,939 | 5,162 |

2. **Significant targets.** It keeps target genes with `adj_p_value < 0.05` (`PADJ_THRESHOLD` at the top of the script) and excludes the gene's effect on itself. NaN p-values count as not significant.
3. **Size cap.** It keeps at most the top 200 targets per perturbation by |log2 fold change| (`TOP_N_TARGETS`). The full count of significant targets is still recorded.
4. **Knockdown strength.** It records the perturbed gene's own `log_fc` in its own row ("reduced to X% of normal"). If the gene isn't in the measured panel, it falls back to the on-target z-score `ontarget_effect_size`.
5. **Output.**
   - `data/index.json` (1 MB): every gene name, plus per-condition stats for each perturbation (number of targets, knockdown strength, cell count) and the reason for each removal.
   - `data/<Condition>/<Letter>.json`: measured targets for perturbed genes starting with that letter. Each shard is under 0.4 MB, and `data/` totals about 9 MB. The website loads only the shards it needs.

About 420,000 gene→gene edges are stored in total (132k / 145k / 144k for Rest / Stim8hr / Stim48hr).

Reproduce it:

```bash
pip install h5py s3fs numpy pandas
python preprocess.py           # writes ./data
node scripts/model_check.js    # writes data/model_check.json
```

## The chain-reaction model (`js/model.js`)

- **Step 1** uses the measured effects of knocking down gene A.
- **Step 2+**: if an affected gene B was itself knocked down in the same condition, B passes on its own measured effects, scaled by how much B changed relative to its own knockdown:

  ```
  scale(B)      = change(B) / knockdown(B)        (both log2 fold changes)
  change(C)    += scale(B) × measured effect of B on C
  ```

- Contributions from multiple paths are **added up** (a linear model).
- Safety limits:
  - |scale| is capped at 2.
  - Each gene's predicted total is capped at ±5 log2FC.
  - B passes effects on only if its own knockdown was real (log2FC ≤ −0.3).
- Genes never perturbed in the screen cannot pass anything on.
- Propagation stops at the chosen depth (1–3) or when changes fall below the minimum-effect cutoff.
- **By default, only genes that went DOWN pass effects on.** The screen measured what happens when a gene goes *down*. Assuming that a gene going *up* does the exact opposite turned out to be wrong in this data (see below). A checkbox turns it on for comparison.

**Model check.** Step 1 *is* the measurement, so the check compares only the **indirect (step 2+) prediction** with the measured effect, for genes that have both. `scripts/model_check.js` runs the same `model.js` across every knockdown with ≥20 targets. As a **baseline**, it also compares each prediction with a *random other* knockdown's measurements.

| Condition (depth 2, cutoff 0.1) | Down-only: median r | Random baseline | Up + down: median r |
|---|---|---|---|
| Rest | 0.46 (83% positive) | 0.41 | −0.14 |
| Stim8hr | 0.38 (81% positive) | 0.39 | −0.22 |
| Stim48hr | 0.47 (83% positive) | 0.40 | 0.02 |

What this shows:

- Effects passed through genes that went down match the measurements reasonably well.
- Passing effects through genes that went up gives *anti*-correlated predictions.
- Most of the down-only match is also reached by the random baseline. So the model mostly captures a **generic response shared by many knockdowns**, not regulation specific to the chosen gene.

## Limitations (please read)

- **Double counting:** measured effects already include indirect effects, so adding predicted chains on top can count the same effect twice.
- **Linear and additive:** real gene regulation saturates, has thresholds, and involves combinations of factors.
- **Incomplete network:** only significant targets (adj. p < 0.05), only the top 200 per gene, and only genes that were perturbed *and* passed QC can propagate effects. About 17% of perturbed genes aren't in the measured panel, so their knockdown strength can't be used for scaling.
- **"Up" ≠ opposite of "down":** CRISPRi only tells us what happens when a gene goes down.
- **Weak gene-specific signal:** the random-pairing baseline is nearly as good as the real prediction (see above).
- **Compare view:** "affected in only one condition" means "not among the stored significant targets of the other condition." For genes with more than 200 targets, a gene might be a weaker, unstored hit there.
- Educational tool, not for clinical or research conclusions.

## Running locally

Browsers block `fetch()` of local files, so use a tiny web server:

```bash
python3 -m http.server 8000
```

Then open <http://localhost:8000>.

Useful links: the URL updates as you explore (for example `#g=CD3E&c=Stim8hr&m=compare&c2=Rest`), so you can bookmark specific views for a demo.

## Files

| Path | What it is |
|---|---|
| `index.html`, `css/style.css` | Page layout and styles (light + dark mode) |
| `js/app.js` | User interface: search, network (Cytoscape.js 3.30.2 from jsDelivr), table, charts |
| `js/model.js` | Chain-reaction model + model check (shared by the website and Node) |
| `preprocess.py` | Streams the raw `.h5ad` from S3 and writes `data/` |
| `scripts/model_check.js` | Runs the model across all genes, writes `data/model_check.json` |
| `data/` | Small JSON files used by the website (~9 MB total) |
| `presentation/` | Model-check summary chart (HTML + PNG, light and dark) for slides |

The raw `.h5ad` file is never committed (see `.gitignore`).
