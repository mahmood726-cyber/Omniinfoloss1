# 786-M21 — Multi-method meta-analysis workbench (with information-recovery panel)

A single-file, **fully offline** dashboard for evidence synthesis. It covers
pairwise meta-analysis (fixed / DerSimonian–Laird / Knapp–Hartung, plus
multilevel and robust-variance options), network meta-analysis (fixed / DL /
REML / multilevel), component NMA, diagnostic-test-accuracy (bivariate SROC),
dose-response, GOSH, a RoBMA-style ensemble, a small Bayesian MCMC sampler,
RoB-2 and GRADE panels, and an **"information-recovery" panel** that combines a
traditional pooled estimate with median/IQR distribution reconstruction,
quantile treatment effects and baseline-risk meta-regression.

**Live app:** open `index.html` (or the GitHub Pages link). No build step, no
network, no external CDN.

## Layout

```
index.html   single-file UI (loads engine.js + plotly.min.js)
engine.js    pure statistical core — runs in Node and the browser
tests.js     Node test harness, 48 assertions
plotly.min.js  vendored Plotly 2.24.1 (offline)
LICENSE      Apache-2.0
```

## Statistical core (`engine.js`)

All pure (no DOM / Plotly) and exported for Node:

| Object | What it does |
|---|---|
| `Stat` | `pnorm` (normal CDF `½(1+erf(x/√2))`), `qnorm` (inverse), `qt` (t-quantile approx), `invLogit` |
| `Matrix` | dense linear algebra (zeros / transpose / dot / Gauss–Jordan inverse) |
| `Optim.nelderMead` | derivative-free simplex minimiser |
| `Pooling` | log-scale pairwise pooler: τ² via DL `(Q−(k−1))/C`, I² `(Q−df)/Q`, 95% CI, prediction interval, Knapp–Hartung variance inflation, and Egger's small-study test |
| `Multilevel` / `RVE` | 3-level random-effects fit; cluster-robust (sandwich) variance with `t_{m−1}` |
| `NMA_RE` / `NMA_Multilevel` | network meta-analysis via GLS on the contrast design (fixed / DL / REML / multilevel) |
| `InfoMetrics` | traditional-vs-recovered information accounting (the "information loss" view) |
| `DistributionRecovery` | median/IQR → mean/SD reconstruction |
| `QTE` / `BaselineRisk` | quantile treatment effects; effect-on-logit(control-risk) meta-regression |
| `EnhancedMA.run` | combines the traditional pooled estimate with the recovered-distribution components |

## Fixes applied during revival (2026-06-05)

- **Made fully offline.** Vendored Plotly 2.24.1 into `plotly.min.js` and pointed
  the script tag at the local file; removed the Google Fonts `<link>` (system
  fonts fall back). `grep` confirms **zero** external `http(s)` references remain.
- **Single source of truth.** Extracted the pure statistical objects verbatim
  into `engine.js`; the inline duplicates were deleted and the page now loads
  `engine.js`. The DOM-coupled `Pairwise.pool`/`Pairwise.egger` are now thin
  delegates to the pure `Pooling.*`.
- **Fixed an unambiguous correctness bug (k=1 → NaN).** The original inline
  pooler computed `tau2 = (Q−df)/C` with no guard; for a single study `df=0`,
  `C=0` and `Q=0`, so `tau2 = 0/0 = NaN`, which poisoned the estimate and CI.
  The pure `Pooling.pool` now degrades to the fixed-effect result (`tau2=0`) when
  `df<1` or `C≤0`, and skips the Knapp–Hartung `t`-inflation when `df<1`. Locked
  in by tests.
- Added `tests.js` (48 assertions, all passing) with an independently
  hand-derived two-study DL example and edge cases.
- Added Pages scaffold (`.nojekyll`, `.gitignore`); renamed
  `omniinfo1.html` → `index.html`.

## Tests

```
node tests.js
# 48 passed, 0 failed
```

Checks include normal-CDF reference points (Φ(0)=0.5, Φ(1.96)≈0.975), a
hand-worked two-study DL case (τ²=0.25, pooled=0.5, I²=50%, PI=0.5±1.96√0.5),
the k=1 no-NaN degradation, two-identical-studies (Q=0 ⇒ τ²=0, I²=0), an empty
guard, a matrix `A·A⁻¹=I` check, and information-recovery property checks
(no-extra-info ⇒ gain=1; distributional info strictly increases recovered %).

## Caveats

DerSimonian–Laird under-estimates τ² for small *k* (REML/Paule–Mandel are
preferred for k<10); the dashboard preserves the original methods for continuity
and always reports τ² and I² alongside each estimate. The "information-recovery"
panel is an exploratory heuristic for visualising how much extra signal
median/IQR and baseline-risk inputs add to a pooled estimate — it is not a
validated inferential procedure, and the QTE simulation is stochastic. Treat all
outputs as hypothesis-generating. Apache-2.0 licensed.
