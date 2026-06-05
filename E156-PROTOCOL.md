# E156-PROTOCOL — 786-M21 (multi-method meta-analysis workbench)

- **Project:** Omniinfoloss1 (GitHub repo `Omniinfoloss1`, user `mahmood726-cyber`)
- **Revived:** 2026-06-05 (from a single-file `omniinfo1.html` dump, title "786-M21 Robust Visuals")
- **Type:** single-file offline browser tool + Node-testable engine
- **Dashboard:** GitHub Pages (`index.html`)

## What changed in the revival

- Made **fully offline**: vendored Plotly 2.24.1 into `plotly.min.js` (local
  `<script src>`); removed the Google Fonts CDN `<link>`. No external `http(s)`
  reference remains.
- Extracted the pure statistical core verbatim into `engine.js` (single source
  of truth); deleted the inline duplicates and pointed the page at `engine.js`.
  `Pairwise.pool`/`egger` now delegate to the pure `Pooling.*`.
- **Fixed a correctness bug**: the pooler returned `NaN` for a single study
  (`tau2=(Q−df)/C=0/0`); it now degrades to the fixed-effect result (`tau2=0`)
  when `df<1` or `C≤0`, with the Knapp–Hartung `t`-inflation skipped at `df<1`.
- Added `tests.js` (48 assertions, all passing) with an independently
  hand-derived two-study DL example and edge cases.
- Added Pages scaffold (`.nojekyll`, `.gitignore`); renamed
  `omniinfo1.html` → `index.html`.

## Body (E156 draft — CURRENT BODY)

Can one offline browser tool carry a full evidence-synthesis workflow without
trusting a remote chart server or a buried, untested statistical core? This
revival takes a single-file dashboard covering pairwise, network, component,
diagnostic-accuracy and dose-response meta-analysis, plus an information-recovery
panel that adds median/IQR reconstruction and baseline-risk regression. The pure
mathematics, log-scale DerSimonian–Laird and Knapp–Hartung pooling with τ², I²,
confidence and prediction intervals, were extracted verbatim into a Node-testable
engine and locked behind forty-eight hand-derived assertions. An audit found that
single-study pooling returned NaN because the heterogeneity term divided zero by
zero, so it now degrades correctly to the fixed-effect estimate. Plotly was
vendored locally and the font CDN removed, leaving zero external references. The
result is a reproducible, offline workbench whose numbers match independent hand
computation rather than an opaque inline script. It remains an exploratory
synthesis aid, not a validated clinical or inferential instrument.

SUBMITTED: [ ]
