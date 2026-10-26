# Two-Compartment Vancomycin Pharmacokinetics

A reproducible research pipeline for synthetic infusion histories and individualized exposure modeling.

## Methods

- Apply correction chains to administration and serum-level records, using event-effective and specimen-collection times.
- Reconstruct time-varying drug delivery and serum creatinine.
- Propagate a two-compartment amount model using exact matrix exponentials.
- Fit clearance and distribution parameters using constrained penalized nonlinear least squares.
- Integrate 24–48 h exposure; project seven candidate regimens over 48–72 h, then rank them using the protocol's therapeutic-range rules.

## Run

```bash
python -m pip install -r requirements.txt
python full_pipeline.py --input data --output results
python full_pipeline.py --input data_full --output results_full
python -m unittest discover -s tests -v
```

`data/` is a small example; `data_full/` contains 240 synthetic courses.

## Validation status

All 240 courses are processed. Current code matches the separate reference for course completeness, renal-adjusted clearance, central volume, observed 24–48 h exposure, and QC totals. **The 48–72 h projections are not yet reference-equivalent** (1,466 values outside the strict 0.6 mg·h/L comparison tolerance), and 128 rankings differ. The software is a research prototype, not a verified dosing tool.

**No clinical use.** This model is simplified, and must not guide actual medication decisions.
