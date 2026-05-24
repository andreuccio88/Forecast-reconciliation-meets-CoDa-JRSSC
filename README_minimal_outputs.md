# JRSSC replication code - minimal outputs

This version keeps the computational workflow unchanged, but prunes the final `outputs/` folder so that only the CSV and LaTeX files needed to reproduce the manuscript tables/results remain.

After running:

```r
source("run_all.R")
```

the final files are in `outputs/tables/`:

- `RMSE_all_structures.csv` and `point_RMSE_table.tex` for the main point-forecast RMSE table.
- `CRPS_all_structures.csv` and `interval_CRPS_table.tex` for the model-based CoDA bootstrap CRPS table.
- `BB_CRPS_all_structures.csv` and `BB_interval_CRPS_table.tex` for the residual block-bootstrap robustness CRPS table.
- `MAE_all_structures.csv` and `point_MAE_table.tex` for the appendix MAE table.

Intermediate CSV diagnostics, per-structure duplicate CSVs, figures, and RDS/RData cache files are removed at the end of the run.

To restore figure exports, set `CREATE_FIGURES <- TRUE` in `run_all.R`. To keep all intermediate outputs, set `PRUNE_TO_MANUSCRIPT_TABLES <- FALSE` in `run_all.R`.
