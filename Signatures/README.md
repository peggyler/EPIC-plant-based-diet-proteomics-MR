# ProPlant: plant-based diet protein signatures

R analysis scripts by Peggy Ler for hPDI and uPDI protein discovery in the EPIC
subcohort. This repository currently contains LIMMA, elastic-net and XGBoost
signature analyses. No participant data or fitted outputs are included.

## Files

- `signatures/limma`: full cohort LIMMA scripts.
- `signatures/enet`: subcohort elastic-net scripts.
- `signatures/xgb/hpdi/`: subcohort hPDI XGBoost scripts.
- `signatures/xgb/updi/`: subcohort uPDI XGBoost scripts.

The XGB folders use the `xgboost` R package. Keep each `00_data_prep*.R` file in
its corresponding XGB folder: scripts 01–05 source that helper automatically.
The package needs no separate configuration folder or launcher.

## Before running

1. Set the R working directory to the repository root containing this README
   and the `signatures` folder. Start each script in a fresh R session.
2. Edit `DATA_FILE` and `OUTPUT_ROOT` near the top of each script you run.
   Use the same input file and output root throughout a workflow. Forward
   slashes work in Windows paths. Defaults are `data/subcohort_20251106.RData`
   and `outputs`; no such data file is provided.
3. The RData file must contain a data frame named `subcohort`, including the
   appropriate numeric diet score (`Hpdi_NoVegoils` or `Updi_NoVegoils`), numeric
   protein columns beginning with `seq`, and the required covariates below.
4. Install the required packages, preferably matching the original analysis
   environment: `limma`, `statmod`, `glmnet`, `dplyr`, `tibble`, `readr`,
   `data.table`, `Matrix`, `xgboost`, `caret` and R's bundled `splines`.
   Original R/package versions were not supplied with these scripts.

Covariates are `Sex`, `Age_Recr`, `L_School`, `alc_freq`, `smoke`, `Center` and
physical activity. LIMMA/elastic net use `PA_total_MET_quint`; XGB uses
`PA_total_MET`. Smoking and alcohol variables must include the reference value
`Never`. Education is recoded from `L_School` (0/1; 2/3/4; 5/missing).
`Idepic` is an optional ID column used in elastic-net prediction exports.

For example, open a script in RStudio and use Source, or run from the repository root:

```r
source("signatures/Limma.R")
```

Output folders are created automatically. Keep your local data, generated
outputs and private file paths out of the public repository. Outputs include
participant-level predictions, prep bundles and SHAP matrices.

## Running order

### Full-subcohort LIMMA

Run `signatures/Limma.R` for hPDI and `signatures/Limma_uPDI.R` for uPDI.
These are independent association analyses and use linear age adjustment.

### Elastic net

Run 1.1 → 1.2 → 1.3 separately for each outcome:

| Step | hPDI file in `signatures/` | uPDI file in `signatures/` |
|---|---|---|
| Split; training-only LIMMA filter | `1.1_Preprocess_ENET_hPDI_with_filter.R` | `1.1_Preprocess_ENET_uPDI_with_filter.R` |
| Tune; fit; bootstrap | `1.2_Train_Lasso_hPDI.R` | `1.2_Train_Lasso_uPDI.R` |
| Test-set evaluation | `1.3_Test_Lasso_hPDI.R` | `1.3_Test_Lasso_uPDI.R` |

Historical filenames contain `Lasso`, but the main models search elastic-net
alpha values 0.1, 0.3, 0.5, 0.7 and 0.9. Outputs are saved under `hpdi/` and
`updi/` within your chosen output root. Training uses 1,000 bootstrap replicates.
The training/test scripts require exactly one matching prep/model RDS in the
outcome folder. Use a separate output root for a separate analysis run to avoid
overwriting files or mixing model versions.

`1.2_Train_Lasso_uPDI_alpha1.R` is an optional wider alpha grid (0.1–1.0 in
steps of 0.1). It is not part of the default sequence and does not fix alpha at
1. It reads the main uPDI prep file and writes separately under
`updi/alpha_grid_01_to_10/`. The supplied 1.3 script evaluates the main model.

### XGBoost / XGB

Within each outcome's GBM folder, run the numbered scripts in order:

| Prefix | Step |
|---|---|
| `01_` | Tune tree structure; create and save the split/CV folds |
| `02_` | Tune regularisation |
| `03_` | Tune column sampling |
| `04_` | Fit the final model and evaluate the test set |
| `05_` | Compute SHAP values and summaries |

The `00_` files are sourced helpers and should not be run independently.
Stage A creates a timestamped run directory under `<outcome>/xgb/` and saves
its location in `current_run_dir.txt`. Subsequent stages use that same run.
Re-running Stage A selects a new run for subsequent stages. Keep the same data,
row order and output root; if you move saved outputs, update the pointer file.
The grids refine earlier broad searches that are not included in this upload.
The supplied setting of 16 XGBoost threads is retained.

## Reproducibility notes

- The duplicate post-filter scaling of diet score and age was removed from
  both 1.3 test scripts. Earlier scaling remains, matching 1.1. Peggy reported
  unchanged results after rerunning 1.3; numerical outputs were not provided
  for independent verification. Existing fitted models can be used for this check.
- The elastic-net split uses seed 258 and the bootstrap uses seed 42. The
  fitting scripts do not fix a CV seed or use a shared `foldid` across alphas;
  rerunning training can therefore choose different tuning parameters.
- Elastic-net test R2 is squared Pearson correlation; XGB test R2 is
  `1 - SSE/SST`. They are different definitions and need explicit labels when
  compared. Both calculations remain as supplied.
- Predictive preprocessing scales age and the outcome before splitting.
  Protein preprocessing is assumed to have occurred upstream; ENET uses
  `standardize = FALSE`, while XGB replaces missing protein values with zero.
  Confirm these choices and the differing physical-activity covariates against
  the manuscript. ENET's protein filter is fitted on the full training split,
  outside the subsequent internal CV folds.
- The SHAP signed-weight calculation uses the sign of the mean SHAP
  contribution. This alone does not establish the direction of the
  protein-abundance relationship. Its numerical calculation is retained.
- In XGB tuning, the column `proteins_used` counts all used features, including
  covariates. Stage 04 separately extracts protein-only features. Historical
  commented results are annotations, not outputs verified during cleanup.
- The five alternative hPDI scripts originally nested inside `files.zip` are
  omitted from this simpler upload. The surrounding hPDI scripts are retained;
  confirm that these are the manuscript-producing versions before final release.

The scripts have been checked with an R syntax grammar parser. A native R
runtime and study data were unavailable during cleanup, so the full workflow
has not been executed or numerically validated here. Record `sessionInfo()`
when running in the original analysis environment.
