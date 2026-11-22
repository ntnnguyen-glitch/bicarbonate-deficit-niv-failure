# Bicarbonate deficit and NIV failure in acute hypercapnic acidotic respiratory failure

R code for the cohort study "Bicarbonate below the level expected for PaCO2 and failure of noninvasive ventilation in acute hypercapnic acidotic respiratory failure: a two-database cohort study".

## Data

The code uses MIMIC-IV v3.1 (https://doi.org/10.13026/kpb9-mt58) and eICU-CRD v2.0 (https://doi.org/10.13026/C2WM1R). Both are available on PhysioNet to credentialed users who have completed the required training and signed the data use agreements. No patient-level data are included in this repository.

## Requirements

R 4.6.1 with the packages data.table, DBI, duckdb, survival, cmprsk, splines, WeightIt, MatchIt, sandwich, mice, pROC and ggplot2.

## Run order

Set the paths to the local PhysioNet files, create a `fig/` folder in the working directory, then run the scripts from the `R/` folder in this order:

```sh
export MIMIC_DIR=/path/to/mimiciv/3.1
export EICU_DIR=/path/to/eicu-crd/2.0
export SOURCE_DB=hco3_source.duckdb
Rscript 00_extract_source_tables.R   # source tables from the PhysioNet CSV files
Rscript s1.R    # cohort selection, exposure, outcomes and covariates -> D.rds
Rscript s3.R    # primary, secondary, confounding and sensitivity analyses -> A.rds, R.rds
Rscript s5.R    # COPD strata, alternative slope, early failures, spline, ROC -> X.rds
Rscript s5b.R   # multiple imputation (continuous exposure, m = 10), hospital counts -> X.rds
Rscript s9.R    # follow-up checks, pre-NIV blood gas, DNI timing, DeltaHCO3 with pH -> O.rds
Rscript s10.R   # Model 3 with multiple imputation (m = 50), complete-case analysis, E-value -> MI.rds
Rscript s6.R    # Table 1, Table S3, Figures 1-4 and Figure S1
```

`s2.R` holds helper functions read by `s3.R` and `s9.R`.

## Notes

All pre-specified analyses are run and printed regardless of their P values. Some printed table labels are in Vietnamese; they do not affect the results.

## License

MIT
