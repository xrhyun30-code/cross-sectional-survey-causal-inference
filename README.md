# Cross-sectional survey causal-inference analyses

R code for the simulation study and KNHANES 2020 analysis in “What do causal
inference methods change in a cross-sectional survey? A simulation study and
KNHANES 2020 example.”

## Contents

- `R/simulation_baseline.R`: baseline simulation scenarios S1, S2, S4, and S5.
- `R/run_additional_analyses.R`: revised scenarios S3a, S3b, and S6, plus the
  unweighted propensity-score, multiple-imputation, and S1 truth-diagnostic
  analyses.
- `R/empirical_knhanes2020.R`: primary KNHANES analysis and weight diagnostics.
- `R/sensitivity_knhanes2020.R`: alternative exposure/outcome definitions and
  age-stratified analyses.

## Data and privacy

No KNHANES data or analysis outputs are included. Obtain the KNHANES 2020
public-use file from the Korea Disease Control and Prevention Agency and follow
its registration and use conditions. Do not commit the SAS data file or
participant-level data to this project.

By default, the scripts look for `data/hn20_all.sas7bdat`. Set
`KNHANES2020_SAS` to another local path if needed. The data directory and
generated output directories are excluded by `.gitignore`.

## Requirements

The analyses were run with R 4.2.2. Install the packages used by the relevant
script:

- `haven` and `survey` for the KNHANES analyses;
- `mice` for multiple imputation in `run_additional_analyses.R`;
- base R's `parallel` package for the simulation.

## Running the analyses

Run commands from the project directory.

```sh
Rscript R/simulation_baseline.R 1000 output/simulation_baseline.csv
Rscript R/run_additional_analyses.R --sim-reps 1000 --scenarios S3a,S3b,S6 --skip-knhanes
KNHANES2020_SAS=data/hn20_all.sas7bdat Rscript R/empirical_knhanes2020.R
KNHANES2020_SAS=data/hn20_all.sas7bdat Rscript R/sensitivity_knhanes2020.R
KNHANES2020_SAS=data/hn20_all.sas7bdat Rscript R/run_additional_analyses.R --skip-simulations
```

The simulation and KNHANES analyses can take several minutes or longer,
depending on the number of replicates, imputations, and available CPU cores.
The additional-analysis runner writes its output under `R/analysis_outputs`;
the other KNHANES scripts write under `output/`.

## License

No license has been selected. All rights are reserved by the authors; public
availability alone does not grant permission to reuse or redistribute the code.
