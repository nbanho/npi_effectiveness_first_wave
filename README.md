This repository provides data and code for the paper on ``Estimating the effects of non-pharmaceutical interventions on the number of new infections with COVID-19 during the first epidemic wave''.

- A preprocessed data file is in the `data` folder.
- The model code can be found in `models/sim-g.stan`.
- The default model can be run with `run_base.r`.
- All models for the sensitivity analysis are run with `robustness-checsk\run_rc.r`.
- After running the model(s), the results can be analyzed with `analysis.Rmd`. 
- A comparison with a similar study is included and can be run and analyzed with `robustness-checks/comparison_with_brauner.rmd`. 

