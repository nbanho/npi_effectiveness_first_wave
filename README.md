Data and code from the paper "Estimating the effects of non-pharmaceutical interventions on the number of new infections with COVID-19 during the first epidemic wave" by Banholzer et al.

*Structure of the repository*:
- A preprocessed data file is in the `data/data_preprocessed.csv`.
- The Stan model code is in `models/sim-g.stan`.
- The default model can be run with `run_base.r`.
- All models for the sensitivity analysis can be run with `robustness-checsk\run_rc.r`.
- After running the model(s), the results presented in the paper can be reproduced with `analysis.Rmd`. 
- A comparison with a similar study is included and can be run and reproduced with `robustness-checks/comparison_with_brauner.rmd`. 

