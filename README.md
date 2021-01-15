Data and code from the paper "Estimating the effects of non-pharmaceutical interventions on the number of new infections with COVID-19 during the first epidemic wave" by Banholzer et al.

*Structure of the repository*:
- **A preprocessed data file is in the `data/data_preprocessed.csv`.**
- **The Stan model code is in `models/sim-g.stan`.**
- The default model can be run with `run_base.r`.
- All models for the sensitivity analysis can be run with `robustness-checsk\run_rc.r`.
- After running the model(s), the results presented in the paper can be reproduced with `analysis.rmd`. 
- A comparison with a similar study is included and can be run and reproduced with `robustness-checks/comparison_with_brauner.rmd`. 

*Run and analysis the default model:*
Open the R-Project `npi_effectiveness_first_wave.Rproj` and ensure that all required libraries and their dependencies are installed. After that, run `run_base.r` to fit the default model. It should take approximately one hour to run the default model with the specified settings. After the sampling is finished, the sampled quantities from the default model are stored in `fitted-models/base.rds`. Open `analysis.rmd` to load the sampled quantities and run the code up to the sensitivity checks in order to show the results from the default model as presented in the paper.
