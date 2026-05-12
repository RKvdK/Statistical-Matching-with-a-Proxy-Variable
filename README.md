# Statistical Matching with a Proxy Variable

## Project description

This repository contains the research archive accompanying the MSc thesis "Statistical Matching with a Proxy Variable: The Role of Sample Size, Unit Overlap, Selectivity, Proxy Quality and the Validity of the Conditional Independence Assumption." for the Research Master's programme Methodology and Statistics for the Biomedical, Behavioural and Social Sciences. This project examines statistical matching using doubly robust estimation methods and simulation-based analyses.

## Repository structure

### Code

The `Code/` folder contains the following scripts, which should be run in the indicated numerical order:

- `01_packages.R`: Install and load the required packages.
- `02_utils.R`: Helper functions. 
- `03_parametergrid.R`: Definition of the simulation parameter grid.
- `04_selectivity.R`: Operationalization of the selectivity mechanisms.
- `05_generate_population.R`: Function to generate simulation populations.
- `06_monte_carlo.R`: Function to perform Monte Carlo simulation. 
- `07_run_simulation.R`: Function to parallelise the full grid simulation.
- `08_fullgrid_results.R`: Performane of the full grid simulation.
- `09_homogeneity_within_joints.R`: Homogeneity evaluation within joint distributions.
- `10_biasvariance_convergence.R`: Analysis of the bias-variance estimate convergence.
- `11_fullgrid_visualizations.R`: Visualization for the full-grid results.
- `12_scenarios.R`: Definition of the simulation scenarios.
- `13_scenario_results.R`: Calculation and visualization of the scenario-based analysis.

Additionally, this folder contains `SMPV.Rproj`, which is the corresponding Rstudio project file.

### Data

The `Data/` folder contains the following objects:

- `MC50.rds`: RDS-file containing the full grid simulation results.
- `MC50.txt`: Text file containing all parallized simulation logs.

### Manuscript

The `Manuscript/` folder contains the Quarto source file `Thesis.qmd` used to render the thesis manuscript, as well as the corresponding PDF version: `Thesis.pdf`. 
Additionally, the `Manuscript/` folder contains the `Manuscript/References` folder. This folder contains:

- `references.bib`: Bibliography database to manage the manuscript citations.
- `apa.csl`: Citation style language file.

Finally, the `Manuscript/` folder contains the `Manuscript/Figures` folder. This folder contains all .png-files that are used in the PDF manuscript:

- `METFIG.png`: Overview of the performed simulation study.
- `patdeg.png`: `Overview of the implemented selectivity patterns and degrees under MAR (generated in `04_selectivity.R`).
- `convplot.png`: Convergence plot for the bias-variance estimates (generated in `10_biasvariance_convergence.R`).
- `jointcell.png`: Evaluation of the homogeneity within the joint distributions (generated in `09_homogeneity_within_joints.R`).
- `sensmar.png`: Evaluation of the sensitivity of the estimates with respect to the MAR patterns and degrees (generated in `11_fullgrid_visualizations.R`).
- `sensmnar.png`: Evaluation of the sensitivity of the estimates with respect to the MNAR main and interaction effects scenarios (generated in `11_fullgrid_visualizations.R`).
- `mardiff.png`: Full-grid driver plot under MAR (generated in `11_fullgrid_visualizations.R`).
- `marstep.png`: More detailed driver plot regarding overlap and external sample size under MAR (generated in `11_fullgrid_visualizations.R`).
- `mnardiff.png`: Full-grid driver plot under MNAR (generated in `11_fullgrid_visualizations.R`).
- `mnarstep.png`: More detailed driver plot regarding overlap and external sample size under MNAR (generated in `11_fullgrid_visualizations.R`).
- `marscenres.png`: Scenario-based results under MAR (generated in `13_scenario_results.R`).
- `mnarscenres.png`: Scenario-based results under MNAR (generated in `13_scenario_results.R`).

### Output

The `Output/` folder contains the following objects:

- `A1_benchmark`: Monte Carlo convergence results for scenario A1 (generated in `10_biasvariance_convergence.R`).
- `D7_benchmark`: Monte Carlo convergence results for scenario D7 (generated in `10_biasvariance_convergence.R`).
- `H16_benchmark`: Monte Carlo convergence results for scenario H16 (generated in `10_biasvariance_convergence.R`).
- `Scenario Calculation Results`: Raw Monte Carlo simulation results for all scenarios (generated in `13_scenario_results.R`).
- `Final Estimates`: Aggregated scenario-level performance estimates for all three estimators (generated in `13_scenario_results.R`).

### Root files



## Software requirements

This project was developed using R (version 4.4.3), Quarto and renv. The latter was used to manage all package dependencies.

## Reproducibility



## Data



## Ethics and privacy

Ethical consent for this study was granted (FETC File Number: 25-1983) by the Ethics Review Board of the Faculty of Social and Behavioural Sciences of Utrecht University. Moreover, this project does not contain personally identifiable information.

## Contact

Author: Robert van der Kaap
Research Master's programme: Methodology and Statistics for the Biomedical, Behavioural and Social Sciences
Utrecht University
