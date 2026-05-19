# Statistical Matching with a Proxy Variable

## Project description

This repository contains the research archive accompanying the MSc thesis "Statistical Matching with a Proxy Variable: The Role of Sample Size, Unit Overlap, Selectivity, Proxy Quality and the Validity of the Conditional Independence Assumption." for the Research Master's programme Methodology and Statistics for the Biomedical, Behavioural and Social Sciences. This project examines statistical matching using doubly robust estimation methods and simulation-based analyses. This archive is intended to facilitate transparency and reproducibility of the simulation study and thesis results.

## Repository structure

### Root files

The root folder contains the following files on the project level:

- `README.md`: Overview of the research archive and reproduction instructions.
- `LICENSE`: License file for this research repository.
- `.gitignore`: Specification of files and folders that should not be tracked by Git.
- `renv.lock`: Lockfile containing the exact R package versions used in this project.
- `renv/`: Folder containing the files required to activate the project-specific `renv` environment.

### Code

The `Code/` folder contains the following scripts, which should be run in the indicated numerical order:

- `01_packages.R`: Install and load the required packages.
- `02_utils.R`: Helper functions. 
- `03_parametergrid.R`: Definition of the simulation parameter grid.
- `04_selectivity.R`: Operationalization of the selectivity mechanisms.
- `05_generate_population.R`: Function to generate simulation populations.
- `06_monte_carlo.R`: Function to perform Monte Carlo simulation. 
- `07_run_simulation.R`: Function to parallelise the full grid simulation.
- `08_fullgrid_results.R`: Performance of the full grid simulation.
- `09_homogeneity_within_joints.R`: Homogeneity evaluation within joint distributions.
- `10_fullgrid_visualizations.R`: Visualization for the full-grid results.
- `11_scenarios.R`: Definition of the simulation scenarios.
- `12_biasvariance_convergence.R`: Analysis of the bias-variance estimate convergence.
- `13_scenario_results.R`: Calculation and visualization of the scenario-based analysis.

Additionally, this folder contains `SMPV.Rproj`, which is the corresponding RStudio project file.

### Data

The `Data/` folder contains the following objects:

- `MC50.rds`: RDS-file containing the full grid simulation results.
- `MC50.txt`: Text file containing all parallelized simulation logs.

### Manuscript

The `Manuscript/` folder contains the Quarto source file `Thesis.qmd` used to render the thesis manuscript, as well as the corresponding PDF version: `Thesis.pdf`. 
Additionally, the `Manuscript/` folder contains the `Manuscript/References` folder. This folder contains:

- `references.bib`: Bibliography database to manage the manuscript citations.
- `apa.csl`: Citation style language file.

Finally, the `Manuscript/` folder contains the `Manuscript/Figures` folder. This folder contains all .png-files that are used in the PDF manuscript:

- `METFIG.png`: Overview of the performed simulation study.
- `patdeg.png`: Overview of the implemented selectivity patterns and degrees under MAR (generated in `04_selectivity.R`).
- `convplot.png`: Convergence plot for the bias-variance estimates (generated in `12_biasvariance_convergence.R`).
- `jointcell.png`: Evaluation of the homogeneity within the joint distributions (generated in `09_homogeneity_within_joints.R`).
- `sensmar.png`: Evaluation of the sensitivity of the estimates with respect to the MAR patterns and degrees (generated in `10_fullgrid_visualizations.R`).
- `sensmnar.png`: Evaluation of the sensitivity of the estimates with respect to the MNAR main and interaction effects scenarios (generated in `10_fullgrid_visualizations.R`).
- `mardiff.png`: Full-grid driver plot under MAR (generated in `10_fullgrid_visualizations.R`).
- `marstep.png`: More detailed driver plot regarding overlap and external sample size under MAR (generated in `10_fullgrid_visualizations.R`).
- `mnardiff.png`: Full-grid driver plot under MNAR (generated in `10_fullgrid_visualizations.R`).
- `mnarstep.png`: More detailed driver plot regarding overlap and external sample size under MNAR (generated in `10_fullgrid_visualizations.R`).
- `marscenres.png`: Scenario-based results under MAR (generated in `13_scenario_results.R`).
- `mnarscenres.png`: Scenario-based results under MNAR (generated in `13_scenario_results.R`).

### Output

The `Output/` folder contains the following objects:

- `A1_benchmark`: Monte Carlo convergence results for scenario A1 (generated in `12_biasvariance_convergence.R`).
- `D7_benchmark`: Monte Carlo convergence results for scenario D7 (generated in `12_biasvariance_convergence.R`).
- `H16_benchmark`: Monte Carlo convergence results for scenario H16 (generated in `12_biasvariance_convergence.R`).
- `Scenario Calculation Results`: Raw Monte Carlo simulation results for all scenarios (generated in `13_scenario_results.R`).
- `Final Estimates`: Aggregated scenario-level performance estimates for all three estimators (generated in `13_scenario_results.R`).

## Software requirements

This project was developed using R (version 4.4.3), Quarto and renv. The latter was used to manage all package dependencies.

## Reproducibility

To reproduce the results, download or clone this repository and open `Code/SMPV.Rproj` in RStudio. First, restore the project-specific package environment by running `renv::restore()` in the R console. Then run the scripts in the `Code/` folder in numerical order, starting from `01_packages.R` and ending with `13_scenario_results.R`. The main full-grid simulation results are stored in `Data/MC50.rds`. Additional simulation outputs required for the manuscript figures are stored in the `Output/` folder. Finally, render the thesis manuscript by opening `mMnuscript/Thesis.qmd` in Quarto or RStudio and selecting `Render`. 

## Data, ethics and privacy

This project is fully based on simulated data. No empirical or personally identifiable data are included in this repository. Ethical consent for this study was granted (FETC File Number: 25-1983) by the Ethics Review Board of the Faculty of Social and Behavioural Sciences of Utrecht University. 

## Contact

- Author: Robert van der Kaap (responsible for this archive)
- Research Master's programme: Methodology and Statistics for the Biomedical, Behavioural and Social Sciences
- Utrecht University
- DOI: 10.5281/zenodo.xxxxxxx
