# Statistical Matching with a Proxy Variable

## Project description

This repository contains the research archive accompanying the MSc thesis "Statistical Matching with a Proxy Variable: The Role of Sample Size, Unit Overlap, Selectivity, Proxy Quality and the Validity of the Conditional Independence Assumption." for the Research Master's programme Methodology and Statistics for the Biomedical, Behavioural and Social Sciences. 

This project examines statistical matching using doubly robust estimation methods and simulation-based analyses.

## Repository structure

### Code

The code folder contains the following scripts, which should be run in the indicated numerical order:

01_packages.R: Install and load the required packages.
02_utils.R: Helper functions. 
03_parametergrid.R: Definition of the simulation parameter grid.
04_selectivity.R: Operationalization of the selectivity mechanisms.
05_generate_population.R: Function to generate simulation populations.
06_monte_carlo.R: Function to perform Monte Carlo simulation. 
07_run_simulation.R: Function to parallelise the full grid simulation.
08_fullgrid_results.R: Performane of the full grid simulation.
09_homogeneity_within_joints.R: Homogeneity evaluation within joint distributions.
10_biasvariance_convergence.R: Analysis of the bias-variance estimate convergence.
11_fullgrid_visualizations.R: Visualization for the full-grid results.
12_scenarios.R: Definition of the simulation scenarios.
13_scenario_results.R: Calculation and visualization of the scenario-based analysis.

### Data

This folder contains the following objects:

MC50.rds: RDS-file containing the full grid simulation results.
MC50.text: Text file containing all parallized simulation logs.

### Manuscript


### Output

## Software requirements

This project was developed using R (version 4.4.3), quarto and renv. The latter was used to manage all package dependencies.

## Reproducibility



## Data



## Ethics and privacy

Ethical consent for this study was granted (FETC File Number: 25-1983) by the Ethics Review Board of the Faculty of Social and Behavioural Sciences of Utrecht University. Moreover, this project does not contain personally identifiable information.

## Contact

Author: Robert van der Kaap
Research Master's programme: Methodology and Statistics for the Biomedical, Behavioural and Social Sciences
Utrecht University
