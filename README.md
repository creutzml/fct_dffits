# fct_dffits
A collection of code and results Chapter 4, Creutzinger (2024): "Identifying Influential Observations in a Concurrent Functional Regression Model". 

## General Structure of Repository
The contents of this repository are organized into three main folders: R Code, Data, and Figures. A short description of each folder is given below.

- **R Code**: folder that contains all R scripts for implementing Theoretical functional DFFITS, creating a random sample of functional data with influential observations present, reproducing simulations, and reproducing the case study.
- **Data**: the saved simulation .RData files are provided here for easy reproducibility of the simulation figures and tables. If any new simulation iterations are run, a new folder is created in this directory and the new simulation results are saved there. 
- **Figures**: copies of the plots and figures presented in the manuscript.

## R Code: Details
The primary R function of interest is `fRegress_concurrent()`, which comes from the `creutzml/ffscbExtra` package. When `fRegress_concurrent()` is used, it returns back a list of objects which includes the functional DFFITS.

**An Example**:
```
# Load the Ojo et al. (2021) library for easy data generation
library(fdaoutlier)

# Create an example data set
test_data_obj <- simulation_model7()
test_data <- test_data_obj$data

# Implement POD
pod_test <- pod_fda(test_data, cutoff = "classical1.5")
```

The remaining scripts in the folder allow one to recreate the simulation results and case study results found in the manuscript. Specifically, a description of each script is given below:
- `pittman_data_gen.R`: this allows one to use Pittman's (2022) parametric bootstrapping algorithm and to produce random samples of functional data for use in the simulations.
- `pittman_data_gen.R`: this allows one to use Pittman's (2022) parametric bootstrapping algorithm, without the use of any regularization or smoothing
- `pittman_case_study.R`: this will reproduce the case study application presented in Chapter 3, Pittman (2022) and Chapter 4, Creutzinger (2024)
- `run_dffits_pointwise_sim.R`: this reproduces the first, pointwise, simulation in Chapter 4, Creutzinger (2024)
- `results_dffits_pointwise_sim.R`: this script summarizes the simulation results acquired from `run_dffits_pointwise_sim.R`
- `run_dffits_creutz_v_pitt_sim.R`: this script reproduces the simulation used to compare Theoretical and Bootstrapped methodologies in Chapter 4, Creutzinger (2024)
- `results_dffits_creutz_v_pitt_sim.R`: this script summarizes the simulation results acquired from `run_dffits_creutz_v_pitt_sim.R`

## Data: Details
Simulation results for the pointwise simulation are contained in the following .RData files: `dffits_exp_mod1_df_11_23_23.RData`, `dffits_exp_mod2_df_11_23_23.RData`, and `dffits_exp_mod3_df_11_23_23.RData`. The summary results of those simulation results are found in `dffits_exp_mod1_sum_11_23_23.RData`, `dffits_exp_mod2_sum_11_23_23.RData`, and `dffits_exp_mod3_sum_11_23_23.RData`.

The simulation results for the functional simulation, which compares Creutzinger and Pittman, are contained in: `dffits_creutz_v_pitt_sim_a_0.RData`, `dffits_creutz_v_pitt_sim_a_0_T_100.RData`, `dffits_creutz_v_pitt_sim_a_0.25.RData`, `dffits_creutz_v_pitt_sim_a_0.25_T_100.RData`, `dffits_creutz_v_pitt_sim_a_0.5.RData`, and `dffits_creutz_v_pitt_sim_a_0.5_T_100.RData`.

Lastly, the data for the case study is stored in `FinalXtL1star.csv`, `FinalYtL1star.csv`, and `Oct15CongHt.txt`. These files were obtained from [Pittman's GitHub](https://github.com/rpittman188/fdaconcur).


## Figures: Details
This folder contains pdf copies of figures created in the manuscript. The file names and corresponding figure number are given below:
- `x_data_ex.pdf`: Figure 4.1
- `dffits_sim_models_dffits.pdf1`: Figure 4.2
- `dffits_creutz_v_pitt_mcc_n.pdf`: Figure 4.3
- `dffits_creutz_v_pitt_mcc_T.pdf`: Figure 4.4
- `congaree_river_realigned.pdf`: Figure 4.5
- `pitt_dffits.pdf`: Figure 4.6
- `creutz_dffits.pdf`: Figure 4.7
- `dffits_sim_mod1_n=10_mcc.pdf`: Figure C.1
- `dffits_sim_mod1_n=50_mcc.pdf`: Figure C.2
- `dffits_sim_mod1_n=100_mcc.pdf`: Figure C.3
- `dffits_sim_mod2_n=10_mcc.pdf`: Figure C.4
- `dffits_sim_mod2_n=50_mcc.pdf`: Figure C.5
- `dffits_sim_mod2_n=100_mcc.pdf`: Figure C.6
- `dffits_sim_mod3_n=10_mcc.pdf`: Figure C.7
- `dffits_sim_mod3_n=50_mcc.pdf`: Figure C.8
- `dffits_sim_mod3_n=100_mcc.pdf`: Figure C.9
- `dffits_creutz_v_pitt_specificity_n.pdf`: Figure C.10
- `dffits_creutz_v_pitt_specificity_T.pdf`: Figure C.11
- `dffits_creutz_v_pitt_sensitivity_n.pdf`: Figure C.12
- `dffits_creutz_v_pitt_sensitivity_T.pdf`: Figure C.13

## References
- Creutzinger, Michael L. (2024, May). Chapter 3

- Liebl, Dominik, and Matthew Reimherr. (2023, July). “Fast and Fair Simultaneous Confidence Bands for Functional Parameters.” Journal of the Royal Statistical Society Series B: Statistical Methodology 85(3):842–68. doi: 10.1093/jrsssb/qkad026.

- Pittman, R. (2022), “Using Concurrent Functional Regression to Reconstruct River Stage Data During Flood Events and Identify Influential Functional Measurements.”
  
- R Development Core Team. (2021). R: A Language and Environment for Statistical Computing. Retrieved from http://www.r-project.org

## Contributors:
- Michael L. Creutzinger
  - Doctoral Candidate (May 2024)
  - Department of Statistics
  - Colorado State University

- Julia L. Sharp, PhD
  - Mathemathical Statistican
  - National Institute of Science and Technology
