# Semester-Project-Multicarving
This project is comparing the carving estimator used in "Multicarving for high-dimensional post-selection inference"(https://arxiv.org/pdf/2006.04613.pdf) with the carving estimator presented in "A parametric distribution for exact post-selection inference with data carving"(https://arxiv.org/pdf/2305.12581.pdf).

## Content
* Multicarving_paper: Whole repository obtained from https://github.com/cschultheiss/Multicarving, implementing the multicarving procedure
* carve_combined: Implementation of the carving estimator presented by Erik Drysdale in his paper, as well as PoSI estimator from the paper written by Lee et al. (https://arxiv.org/pdf/1311.6238) and data splitting.
* SNTN_distribution: CDF and PDF of a sum of normal and truncated normal distribution following again Erik Drysdale
* SNTN_visualization: Plots an example CDF and PDF of the sntn distribution
* split_select: This is exactly the "one.split" function extracted from "multi.carve" in "Multicarving Christoph". Performs Lasso selection while ensuring correct conditions for carving procedure.
* parallelized_power_study: Main simulation file comparing the power and FWER of both carving estimators. Here we included also the carving estimator in the saturated view, data splitting estimator and pure PoSI estimator. It is possible to specify how many selection attempts are given to data splitting and combined carving approaches with the flexible_selection_count variable.
* parallelized_power_study_diff_fraq: Simulation for the same estimators as in parallelized_power_study, but now reducing fractions of data used for selection in the case that data splitting and combined carving are not well-defined
* p_vals_screened_distribution: Calculation and visualization of the distribution of p-values obtained from carve_linear under given screening (all truly active predictor variables where selected by "split_select") and given that the underlying coefficients were truly inactive
* Plots: Contains plots obtained from p_vals_screened_distribution.R, parallelized_power_study.R and parallelized_power_study_diff_fraq.R
* Environments: the simulation environments that generate the plots
