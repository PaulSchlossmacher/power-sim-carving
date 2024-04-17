# Semester-Project-Multicarving
This project is comparing the carving estimator used in "Multicarving for high-dimensional post-selection inference"(https://arxiv.org/pdf/2006.04613.pdf) with the carving estimator presented in "A parametric distribution for exact post-selection inference with data carving"(https://arxiv.org/pdf/2305.12581.pdf).

## Content
* Multicarving Christoph: Whole repository obtained from https://github.com/cschultheiss/Multicarving, implementing the multicarving procedure
* carve_linear: Implementation of the carving estimator presented by Erik Drysdale in his paper
* SNTN_distribution: CDF and PDF of a sum of normal and truncated normal distribution following again Erik Drysdale
* SNTN_visualization: Plots an example CDF and PDF of the sntn distribution
* split_select: This is exactly the "one.split" function extracted from "multi.carve" in "Multicarving Christoph". Performs Lasso selection while ensuring correct conditions for carving procedure.
* Power Study Toeplitz: Main simulation file comparing the power and Type I error rate of both carving estimators
* p_vals_screened_distribution: Calculation and visualization of the distribution of p-values obtained from carve_linear under given screening (all truly active predictor variables where selected by "split_select") and given that the underlying coefficients where truly inactive
* Manual_Calc: For debugging purposes on carve_linear
* Plots: Contains plots obtained from "p_vals_screened_distribution" and "SNTN_visualization"
* Multicarving.pdf: Blog-like history of the project
