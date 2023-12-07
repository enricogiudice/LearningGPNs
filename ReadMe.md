Results/ contains the collected data from the GP_sims.R and Sachs.R;
Plots/ contains all generated plots and the code used to create Figures 1, A2, A3, A4, A6, A7 & A8.

The .Stan files contain code for the Gaussian (Gauss.stan), GP (Add.stan) & GP^2 (Add_interact.stan) models.

BayesStanFns.R and BayesStanFns_interact.R contain functions to compute the score according to the GP and GP^2 models respectively.

Fourier_fns.R contains functions to generate data according to Equation (14).

GP_sims.R contains the code for running the main simulation study (the data for Plots/Plots.R).

Sachs.R generates the data for Table 1.

comparison_algs.R contains functions for implementing k_PC, DiBS+ and comparing the results.

posterior.R generates Figures 2 and A5.

sampling_fns.R contains the main functions for BGe & GP score-based structure learning.

score_equivalence.R generates Figure A1.
