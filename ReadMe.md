This repository contains all the code to reproduce the results in [A Bayesian Take on Gaussian Process Networks](https://papers.nips.cc/paper_files/paper/2023/hash/b146e7c87685fa208bd95ce4b08e330c-Abstract-Conference.html).

- Results/ contains the collected data from the GP_sims.R and Sachs.R.
- Plots/ contains all generated plots and the code used to create Figures 1, 2, A2, A3, A7, A8, & A9.

The .Stan files contain code for the Gaussian (Gauss.stan), GP (Add.stan) and GP<sup>2</sup> (Add_interact.stan) models.

- BayesStanFns.R and BayesStanFns_interact.R contain functions to compute the score according to the GP and GP^2 models respectively.
- Fourier_fns.R contains functions to generate data according to Equation (15).
- GP_sims.R contains the code for running the main simulation study (the data for Plots/Plots.R).
- Sachs.R generates the data for Table 1.
- comparison_algs.R contains functions for implementing k_PC, DiBS+ and comparing the results.
- posterior.R generates Figures 3 and A4.
- sampling_fns.R contains the main functions for BGe & GP score-based structure learning.
- score_equivalence.R generates Figure A1.

Reference
---------

```
@inproceedings{NEURIPS2023_b146e7c8,
 author = {Giudice, Enrico and Kuipers, Jack and Moffa, Giusi},
 booktitle = {Advances in Neural Information Processing Systems},
 editor = {A. Oh and T. Neumann and A. Globerson and K. Saenko and M. Hardt and S. Levine},
 pages = {56602--56614},
 publisher = {Curran Associates, Inc.},
 title = {A Bayesian Take on Gaussian Process Networks},
 url = {https://proceedings.neurips.cc/paper_files/paper/2023/file/b146e7c87685fa208bd95ce4b08e330c-Paper-Conference.pdf},
 volume = {36},
 year = {2023}
}
```
