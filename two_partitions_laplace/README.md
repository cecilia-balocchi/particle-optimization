## Overview

The directory `src` contains an Rcpp implementation of the approximate particle optimization procedure using Laplace approximation, where count data are modeled as Poisson random variables and the marginal likelihood is computed using the Laplace approximation.

The R scripts to replicate the results shown in Section S6 of the Supplementary Materials can be found in the directory `scripts`.

### Details

To replicate the small illustration, first generate the data using `scripts/generate_simulation_data/R`. This will create several .RData objects in the `data` directory (`alphabetas.RData` and `partitions.RData`), which contain a pair of partitions of the 10 x 10 grid and different vectors of means and time trends.

Then run the illustrations with `scripts/simulation_counts.R` (this was run with parameter `sim_number` equal to 1, 2, and 3 to generate the high, medium and low cluster separation configurations).
This saves the results in an .RData object in the `results` directory.

The script `scripts/plot_results.R` will create the figures shown in the Supplementary Materials.