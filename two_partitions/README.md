## Overview

The directory `src` contains an Rcpp implementation of the particle optimization procedure, when each particle corresponds to a pair of partitions (one for the mean level and one for the temporal trend).
It also contains the code to implement the method by Anderson et al (2017) "Spatial clustering of average risks and risk trends in bayesian disease mapping": the folder `src/anderson2017_code/` contains the scripts downloaded from https://bitbucket.org/craiganderson1989/model-code, while `src/load_anderson2017.R` is a helper function I wrote. 
Lastly, the folder `src/equalpart/` contains a different implementation of the particle optimization procedure, when each particle corresponds to a pair of *equal* partitions, i.e. when the partition of the mean levels is constrained to be equal to the partition of the time trends.

The R scripts used in Section 4 and 5 of the main manuscript and Section S3 and S4 of the Supplementary Materials may be found in the directory `scripts`.

### Section 4, or Simulation analysis

To replicate the simulation analysis, run the script `scripts/generate_simulation_data.R`. This will create several .RData objects in the `data` directory (`alphabetas.RData` and `partitions.RData`), which contain a pair of partitions of the 20 x 20 grid and different vectors of means and time trends.

The script `scripts/figure3.R` uses these .RData objects to produce Figue 3.

The script `scripts/simulation.R` runs the particle optimization method together with other competing methods on a batch of simulated datasets for different values of the cluster separation, and saves the results in an .RData object in the `results` directory.

To perform the full simulation analysis we run the simulation script for several cluster separation settings and for each of them we evaluated our method on 100 datasets, generated with 100 different draws of `alpha_i` and `beta_i`, as computed in `scripts/generated_simulation_data.R` and saved in `data/alphabetas.RData`. These simulations were run on a high-performance computing cluster. Moreover, the simulations were split into several batches and called directly from the command line. Each of these simulation scripts contains code of the form
```{r}
args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])
batch <- as.numeric(args[2])
```
which take in command line arguments. The argument `sim_number` was varied from 1 (high cluster separation) to 3 (low cluster separation), and `batch` from 1 to 25 when `batch_size` is set equal to 4.
To run the scripts within an IDE like RStudio, you will need to remove these lines and set the simulation number (which corresponds to the specification of the cluster separation setting) and the batch manually. Each of the simulation scripts saves a set of results to the `results` directory.

The main simulation study reported in Section 4 was performed by running `scripts/simulation.R` with `sim_number` equal to 2. Section S3 of the Supplementary Materials reports additional results the same cluster separation setting, but also for high and low cluster separations (`sim_number` equal respectively 1 and 3). 
A separate script was used to run the simulation with the SCC method, `scripts/simulation_scc/R`. Before running this script, the script `scripts/prepare_scc_design.R` should be run. 
Moreover, `scripts/simulation_counts.R` is the script used to perform the analysis when the data was generated from a Poisson distribution, `scripts/simulation_sensitivity.R` studies the performance of our method for different values of `rho`, and `scripts/simulation_equalpart.R` analyzes how the particle optimization procedure would perform if we constrained the two partitions in each particle to coincide. All of these analyses are reported in Section S3 of the Supplementary Materials.

We combine the all outputs from the simulation batches into a unique .RData file using `scripts/compile_simulation_results.R`. For convenience, the .RData files containing combined results are reported in `results/summaries/`.

The script `scripts/figure4.R` takes these combined results and produces Figure 4 together with Figures S1-S3. Figures S4-S6 can be reproduced using a similar script.

### Section 5, or Analysis of Crime in Philadelphia

The directory `get_data` contains instructions and code to download and aggregate the data used for the analysis of crime in Philadelphia. See detailed instructions in `get_data/readme.md`. For convenience, the final aggregated datasets are already provided in the `data` folder.

The script `scripts/figure1.R` uses the crime density data to produce Figue 1.

The script `scripts/tracts_partopt.R` was used to analyze crime density in Philadelphia's census tracts using the Particle Optimization method. The output from this analysis is saved in the `results` folder for convenience.

The maps in Figure 5 are produced from these results using `scripts/figure5.R`.

The scripts `scripts/tracts_KmSc.R`, `scripts/tracts_anderson2017.R` and `scripts/tracts_scc.R` implement the competing methods on the same crime data. Before running `scripts/tracts_scc.R`, make sure to run `scripts/compile_simulation_results_tracts.R`.

With `scripts/table1_predictions.R` we compare the performance of all these methods on out-of-sample crime density data (corresponding to crimes for the year 2018). The same script reproduces the results reported in Table 1.
