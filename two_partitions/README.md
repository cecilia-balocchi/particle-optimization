## Overview

The directory `get_data` contains instructions and code to download and aggregate the data used. For convenience, the final aggregated datasets are already provided in the `data` folder.

The directory `src` contains an Rcpp implementation of the particle optimization procedure, specifically written for the case where each particle corresponds to a pair of partititions (one partition for the mean level of crime and one for the temporal trend). 
The directory `scripts` contains the R scripts for running this code and analyzing the results, as reported in Section 5 and S3 of the Supplementary Materials.

To replicate the analysis, you can use the script `particle_opt_tracts.R` in the `script` directory. 
Our simulations were run on a high-performance computing cluster and script was slightly changed to set the specific prior distribution. The output will be saved in the `results` folder. For convenience I have saved in `results` the output from the simulations that I ran.

To replicate the plots, you can use the script `plots.R`, which calls the functions in `plot_functions.R`. 
The script `predictions.R` analyzed the out of sample RMSE for the different prior choices and models.

The folder `data` contains the crime data, which was downloaded from [Open Data Philly](https://www.opendataphilly.org/dataset/crime-incidents), aggregated and used for the analysis. The folder `data` also contains some additional datasets to compute the out of sample predictions and plots.

### How to run the algorithm with the different priors

In the paper we use our algorithm with different priors. To replicate that analysis you can use the script `particle_opt_tracts.R` and tune prior parameters: 
- for the EP-EP prior, set `priorA_input = 0, priorB_input = 0`
- for the EP-Unif prior, set `priorA_input = 0, priorB_input = 3`
- for the Unif-EP prior, set `priorA_input = 3, priorB_input = 0`
- for the Unif-Unif prior, set `priorA_input = 3, priorB_input = 3`
