## Overview

The directory `src` contains an Rcpp implementation of the particle optimization procedure, specifically written for the case where each particle corresponds to a pair of partititions (one partition for the mean level of crime and one for the temporal trend). 
The directory `scripts` contains the R scripts for running this code and analyzing the results, as reported in Section 5  and S3.

To replicate the analysis, you can use the script `particle_opt_tracts.R` in the `script` directory. 
Our simulations were run on a high-performance computing cluster and script was slightly changed to set the specific prior distribution. The output will be saved in the `results` folder.

To replicate the plots, you can use the script `plots.R`, which calls the functions in `plot_functions.R`. 
The script `predictions.R` analyzed the out of sample RMSE for the different prior choices and models.

The folder `data` contains the crime data, which was downloaded from [Open Data Philly](https://www.opendataphilly.org/dataset/crime-incidents) and used for the analysis. Moreover, it contains some additional datasets to compute the out of sample predictions and plots.

The data has been aggregated within census tracts and years using the code in the repo [Urban-project](https://github.com/cecilia-balocchi/Urban-project).
