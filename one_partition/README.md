## Overview

The directory `src` contains an Rcpp implementation of the particle optimization procedure.
The R scripts used in Section 4 of the main manuscript and Section S2 of the Supplementary Materials may be found in the directory `scripts`.
To replicate these analyses, run the script "generated_data.R".
This will create several .RData objects in the `data` directory that contain (i) different partitions of the 20 x 20 grid, (ii) different vectors of cluster means, (iii) several different draws of the parameters alpha_i.

The script "illustration.R" runs the particle optimization method on a single simulated dataset for each of the six different settings of the cluster means and saves the results in several .RData objects in the `results` directory.
The scripts "figure4.R", "figure5.R", and  "figureS1.R" take these results and respectively produce Figures 4, 5, and S1 in the main manuscript and supplementary materials.

The main simulation was divided into several scripts: 
* particle_simulation.R: runs the proposed particle optimization method
* map_simulation.R: finds the MAP particle
* fixed_km_sc_simulation.R: runs k-means, spectral clustering, and other non-adaptive methods

These simulations were run on a high-performance computing cluster.
Moreover, the simulations were split into several batches and called directly from the command line.
Each of these simulation scripts contains code of the form
```{r}
args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])
batch <- as.numeric(args[2])
```
which take in command line arguments.
To run the scripts within an IDE like RStudio, you will need to remove these lines and set the simulation number (which corresponds to the specification of the cluster means) and the batch manually.
Each of the simulation scripts saves a set of results to the `results` directory.
The script "compile_results.R" combines all of the results. 

The script "one_partition_rmse_L.R" explores how the RMSE of the BMA estimator varies as a function of the number of particles used (L). 

The remaining scripts are utility scripts that are calld by the others.  
