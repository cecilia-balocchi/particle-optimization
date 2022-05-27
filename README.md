# Crime in Philadelphia: Bayesian Clustering with Particle Optimization

This repository contains the code for replication of the results in the paper "Crime in Philadelphia: Bayesian Clustering with Particle Optimization" by Balocchi, Deshpande, George and Jensen ([arxiv link](https://arxiv.org/abs/1912.00111)).
It also contains an R package implementing the particle optimization (PARTOPT) algorithm.

If you want to use PARTOPT, see below for installation instructions of `PARTOPT`.

The sub-directory `two_partitions` is the main folder and contains the code for a model where there are two latent partitition, such as the partition of the mean level of crime and of the temporal trend, as presented in the paper.

The sub-directory `one_partition` contains the code for a model where there is only one latent partitition. This was used in a previous version of the paper, and will be deprecated soon.

The sub-directory `two_partitions_laplace` contains the code for a model where there are two latent partitition, but the outcome is count data modeled as a Poisson random variable. The marginal likelihood is approximated using a Laplace approximation.

## R Package - Installation

The package source files are contained in the sub-directory `PARTOPT/` . To install, you can either download that directory and then build and install the package from the command line (e.g. `R CMD BUILD ...` and `R CMD INSTALL ...`). You can also install the package using `devtools::install_github` as follows:

```r
library(devtools)
devtools::install_github(repo = "cecilia-balocchi/particle-optimization/PARTOPT")
```

A demo script can be found in `PARTOPT/R/demo.R`. This generates some data and shows some examples of usages of `PARTOPT`.