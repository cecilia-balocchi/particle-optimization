# Crime in Philadelphia: Bayesian Clustering with Particle Optimization

This repository contains the code for replication of the results in the paper "Crime in Philadelphia: Bayesian Clustering with Particle Optimization" by Balocchi, Deshpande, George and Jensen.
It also contains an Rcpp implementation of the particle optimization algorithm.

The sub-directory `one_partition` contains the code for a model where there is only one latent partitition, and corresponds to the simulations reported in section 4.
The sub-directory `two_partitions` contains the code for a model where there are two latent partitition, such as the partition of the mean level of crime and of the temporal trend, and corresponds to the results reported in section 5.
