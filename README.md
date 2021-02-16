# Crime in Philadelphia: Bayesian Clustering with Particle Optimization

This repository contains the code for replication of the results in the paper "Crime in Philadelphia: Bayesian Clustering with Particle Optimization" by Balocchi, Deshpande, George and Jensen.
It also contains an Rcpp implementation of the particle optimization algorithm.

The sub-directory `two_partitions` is the main folder and contains the code for a model where there are two latent partitition, such as the partition of the mean level of crime and of the temporal trend, as presented in the paper.

The sub-directory `one_partition` contains the code for a model where there is only one latent partitition. This was used in a previous version of the paper, and will be deprecated soon.