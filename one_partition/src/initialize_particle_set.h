//
//  initialize_particle_set.hpp


#ifndef GUARD_initialize_particle_set_h
#define GUARD_initialize_particle_set_h

#include<RcppArmadillo.h>
#include "partition.h"
#include "various_functions.h"
#include "partition_functions.h"
#include "rng.h"
#include <stdio.h>


void initialize_particle_set(std::vector<LPPartition> &particle_set, const int L, const arma::vec ybar, const double total_ss, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double eta, const double nu_sigma, const double lambda_sigma, const int init_id, RNG &gen);


void initialize_particle(LPPartition &gamma_init, const arma::vec ybar, const double total_ss, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double eta, const double nu_sigma, const double lambda_sigma, RNG &gen);


#endif
