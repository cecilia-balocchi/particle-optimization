//
//  partition_summary.cpp
//  
//
//

#include <stdio.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "partition.h"
#include "partition_functions.h"
#include <vector>
#include <ctime>

using namespace arma;
using namespace std;

// [[Rcpp::export]]
Rcpp::List partition_summary(arma::mat Y,
                             const arma::mat A_block,
                             Rcpp::List gamma_init,
                             const double a1 = 1.0,
                             const double a2 = 1.0,
                             const double nu_sigma = 3,
                             const double lambda_sigma = 1,
                             const double eta = 1.0,
                             const double rho = 0.99)
{
  int n = Y.n_rows;
  int T = Y.n_cols;
  
  arma::vec ybar(n);
  double total_ss = 0;
  for(int i = 0; i < n; i++){
    ybar(i) = arma::mean(Y.row(i));
    total_ss += (T-1) * arma::var(Y.row(i));
  }
  
  LPPartition gamma_0 = new Partition(n, gamma_init, ybar, T, A_block, rho, a1, a2, eta);
  Rcpp::List particle_out(1);
  particle_out[0] = gamma_init;
  
  Rcpp::List results;
  results["particles"] = particle_out;
  results["log_like"] = total_log_like(gamma_0, total_ss, T, nu_sigma, lambda_sigma);
  results["log_prior"] = total_log_prior(gamma_0);
  results["log_post"] = total_log_post(gamma_0, total_ss, T, nu_sigma, lambda_sigma);
  results["alpha"] = gamma_0->alpha_hat;
  return(results);
}
