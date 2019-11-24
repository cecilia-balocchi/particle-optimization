//
//  partition.h


#ifndef GUARD_PARTITION_H_
#define GUARD_PARTITION_H_

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <vector>

#define LTPI 1.83787706640934536

using Rcpp::Rcout;
using arma::endl;



typedef class Partition* LPPartition;
class Partition
{
public:
  int nObs; // number of indices
  int K; // number of clusters
  std::vector<int> cluster_config; // sizes of each cluster
  std::vector<int> cluster_assignment; // size nObs. Gives us the cluster label
  std::vector<std::vector<int> > clusters; // the actual clusters
  arma::mat pairwise_assignment; // array of pairwise assignments
  
  //std::vector<double> log_like;
  std::vector<double> log_prior;
  std::vector<double> alpha_hat;
  std::vector<double> alpha_bar;
  //std::vector<arma::mat> Omegay; // holds the Omega_y
  std::vector<double> log_det_Omegay; // holds log determinant of Omegay
  std::vector<double> y_Omegay_y; // holds
  
  
public:
  Partition();// constructor
  Partition(LPPartition initial_partition);
  Partition(int n, Rcpp::List gamma_init, const arma::vec &ybar, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double eta);
  Partition(int n, std::vector<std::vector<int> > init_clusters, const arma::vec &ybar, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double eta);
  
  ~Partition(); // destructor
public:
  void Print_Partition(const double total_ss, const int T, double nu_sigma, double lambda_sigma);
  void Copy_Partition(LPPartition initial_partition); // overwrite attributes with the ones from initial_partition
public:
  void get_pairwise();
  void get_Omega(int cluster_id, const arma::vec &ybar, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2);
  void log_pi_ep(int cluster_id, const double eta);
  void alpha_postmean(int cluster_id, const arma::vec &ybar, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2);
  void Split(int split_k, std::vector<std::vector<int> > &new_clusters, const arma::vec &ybar, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double eta);
  void Merge(int k_1, int k_2, const arma::vec &ybar, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double eta);
  void Split_Merge(int split_k, std::vector<std::vector<int> > &new_clusters, std::vector<int> k_star, const arma::vec &ybar, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double eta);

};
#endif /* partition_hpp */
