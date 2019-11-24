//
//  partition.cpp
//


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
#include "partition.h"
#include "various_functions.h"



using Rcpp::Rcout;
using namespace arma;
Partition::Partition()
{
  nObs = 1;
  K = 1;
  cluster_config = std::vector<int>(1,0);
  cluster_assignment = std::vector<int>(1,0); // size nObs. Gives us the cluster label
  clusters = std::vector<std::vector<int> >(1, std::vector<int>(1,0));
  pairwise_assignment = arma::zeros<arma::mat>(1,1);
  //log_like = std::vector<double>(1,0.0);
  //Omegay = std::vector<arma::matrix>(1, arma::matrix(1,1));
  y_Omegay_y = std::vector<double>(1,0.0);
  log_det_Omegay = std::vector<double>(1,0.0);
  log_prior = std::vector<double>(1, 0.0);
  alpha_hat = std::vector<double>(1, 0.0);
  alpha_bar = std::vector<double>(1, 0.0);
}

Partition::Partition(LPPartition initial_partition)
{
  // need to make a deep copy
  nObs = initial_partition->nObs;
  K = initial_partition->K;
  cluster_config = std::vector<int>(K,0);
  cluster_assignment = std::vector<int>(nObs,0);
  clusters = std::vector<std::vector<int> >(K, std::vector<int>(1,0));
  pairwise_assignment = arma::zeros<arma::mat>(nObs, nObs);
  //log_like = std::vector<double>(K,0.0);
  //Omegay = std::vector<arma::mat>(K, arma::matrix(1,1));
  log_det_Omegay = std::vector<double>(K,0.0);
  y_Omegay_y = std::vector<double>(K, 0.0);
  log_prior = std::vector<double>(K,0.0);
  alpha_hat = std::vector<double>(nObs,0.0);
  alpha_bar = std::vector<double>(K,0.0);
  //Rcpp::Rcout << "Created containers of right size" << endl;
  
  // update all of the attributes
  for(int k = 0; k < K; k++){
    cluster_config[k] = initial_partition->cluster_config[k];
    //Rcpp::Rcout << "updated cluster config for k = " << k << endl;
    clusters[k].resize(cluster_config[k],0);
    //Rcpp::Rcout << "re-sized clusters[k] for k = " << k << endl;
    for(int i = 0; i < cluster_config[k]; i++){
      clusters[k][i] = initial_partition->clusters[k][i];
    }
    //Rcpp::Rcout << "Successfully wrote clusters" << endl;
    //Omegay[k] = initial_partition->Omegay[k];
    y_Omegay_y[k] = initial_partition->y_Omegay_y[k];
    log_det_Omegay[k] = initial_partition->log_det_Omegay[k];
    
    //Rcpp::Rcout << "Successfully wrote log_like" << endl;
    log_prior[k] = initial_partition->log_prior[k];
    // Rcpp::Rcout << "Successfully wrote log_prior" << endl;
    alpha_bar[k] = initial_partition->alpha_bar[k];
    //Rcpp::Rcout << "Successfull wrote alpha_bar" << endl;
  }
  for(int i = 0; i < nObs; i++){
    alpha_hat[i] = initial_partition->alpha_hat[i];
    cluster_assignment[i] = initial_partition->cluster_assignment[i];
    for(int j = i; j < nObs;j++){
      pairwise_assignment(i,j) = initial_partition->pairwise_assignment(i,j);
      pairwise_assignment(j,i) = initial_partition->pairwise_assignment(j,i);
    }
  }
}

Partition::Partition(int n, Rcpp::List gamma_init, const arma::vec &ybar, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double eta)
{
  nObs = n;
  K = gamma_init.size();
  cluster_config = std::vector<int>(K, 0);
  cluster_assignment = std::vector<int>(nObs,0);
  clusters = std::vector<std::vector<int> >(K, std::vector<int>(1,0));
  pairwise_assignment = arma::zeros<arma::mat>(nObs,nObs);
  //log_like = std::vector<double>(K, 0.0);
  //Omegay = std::vector<arma::matrix>(K);
  log_det_Omegay = std::vector<double>(K,0.0);
  y_Omegay_y = std::vector<double>(K,0.0);
  
  log_prior = std::vector<double>(K, 0.0);
  alpha_hat = std::vector<double>(nObs,0.0);
  alpha_bar = std::vector<double>(K,0.0);
  
  //Rcpp::Rcout << "created tmp_vec" << endl;
  Rcpp::NumericVector tmp_vec;
  for(int k = 0; k < K; k++){
    //Rcpp::Rcout << "[partition]: k = " << k << endl;
    tmp_vec = gamma_init[k];
    cluster_config[k] = tmp_vec.size();
    clusters[k].resize(cluster_config[k], 0);
    for(int i = 0; i < cluster_config[k]; i++){
      clusters[k][i] = tmp_vec[i] - 1; // R is 1-indexed and C++ is 0-indexed
      cluster_assignment[clusters[k][i]] = k;
    }
    // now is the time to compute log-likelihood, log-prior, alpha_hat, and alpha_bar
    //log_likelihood(k, ybar, T, A_block, rho, a1, a2);
    get_Omega(k, ybar, T, A_block, rho, a1, a2);
    alpha_postmean(k, ybar, T, A_block, rho, a1, a2);
    log_pi_ep(k, eta);
    
  }
  // also run get_pairwise
  get_pairwise();
}

Partition::Partition(int n, std::vector<std::vector<int> > init_clusters, const arma::vec &ybar, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double eta)
{
  nObs = n;
  K = init_clusters.size();
  cluster_config = std::vector<int>(K,0);
  cluster_assignment = std::vector<int>(nObs,0);
  clusters = init_clusters;
  pairwise_assignment = arma::zeros<arma::mat>(nObs, nObs);
  log_det_Omegay = std::vector<double>(K, 0.0);
  y_Omegay_y = std::vector<double>(K, 0.0);
  
  log_prior = std::vector<double>(K, 0.0);
  alpha_hat = std::vector<double>(nObs,0.0);
  alpha_bar = std::vector<double>(K,0.0);
  
  for(int k = 0; k < K; k++){
    cluster_config[k] = init_clusters[k].size();
    for(int i = 0; i < cluster_config[k]; i++){
      cluster_assignment[clusters[k][i]] = k;
    }
    get_Omega(k, ybar, T, A_block, rho, a1, a2);
    alpha_postmean(k, ybar, T, A_block, rho, a1, a2);
    log_pi_ep(k, eta);
  }
  get_pairwise();
}



Partition::~Partition(){}

void Partition::Print_Partition(const double total_ss, const int T, double nu_sigma, double lambda_sigma)
{
  Rcpp::Rcout << "nObs = " << nObs << endl;
  Rcpp::Rcout << "K = " << K << endl;
  
  
  Rcpp::Rcout << "Size of clusters:";
  for(int k = 0; k < K; k++){
    Rcpp::Rcout << cluster_config[k] << " ";
  }
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "Clusters:" << std::endl;
  for(int k = 0; k < K; k++){
    Rcpp::Rcout << "Cluster " << k  << " : ";
    for(int j = 0; j < cluster_config[k]; j++){
      Rcpp::Rcout << clusters[k][j] << " ";
    }
    Rcpp::Rcout << std::endl;
  }
  
  // compute the total log-likelihood
  
  Rcpp::Rcout << "y_Omega_y:" ;
  for(int k = 0; k < K; k++){
    Rcpp::Rcout << y_Omegay_y[k] << " ";
  }
  Rcpp::Rcout << endl;
  
  Rcpp::Rcout << "log_det:" ;
  for(int k = 0; k < K; k++){
    Rcpp::Rcout << log_det_Omegay[k] << " ";
  }
  Rcpp::Rcout << endl;
  
  
  double log_like = 0.0;
  double log_det = 0.0;
  double quad_form = 0.0;
  for(int k = 0; k < K; k++){
    log_det += log_det_Omegay[k];
    quad_form += y_Omegay_y[k];
  }
  //log_like = 0.5 * log_det - ( (nu_sigma + nObs)/2) * log( (nu_sigma * lambda_sigma + quad_form)/2);
  log_like = 0.5 * log_det - ( (nu_sigma + ((double) T) * ((double) nObs))/2 ) * log( (nu_sigma + lambda_sigma + quad_form + total_ss) / 2);
  Rcpp::Rcout << "Log-likelihood: " << log_like;
  
  
  //double log_post = 0.0;
  //for(int k = 0; k < K; k++){
  //  log_post += log_like[k] + log_prior[k];
  //}
  //Rcpp::Rcout << "Log-posterior : " << log_post << std::endl;
  
  //Rcpp::Rcout << std::endl;
  //Rcpp::Rcout << "Log-likelihood:" ;
  //for(int k = 0; k < K; k++){
  //  Rcpp::Rcout << log_like[k] << " ";
  //}
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "Log-prior:" ;
  for(int k = 0; k < K; k++){
    Rcpp::Rcout <<  log_prior[k] << " ";
  }
  Rcpp::Rcout << std::endl;
  //Rcpp::Rcout << "Log-post:";
  //for(int k = 0; k < K; k++){
  //  Rcpp::Rcout << log_like[k] + log_prior[k] << " ";
  //}
  Rcpp::Rcout << std::endl;
  
  Rcpp::Rcout << "Alpha-bar: ";
  for(int k = 0; k < K; k++){
    Rcpp::Rcout << alpha_bar[k] << " ";
  }
  Rcpp::Rcout << endl;
  
  //Rcpp::Rcout << "Cluster assignment:";
  //for(int i = 0; i < nObs; i++){
  //  Rcpp::Rcout << cluster_assignment[i] << " ";
  //}
  //Rcpp::Rcout << endl;
  
}
void Partition::Copy_Partition(LPPartition initial_partition){
  // need to make a deep copy
  nObs = initial_partition->nObs;
  K = initial_partition->K;
  
  // need to clear and resize existing attributes
  cluster_config.clear();
  cluster_config.resize(K,0);
  for(int kk = 0; kk < clusters.size(); kk++){
    clusters[kk].clear();
  }
  clusters.resize(K, std::vector<int>(1,0));
  cluster_assignment.clear();
  cluster_assignment.resize(nObs, 0);
  
  pairwise_assignment.resize(nObs, nObs);
  
  alpha_hat.clear();
  alpha_hat.resize(nObs,0.0);
  
  //log_like.clear();
  //log_like.resize(K, 0.0);
  
  //Omegay.clear();
  //Omegay.resize(K);
  log_det_Omegay.clear();
  log_det_Omegay.resize(K);
  y_Omegay_y.clear();
  y_Omegay_y.resize(K);
  
  
  log_prior.clear();
  log_prior.resize(K,0.0);
  alpha_bar.clear();
  alpha_bar.resize(K,0.0);
  
  // update the value of cluster_config
  for(int k = 0; k < K; k++){
    cluster_config[k] = initial_partition->cluster_config[k];
    clusters[k].resize(cluster_config[k], 0);
    for(int i = 0; i < cluster_config[k]; i++){
      clusters[k][i] = initial_partition->clusters[k][i];
    }
    //log_like[k] = initial_partition->log_like[k];
    //Omegay[k] = initial_partition->Omegay[k];
    log_det_Omegay[k] = initial_partition->log_det_Omegay[k];
    y_Omegay_y[k] = initial_partition->y_Omegay_y[k];
    log_prior[k] = initial_partition->log_prior[k];
    alpha_bar[k] = initial_partition->alpha_bar[k];
  }
  for(int i = 0; i < nObs; i++){
    alpha_hat[i] = initial_partition->alpha_hat[i];
    cluster_assignment[i] = initial_partition->cluster_assignment[i];
    for(int j = i; j < nObs; j++){
      pairwise_assignment(i,j) = initial_partition->pairwise_assignment(i,j);
      pairwise_assignment(j,i) = initial_partition->pairwise_assignment(j,i);
    }
  }
  
}



// methods

void Partition::get_pairwise(){
  for(int i = 0; i < nObs; i++){
    for(int j = i; j < nObs; j++){
      if(cluster_assignment[i] == cluster_assignment[j]){
        pairwise_assignment(i,j) = 1;
        pairwise_assignment(j,i) = 1;
      } else{
        pairwise_assignment(i,j) = 0;
        pairwise_assignment(j,i) = 0;
      }
    }
  }
}

void Partition::get_Omega(int cluster_id, const arma::vec &ybar, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2)
{
  int n_k = cluster_config[cluster_id];
  
  arma::mat tmp_Omegay = arma::zeros<mat>(n_k,n_k);
  if(n_k == 1){
    //tmp_Omegay(0,0) = 1.0/(1.0/( (double) T) + a2 + a1/(1.0 - rho));
    //Omegay[cluster_id] = tmp_Omega_y;
    //log_det_Omegay[cluster_id] = log(tmp_Omegay(0,0));
    //y_Omegay_y[cluster_id] = ybar(clusters[cluster_id][0]) * ybar(clusters[cluster_id][0]) * tmp_Omegay(0,0);
    tmp_Omegay(0,0) = 1.0/( 1.0/( (double) T) + a2 + a1/(1.0 - rho));
    log_det_Omegay[cluster_id] = log(tmp_Omegay(0,0));
    y_Omegay_y[cluster_id] = ybar(clusters[cluster_id][0]) * ybar(clusters[cluster_id][0]) * tmp_Omegay(0,0);
    
  } else{
    arma::mat A_block_k = Submatrix(A_block, n_k, n_k, clusters[cluster_id], clusters[cluster_id]);
    arma::vec row_sums = arma::zeros<vec>(n_k);
    arma::vec y_vec = arma::zeros<vec>(n_k);
    for(int i = 0; i < n_k; i++){
      row_sums(i) = accu(A_block_k.row(i));
      y_vec(i) = ybar[clusters[cluster_id][i]];
    }
    arma::mat D = diagmat(row_sums);
    arma::mat A_star_k = D - A_block_k; // the unweighted Laplacian
    arma::mat Omega_alpha = rho * A_star_k;
    Omega_alpha.diag() += 1 - rho;
    arma::mat Sigma_alpha = arma::inv_sympd(Omega_alpha);
    arma::mat Sigma_y = a1 * Sigma_alpha;
    Sigma_y += a2;
    Sigma_y.diag() += 1/( (double) T);
    // need to invert a1 * Sigma_alpha + a2 * 11'
    // [a1 * Sigma_alpha + a2 * 11']^-1 = a1^-1 * Omega_alpha - a1^-2 * a2 * [Omega_alpha 1 1' Omega_alpha]/(1 + a1^-1a2 * 1'Omega_alpha1]
    // Note: Omega_alpha 1 is constant with entries 1 - rho
    //
    //arma::mat Sigma_y = 1.0/a1 * Omega_alpha;
    //Sigma_y -= 1.0/(a1 * a1) * a2 * (1 - rho) * (1 - rho) / (1 + n_k/a1 * a2 * (1 - rho));
    //Sigma_y.diag() += 1.0/( (double) T);
    tmp_Omegay = inv_sympd(Sigma_y);
    double tmp_log_det = 0.0;
    double tmp_sgn = 1.0;
    arma::log_det(tmp_log_det, tmp_sgn, tmp_Omegay);
    //Omega_y[cluster_id] = tmp_Omegay;
    log_det_Omegay[cluster_id] = tmp_log_det;
    y_Omegay_y[cluster_id] = as_scalar(y_vec.t() * tmp_Omegay * y_vec);
  }
}

void Partition::log_pi_ep(int cluster_id, const double eta)
{
  log_prior[cluster_id] = log(eta) + lgamma(cluster_config[cluster_id]);
}

void Partition::alpha_postmean(int cluster_id, const arma::vec &ybar, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2)
{
  int n_k = cluster_config[cluster_id];
  if(n_k == 1){
    //alpha_hat[clusters[cluster_id][0]] = T * ybar[clusters[cluster_id][0]]/(T + 1.0/(a1/(1-rho) + a2));
    //alpha_hat[clusters[cluster_id][0]] = T * ybar[clusters[cluster_id][0]]/(T + 1.0/( (a1/(1-rho) + a2)));
    //alpha_bar[cluster_id] = (1.0/a1) * (1.0 - rho) * alpha_hat[clusters[cluster_id][0]]/( 1.0/a1 * (1.0 - rho) + 1.0/a2);
    alpha_hat[clusters[cluster_id][0]] = (double) T * ybar[clusters[cluster_id][0]]/( (double) T + (1.0 - rho)/(a1 + a2 * (1.0 - rho)));
    alpha_bar[cluster_id] = (1.0/a1) * (1.0 - rho) * alpha_hat[clusters[cluster_id][0]]/(1.0/a2 + (1.0 - rho)/a1);
  } else {
    arma::mat A_block_k = Submatrix(A_block, n_k, n_k, clusters[cluster_id], clusters[cluster_id]); // pulls out the relevant submatrix
    arma::vec row_sums = arma::zeros<vec>(n_k);
    arma::vec y_vec = arma::zeros<vec>(n_k);
    for(int i = 0; i < n_k; i++){
      row_sums(i) = accu(A_block_k.row(i));
      y_vec(i) = ybar[clusters[cluster_id][i]];
    }
    //Rcpp::Rcout << "Got row_sums" << endl;
    
    mat D = diagmat(row_sums);
    mat A_star_k = D - A_block_k; // this is the graph laplacian
    
    arma::mat Omega_alpha = rho * A_star_k;
    Omega_alpha.diag() += 1.0 - rho;
    
    arma::mat V_inv = 1.0/a1 * Omega_alpha;
    V_inv -= (1.0 - rho) * (1.0 - rho)/a1 * a2/(a1 + a2 * (1.0 - rho) * (double) n_k);
    V_inv.diag() += (double) T;
    arma::mat V = arma::inv_sympd(V_inv);
    arma::vec tmp_alpha = ( (double) T) * V * y_vec;
    double tmp_alpha_sum = 0.0;
    for(int i = 0; i < n_k; i++){
      alpha_hat[clusters[cluster_id][i]] = tmp_alpha(i);
      tmp_alpha_sum += tmp_alpha(i);
    }
    alpha_bar[cluster_id] = (1.0/a1) * (1.0 - rho) * tmp_alpha_sum/(1.0/a2 + n_k * (1.0 - rho)/a1);
    
    /*
     arma::mat V_inv = 1.0/a1 * Omega_alpha;
     V_inv -= 1.0/(a1 * a1) * (1 - rho) * (1 - rho)/(1 + 1.0/a1 * a2 * n_k * (1 - rho));
     V_inv.diag()+=T;
     arma::mat V = inv_sympd(V_inv);
     arma::vec tmp_alpha = ( (double) T)* V * y_vec;
     double tmp_alpha_sum = 0.0;
     for(int i = 0; i < n_k; i++){
     alpha_hat[clusters[cluster_id][i]] = tmp_alpha(i);
     tmp_alpha_sum += tmp_alpha(i);
     }
     alpha_bar[cluster_id] = (1.0/a1 * tmp_alpha_sum * (1 - rho))/(1.0/a2 + n_k * (1-rho)/a1);
     */
  }
}

void Partition::Split(int split_k, std::vector<std::vector<int> > &new_clusters, const arma::vec &ybar, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double eta)
{
  // create a temporary copy of the main attributes
  int orig_K = K;
  int k = 0;
  std::vector<int> orig_cluster_config(orig_K,0);
  std::vector<std::vector<int> > orig_clusters(orig_K, std::vector<int>(1,0));
  std::vector<int> orig_cluster_assignment(nObs, 0);
  //std::vector<double> orig_log_like(orig_K,0.0);
  //std::vector<arma::mat> orig_Omegay;
  std::vector<double> orig_log_det_Omegay;
  std::vector<double> orig_y_Omegay_y;
  std::vector<double> orig_log_prior(orig_K,0.0);
  std::vector<double> orig_alpha_hat(nObs, 0.0);
  std::vector<double> orig_alpha_bar(orig_K, 0.0);
  arma::mat orig_pairwise_assignment = arma::zeros<mat>(nObs, nObs);
  
  for(int k = 0; k < orig_K; k++){
    orig_cluster_config[k] = cluster_config[k];
    orig_clusters[k].resize(cluster_config[k],0);
    for(int i = 0; i < cluster_config[k]; i++){
      orig_clusters[k][i] = clusters[k][i];
    }
    //orig_log_like[k] = log_like[k];
    //orig_Omegay[k] = Omegay[k];
    orig_log_det_Omegay[k] = log_det_Omegay[k];
    orig_y_Omegay_y[k] = y_Omegay_y[k];
    orig_log_prior[k] = log_prior[k];
    orig_alpha_bar[k] = alpha_bar[k];
    
  }
  for(int i = 0; i < nObs; i++){
    orig_cluster_assignment[i] = cluster_assignment[i];
    orig_alpha_hat[i] = alpha_hat[i];
    for(int j = i; j < nObs; j++){
      orig_pairwise_assignment(i,j) = pairwise_assignment(i,j);
      orig_pairwise_assignment(j,i) = pairwise_assignment(j,i);
    }
  }
  
  // delete and re-size
  int num_splits = new_clusters.size();
  K = orig_K + num_splits - 1;
  cluster_config.clear();
  cluster_config.resize(K, 0);
  
  for(int k = 0; k < orig_K; k++){
    clusters[k].clear();
  }
  clusters.clear();
  clusters.resize(K, std::vector<int>(0,1));
  
  //log_like.clear();
  //log_like.resize(K, 0.0);
  //Omegay.clear();
  //Omegay.resize(K);
  log_det_Omegay.clear();
  log_det_Omegay.resize(K);
  y_Omegay_y.clear();
  y_Omegay_y.resize(K);
  log_prior.clear();
  log_prior.resize(K,0.0);
  alpha_bar.clear();
  alpha_bar.resize(K,0.0);
  
  // loop over old clusters, copying over everything except for split_k
  for(int k = 0; k < orig_K; k++){
    if(k != split_k){
      cluster_config[k] = orig_cluster_config[k];
      clusters[k].resize(cluster_config[k], 0);
      for(int i = 0; i < cluster_config[k]; i++){
        clusters[k][i] = orig_clusters[k][i];
        alpha_hat[clusters[k][i]] = orig_alpha_hat[clusters[k][i]];
        cluster_assignment[clusters[k][i]] = orig_cluster_assignment[clusters[k][i]];
      }
      //log_like[k] = orig_log_like[k];
      //Omegay[k] = orig_Omegay[k];
      log_det_Omegay[k] = orig_log_det_Omegay[k];
      y_Omegay_y[k] = orig_y_Omegay_y[k];
      log_prior[k] = orig_log_prior[k];
      alpha_bar[k] = orig_alpha_bar[k];
    }
  }
  // clusters labeled split_k needs to be updated to contain only those elements in new_clusters[0]
  cluster_config[split_k] = new_clusters[0].size();
  clusters[split_k].resize(cluster_config[split_k], 0);
  for(int i = 0; i < cluster_config[split_k]; i++){
    clusters[split_k][i] = new_clusters[0][i];
    cluster_assignment[clusters[split_k][i]] = split_k;
  }
  //log_likelihood(split_k, ybar, T, A_block, rho, a);
  get_Omega(split_k, ybar, T, A_block, rho, a1, a2);
  log_pi_ep(split_k, eta);
  alpha_postmean(split_k, ybar, T, A_block, rho, a1, a2);
  // now loop over all of the other splits: start index at 1 because we already have dealt with new_clusters[0]
  for(int split_ix = 1; split_ix < num_splits; split_ix++){
    k = orig_K + split_ix - 1; // new cluster label
    cluster_config[k] = new_clusters[split_ix].size();
    clusters[k].resize(cluster_config[k], 0);
    for(int i = 0; i < cluster_config[k]; i++){
      clusters[k][i] = new_clusters[split_ix][i];
    }
    //log_likelihood(k, ybar, T, A_block, rho, a);
    get_Omega(k, ybar, T, A_block, rho, a1,a2);
    log_pi_ep(k, eta);
    alpha_postmean(k, ybar, T, A_block, rho, a1, a2);
  }
  get_pairwise();
}

void Partition::Merge(int k_1, int k_2, const arma::vec &ybar, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double eta)
{
  
  // make temporary copy of original attributes
  int orig_K = K;
  std::vector<int> orig_cluster_config(orig_K, 0);
  std::vector<std::vector<int> > orig_clusters(orig_K, std::vector<int>(1,0));
  std::vector<int> orig_cluster_assignment(nObs,0);
  //std::vector<double> orig_log_like(orig_K,0.0);
  //std::vector<arma::mat> orig_Omegay(K);
  std::vector<double> orig_log_det_Omegay(K);
  std::vector<double> orig_y_Omegay_y(K);
  std::vector<double> orig_log_prior(orig_K, 0.0);
  std::vector<double> orig_alpha_hat(nObs, 0.0);
  std::vector<double> orig_alpha_bar(K, 0.0);
  
  for(int k = 0; k < orig_K; k++){
    orig_cluster_config[k] = cluster_config[k];
    orig_clusters[k].clear();
    orig_clusters[k].resize(orig_cluster_config[k]);
    //Rcpp::Rcout << "orig_clusters[" << k << "]: " ;
    for(int i = 0; i < orig_cluster_config[k]; i++){
      orig_clusters[k][i] = clusters[k][i];
      orig_cluster_assignment[orig_clusters[k][i]] = cluster_assignment[orig_clusters[k][i]];
      //Rcpp::Rcout << orig_clusters[k][i] << "  " ;
      orig_alpha_hat[orig_clusters[k][i]] = alpha_hat[clusters[k][i]];
    }
    //Rcpp::Rcout << endl;
    //orig_log_like[k] = log_like[k];
    //orig_Omegay[k] = Omegay[k];;
    orig_log_det_Omegay[k] = log_det_Omegay[k];
    orig_y_Omegay_y[k] = y_Omegay_y[k];
    orig_log_prior[k] = log_prior[k];
    orig_alpha_bar[k] = alpha_bar[k];
    clusters[k].clear();
  }
  
  int k_max = std::max(k_1, k_2);
  int k_min = std::min(k_1, k_2);
  
  // re-size everything
  K = orig_K - 1;
  cluster_config.clear();
  cluster_config.resize(K,0);
  clusters.resize(K, std::vector<int>(1,0));
  //log_like.clear();
  //log_like.resize(K, 0.0);
  //Omegay.clear();
  //Omegay.resize(K);
  log_det_Omegay.clear();
  log_det_Omegay.resize(K);
  y_Omegay_y.clear();
  y_Omegay_y.resize(K);
  log_prior.clear();
  log_prior.resize(K, 0.0);
  alpha_bar.clear();
  alpha_bar.resize(K, 0.0);
  // loop over the OLD labels
  for(int k = 0; k < orig_K; k++){
    if(k == k_min){ // perform the merge
      //Rcpp::Rcout << "Performing the merge" << endl;
      cluster_config[k] = orig_cluster_config[k_min] + orig_cluster_config[k_max];
      clusters[k].resize(cluster_config[k]);
      for(int i = 0; i < orig_cluster_config[k_min]; i++){
        clusters[k][i] = orig_clusters[k_min][i];
      }
      for(int i = 0; i < orig_cluster_config[k_max]; i++){
        clusters[k][orig_cluster_config[k_min] + i] = orig_clusters[k_max][i];
      }
      //clusters[k].insert(clusters[k].end(), orig_clusters[k_min].begin(), orig_clusters[k_min].end());
      //clusters[k].insert(clusters[k].end(), orig_clusters[k_max].begin(), orig_clusters[k_max].end());
      //Rcpp::Rcout << orig_clusters[k_min].begin() << "   " << orig_clusters[k_min].end() << endl;
      
      
      for(int i = 0; i < cluster_config[k]; i++){
        cluster_assignment[clusters[k][i]] = k;
      }
      //log_likelihood(k, ybar, T, A_block, rho, a);
      get_Omega(k, ybar, T, A_block, rho, a1, a2);
      log_pi_ep(k, eta);
      alpha_postmean(k, ybar, T, A_block, rho, a1, a2);
      
    } else if(k < k_max){ // label does not change
      cluster_config[k] = orig_cluster_config[k];
      clusters[k].resize(cluster_config[k],0);
      for(int i = 0; i < cluster_config[k]; i++){
        clusters[k][i] = orig_clusters[k][i];
        alpha_hat[clusters[k][i]] = orig_alpha_hat[clusters[k][i]];
        cluster_assignment[clusters[k][i]] = orig_cluster_assignment[clusters[k][i]];
      }
      //log_like[k] = orig_log_like[k];
      //Omegay[k] = orig_Omegay[k];
      log_det_Omegay[k] = orig_log_det_Omegay[k];
      y_Omegay_y[k] = orig_y_Omegay_y[k];
      log_prior[k] = orig_log_prior[k];
      alpha_bar[k] = orig_alpha_bar[k];
    } else if(k > k_max){ // new label is k-1
      cluster_config[k-1] = orig_cluster_config[k];
      clusters[k-1].resize(cluster_config[k-1], 0);
      for(int i = 0; i < cluster_config[k-1]; i++){
        clusters[k-1][i] = orig_clusters[k][i];
        alpha_hat[clusters[k-1][i]] = orig_alpha_hat[clusters[k-1][i]];
        cluster_assignment[clusters[k-1][i]] = k-1;
      }
      //log_like[k-1] = orig_log_like[k];
      //Omegay[k-1] = orig_Omegay[k];
      log_det_Omegay[k-1] = orig_log_det_Omegay[k];
      y_Omegay_y[k-1] = orig_y_Omegay_y[k];
      log_prior[k-1] = orig_log_prior[k];
      alpha_bar[k-1] = orig_alpha_bar[k];
    }
    
  }
  get_pairwise();
}

void Partition::Split_Merge(int split_k, std::vector<std::vector<int> > &new_clusters, std::vector<int> k_star, const arma::vec &ybar, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double eta)
{
  if(new_clusters.size() != k_star.size()){
    Rcpp::Rcout << "[Split_Merge]: new_clusters and k_star must be of the same size!" << endl;
  } else{
    int num_new_clusters = new_clusters.size();
    int remain_ix = -1;
    int n_trailing = 0; // how many new clusters get tacked on to the end
    int k = 0;
    int counter = 0;
    int orig_K = K;
    int unik_flag = 1;
    std::vector<int> unik_k_star;
    std::vector<std::vector<int> > unik_k_star_map; // unik_k_star_map[u_ix] is a vector holding the indices nc_ix that get mapped to unik_k_star[u_ix];
    std::vector<int> unchanged_k(K,1); //unchanged_k[k] = 1 means that it is unchanged by the split+merge moves
    unchanged_k[split_k] = 0; // clearly the cluster being split will be changed by our Split_Merge
    
    for(int nc_ix = 0; nc_ix < num_new_clusters; nc_ix++){ // check to see that all of the k_star's are valid. also find remain_ix
      if(k_star[nc_ix] > 0 & k_star[nc_ix] >= K){
        Rcpp::Rcout << "[Split_Merge]: k_star must be between -1 and K-1 (inclusive)" << endl;
      }
      if(remain_ix == -1 & k_star[nc_ix] == -1){
        remain_ix = nc_ix;
        //Rcpp::Rcout << "[Split_Merge]: Found remain_ix = " << remain_ix << endl;
      }
      // scan unique k_star to see if nearest neighbor of new_cluster[nc_ix] is someone else's nearest neighbor
      unik_flag = 1;
      for(int u_ix = 0; u_ix < unik_k_star.size(); u_ix++){
        if(k_star[nc_ix] == unik_k_star[u_ix]){ // nearest neighbor of new cluster nc_ix is a nearest neighbor of another new cluster!
          unik_flag = 0;
          if(k_star[nc_ix] != -1){
            unik_k_star_map[u_ix].push_back(nc_ix);
          } else if( (k_star[nc_ix] == -1) & (nc_ix != remain_ix)){
            unik_k_star_map[u_ix].push_back(nc_ix);
            n_trailing++;
          }
          break;
        }
      }
      if(unik_flag == 1){ // we have a new k_star!
        if(k_star[nc_ix] == -1){
          // if nc_ix = remain_ix, we don't want to add -1 to the set of unique k_star's yet.
          // This is to handle the situation where there are no trailing clusters
          if(nc_ix != remain_ix){
            unik_k_star.push_back(k_star[nc_ix]); // add k_star[nc_ix] to running list of unique k_stars
            unik_k_star_map.push_back(std::vector<int>(1, nc_ix));
            n_trailing++;
          }
        } else{ // we have a unique k_star that is non-negative
          unik_k_star.push_back(k_star[nc_ix]); // add k_star[nc_ix] to running list of unique k_stars
          unik_k_star_map.push_back(std::vector<int>(1,nc_ix)); // add nc_ix to the running list of new cluster indices that get merged to unik_k_star[u_ix]
          // if k_star[nc_ix] > -1, original cluster with that label will be changed by the split-merge
          unchanged_k[k_star[nc_ix]] = 0;
        }
      }
    } // closes loop over the new cluster indices
    
    
    // some print statements here to check our progress
    /*
     Rcpp::Rcout << "remain_ix = " << remain_ix << endl;
     Rcpp::Rcout << "num unique k_star = " << unik_k_star.size() << endl;
     for(int u_ix = 0; u_ix < unik_k_star.size(); u_ix++){
     Rcpp::Rcout << "new_clusters being added to cluster labelled " << unik_k_star[u_ix] << " : " ;
     for(int nc_ix = 0; nc_ix < unik_k_star_map[u_ix].size(); nc_ix++){
     Rcpp::Rcout << unik_k_star_map[u_ix][nc_ix] << " ";
     }
     Rcpp::Rcout << endl;
     }
     Rcpp::Rcout << "n_trailing = "<< n_trailing << endl;
     Rcpp::Rcout << "Existing clusters that are unchanged: " ;
     for(int kk = 0; kk < orig_K; kk ++){
     if(unchanged_k[kk] == 1) Rcpp::Rcout << kk << " ";
     }
     Rcpp::Rcout << endl;
     */
    // Now ready to perform the merges //
    // first create copies
    std::vector<int> orig_cluster_config(orig_K,0);
    std::vector<std::vector<int> > orig_clusters(orig_K, std::vector<int>(1,0));
    std::vector<int> orig_cluster_assignment(nObs, 0);
    std::vector<double> orig_alpha_hat(nObs,0);
    //std::vector<double> orig_log_like(orig_K,0.0);
    //std::vector<arma::mat> orig_Omegay(K);
    std::vector<double> orig_log_det_Omegay(orig_K);
    std::vector<double> orig_y_Omegay_y(orig_K);
    
    std::vector<double> orig_log_prior(orig_K, 0.0);
    std::vector<double> orig_alpha_bar(orig_K, 0.0);
    for(int kk = 0; kk < orig_K; kk++){
      orig_cluster_config[kk] = cluster_config[kk];
      orig_clusters[kk].resize(orig_cluster_config[kk],0);
      for(int i = 0; i < orig_cluster_config[kk]; i++){
        orig_clusters[kk][i] = clusters[kk][i];
      }
      orig_log_det_Omegay[kk] = log_det_Omegay[kk];
      orig_y_Omegay_y[kk] = y_Omegay_y[kk];
      orig_log_prior[kk] = log_prior[kk];
      orig_alpha_bar[kk] = alpha_bar[kk];
    }
    for(int i = 0; i < nObs; i++){
      orig_alpha_hat[i] = alpha_hat[i];
      //orig_cluster_assignment[i] = orig_cluster_assignment[i];
      orig_cluster_assignment[i] = cluster_assignment[i];
    }
    //Rcpp::Rcout << "[Split_Merge]: Made copies of original attributes" << endl;
    // clearing all of the vectors here would be a good idea probably
    
    
    if(remain_ix != -1){ // one of the new clusters will remain labelled split_k
      // Re-initialize all of the attributes //
      //Rcpp::Rcout << "[Split_Merge]: remain_ix != -1. One new cluster will remain labelled split_k" << endl;
      
      K = orig_K + n_trailing;
      cluster_config.clear();
      cluster_config.resize(K, 0);
      clusters.clear();
      clusters.resize(K, std::vector<int>(1,0));
      cluster_assignment.clear();
      cluster_assignment.resize(nObs, -1);
      alpha_hat.clear();
      alpha_hat.resize(nObs, 0.0);
      
      //Omegay.clear();
      //Omegay.resize(K);
      log_det_Omegay.clear();
      log_det_Omegay.resize(K);
      y_Omegay_y.clear();
      y_Omegay_y.resize(K);
      //log_like.resize(K, 0.0);
      log_prior.clear();
      log_prior.resize(K,0.0);
      alpha_bar.clear();
      alpha_bar.resize(K,0.0);
      // update all of the unchanged clusters //
      //Rcpp::Rcout << "    Updating unchanged clusters: " ;
      
      for(int kk = 0; kk < orig_K; kk++){
        if(unchanged_k[kk] == 1){
          //Rcpp::Rcout << kk << " " ;
          cluster_config[kk] = orig_cluster_config[kk];
          clusters[kk].resize(cluster_config[kk],0);
          for(int i = 0; i < cluster_config[kk]; i++){
            clusters[kk][i] = orig_clusters[kk][i];
            cluster_assignment[clusters[kk][i]] = orig_cluster_assignment[clusters[kk][i]];
            alpha_hat[clusters[kk][i]] = orig_alpha_hat[clusters[kk][i]];
          }
          //log_like[kk] = orig_log_like[kk];
          //Omegay[kk] = orig_Omegay[kk];
          log_det_Omegay[kk] = orig_log_det_Omegay[kk];
          y_Omegay_y[kk] = orig_y_Omegay_y[kk];
          log_prior[kk] = orig_log_prior[kk];
          alpha_bar[kk] = orig_alpha_bar[kk];
        }
      } // closes loop that updates the unchanged clusters
      //Rcpp::Rcout << endl;
      // now update split_k //
      
      cluster_config[split_k] = new_clusters[remain_ix].size(); // remember remain_ix is the index of the new cluster that remains at k-1
      clusters[split_k].resize(cluster_config[split_k]);
      for(int i = 0; i < new_clusters[remain_ix].size(); i++){
        clusters[split_k][i] = new_clusters[remain_ix][i];
        cluster_assignment[clusters[split_k][i]] = split_k;
      }
      //log_likelihood(split_k, ybar, T, A_block, rho, a);
      get_Omega(split_k, ybar, T, A_block, rho, a1,a2);
      log_pi_ep(split_k, eta);
      alpha_postmean(split_k, ybar, T, A_block, rho, a1, a2);
      //Rcpp::Rcout << "    Updating cluster originally labelled split_k" << endl;
      
      // now go through unik_k_star's and update these //
      for(int u_ix = 0; u_ix < unik_k_star.size(); u_ix++){ // loop over the indices of the unique k-star's
        k = unik_k_star[u_ix];
        if(k == -1){// this designates a trailing cluster
          //Rcpp::Rcout << "    Creating trailing clusters " ;
          // all trailing clusters are updated now.
          for(int nc_ix = 0; nc_ix < n_trailing; nc_ix++){
            k = orig_K + nc_ix; // new label. this cluster label did not exist before the split/merge
            //Rcpp::Rcout << k << " ";
            cluster_config[k] = new_clusters[unik_k_star_map[u_ix][nc_ix]].size();
            clusters[k].clear();
            clusters[k].resize(cluster_config[k]);
            for(int i = 0; i < cluster_config[k]; i++){
              clusters[k][i] = new_clusters[unik_k_star_map[u_ix][nc_ix]][i];
              cluster_assignment[clusters[k][i]] = k;
            }
            // need to update the log-likelihood
            //log_likelihood(k, ybar, T, A_block, rho, a);
            get_Omega(k, ybar, T, A_block, rho, a1, a2);
            log_pi_ep(k, eta);
            alpha_postmean(k, ybar, T, A_block, rho, a1, a2);
          } // closes loop over the trailing clusters
          //Rcpp::Rcout << endl;
        } else{ // k is an existing label
          //Rcpp::Rcout << "    Updating cluster originally labelled" << k << endl;
          cluster_config[k] = orig_cluster_config[k];
          // loop over the indices of the new clusters being added to k, incrementing cluster_config
          for(int nc_ix = 0; nc_ix < unik_k_star_map[u_ix].size(); nc_ix++){ // loop over the indices of new clusters being added to k
            cluster_config[k] += new_clusters[unik_k_star_map[u_ix][nc_ix]].size();
          }
          clusters[k].clear();
          clusters[k].resize(cluster_config[k],0);
          // add the elements from the original cluster k back to it
          counter = 0;
          for(int i = 0; i < orig_cluster_config[k]; i++){
            clusters[k][counter] = orig_clusters[k][i];
            cluster_assignment[clusters[k][counter]] = k;
            counter++;
          }
          // loop over indices of new clusters being added to k and then add the elements of each new_cluster to k
          //Rcpp::Rcout << "    Adding new_clusters " ;
          for(int nc_ix = 0; nc_ix < unik_k_star_map[u_ix].size(); nc_ix++){
            //Rcpp::Rcout << unik_k_star_map[u_ix][nc_ix] << " " ;
            for(int i = 0; i < new_clusters[unik_k_star_map[u_ix][nc_ix]].size(); i++){
              clusters[k][counter] = new_clusters[unik_k_star_map[u_ix][nc_ix]][i];
              cluster_assignment[clusters[k][counter]] = k;
              counter++;
            }
          } // closes loop that adds elements of new_clusters to cluster originally labelled k
          //Rcpp::Rcout << " to cluster originally labelled " << k << endl;
          //log_likelihood(k, ybar, T, A_block, rho, a);
          get_Omega(k, ybar, T, A_block, rho, a1, a2);
          log_pi_ep(k, eta);
          alpha_postmean(k, ybar, T, A_block, rho, a1, a2);
        } // closes if/else checking whether unik_k_star[u_ix] == -1 (i.e. make the trailing clusters) or not
      } // closes loop over the unik_k_stars (i.e. the labels of the original clusters that are changed by this Split/Merge
      //Rcpp::Rcout << "    Updated the clusters changed by the Split_Merge" << endl;
      
    } else{ // all new clusters will be merged to an existing cluster
      //Rcpp::Rcout << "[Split_Merge]: remain_ix == -1." << endl;
      
      K = orig_K - 1;
      // re-initialize all of the attributes //
      cluster_config.clear();
      cluster_config.resize(K, 0);
      clusters.clear();
      clusters.resize(K, std::vector<int>(1,0));
      cluster_assignment.clear();
      cluster_assignment.resize(nObs,-1);
      alpha_hat.clear();
      alpha_hat.resize(nObs, 0.0);
      
      //log_like.resize(K,0.0);
      //Omegay.clear();
      //Omegay.resize(K);
      log_det_Omegay.clear();
      log_det_Omegay.resize(K);
      
      y_Omegay_y.clear();
      y_Omegay_y.resize(K);
      log_prior.clear();
      log_prior.resize(K,0.0);
      alpha_bar.clear();
      alpha_bar.resize(K,0.0);
      //Rcpp::Rcout << "    attributes re-initialized" << endl;
      
      // loop over the original labels
      //Rcpp::Rcout << "   Updating clusters left unchanged by Split_Merge: " << endl;
      for(int kk = 0; kk < orig_K; kk++){
        //Rcpp::Rcout << "    kk = " << kk ;
        if(unchanged_k[kk] == 1){ // cluster contents not changed by split+merge
          if(kk < split_k){ // cluster label not changed
            //Rcpp::Rcout << "    kk = " << kk << "  original label left unchanged" << endl;
            cluster_config[kk] = orig_cluster_config[kk];
            clusters[kk].clear();
            clusters[kk].resize(cluster_config[kk], 0);
            for(int i = 0; i < cluster_config[kk]; i++){
              clusters[kk][i] = orig_clusters[kk][i];
              cluster_assignment[clusters[kk][i]] = kk;
              alpha_hat[clusters[kk][i]] = orig_alpha_hat[clusters[kk][i]];
            }
            //log_like[kk]  = orig_log_like[kk];
            //Omegay[kk] = orig_Omegay[kk];
            log_det_Omegay[kk] = orig_log_det_Omegay[kk];
            y_Omegay_y[kk] = orig_y_Omegay_y[kk];
            log_prior[kk] = orig_log_prior[kk];
            alpha_bar[kk] = orig_alpha_bar[kk];
          } else if(kk > split_k){ // we need to subtract 1 from the cluster label!
            //Rcpp::Rcout << "    kk =  "<< kk << "    need to adjust the label" << endl;
            cluster_config[kk-1] = orig_cluster_config[kk];
            clusters[kk-1].clear();
            clusters[kk-1].resize(cluster_config[kk-1],0);
            // if I were to guess, this loop causes some problems:
            //   we should loop over i until cluster_config[kk-1] not cluster_config[kk], which presumably is still 0 since
            // the outer loop hasn't hit kk yet, and cluster_config[kk]
            // confirmation: when I tried running it again, I noticed that the loop never entered i.
            /*
             for(int i = 0; i < cluster_config[kk]; i++){
             Rcpp::Rcout << "        i = " << i << endl;
             clusters[kk-1][i] = orig_clusters[kk][i];
             cluster_assignment[clusters[kk-1][i]] = kk-1;
             alpha_hat[clusters[kk-1][i]] = orig_alpha_hat[clusters[kk-1][i]];
             }
             */
            for(int i = 0; i < cluster_config[kk-1];i++){
              clusters[kk-1][i] = orig_clusters[kk][i];
              cluster_assignment[clusters[kk-1][i]] = kk-1;
              alpha_hat[clusters[kk-1][i]] = orig_alpha_hat[clusters[kk-1][i]];
            }
            //log_like[kk-1] = orig_log_like[kk];
            //Omegay[kk-1] = orig_Omegay[kk];
            log_det_Omegay[kk-1] = orig_log_det_Omegay[kk];
            y_Omegay_y[kk-1] = orig_y_Omegay_y[kk];
            log_prior[kk-1] = orig_log_prior[kk];
            alpha_bar[kk-1] = orig_alpha_bar[kk];
          }
        } // closes if statement checking whether cluster labelled kk is being modified
      } // closes loop over the cluster labels. only the unchanged ones have been copied
      //Rcpp::Rcout << "    Updating clusters changed by Split_Merge" << endl;
      
      for(int u_ix = 0; u_ix < unik_k_star.size(); u_ix++){
        k = unik_k_star[u_ix];
        
        if(k < split_k){ // label does not change but contents do!
          //Rcpp::Rcout << "    k = " << k << "  label does not change" << endl;
          cluster_config[k] = orig_cluster_config[k];
          for(int nc_ix = 0; nc_ix < unik_k_star_map[u_ix].size(); nc_ix++){ // loop over the indices of new clusters being added to k
            cluster_config[k] += new_clusters[unik_k_star_map[u_ix][nc_ix]].size();
          }
          clusters[k].clear();
          clusters[k].resize(cluster_config[k],0);
          // add the elements from the original cluster k back to it
          counter = 0;
          for(int i = 0; i < orig_cluster_config[k]; i++){
            clusters[k][counter] = orig_clusters[k][i];
            cluster_assignment[clusters[k][counter]] = k;
            counter++;
          }
          // loop over indices of new clusters being added to k
          for(int nc_ix = 0; nc_ix < unik_k_star_map[u_ix].size(); nc_ix++){
            for(int i = 0; i < new_clusters[unik_k_star_map[u_ix][nc_ix]].size(); i++){
              clusters[k][counter] = new_clusters[unik_k_star_map[u_ix][nc_ix]][i];
              cluster_assignment[clusters[k][counter]] = k;
              counter++;
            }
          } // closes loop that adds elements of new_clusters to cluster originally labelled k
          //log_likelihood(k, ybar, T, A_block, rho, a);
          get_Omega(k, ybar, T, A_block, rho, a1, a2);
          log_pi_ep(k, eta);
          alpha_postmean(k, ybar, T, A_block, rho, a1, a2);
        } else if(k > split_k){
          //Rcpp::Rcout << "     k = " << k << " label changes" << endl;
          cluster_config[k-1] = orig_cluster_config[k];
          for(int nc_ix = 0; nc_ix < unik_k_star_map[u_ix].size(); nc_ix++){// loop over indices of new clusters being added to what was labelled k
            cluster_config[k-1] += new_clusters[unik_k_star_map[u_ix][nc_ix]].size();
          }
          clusters[k-1].resize(cluster_config[k-1],0);
          counter = 0;
          // add the elements from the original cluster labelled k (remember the new label is now k-1)
          for(int i = 0; i < orig_cluster_config[k]; i++){
            clusters[k-1][counter] = orig_clusters[k][i];
            cluster_assignment[clusters[k-1][counter]] = k-1;
            counter++;
          }
          // add in the elements from the new clusters
          for(int nc_ix = 0; nc_ix < unik_k_star_map[u_ix].size(); nc_ix++){
            for(int i = 0; i < new_clusters[unik_k_star_map[u_ix][nc_ix]].size(); i++){
              clusters[k-1][counter] = new_clusters[unik_k_star_map[u_ix][nc_ix]][i];
              cluster_assignment[clusters[k-1][counter]] = k-1;
              counter++;
            }
          } // closes loop that adds elements of new_clusters to cluster originally labelled k (but is now labelled k-1).
          
          //log_likelihood(k-1, ybar, T, A_block, rho, a);
          get_Omega(k-1, ybar, T, A_block, rho, a1, a2);
          log_pi_ep(k-1, eta);
          alpha_postmean(k-1, ybar, T, A_block, rho, a1, a2);
          
        } // closes if/else checking that k < split_k or k > split_k
      } // closes loop over the unique k_stars
    } // closes if/else checking whether any new cluster gets to be labelled split_k
    
    
    // [SKD]: 14 January update. For some reason the y_Omegay_y and log_det_Omegay attributes are not updating correctly
    // This is especially the case with clusters that are supposedly left unchanged by the split/merge operations
    for(int k = 0; k < K; k++){
      if( (y_Omegay_y[k] == 0) || (log_det_Omegay[k] == 0)){
        Rcpp::Rcout << "log_det_Omegay or y_Omega_y is zero" << endl;
        get_Omega(k, ybar, T, A_block, rho, a1, a2);
      }
    }
    
    // need to update the pairwise_allocation matrix!
    get_pairwise();
    
    
  } // closes if/else checking whether new_clusters and k_star have the same size
}


