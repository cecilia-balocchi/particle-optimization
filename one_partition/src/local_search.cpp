//
//  local_search.cpp
//  Given an initial partition, form each island candidate and return the log-posterior
//  of each.
//

#include <stdio.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "partition.h"
#include "partition_functions.h"
#include <vector>

using namespace std;

// [[Rcpp::export]]
Rcpp::List local_search(arma::mat Y,
                        const arma::mat A_block,
                        Rcpp::List gamma_init,
                        const double a1 = 1.0,
                        const double a2 = 1.0,
                        const double nu_sigma = 3,
                        const double lambda_sigma = 1,
                        const double rho = 0.99,
                        const double eta = 1.0,
                        const double eps = 1e-3)
{
  int n = Y.n_rows;
  int T = Y.n_cols;
  
  arma::vec ybar(n);
  double total_ss = 0;
  for(int i = 0; i < n; i++){
    ybar(i) = arma::mean(Y.row(i));
    total_ss += (T-1) * arma::var(Y.row(i));
  }
  
  Rcpp::Rcout << "n = " << n << endl;
  LPPartition gamma_0 = new Partition(n, gamma_init, ybar, T, A_block, rho, a1, a2, eta);
  Rcpp::Rcout << "Created gamma_0" << endl;
  
  split_info local_si;
  get_local(local_si, gamma_0, T, A_block, rho, a1, a2);
  /*
  Rcpp::Rcout << "got " << local_si.num_splits << " local candidates" << endl;
  for(int split_ix = 0; split_ix < local_si.num_splits; split_ix++){
    Rcpp::Rcout << "creating " << local_si.new_clusters[split_ix].size() << " new clusters" << endl;
    for(int nc_ix = 0; nc_ix < local_si.new_clusters[split_ix].size(); nc_ix++) Rcpp::Rcout << " size " << local_si.new_clusters[split_ix][nc_ix].size() << " and neighbor = " << local_si.nearest_neighbor[split_ix][nc_ix] ;
    Rcpp::Rcout << endl;
  }
  */

  std::vector<int> k_star(n,-1);
  int L = n+1;
  std::vector<LPPartition> particle_set(L);
  std::vector<double> w(L);
  particle_set[0] = new Partition(gamma_0);
  for(int i = 1; i < L; i++){
    //Rcpp::Rcout << "i = " << i;
    particle_set[i] = new Partition(gamma_0);
    k_star.clear();
    k_star.resize(local_si.new_clusters[i-1].size());
    for(int nc_ix = 0; nc_ix < local_si.new_clusters[i-1].size(); nc_ix++) k_star[nc_ix] = -1;
    //Rcpp::Rcout << local_si.new_clusters[i-1].size() << " " << k_star.size() << endl;
    
    particle_set[i]->Split_Merge(local_si.split_k[i-1], local_si.new_clusters[i-1], k_star, ybar, T, A_block, rho, a1, a2, eta);
    //Rcpp::Rcout << "  log-post = " << total_log_post(particle_set[i], nu_sigma, lambda_sigma) << endl;
    
    w[i] = 1/(n+1);
  }
  Rcpp::Rcout << "Finished main loop" << endl;
  // Run update w.
  update_w(particle_set, w, L, total_ss, T, nu_sigma, lambda_sigma, 1.0); // use lambda = 1 so that w is the re-weighted log-posterior

  // Find unique particles. For this particular function, this is slightly overkill since we know there are n+1 unique particles
  // However, it will sort everything, which is nice
  std::vector<std::vector<int> > particle_map;
  std::vector<double> pstar;
  std::vector<int> counts;
  get_unik_particles(particle_map, pstar, counts, particle_set, w);
  
  int L_star = particle_map.size();
  std::vector<LPPartition> unik_particles(L_star);
  std::vector<double> log_like(L_star);
  std::vector<double> log_prior(L_star);
  std::vector<double> log_post(L_star);
  arma::mat alpha_out(n, L_star);
  for(int l = 0; l < L_star; l++){
    unik_particles[l] = new Partition(particle_set[particle_map[l][0]]);
    log_like[l] = total_log_like(unik_particles[l], total_ss, T, nu_sigma, lambda_sigma);
    log_prior[l] = total_log_prior(unik_particles[l]);
    log_post[l] = total_log_post(unik_particles[l], total_ss, T, nu_sigma, lambda_sigma);
    for(int i = 0; i < n; i++) alpha_out(i,l) = unik_particles[l]->alpha_hat[i];
  }
  Rcpp::List unik_particles_out;
  format_particle_set(unik_particles, unik_particles_out); // format the particle set so that it can returned as an R list
  Rcpp::List results;
  results["particles"] = unik_particles_out;
  results["pstar"] = pstar;
  results["counts"] = counts;
  results["log_like"] = log_like;
  results["log_prior"] = log_prior;
  results["log_post"] = log_post;
  results["alpha"] = alpha_out;
  return results;
}
