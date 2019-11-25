//
//  map_cluster.cpp
//  
//
//

#include <stdio.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "partition.h"
#include "partition_functions.h"
#include "initialize_particle_set.h"
#include "rng.h"
#include <vector>
#include <ctime>

using namespace std;

// [[Rcpp::export]]
Rcpp::List map_partition(arma::mat Y,
                         const arma::mat A_block,
                         const double a1 = 1.0,
                         const double a2 = 1.0,
                         const double nu_sigma = 3,
                         const double lambda_sigma = 1,
                         const double rho = 0.99,
                         const double eta = 1.0,
                         const int max_iter = 10, const double eps = 1e-3,
                         const double split_frac = 0.1,
                         bool verbose = false)
{
  Rcpp::RNGScope scope;
  RNG gen;
  int n = Y.n_rows;
  int T = Y.n_cols;
  
  arma::vec ybar(n);
  double total_ss = 0;
  for(int i = 0; i < n; i++){
    ybar(i) = arma::mean(Y.row(i));
    total_ss += (T-1) * arma::var(Y.row(i));
  }
  LPPartition gamma_0 = new Partition();
  initialize_particle(gamma_0, ybar, total_ss, T, A_block, rho, a1, a2, eta, nu_sigma, lambda_sigma, gen);
  
  if(verbose == true) Rcpp::Rcout << "n = " << n << endl;
  if(verbose == true) Rcpp::Rcout << "Created gamma_0" << endl;
  
  // for the main loop, we need the following quantities
  LPPartition spec_split_candidate = new Partition(gamma_0); // for spectral splits
  LPPartition tail_split_candidate = new Partition(gamma_0); // for tail splits
  LPPartition km_split_candidate = new Partition(gamma_0); // for km splits
  LPPartition merge_candidate = new Partition(gamma_0); // for merges
  LPPartition border_candidate = new Partition(gamma_0); // for border candidate
  LPPartition island_candidate = new Partition(gamma_0); // for island candidates
  LPPartition local_candidate = new Partition(gamma_0); // for local candidates
  
  split_info spec_si;
  split_info tail_si;
  split_info km_si;
  split_info bi;
  split_info isl_i;
  merge_info mi;
  split_info local_si;
  
  double spec_split_obj = 0.0;
  double tail_split_obj = 0.0;
  double km_split_obj = 0.0;
  double merge_obj = 0.0;
  double border_obj = 0.0;
  double island_obj = 0.0;
  double local_obj = 0.0;
  double accepted_obj = 0.0;
  
  int spec_split_flag = 1;
  int tail_split_flag = 1;
  int km_split_flag = 1;
  int merge_flag = 1;
  int border_flag = 1;
  int island_flag = 1;
  int local_flag = 1;
  
  
  std::vector<LPPartition> particle_set(1);
  particle_set[0] = new Partition(gamma_0);
  
  int iter = 0;
  double old_objective = 0.0;
  double objective = total_log_post(particle_set[0], total_ss, T, nu_sigma, lambda_sigma);
  int flag = 0;

  bool try_local = false;
  
  std::vector<Rcpp::List> particle_set_trajectory;
  Rcpp::List tmp_list;
  format_particle_set(particle_set, tmp_list);
  particle_set_trajectory.push_back(tmp_list);
  
  std::vector<double> objective_trajectory;
  std::vector<double> log_like_trajectory;
  std::vector<double> log_prior_trajectory;
  std::vector<double> log_post_trajectory;
  std::vector<double> tmp_alpha_hat(n);
  std::vector<std::vector<double> > alpha_trajectory;
  
  objective_trajectory.push_back(objective);
  log_like_trajectory.push_back(total_log_like(particle_set[0], total_ss, T, nu_sigma, lambda_sigma));
  log_prior_trajectory.push_back(total_log_prior(particle_set[0]));
  log_post_trajectory.push_back(total_log_post(particle_set[0], total_ss, T, nu_sigma, lambda_sigma));
  for(int i = 0; i < n; i++) tmp_alpha_hat[i] = particle_set[0]->alpha_hat[i];
  alpha_trajectory.push_back(tmp_alpha_hat);

  time_t tp;
  int time1 = time(&tp);
  while( (iter < max_iter) & (flag == 0)){
    if(verbose == true) Rcpp::Rcout << "Starting iter = " << iter << endl;
    old_objective = objective;
    Rcpp::checkUserInterrupt();
    spec_split_flag = 1;
    tail_split_flag = 1;
    km_split_flag = 1;
    border_flag = 1;
    merge_flag = 1;
    island_flag = 1;
    local_flag = 1;
    try_local = true;
    
    // Spectral Split
    get_spectral_split(spec_si, particle_set[0], T, A_block, rho, a1, a2, 500);
    delete spec_split_candidate;
    spec_split_candidate = new Partition(particle_set[0]);
    best_split_map(spec_si, spec_split_candidate, particle_set[0], ybar, total_ss, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, true);
    spec_split_obj = total_log_post(spec_split_candidate, total_ss, T, nu_sigma, lambda_sigma);
    spec_split_flag = Partition_Equal(spec_split_candidate, particle_set[0]);
    if(spec_split_flag == 0) try_local = false; // there is a spectral split that improves objective. don't try all local moves
    accepted_obj = spec_split_obj;
    
    // Tail splits
    get_tail_split(tail_si, particle_set[0], T, A_block, rho, a1, a2, 0.025);
    delete tail_split_candidate;
    tail_split_candidate = new Partition(particle_set[0]);
    best_split_map(tail_si, tail_split_candidate, particle_set[0], ybar, total_ss, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, true);
    tail_split_obj = total_log_post(tail_split_candidate, total_ss, T, nu_sigma, lambda_sigma);
    tail_split_flag = Partition_Equal(tail_split_candidate, particle_set[0]);
    if(tail_split_flag == 0) try_local = false; // there is a tail split that improves objective. don't try all local moves
    if(tail_split_obj > accepted_obj) accepted_obj = tail_split_obj;
    
    // KM splits
    get_km_split(km_si, particle_set[0], T, A_block, rho, a1, a2, 1000);
    delete km_split_candidate;
    km_split_candidate = new Partition(particle_set[0]);
    best_split_map(km_si, km_split_candidate, particle_set[0], ybar, total_ss, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, true);
    km_split_obj = total_log_post(km_split_candidate, total_ss, T, nu_sigma, lambda_sigma);
    km_split_flag = Partition_Equal(km_split_candidate, particle_set[0]);
    if(km_split_flag == 0) try_local = false;
    if(km_split_obj > accepted_obj) accepted_obj = km_split_obj;
    
    // Merges
    get_merge(mi, particle_set[0], A_block);
    delete merge_candidate;
    merge_candidate = new Partition(particle_set[0]);
    best_merge_map(mi, merge_candidate, particle_set[0], ybar, total_ss, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta);
    merge_obj = total_log_post(merge_candidate, total_ss, T, nu_sigma, lambda_sigma);
    merge_flag = Partition_Equal(merge_candidate, particle_set[0]);
    if(merge_flag == 0) try_local = false;
    if(merge_obj > accepted_obj) accepted_obj = merge_obj;
    
    // Border
    get_border(bi, particle_set[0], T, A_block, rho, a1, a2);
    delete border_candidate;
    border_candidate = new Partition(particle_set[0]);
    best_split_map(bi, border_candidate, particle_set[0], ybar, total_ss, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, false);
    border_obj = total_log_post(border_candidate, total_ss, T, nu_sigma, lambda_sigma);
    border_flag = Partition_Equal(border_candidate, particle_set[0]);
    if(border_flag == 0) try_local = false;
    if(border_obj > accepted_obj) accepted_obj = border_obj;
    
    // Island
    get_island(isl_i, particle_set[0], T, A_block, rho, a1, a2, 0.05);
    delete island_candidate;
    island_candidate = new Partition(particle_set[0]);
    best_split_map(isl_i, island_candidate, particle_set[0], ybar, total_ss, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, false);
    island_obj = total_log_post(island_candidate, total_ss, T, nu_sigma, lambda_sigma);
    island_flag = Partition_Equal(island_candidate, particle_set[0]);
    if(island_flag == 0) try_local = false;
    if(island_obj > accepted_obj) accepted_obj = island_obj;
    
    if(try_local == true){
      if(verbose == true) Rcpp::Rcout << "Trying local moves" << endl;
      get_local(local_si, particle_set[0], T, A_block, rho, a1, a2);
      delete local_candidate;
      local_candidate = new Partition(particle_set[0]);
      best_split_map(local_si, local_candidate, particle_set[0], ybar, total_ss, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, false);
      local_obj = total_log_post(local_candidate, total_ss, T, nu_sigma, lambda_sigma);
      local_flag = Partition_Equal(local_candidate, particle_set[0]);
      if(local_obj > accepted_obj) accepted_obj = local_obj;
    }
    
    
    /*
    Rcpp::Rcout << "spec_split_obj = " << spec_split_obj << " spec_split_flag = " << spec_split_flag << endl;
    Rcpp::Rcout << "tail_split_obj = " << tail_split_obj << " tail_split_flag = " << tail_split_flag << endl;
    Rcpp::Rcout << "km_split_obj = " << km_split_obj << " km_split_flag = " << km_split_flag << endl;
    Rcpp::Rcout << "merge_obj = " << merge_obj << " merge_flag = " << merge_flag << endl;
    Rcpp::Rcout << "border_obj = " << border_obj << " border_flag = " << border_flag << endl;
    Rcpp::Rcout << "island_obj = " << island_obj << " island_flag = " << island_flag << endl;
    Rcpp::Rcout << "accepted_obj = " << accepted_obj << endl;
    */
    
    if(accepted_obj == spec_split_obj){
      //Rcpp::Rcout << "accepted spec split" << endl;
      if(spec_split_flag == 0) particle_set[0]->Copy_Partition(spec_split_candidate);
      else flag = 1; // best move is to not move at all
    } else if(accepted_obj == tail_split_obj){
      //Rcpp::Rcout << "accepted tail split" << endl;
      if(tail_split_flag == 0) particle_set[0]->Copy_Partition(tail_split_candidate);
      else flag = 1;
    } else if(accepted_obj == km_split_obj){
      //Rcpp::Rcout << "accepted km split" << endl;
      if(km_split_flag == 0) particle_set[0]->Copy_Partition(km_split_candidate);
      else flag = 1;
    } else if(accepted_obj == merge_obj){
      //Rcpp::Rcout << "accepted merge" << endl;
      if(merge_flag == 0) particle_set[0]->Copy_Partition(merge_candidate);
      else flag = 1;
    } else if(accepted_obj == border_obj){
      //Rcpp::Rcout << "accepted border" << endl;
      if(border_flag == 0) particle_set[0]->Copy_Partition(border_candidate);
      else flag = 1;
    } else if(accepted_obj == island_obj){
      //Rcpp::Rcout << "accepted island" << endl;
      if(island_flag == 0) particle_set[0]->Copy_Partition(island_candidate);
      else flag = 1;
    } else if( (try_local == true) && (accepted_obj == local_obj)){
      if(local_flag == 0) particle_set[0]->Copy_Partition(local_candidate);
      else flag = 1;
    }
    //Rcpp::Rcout << "flag = " << flag << endl;
    objective = total_log_post(particle_set[0], total_ss, T, nu_sigma, lambda_sigma);
    if(verbose == true) Rcpp::Rcout << "objective = " << objective << "  old_objective = " << old_objective << " % diff = " << 100.0 * abs( (objective - old_objective)/objective) << endl;
    iter++;
    
    
    // now we can update the trajectories
    format_particle_set(particle_set, tmp_list);
    particle_set_trajectory.push_back(tmp_list);
    objective_trajectory.push_back(objective);
    log_like_trajectory.push_back(total_log_like(particle_set[0], total_ss, T, nu_sigma, lambda_sigma));
    log_prior_trajectory.push_back(total_log_prior(particle_set[0]));
    log_post_trajectory.push_back(total_log_post(particle_set[0], total_ss, T, nu_sigma, lambda_sigma));
    for(int i = 0; i < n; i++) tmp_alpha_hat[i] = particle_set[0]->alpha_hat[i];
    alpha_trajectory.push_back(tmp_alpha_hat);
  } // closes while loop
  int time2 = time(&tp);
  
  Rcpp::List particle_set_trajectory_out(particle_set_trajectory.size());
  arma::mat alpha_trajectory_out(n, alpha_trajectory.size());
  for(int ix = 0; ix < particle_set_trajectory.size(); ix++){
    particle_set_trajectory_out[ix] = particle_set_trajectory[ix];
    for(int i = 0; i < n; i++) alpha_trajectory_out(i,ix) = alpha_trajectory[ix][i];
  }
  
  Rcpp::List results;
  results["particles"] = particle_set_trajectory[particle_set_trajectory.size()-1]; // the map estimate
  results["pstar"] = 1.0;
  results["counts"] = 1;
  results["log_like"] = total_log_like(particle_set[0], total_ss, T, nu_sigma, lambda_sigma);
  results["log_prior"] = total_log_prior(particle_set[0]);
  results["log_post"] = total_log_post(particle_set[0], total_ss, T, nu_sigma, lambda_sigma);
  results["alpha"] = particle_set[0]->alpha_hat;
  results["time"] = time2- time1;
  //results["particle_trajectory"] = particle_set_trajectory_out;
  //results["log_like_trajectory"] = log_like_trajectory;
  //results["log_prior_trajectory"] = log_prior_trajectory;
  //results["log_post_trajectory"] = log_post_trajectory;
  //results["objective_trajectory"] = objective_trajectory;
  //results["alpha_trajectory"] = alpha_trajectory_out;
  
  return(results);
}
