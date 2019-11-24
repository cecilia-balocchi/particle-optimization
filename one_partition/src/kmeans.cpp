//
//  kmeans.cpp
//  


#include <stdio.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "partition.h"
#include "various_functions.h"
#include "partition_functions.h"
#include <vector>
#include <ctime>

using namespace arma;
using namespace std;


using Rcpp::List;
using Rcpp::Rcout;
using Rcpp::NumericVector;

// [[Rcpp::export]]
Rcpp::List kmeans_particle(arma::mat Y,
                  const arma::mat A_block,
                  const int max_splits = 5,
                  const double a1 = 1.0,
                  const double a2 = 1.0,
                  const double nu_sigma = 3,
                  const double lambda_sigma = 1,
                  const double rho = 0.99,
                  const double eta = 1.0)
{
  
  int n = Y.n_rows;
  int T = Y.n_cols;
  
  arma::vec ybar(n);
  double total_ss = 0;
  for(int i = 0; i < n; i++){
    ybar(i) = arma::mean(Y.row(i));
    total_ss += (T-1) * arma::var(Y.row(i));
  }
  
  std::vector<std::vector<int> > tmp_clusters(1);
  for(int i = 0; i < n; i++) tmp_clusters[0].push_back(i);
  
  LPPartition gamma_0 = new Partition(n, tmp_clusters, ybar, T, A_block, rho, a1, a2, eta);

  LPPartition gamma_new = new Partition(gamma_0);
  
  std::vector<LPPartition> particle_set(max_splits); // holds the results of kmeans after finding connected components
  std::vector<LPPartition> init_particle_set(max_splits); // hold the actual results of kmeans
  for(int l = 0; l < max_splits; l++){
    init_particle_set[l] = new Partition(gamma_0);
    particle_set[l] = new Partition(gamma_0);
  }
  
  
  // we will run k-means for K = 1, ..., 10
  int split_k = 0; // we will only ever start with gamma_init being the partition with a single cluster
  int n_k = 1;
  
  arma::mat U = arma::zeros<arma::mat>(n_k,1); // holds the data passed to kmeans_repeat
  arma::mat means = arma::zeros<arma::mat>(1,max_splits); // holds the means returned by kmeans_repeat
  double min_score = 0.0; // holds sum-of-squared distances from observations to cluster means
  bool status = true; // tells us whether k-means was successful
  
  arma::vec cluster_dist = arma::zeros<arma::vec>(max_splits); // holds distance from single observaiton to each cluster mean
  arma::uvec cluster_dist_indices(max_splits); // used to sort cluster_dist

  std::vector<std::vector<int> > init_new_clusters; // stores output of k-means. NOTE: these clusters may not be connected.
  arma::mat A_tmp = arma::zeros<arma::mat>(n_k,n_k); // store submatrices of A_block. used to determine connected components
  std::vector<std::vector<int> > connected_components; // connected components of the individual sub clusters found by k-means

  std::vector<std::vector<int> > tmp_new_clusters; // used for figuring out the connected components of the new clusters
  std::vector<std::vector<int> > new_clusters; // holds the final new clusters
  std::vector<int> k_star; // actual nearest neighbors of newly formed clusters
  
  int l = 0;
  
  time_t tp;
  int time1 = time(&tp);
  
  for(int num_splits = 1; num_splits <= max_splits; num_splits++){
    if(num_splits > 1){
      n_k = gamma_0->cluster_config[split_k];
      U.set_size(n_k,1);
      means.set_size(1, num_splits);
      min_score = 0.0;
      status = true;
      cluster_dist.set_size(num_splits);
      cluster_dist_indices.set_size(num_splits);
      
      for(int ii = 0; ii < n_k; ii++) U(ii,0) = ybar(gamma_0->clusters[split_k][ii]); // remember this indexing in U
      
      kmeans_repeat(U, means, status, min_score, num_splits, 5000);
      if(status == true){
        init_new_clusters.clear();
        init_new_clusters.resize(num_splits);
        for(int ii = 0; ii < n_k; ii++){
          cluster_dist.zeros();
          for(int new_k = 0; new_k < num_splits; new_k++) cluster_dist(new_k) = arma::norm(U.row(ii).t() - means.col(new_k));
          cluster_dist_indices = arma::sort_index(cluster_dist, "ascend");
          init_new_clusters[cluster_dist_indices(0)].push_back(gamma_0->clusters[split_k][ii]);
        }
        
        // Rcpp::Rcout << "k = " << num_splits << " tot.withinss = " << min_score << std::endl;

        
        tmp_new_clusters.clear();
        for(int kk = 0; kk < init_new_clusters.size(); kk++){
          connected_components.clear();
          A_tmp = Submatrix(A_block, init_new_clusters[kk].size(), init_new_clusters[kk].size(), init_new_clusters[kk], init_new_clusters[kk]);
          new_Connected_Components(A_tmp, init_new_clusters[kk].size(), init_new_clusters[kk], connected_components);
          for(int new_k = 0; new_k < connected_components.size(); new_k++){
            tmp_new_clusters.push_back(connected_components[new_k]);
          }
        }

        
        new_clusters.clear();
        new_clusters.resize(tmp_new_clusters.size());
        for(int nc_ix = 0; nc_ix < new_clusters.size(); nc_ix++){
          for(int ii = 0; ii < tmp_new_clusters[nc_ix].size(); ii++){
            new_clusters[nc_ix].push_back(tmp_new_clusters[nc_ix][ii]);
          }
        }
        delete gamma_new;
        gamma_new = new Partition(n, init_new_clusters, ybar, T, A_block, rho, a1, a2, eta);
        init_particle_set[l]->Copy_Partition(gamma_new);
        
        delete gamma_new;
        gamma_new = new Partition(n, new_clusters, ybar, T, A_block, rho, a1, a2, eta);
        particle_set[l]->Copy_Partition(gamma_new);
        l++;
        //Rcpp::Rcout << std::endl;
      } // closes if checking that kmean_repeat exited successfully
    } else{
      //Rcpp::Rcout << "k = " << num_splits << " tot.withinss = " << (n - 1) * arma::var(ybar) << std::endl;
      
      delete gamma_new;
      gamma_new = new Partition(gamma_0);
      init_particle_set[l]->Copy_Partition(gamma_new);
      
      delete gamma_new;
      gamma_new = new Partition(gamma_0);
      particle_set[l]->Copy_Partition(gamma_new);
      l++;
    } // closes if/else checking whether num_splits == 1
  } // closes loop over num_splits
  int time2 = time(&tp);
  
  Rcpp::List init_unik_particles_out;
  format_particle_set(init_particle_set, init_unik_particles_out);
  arma::mat init_alpha_hat_particle = arma::zeros<mat>(n, init_particle_set.size());
  for(int l = 0; l < init_particle_set.size();l++){
    for(int i = 0; i < n; i++){
      init_alpha_hat_particle(i,l) = init_particle_set[l]->alpha_hat[i];
    }
  }
  
  Rcpp::List unik_particles_out;
  std::vector<double> log_like(particle_set.size());
  std::vector<double> log_prior(particle_set.size());
  std::vector<double> log_post(particle_set.size());
  format_particle_set(particle_set, unik_particles_out);
  arma::mat alpha_hat_particle = arma::zeros<mat>(n, particle_set.size());
  for(int l = 0; l < particle_set.size();l++){
    for(int i = 0; i < n; i++){
      alpha_hat_particle(i,l) = particle_set[l]->alpha_hat[i];
    }
    log_like[l] = total_log_like(particle_set[l], total_ss, T, nu_sigma, lambda_sigma);
    log_prior[l] = total_log_prior(particle_set[l]);
    log_post[l] = total_log_post(particle_set[l], total_ss, T, nu_sigma, lambda_sigma);
  }
  

  Rcpp::List results;
  results["Y"] = Y;
  results["particles"] = unik_particles_out;
  results["alpha"] = alpha_hat_particle;
  results["log_like"] = log_like;
  results["log_prior"] = log_prior;
  results["log_post"] = log_post;
  results["init_particles"] = init_unik_particles_out;
  results["init_alpha_hat_particle"] = init_alpha_hat_particle;
  results["time"] = time2 - time1;
  return(results);
  
  
  
}
