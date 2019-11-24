
// partition_functions.cpp

// [[Rcpp::depends(RcppArmadillo)]]

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <vector>
#include <algorithm>
#include "partition.h"
#include "various_functions.h"
#include "partition_functions.h"



// Compare two partitions and see if they are equal
int Partition_Equal(Partition *partition1, Partition *partition2){
  int flag = 1;
    // simpler to just compare pairwise allocation
  for(int i = 0 ; i < partition1->nObs;i++){
    for(int j = 0; j < partition1->nObs;j++){
      if(partition1->pairwise_assignment(i,j) != partition2->pairwise_assignment(i,j)){
        //Rcpp::Rcout << "[Partition_Equal]: i = " << i << " j = " << j << endl;
        flag = 0;
        break;
      }
    }
    if(flag == 0) break;
  }
  return flag;
}

void get_unik_particles(std::vector<std::vector<int> > &particle_map, std::vector<double> &p_star, std::vector<int> &counts, std::vector<LPPartition> particle_set, std::vector<double> w)
{
  particle_map.clear();
  p_star.clear();
  counts.clear();
  std::vector<std::vector<int> > tmp_particle_map(1, std::vector<int>(1,0)); // first particle always considered unique
  std::vector<double> tmp_w(1, w[0]);
  int counter = 0;
  for(int l = 1; l < particle_set.size(); l++){
    counter = 0;
    for(int ul = 0; ul < tmp_particle_map.size(); ul++){
      if(Partition_Equal(particle_set[l], particle_set[tmp_particle_map[ul][0]]) == 1){
        tmp_particle_map[ul].push_back(l);
        tmp_w[ul] += w[l];
      } else{
        counter++;
      }
    }
    if(counter == tmp_particle_map.size()){
      tmp_particle_map.push_back(std::vector<int>(1,l));
      tmp_w.push_back(w[l]);
    }
  }
  
  arma::vec w_vec = arma::zeros<vec>(tmp_particle_map.size());
  arma::uvec w_indices(tmp_particle_map.size());
  for(int l = 0; l < tmp_particle_map.size(); l++){
    w_vec(l) = tmp_w[l];
  }
  w_indices = arma::sort_index(w_vec, "descend");
  for(int l = 0; l < tmp_particle_map.size(); l++){
    particle_map.push_back(tmp_particle_map[w_indices(l)]);
    p_star.push_back(tmp_w[w_indices(l)]);
    counts.push_back(tmp_particle_map[w_indices(l)].size());
  }
  
}

// function to compute entropy
// when we replace the (current_l)^th particle with the candidate_clusters
//
double Entropy(unsigned current_l, Partition* candidate_particle, std::vector<LPPartition> particle_set, std::vector<double> w){
  unsigned L = particle_set.size();
  // need to loop over to extract the unique partitions
  std::vector<LPPartition> unik_particles;
  std::vector<double> p_star;

  unik_particles.push_back(candidate_particle);
  p_star.push_back(w[current_l]);

  // in a sense, we are replacing particle_set[current_l] with candidate_particle
  // by adding it to the unik_particles vector

  int num_unik_particles = 1;
  int counter = 0;
  for(unsigned l = 0; l < L; l++){ // loop over current particle set
	counter = 0;
	if(l != current_l){
	  for(unsigned ul = 0; ul < num_unik_particles; ul++){
	    if(Partition_Equal(particle_set[l], unik_particles[ul]) == 1){
	      p_star[ul] += w[l]; // update p_star
	      break;
        } else {
	      counter++;
        }
	  }

	  if(counter == num_unik_particles){
	    // we have found a new unique particle
	    unik_particles.push_back(particle_set[l]);
	    p_star.push_back(w[l]);
	    num_unik_particles++;
      }
    }
  }
  double entropy = 0.0;
  double tmp = 0.0;
  for(unsigned ul = 0; ul < num_unik_particles; ul++){
    tmp = p_star[ul] * log(p_star[ul]);
    if(!std::isnan(tmp)){
      entropy += tmp;
    }
  }
  //std::cout << std::endl;
  return -1.0 * entropy;

}

// a silly function to add up log-likelihoods and log-priors

double total_log_post(LPPartition partition, const double total_ss, const int T, const double nu_sigma, const double lambda_sigma){
	//double log_post = total_log_like(partition, a_sigma, nu_sigma);
  double log_post = total_log_like(partition, total_ss, T, nu_sigma, lambda_sigma);
  for(int k = 0; k < partition->K; k++){
    log_post += partition->log_prior[k];
  }
  
	//for(int k = 0; k < partition->K; k++){
	//	log_post += partition->log_like[k] + partition->log_prior[k];
	//}
	return log_post;
}

//double total_log_like(LPPartition partition, const double a_sigma, const double nu_sigma){
double total_log_like(LPPartition partition, const double total_ss, const int T, const double nu_sigma, const double lambda_sigma){

  double log_like = 0.0;
  double log_det = 0.0;
  double quad_form = 0.0;
  
  for(int k = 0; k < partition->K; k++){
    log_det += partition->log_det_Omegay[k];
    quad_form += partition->y_Omegay_y[k];
  }
  
  log_like = 0.5 * log_det - ( (nu_sigma + ((double) T) * ((double) partition->nObs))/2 ) * log( (nu_sigma + lambda_sigma + quad_form + total_ss) / 2);
  //log_like = 0.5 * log_det - ( (nu_sigma + partition->nObs)/2) * log( (nu_sigma * lambda_sigma + quad_form)/2);

  return log_like;
}
double total_log_prior(LPPartition partition){
  double log_prior = 0.0;
  for(int k = 0; k < partition->K; k++){
	log_prior += partition->log_prior[k];
  }
  return log_prior;
}

double alpha_bar_func(std::vector<int> new_cluster, LPPartition gamma_l, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2)
{
  int n_k = new_cluster.size();
  if(n_k == 1){
    return(gamma_l->alpha_hat[new_cluster[0]]);
  } else{
    double tmp_alpha_bar = 0.0;
    for(int i = 0; i < n_k; i++){
      tmp_alpha_bar += (1 - rho) * gamma_l->alpha_hat[new_cluster[i]];
    }
    return( ((1.0/a1) * tmp_alpha_bar)/(1/a2 + n_k * (1-rho)/a1));
  }
}

void update_w(std::vector<LPPartition> particle_set, std::vector<double> &w, const int L, const double total_ss, const int T, const double nu_sigma, const double lambda_sigma, const double lambda)
{
  
  double max_log_post = 0.0;
  double tmp_log_post = 0.0;
  double tmp_norm = 0.0;
  double tmp_p_star = 0.0;
  //Rcpp::Rcout << "[update_w]: Entering" << endl;

  // First we need to identify the unique particles
  std::vector<LPPartition> unik_particles;
  unik_particles.push_back(particle_set[0]);
  std::vector<int> particle_assignment(L,-1); // tells us to which unique particle each element of particle_set corresponds
  particle_assignment[0] = 0;
  std::vector<int> particle_counts;
  
  std::vector<double> p_star;
  std::vector<double> log_post;

  log_post.push_back(total_log_post(particle_set[0], total_ss, T, nu_sigma, lambda_sigma));
  p_star.push_back(0);
  //max_log_post = total_log_post(particle_set[0], a_sigma, nu_sigma);
  max_log_post = total_log_post(particle_set[0], total_ss, T, nu_sigma, lambda_sigma);
  //Rcpp::Rcout << "[update_w]: Ready to find unik_particles" << endl;
  
  int num_unik_particles = 1;
  int counter = 0;

  for(int l = 1; l < L; l++){ // loop over the particle set
    counter = 0;
    for(int ul = 0; ul < num_unik_particles; ul++){
      if(Partition_Equal(particle_set[l], unik_particles[ul]) == 1){
        // l^th particle is equal to the ul^th unique partition
        particle_assignment[l] = ul;
        break;
      } else{
        counter++;
      }
    } // closes loop over the unique particles
    if(counter == num_unik_particles){
      // we have found a new unique particle
      unik_particles.push_back(particle_set[l]);
      particle_assignment[l] = num_unik_particles;
      p_star.push_back(0.0); // for now we populate p_star with 0's
      //tmp_log_post = total_log_post(particle_set[l], a_sigma, nu_sigma);
      tmp_log_post = total_log_post(particle_set[l], total_ss, T, nu_sigma, lambda_sigma);
      log_post.push_back(tmp_log_post);
      if(tmp_log_post > max_log_post){
        max_log_post = tmp_log_post;
      }
      num_unik_particles++;
    }
  }
  //Rcpp::Rcout << "[update_w] : num_unik_particles = " << num_unik_particles << endl;

  particle_counts.clear();
  particle_counts.resize(num_unik_particles,0);
  for(int l = 0; l < L; l++){
    particle_counts[particle_assignment[l]]++;
  }
  for(int ul = 0; ul < num_unik_particles;ul++){
    tmp_log_post = log_post[ul] - max_log_post;
    tmp_p_star = exp(1/lambda * tmp_log_post);
    tmp_norm += tmp_p_star;
    p_star[ul] = tmp_p_star;
  }
  for(int ul = 0; ul < num_unik_particles;ul++){
    p_star[ul] /= tmp_norm;
  }
  for(int l = 0; l < L; l++){
    w[l] = p_star[particle_assignment[l]]/( (double) particle_counts[particle_assignment[l]]);
  }
}


void get_island(split_info &si, LPPartition gamma_l, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double island_frac)
{
  si.num_splits = 0;
  si.split_k.clear();
  si.nearest_neighbor.clear();
  si.new_clusters.clear();
  
  int orig_K = gamma_l->K;
  int n_k = 1;
  
  double alpha_bar = 0.0; // holds the cluster mean
  arma::vec distance = arma::zeros<arma::vec>(1); // holds distance of each alpha-hat from the overall cluster mean
  arma::uvec indices(1); // used to sort the elements based on how far they are from the cluster mean
  int num_islands = 0; // number of islands we try to create
  
  std::vector<int> island(1);
  std::vector<int> remain(1); // holds indices that remain in cluster split_k (i.e. they are not adjacent to block-group k)
  arma::mat A_tmp = arma::zeros<arma::mat>(1,1); // used to find connected components of remain
  std::vector<std::vector<int> > connected_components; // temporariy hold connected components of each sub-cluster
  std::vector<std::vector<int> > tmp_new_clusters; // will contains connected components of remain and the island
  std::vector<std::vector<int> > new_clusters; // holds the new connected sub-clusters created
  std::vector<int> k_star; // holds the labels of nearest neighbors of new sub-clusters
  
  
  for(int k = 0; k < orig_K; k++){
    n_k = gamma_l->cluster_config[k];
    if(n_k > 1){
      distance.reset();
      distance.set_size(n_k);
      indices.reset();
      indices.set_size(n_k);
      
      alpha_bar = gamma_l->alpha_bar[k];
      for(int i = 0; i < n_k; i++) distance(i) = abs(alpha_bar - gamma_l->alpha_hat[gamma_l->clusters[k][i]]);
      indices = arma::sort_index(distance, "descend");
      num_islands = ceil( (double) n_k * island_frac);
      for(int ix = 0; ix < num_islands; ix++){
        remain.clear();
        island.clear();
        for(int i = 0; i < n_k; i++){
          if(i != ix) remain.push_back(gamma_l->clusters[k][indices(i)]); // includes all but the ix-th furtherst point from cluster mean
          else island.push_back(gamma_l->clusters[k][indices(i)]);
        } // closes loop populating remain and island
        
        // find the connected components of remain
        if(remain.size() > 0 & island.size() > 0){ // this is probably superfluous but it's fine for now
          connected_components.clear();
          tmp_new_clusters.clear();
          A_tmp = Submatrix(A_block, remain.size(), remain.size(), remain, remain);
          new_Connected_Components(A_tmp, remain.size(), remain, connected_components);
          for(int new_k = 0; new_k < connected_components.size(); new_k++){
            tmp_new_clusters.push_back(connected_components[new_k]);
          }
          tmp_new_clusters.push_back(island); // add the island to the mix
          new_clusters.clear();
          new_clusters.resize(tmp_new_clusters.size());
          k_star.clear();
          k_star.resize(tmp_new_clusters.size());
          get_subcluster_neighbor(tmp_new_clusters, new_clusters, k_star, k, gamma_l, T, A_block, rho, a1, a2);
          
          si.split_k.push_back(k);
          si.new_clusters.push_back(new_clusters);
          si.nearest_neighbor.push_back(k_star);
          si.num_splits++;
        } // closes if/else checking that at least 1 elements remains in k.
        
      } // closes loop over the possible islands
    } // closes if/else checking that there is more than 1 element in the cluster
  } // closes loop over all clusters
}

// This is just get_island but it tries to move every group out of its cluster
void get_local(split_info &si, LPPartition gamma_l, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2)
{
  si.num_splits = 0;
  si.split_k.clear();
  si.nearest_neighbor.clear();
  si.new_clusters.clear();
  
  int orig_K = gamma_l->K;
  int n_k = 1;
  
  double alpha_bar = 0.0; // holds the cluster mean
  arma::vec distance = arma::zeros<arma::vec>(1); // holds distance of each alpha-hat from the overall cluster mean
  arma::uvec indices(1); // used to sort the elements based on how far they are from the cluster mean
  int num_islands = 0; // number of islands we try to create
  
  std::vector<int> island(1); // holds index of the island being removed from cluster split_k
  std::vector<int> remain(1); // holds indices that remain in cluster split_k
  arma::mat A_tmp = arma::zeros<arma::mat>(1,1); // used to find connected components of remain
  std::vector<std::vector<int> > connected_components; // temporariy hold connected components of each sub-cluster
  std::vector<std::vector<int> > tmp_new_clusters; // will contains connected components of remain and the island
  std::vector<std::vector<int> > new_clusters; // holds the new connected sub-clusters created
  std::vector<int> k_star; // holds the labels of nearest neighbors of new sub-clusters
  
  
  for(int k = 0; k < orig_K; k++){
    n_k = gamma_l->cluster_config[k];
    if(n_k > 1){
      distance.reset();
      distance.set_size(n_k);
      indices.reset();
      indices.set_size(n_k);

      alpha_bar = gamma_l->alpha_bar[k];
      for(int i = 0; i < n_k; i++) distance(i) = abs(alpha_bar - gamma_l->alpha_hat[gamma_l->clusters[k][i]]);
      indices = arma::sort_index(distance, "descend");
      num_islands = n_k; // will try to remove all elements in sequential order
      for(int ix = 0; ix < num_islands; ix++){
        remain.clear();
        island.clear();
        for(int i = 0; i < n_k; i++){
          if(i != ix) remain.push_back(gamma_l->clusters[k][indices(i)]); // includes all but the ix-th furtherst point from cluster mean
          else island.push_back(gamma_l->clusters[k][indices(i)]);
        } // closes loop populating remain and island
        
        // find the connected components of remain
        if(remain.size() > 0 & island.size() > 0){ // this is probably superfluous but it's fine for now
          connected_components.clear();
          tmp_new_clusters.clear();
          A_tmp = Submatrix(A_block, remain.size(), remain.size(), remain, remain);
          new_Connected_Components(A_tmp, remain.size(), remain, connected_components);
          for(int new_k = 0; new_k < connected_components.size(); new_k++){
            tmp_new_clusters.push_back(connected_components[new_k]);
          }
          tmp_new_clusters.push_back(island); // add the island to the mix
          new_clusters.clear();
          new_clusters.resize(tmp_new_clusters.size());
          k_star.clear();
          k_star.resize(tmp_new_clusters.size());
          get_subcluster_neighbor(tmp_new_clusters, new_clusters, k_star, k, gamma_l, T, A_block, rho, a1, a2);
          
          si.split_k.push_back(k);
          si.new_clusters.push_back(new_clusters);
          si.nearest_neighbor.push_back(k_star);
          si.num_splits++;
        } // closes if/else checking that at least 1 elements remains in k.
        
      } // closes loop over the possible islands
    } // closes if/else checking that there is more than 1 element in the cluster
  } // closes loop over all clusters
}


void get_border(split_info &si, LPPartition gamma_l, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2)
{
  
  si.num_splits = 0;
  si.split_k.clear();
  si.nearest_neighbor.clear();
  si.new_clusters.clear();

  int orig_K = gamma_l->K;
  int split_k = 0;
  
  arma::mat A_sub = arma::zeros<arma::mat>(1,1); // sub-matrix of A_sub used to see if clusters k and kk are adjacenct
  std::vector<int> adj_cluster(1); // temporarily holds the labels of adjacent clusters
  arma::vec adj_cluster_dist = arma::zeros<arma::vec>(1); // holds distance between cluster k and its adjacent clusters
  arma::uvec adj_cluster_indices(1); // for sorting the adjacent clusters
  
  
  std::vector<int> border(1); // holds indices of the border elements that move from split_k to k
  std::vector<int> remain(1); // holds indices that remain in cluster split_k (i.e. they are not adjacent to block-group k)
  arma::mat A_tmp = arma::zeros<arma::mat>(1,1); // used to find connected components of remain
  std::vector<std::vector<int> > connected_components; // temporariy hold connected components of each sub-cluster
  std::vector<std::vector<int> > new_clusters; // holds the new connected sub-clusters created
  std::vector<int> k_star; // holds the labels of nearest neighbors of new sub-clusters
  
  if(orig_K > 1){
    for(int k = 0; k < orig_K; k++){
      adj_cluster.clear();
      for(int kk = 0; kk < orig_K; kk++){
        if(kk != k){
          A_sub = Submatrix(A_block, gamma_l->cluster_config[k], gamma_l->cluster_config[kk], gamma_l->clusters[k], gamma_l->clusters[kk]);
          if(any(vectorise(A_sub) == 1.0)) adj_cluster.push_back(kk); // cluster k and kk are adjacent
        } // closes if checking that kk != k
      } // closes loop over all clusters
      if(adj_cluster.size() > 0){
        adj_cluster_dist.resize(adj_cluster.size());
        adj_cluster_indices.resize(adj_cluster.size());
        for(int kk = 0; kk < adj_cluster.size(); kk++){
          adj_cluster_dist(kk) = abs(gamma_l->alpha_bar[k] - gamma_l->alpha_bar[adj_cluster[kk]]); // distance between cluster k and kk
        }
        adj_cluster_indices = arma::sort_index(adj_cluster_dist, "ascend");
        split_k = adj_cluster[adj_cluster_indices(0)];
        A_sub = Submatrix(A_block, gamma_l->cluster_config[k], gamma_l->cluster_config[split_k], gamma_l->clusters[k], gamma_l->clusters[split_k]); // submatrix of A_block corresponding to clusters k and its neighbor split_k
        
        border.clear();
        remain.clear();
        
        for(int i = 0; i < gamma_l->cluster_config[split_k]; i++){
          // anything in A_sub.col(i) is equal to 1 then element i in cluster split_k is a border
          if(any(vectorise(A_sub.col(i)) == 1.0)) border.push_back(gamma_l->clusters[split_k][i]);
          else remain.push_back(gamma_l->clusters[split_k][i]);
        }
        if(remain.size() > 0 & border.size() > 0){ // condition may be a bit redundant
          // first find the connected components of the elements in remain
          new_clusters.clear();
          k_star.clear();
          connected_components.clear();
          A_tmp = Submatrix(A_block, remain.size(), remain.size(), remain, remain);
          new_Connected_Components(A_tmp, remain.size(), remain, connected_components);
          for(int new_k = 0; new_k < connected_components.size(); new_k++){
            new_clusters.push_back(connected_components[new_k]);
            k_star.push_back(-1);
          }
          // up to now new_clusters and k_star don't include the border elements which are moving to cluster k
          new_clusters.push_back(border);
          k_star.push_back(k);
          
          si.split_k.push_back(split_k);
          si.new_clusters.push_back(new_clusters);
          si.nearest_neighbor.push_back(k_star);
          si.num_splits++;
          
          
          
        } // closes if checking that not all elements in split_k are border (if they are don't do anything!)
      } // closes if checking that there are valid border moves
    } // closes loop over the clusters k
  } // check that there are at least two clusters
}


void get_merge(merge_info &mi, LPPartition gamma_l, const arma::mat &A_block)
{
  mi.num_merges = 0;
  mi.rec_k.clear();
  mi.donor_k.clear();
  
  double alpha_bar = 0.0; // holds the cluster mean
  std::vector<int> tmp_neighbor(1); // holds potential neighbors that are being merged into k
  std::vector<double> tmp_distance(1); // holds distance from k to potential neighbors
  arma::vec distance(1);
  arma::uvec indices(1);
  arma::mat A_sub = arma::zeros<mat>(1,1);
  std::vector<std::vector<int> > new_clusters;
  
  if(gamma_l->K > 1){
    for(int k = 0; k < gamma_l->K; k++){
      alpha_bar = gamma_l->alpha_bar[k];
      tmp_neighbor.clear();
      tmp_distance.clear();
      for(int kk = k+1;kk < gamma_l->K; kk++){
        A_sub = Submatrix(A_block, gamma_l->cluster_config[k], gamma_l->cluster_config[kk], gamma_l->clusters[k], gamma_l->clusters[kk]);
        if(any(vectorise(A_sub) == 1)){ // clusters k and kk are neighbors
          tmp_neighbor.push_back(kk);
          tmp_distance.push_back(abs(alpha_bar - gamma_l->alpha_bar[kk]));
        }
      }
      if(tmp_neighbor.size() > 0){
        distance.reset();
        indices.reset();
        distance.set_size(tmp_neighbor.size());
        indices.set_size(tmp_neighbor.size());
        for(int i = 0; i < tmp_neighbor.size(); i++){
          distance(i) = tmp_distance[i];
        }
        indices = arma::sort_index(distance, "ascend");
        
        mi.rec_k.push_back(k);
        mi.donor_k.push_back(tmp_neighbor[indices(0)]);
        mi.num_merges++;
      }
    }
  }
}

void get_spectral_split(split_info &si, LPPartition gamma_l, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const int reps){

  si.num_splits = 0;
  si.split_k.clear();
  for(int i = 0; i < si.new_clusters.size(); i++){
    si.new_clusters[i].clear();
    si.nearest_neighbor[i].clear();
  }
  si.new_clusters.clear();
  si.nearest_neighbor.clear();
  
  int orig_K = gamma_l->K;
  int max_splits = 5;
  int n_k = 1;
  int split_k = 0;
  
  // initialize all of the stuff for spectral clustering here
  arma::mat A_block_k = arma::zeros<mat>(n_k, n_k);
  arma::vec alpha_hat_cluster = arma::zeros<vec>(n_k);
  arma::mat alpha_sim = arma::zeros<mat>(n_k,n_k);
  arma::mat I_k = arma::zeros<mat>(n_k,n_k);
  I_k.eye();
  arma::mat alpha_dist = arma::zeros<mat>(n_k,n_k);
  arma::mat W_alpha_cl = arma::zeros<mat>(n_k,n_k);
  arma::mat Dinv_sqrt = arma::zeros<mat>(n_k,n_k);
  arma::mat L = arma::zeros<mat>(n_k,n_k);
  arma::mat eigvec = arma::zeros<mat>(n_k,n_k); // holds the eigenvectors of L
  arma::vec eigval = arma::zeros<vec>(n_k); // holds the eigenvalues of L
  arma::mat U = arma::zeros<mat>(n_k,2); // will hold the first several eigenvalues of L
  arma::mat means = arma::zeros<mat>(n_k,2); // the cluster means passed to kmeans
  double min_score = 0.0; // used with kmeans_repeat
  bool status = true; // tells us that kmeans was successful
  
  arma::vec cluster_dist = arma::zeros<vec>(1); // holds distance to cluster means
  arma::uvec cluster_dist_indices(1); // used to sort cluster_dist
  
  std::vector<std::vector<int> > init_new_clusters; // stores output of kmeans. NOTE: these clusters may not be connected
  arma::mat A_tmp = arma::zeros<arma::mat>(n_k,n_k);// stores submatrix of A_block, used to determine connected components
  std::vector<std::vector<int> > connected_components; // connected components of individual sub-clusters discovered
  std::vector<std::vector<int> > tmp_new_clusters;// used for figuring out the connected components of new clusters
  std::vector<std::vector<int> > new_clusters; // holds the final new clusters
  std::vector<int> k_star; // actual nearest neighbors of newly formed clusters
  
  for(int k = 0; k < orig_K; k++){
    n_k = gamma_l->cluster_config[k];
    split_k = k;
    max_splits = 5;
    if(n_k > 1){
      if(sqrt(n_k) < 5) max_splits = ceil(sqrt(n_k));
      else max_splits = 5;
      for(int num_splits = 2; num_splits <= max_splits; num_splits++){
        
        // re-size everything
        A_block_k.set_size(n_k,n_k);
        alpha_hat_cluster.set_size(n_k);
        alpha_sim.set_size(n_k,n_k);
        I_k.set_size(n_k,n_k);
        I_k.eye();
        alpha_dist.set_size(n_k,n_k);
        W_alpha_cl.set_size(n_k,n_k);
        Dinv_sqrt.set_size(n_k,n_k);
        L.set_size(n_k,n_k);
        eigvec.set_size(n_k,n_k);
        eigval.set_size(n_k);
        U.set_size(n_k,num_splits);
        means.set_size(num_splits,num_splits);
        status = true;
        min_score = 0.0;
        cluster_dist.set_size(num_splits);
        cluster_dist_indices.set_size(num_splits);
        
        // Spectral Clustering Begins //
        A_block_k = Submatrix(A_block, n_k, n_k, gamma_l->clusters[split_k], gamma_l->clusters[split_k]);
        for(int i = 0; i  < n_k; i++){
          alpha_hat_cluster(i) = gamma_l->alpha_hat[gamma_l->clusters[split_k][i]];
        }
        alpha_dist = Distance_matrix(alpha_hat_cluster,n_k); // distance matrix
        alpha_sim = exp(-1.0 * square(alpha_dist)/(2 * arma::var(alpha_hat_cluster))); // similarity matrix
        W_alpha_cl = I_k + alpha_sim % A_block_k; // ensure non-adjacent indices have 0 similarity
        Dinv_sqrt = arma::diagmat(1/sqrt(arma::sum(W_alpha_cl,1)));
        L = I_k - Dinv_sqrt * W_alpha_cl * Dinv_sqrt;
        arma::eig_sym(eigval, eigvec, L);
        U = eigvec.cols(0,num_splits-1);
        U = arma::diagmat(1/sqrt(arma::sum(arma::square(U),1))) * U;
        // now run kmeans_repeat
        kmeans_repeat(U, means, status, min_score, num_splits, reps);
        
        if(status == false){
          Rcpp::Rcout << "kmeans failed!" << "k = " << k << "num_splits = " << num_splits << endl;
        } else{
          //Rcpp::Rcout << "kmeans succeeded!" << endl;
          // do more stuff
          init_new_clusters.clear();
          init_new_clusters.resize(num_splits);
          
          // loop over rows of U to figure out to which cluster it belongs
          for(int i = 0; i < n_k; i++){
            cluster_dist.zeros();
            for(int new_k = 0; new_k < num_splits; new_k++){
              cluster_dist(new_k) = arma::norm(U.row(i).t() - means.col(new_k));
            } // closes loop computing distance of U.row(i) to centroid new_k
            cluster_dist_indices = arma::sort_index(cluster_dist, "ascend");
            init_new_clusters[cluster_dist_indices(0)].push_back(gamma_l->clusters[split_k][i]);
          } // closes loop over rows of U
          
          // now loop over init_new_clusters and find the connected components
          tmp_new_clusters.clear();
          for(int kk = 0; kk < init_new_clusters.size(); kk++){
            connected_components.clear();
            A_tmp = Submatrix(A_block, init_new_clusters[kk].size(), init_new_clusters[kk].size(), init_new_clusters[kk], init_new_clusters[kk]); // submatrix of A_block corresponding to newly discovered sub-cluster
            new_Connected_Components(A_tmp, init_new_clusters[kk].size(), init_new_clusters[kk], connected_components);
            for(int new_k = 0; new_k < connected_components.size(); new_k++){
              tmp_new_clusters.push_back(connected_components[new_k]);
            }
          } // closes loop over the number of clusters found by k-means
          // at this point tmp_new_clusters contains all of the new sub-clusters
          
          // get_subcluster_neighbor finds the neighbors of each sub-cluster and arranges them by distance to nearest neighbor
          new_clusters.clear();
          new_clusters.resize(tmp_new_clusters.size());
          k_star.clear();
          k_star.resize(tmp_new_clusters.size());
          get_subcluster_neighbor(tmp_new_clusters, new_clusters, k_star, k, gamma_l, T, A_block, rho, a1, a2);

          
          // write the information about the split to si
          si.split_k.push_back(k);
          si.new_clusters.push_back(new_clusters);
          si.nearest_neighbor.push_back(k_star);
          si.num_splits++;
        } // closes if/else checking that kmeans_repeat succeeded.
      } // closes loop over number of possible splits
    } // closes if checking that there is at least one element in the cluster
  } // closes loop over the clusters
  
}

void get_tail_split(split_info &si, LPPartition gamma_l, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double tail_frac)
{
  si.num_splits = 0;
  si.split_k.clear();
  for(int i = 0; i < si.new_clusters.size(); i++){
    si.new_clusters[i].clear();
    si.nearest_neighbor[i].clear();
  }
  si.new_clusters.clear();
  si.nearest_neighbor.clear();
  
  int orig_K = gamma_l->K;
  int n_k = 1;
  
  arma::vec alpha_hat_cluster = arma::zeros<arma::vec>(1);
  arma::uvec alpha_hat_indices(1);
  
  std::vector<int> left_tail;
  std::vector<int> right_tail;
  std::vector<int> center;
  
  std::vector<int> remain; // contains all of the indicies that are not removed from the cluster
  
  std::vector<std::vector<int> > init_new_clusters; // holds the initial sub-clusters
  arma::mat A_tmp; // submatrix corresponding to an initial subcluster
  std::vector<std::vector<int> > connected_components; // connected components of each initial individual sub-cluster
  std::vector<std::vector<int> > tmp_new_clusters; // temporarily holds the new connected sub-clusters
  std::vector<std::vector<int> > new_clusters; // the final sub-clusters arranged in order
  std::vector<int> k_star; // nearest neighbors of newly formed subclusters
  
  for(int k = 0; k < orig_K; k++){
    n_k = gamma_l->cluster_config[k];
    if(n_k > 1){
      alpha_hat_cluster.clear();
      alpha_hat_indices.clear();
      alpha_hat_cluster.set_size(n_k); // holds the alpha-hat values in the cluster
      alpha_hat_indices.set_size(n_k);
      for(int i = 0; i < n_k; i++){
        alpha_hat_cluster(i) = gamma_l->alpha_hat[gamma_l->clusters[k][i]];
      }
      alpha_hat_indices = arma::sort_index(alpha_hat_cluster, "ascend");
      left_tail.clear();
      right_tail.clear();
      center.clear();
      for(int i = 0; i < ceil(n_k*tail_frac);i++){
        left_tail.push_back(gamma_l->clusters[k][alpha_hat_indices(i)]);
        right_tail.push_back(gamma_l->clusters[k][alpha_hat_indices(n_k - 1 - i)]);
      }
      for(int i = ceil(n_k*tail_frac); i < n_k - ceil(n_k*tail_frac);i++){
        center.push_back(gamma_l->clusters[k][alpha_hat_indices(i)]);
      }
      
      // now try to remove only the left tail
      init_new_clusters.clear();
      init_new_clusters.push_back(left_tail);
      // must figure out which elements are staying in the original cluster
      remain.clear();
      for(int ii = 0; ii < center.size(); ii++) remain.push_back(center[ii]); // everything from the center remains
      for(int ii = 0; ii < right_tail.size(); ii++) remain.push_back(right_tail[ii]); // everything from the right tail remains
      init_new_clusters.push_back(remain);
      
      tmp_new_clusters.clear();
      for(int kk = 0; kk < init_new_clusters.size(); kk++){
        connected_components.clear();
        A_tmp = Submatrix(A_block, init_new_clusters[kk].size(), init_new_clusters[kk].size(), init_new_clusters[kk], init_new_clusters[kk]);
        new_Connected_Components(A_tmp, init_new_clusters[kk].size(), init_new_clusters[kk], connected_components);
        for(int new_k = 0; new_k < connected_components.size(); new_k++) tmp_new_clusters.push_back(connected_components[new_k]);
      } // closes loop over the initial set of sub-clusters
      // at this point tmp_new_clusters contains all of the new sub-clusters and they are each connected
      new_clusters.clear();
      new_clusters.resize(tmp_new_clusters.size());
      k_star.clear();
      k_star.resize(tmp_new_clusters.size());
      get_subcluster_neighbor(tmp_new_clusters, new_clusters, k_star, k, gamma_l, T, A_block, rho, a1, a2);
      
      si.split_k.push_back(k);
      si.new_clusters.push_back(new_clusters);
      si.nearest_neighbor.push_back(k_star);
      si.num_splits++;
      
      // now try to remove the right tail
      init_new_clusters.clear();
      init_new_clusters.push_back(right_tail);
      remain.clear();
      for(int ii = 0; ii < center.size(); ii++) remain.push_back(center[ii]);
      for(int ii = 0; ii < left_tail.size(); ii++) remain.push_back(left_tail[ii]);
      init_new_clusters.push_back(remain);
      tmp_new_clusters.clear();
      for(int kk = 0; kk < init_new_clusters.size(); kk++){
        connected_components.clear();
        A_tmp = Submatrix(A_block, init_new_clusters[kk].size(), init_new_clusters[kk].size(), init_new_clusters[kk], init_new_clusters[kk]);
        new_Connected_Components(A_tmp, init_new_clusters[kk].size(), init_new_clusters[kk], connected_components);
        for(int new_k = 0; new_k < connected_components.size(); new_k++) tmp_new_clusters.push_back(connected_components[new_k]);
      } // closes loop over the initial set of sub-clusters
      // at this point tmp_new_clusters contains all of the new sub-clusters and they are each connected
      new_clusters.clear();
      new_clusters.resize(tmp_new_clusters.size());
      k_star.clear();
      k_star.resize(tmp_new_clusters.size());
      get_subcluster_neighbor(tmp_new_clusters, new_clusters, k_star, k, gamma_l, T, A_block, rho, a1, a2);
      
      si.split_k.push_back(k);
      si.new_clusters.push_back(new_clusters);
      si.nearest_neighbor.push_back(k_star);
      si.num_splits++;
      
      // now try to remove both the left and right tail
      init_new_clusters.clear();
      init_new_clusters.push_back(left_tail);
      for(int ii = 0; ii < right_tail.size(); ii++) init_new_clusters[0].push_back(right_tail[ii]); // everything from the both the left and right tails are together here
      init_new_clusters.push_back(center);
      
      tmp_new_clusters.clear();
      for(int kk = 0; kk < init_new_clusters.size(); kk++){
        connected_components.clear();
        A_tmp = Submatrix(A_block, init_new_clusters[kk].size(), init_new_clusters[kk].size(), init_new_clusters[kk], init_new_clusters[kk]);
        new_Connected_Components(A_tmp, init_new_clusters[kk].size(), init_new_clusters[kk], connected_components);
        for(int new_k = 0; new_k < connected_components.size(); new_k ++) tmp_new_clusters.push_back(connected_components[new_k]);
      }
      new_clusters.clear();
      new_clusters.resize(tmp_new_clusters.size());
      k_star.clear();
      k_star.resize(tmp_new_clusters.size());
      get_subcluster_neighbor(tmp_new_clusters, new_clusters, k_star, k, gamma_l, T, A_block, rho, a1, a2);
      si.split_k.push_back(k);
      si.new_clusters.push_back(new_clusters);
      si.nearest_neighbor.push_back(k_star);
      si.num_splits++;
    } // closes loop checking that there is at least 1 element in the cluster
  } // closes loop over clusters

}

void get_km_split(split_info &si, LPPartition gamma_l, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const int reps)
{
  //Rcpp::Rcout << "[get_km_split]: Entering get_km_split" << endl;
  // re-set everything in si
  si.num_splits = 0;
  si.split_k.clear();
  for(int i = 0; i < si.new_clusters.size(); i++){
    si.new_clusters[i].clear();
    si.nearest_neighbor[i].clear();
  }
  
  si.new_clusters.clear();
  si.nearest_neighbor.clear();
  
  int orig_K = gamma_l->K;
  int max_splits = 5;
  int split_k = 0;
  int n_k = 1;
  
  // initialize the stuff needed for k-means clustering
  arma::mat U = arma::zeros<mat>(n_k, 1); // holds the data passed to k-means
  arma::mat means = arma::zeros<mat>(n_k,1); // holds data passed to k-means
  double min_score = 0.0; // hold sum-of-squared distances from observations to cluster means
  bool status = true; // tells us whether k-means was successful
  
  arma::vec cluster_dist = arma::zeros<vec>(1); // holds distance to cluster means
  arma::uvec cluster_dist_indices(1); // used to sort cluster_dist
  
  std::vector<std::vector<int> > init_new_clusters; // stores output of k-means. NOTE: these clusters may not be connected.
  arma::mat A_tmp = arma::zeros<arma::mat>(n_k,n_k); // store submatrices of A_block. used to determine connected components
  std::vector<std::vector<int> > connected_components; // connected components of the individual sub clusters found by k-means
  
  std::vector<std::vector<int> > tmp_new_clusters; // used for figuring out the connected components of the new clusters
  std::vector<std::vector<int> > new_clusters; // holds the final new clusters
  std::vector<int> k_star; // actual nearest neighbors of newly formed clusters
  
  
  for(int k = 0; k < orig_K; k++){
    //Rcpp::Rcout << "[get_km_split]:  k = " << k << endl;
    n_k = gamma_l->cluster_config[k];
    split_k = k;
    if(n_k > 1){
      if(sqrt(n_k) < 5) max_splits = floor(sqrt(n_k));
      else max_splits = 5;
      for(int num_splits = 2; num_splits <= max_splits; num_splits++){
        U.set_size(n_k, 1);
        means.set_size(1, num_splits); // one thing to do might be to change to means.zeros(1, num_split).
        min_score = 0.0;
        status = true;
        cluster_dist.set_size(num_splits);
        cluster_dist_indices.set_size(num_splits);
        //Rcpp::Rcout << "    set U, means, status " << endl;
        
        // do k-means
        for(int ii = 0; ii < n_k; ii++){
          U(ii,0) = gamma_l->alpha_hat[gamma_l->clusters[split_k][ii]];
          if(U(ii,0) != U(ii,0)){
            Rcpp::Rcout << "[get_km_split]:    possible non-finite value in U" << endl;
            Rcpp::Rcout << "                   split_k = " << split_k << "  ii = " << ii << endl;
            Rcpp::Rcout << gamma_l->clusters[split_k][ii] << endl;
          }
        }
        //Rcpp::Rcout << "Ready to start k-means" << endl;
        kmeans_repeat(U, means, status, min_score, num_splits, reps);
        if(status == false){
          Rcpp::Rcout << "kmeans failed!! n_k = " << n_k << " num_splits = " << num_splits << endl;
        } else{
          //Rcpp::Rcout << "[get_km_split]:     min_score = " << min_score << endl;
          init_new_clusters.clear();
          init_new_clusters.resize(num_splits);
          for(int i = 0; i < n_k; i++){
            cluster_dist.zeros();
            for(int new_k = 0; new_k < num_splits; new_k++){
              cluster_dist(new_k) = arma::norm(U.row(i).t() - means.col(new_k)); // distance from newly created cluster kk to point i
            }
            cluster_dist_indices = arma::sort_index(cluster_dist, "ascend");
            init_new_clusters[cluster_dist_indices(0)].push_back(gamma_l->clusters[split_k][i]);
          } // closes loop over the elements of the original cluster to find new subcluster assignment
          tmp_new_clusters.clear();
          for(int kk = 0; kk < init_new_clusters.size(); kk++){
            connected_components.clear();
            A_tmp = Submatrix(A_block, init_new_clusters[kk].size(), init_new_clusters[kk].size(), init_new_clusters[kk], init_new_clusters[kk]); // submatrix of A_block corresponds to newly discovered sub-cluster
            new_Connected_Components(A_tmp, init_new_clusters[kk].size(), init_new_clusters[kk], connected_components);
            for(int new_k = 0; new_k < connected_components.size(); new_k++){
              tmp_new_clusters.push_back(connected_components[new_k]);
            }
          } // closes loop over the new clusters discovered
          new_clusters.clear();
          new_clusters.resize(tmp_new_clusters.size());
          k_star.clear();
          k_star.resize(tmp_new_clusters.size());
          get_subcluster_neighbor(tmp_new_clusters, new_clusters, k_star, k, gamma_l, T, A_block, rho, a1, a2);
          //Rcpp::Rcout << "[get_km_split]:    num_splits = " << num_splits << " num sub-clusters = " << new_clusters.size() << " : " ;
          //for(int new_k = 0; new_k < new_clusters.size(); new_k++) Rcpp::Rcout << " " << new_clusters[new_k].size() ;
          //Rcpp::Rcout << endl;
          // write split information into si.
          si.split_k.push_back(k);
          si.new_clusters.push_back(new_clusters);
          si.nearest_neighbor.push_back(k_star);
          si.num_splits++;
        } // closes if/else checking whether kmeans_repeat failed
      } // closes loop over num_splits
    } // closes if checking that n_k > 1
  } // closes loop over the clusters of gamma_l
}



// A new version of best_split: the last argument, merge_flag, is used to determine whether we should try merging new sub-clusters with their nearest neighbors
void best_split(split_info &si, LPPartition candidate, const int current_l, const std::vector<LPPartition> particle_set, const std::vector<double> w, const arma::vec &ybar, const double total_ss, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double nu_sigma, const double lambda_sigma, const double eta, const double lambda, const bool merge_flag)
{
  if(si.num_splits > 0){
    LPPartition max_candidate = new Partition(particle_set[current_l]);
    double max_objective = w[current_l]*total_log_post(max_candidate, total_ss, T, nu_sigma, lambda_sigma) + lambda*Entropy(current_l, max_candidate, particle_set, w);
    LPPartition tmp_candidate = new Partition(particle_set[current_l]);
    double tmp_objective = 0.0;
    
    
    
    std::vector<int> k_star; // holds indices of nearest neighbors
    int num_new_clusters = 0; // number of new sub-clusters created by the split
    bool sanity_flag = true;
    
    
    if(merge_flag == false){
      // we will never perform progressive merges. instead we split and merge according to what is stored in si.nearest_neighbor
      for(int split_ix = 0; split_ix < si.num_splits; split_ix++){
        delete tmp_candidate;
        tmp_candidate = new Partition(particle_set[current_l]);
        tmp_candidate->Split_Merge(si.split_k[split_ix], si.new_clusters[split_ix], si.nearest_neighbor[split_ix], ybar, T, A_block, rho, a1, a2, eta);
        sanity_flag = sanity_check(tmp_candidate);
        if(sanity_flag == false){
          Rcpp::Rcout << "[best_split]: Before even trying merges, something is wrong about this split!" << endl;
          Rcpp::Rcout << "   Attempted to split cluster " << si.split_k[split_ix] << " into " << si.new_clusters[split_ix].size() << " parts " << endl;
          for(int nc_ix = 0; nc_ix < si.new_clusters[split_ix].size(); nc_ix++){
            Rcpp::Rcout << "    new cluster " << nc_ix << " of size " << si.new_clusters[split_ix][nc_ix].size() << " and neighbor " << si.nearest_neighbor[split_ix][nc_ix] << " : " << endl;
            for(int ii = 0; ii < si.new_clusters[split_ix][nc_ix].size(); ii++){
              Rcpp::Rcout << si.new_clusters[split_ix][nc_ix][ii] << " " ;
            }
            Rcpp::Rcout << endl;
          }
          Rcpp::Rcout << "Before split candidate is" << endl;
          particle_set[current_l]->Print_Partition(total_ss, T, nu_sigma, lambda_sigma);
          Rcpp::Rcout << "The candidate is now" << endl;
          tmp_candidate->Print_Partition(total_ss, T, nu_sigma, lambda_sigma);
          Rcpp::stop("Terminating");
        } // closes if checking whether there is something wrong with the partition
        tmp_objective = w[current_l]*total_log_post(tmp_candidate, total_ss, T, nu_sigma, lambda_sigma) + lambda*Entropy(current_l, tmp_candidate, particle_set, w);
        if(tmp_objective > max_objective){
          delete max_candidate;
          max_candidate = new Partition(tmp_candidate);
          max_objective = tmp_objective;
        }
      } // closes loop over all of the possible splits
    } else{
      std::vector<int> k_star; // holds the progessively changing indices of the newly formed clusters
      for(int split_ix = 0; split_ix < si.num_splits; split_ix++){
        num_new_clusters = si.new_clusters[split_ix].size();
        k_star.clear();
        k_star.resize(num_new_clusters,-1);
        
        // first try to split all of the candidates but don't do any merges
        delete tmp_candidate;
        tmp_candidate = new Partition(particle_set[current_l]);
        tmp_candidate->Split_Merge(si.split_k[split_ix], si.new_clusters[split_ix], k_star, ybar, T, A_block, rho, a1, a2, eta);
        sanity_flag = sanity_check(tmp_candidate);
        if(sanity_flag == false){
          Rcpp::Rcout << "[best_split]: Before even trying merges, something is wrong about this split!" << endl;
          Rcpp::Rcout << "   Attempted to split cluster " << si.split_k[split_ix] << " into " << si.new_clusters[split_ix].size() << " parts " << endl;
          for(int nc_ix = 0; nc_ix < si.new_clusters[split_ix].size(); nc_ix++){
            Rcpp::Rcout << "    new cluster " << nc_ix << " of size " << si.new_clusters[split_ix][nc_ix].size() << " and neighbor " << si.nearest_neighbor[split_ix][nc_ix] << " : " << endl;
            for(int ii = 0; ii < si.new_clusters[split_ix][nc_ix].size(); ii++){
              Rcpp::Rcout << si.new_clusters[split_ix][nc_ix][ii] << " " ;
            }
            Rcpp::Rcout << endl;
          }
          Rcpp::Rcout << "Before split candidate is" << endl;
          particle_set[current_l]->Print_Partition(total_ss, T, nu_sigma, lambda_sigma);
          Rcpp::Rcout << "The candidate is now" << endl;
          tmp_candidate->Print_Partition(total_ss, T, nu_sigma, lambda_sigma);
          Rcpp::stop("Terminating");
        } // closes if checking whether there is something wrong with the partition
        tmp_objective = w[current_l]*total_log_post(tmp_candidate, total_ss, T, nu_sigma, lambda_sigma) + lambda*Entropy(current_l, tmp_candidate, particle_set, w);
        if(tmp_objective > max_objective){
          delete max_candidate;
          max_candidate = new Partition(tmp_candidate);
          max_objective = tmp_objective;
        }
        // now try to do all of the proposed merges
        for(int nc_ix = 0; nc_ix < num_new_clusters; nc_ix++){
          if(si.nearest_neighbor[split_ix][nc_ix] != -1){
            k_star[nc_ix] = si.nearest_neighbor[split_ix][nc_ix];
            delete tmp_candidate;
            tmp_candidate = new Partition(particle_set[current_l]);
            tmp_candidate->Split_Merge(si.split_k[split_ix], si.new_clusters[split_ix], k_star, ybar, T, A_block, rho, a1, a2, eta);
            sanity_flag = sanity_check(tmp_candidate);
            if(sanity_flag == false){
              Rcpp::Rcout << "[best_split]: Before even trying merges, something is wrong about this split!" << endl;
              Rcpp::Rcout << "   Attempted to split cluster " << si.split_k[split_ix] << " into " << si.new_clusters[split_ix].size() << " parts " << endl;
              for(int nc_ix = 0; nc_ix < si.new_clusters[split_ix].size(); nc_ix++){
                Rcpp::Rcout << "    new cluster " << nc_ix << " of size " << si.new_clusters[split_ix][nc_ix].size() << " and neighbor " << si.nearest_neighbor[split_ix][nc_ix] << " : " << endl;
                for(int ii = 0; ii < si.new_clusters[split_ix][nc_ix].size(); ii++){
                  Rcpp::Rcout << si.new_clusters[split_ix][nc_ix][ii] << " " ;
                }
                Rcpp::Rcout << endl;
              }
              Rcpp::Rcout << "Before split candidate is" << endl;
              particle_set[current_l]->Print_Partition(total_ss, T, nu_sigma, lambda_sigma);
              Rcpp::Rcout << "The candidate is now" << endl;
              tmp_candidate->Print_Partition(total_ss, T, nu_sigma, lambda_sigma);
              Rcpp::stop("Terminating");
            } // closes if checking whether there is something wrong with the partition
            tmp_objective = w[current_l]*total_log_post(tmp_candidate, total_ss, T, nu_sigma, lambda_sigma) + lambda*Entropy(current_l, tmp_candidate, particle_set, w);
            if(tmp_objective > max_objective){
              delete max_candidate;
              max_candidate = new Partition(tmp_candidate);
              max_objective = tmp_objective;
            }
          } // closes if checking to see if we need to update k_star
        } // closes loop over all of the new sub-clusters
      } // closes loop over each potential split
    } // close else part checking whether we do progressive merges or not
    candidate->Copy_Partition(max_candidate);
  } else{
    candidate->Copy_Partition(particle_set[current_l]);
  } // closes if/else checking that we actually try splits at all
  
}

// a version of best_split for MAP estimation
void best_split_map(split_info &si, LPPartition candidate, LPPartition init_particle, const arma::vec &ybar, const double total_ss, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double nu_sigma, const double lambda_sigma, const double eta, const bool merge_flag)
{
  if(si.num_splits > 0){
    LPPartition max_candidate = new Partition(init_particle);
    double max_objective = total_log_post(max_candidate, total_ss, T, nu_sigma, lambda_sigma);
    LPPartition tmp_candidate = new Partition(init_particle);
    double tmp_objective = 0.0;
    
    std::vector<int> k_star; // holds indices of nearest neighbors
    int num_new_clusters = 0; // number of new sub-clusters created by the split
    bool sanity_flag = true;
    
    
    if(merge_flag == false){
      // we will never perform progressive merges. instead we split and merge according to what is stored in si.nearest_neighbor
      for(int split_ix = 0; split_ix < si.num_splits; split_ix++){
        delete tmp_candidate;
        tmp_candidate = new Partition(init_particle);
        tmp_candidate->Split_Merge(si.split_k[split_ix], si.new_clusters[split_ix], si.nearest_neighbor[split_ix], ybar, T, A_block, rho, a1, a2, eta);
        sanity_flag = sanity_check(tmp_candidate);
        if(sanity_flag == false){
          Rcpp::Rcout << "[best_split]: Before even trying merges, something is wrong about this split!" << endl;
          Rcpp::Rcout << "   Attempted to split cluster " << si.split_k[split_ix] << " into " << si.new_clusters[split_ix].size() << " parts " << endl;
          for(int nc_ix = 0; nc_ix < si.new_clusters[split_ix].size(); nc_ix++){
            Rcpp::Rcout << "    new cluster " << nc_ix << " of size " << si.new_clusters[split_ix][nc_ix].size() << " and neighbor " << si.nearest_neighbor[split_ix][nc_ix] << " : " << endl;
            for(int ii = 0; ii < si.new_clusters[split_ix][nc_ix].size(); ii++){
              Rcpp::Rcout << si.new_clusters[split_ix][nc_ix][ii] << " " ;
            }
            Rcpp::Rcout << endl;
          }
          Rcpp::Rcout << "Before split candidate is" << endl;
          init_particle->Print_Partition(total_ss, T, nu_sigma, lambda_sigma);
          Rcpp::Rcout << "The candidate is now" << endl;
          tmp_candidate->Print_Partition(total_ss, T, nu_sigma, lambda_sigma);
          Rcpp::stop("Terminating");
        } // closes if checking whether there is something wrong with the partition
        tmp_objective = total_log_post(tmp_candidate, total_ss, T, nu_sigma, lambda_sigma);
        if(tmp_objective > max_objective){
          delete max_candidate;
          max_candidate = new Partition(tmp_candidate);
          max_objective = tmp_objective;
        }
      } // closes loop over all of the possible splits
    } else{
      std::vector<int> k_star; // holds the progessively changing indices of the newly formed clusters
      for(int split_ix = 0; split_ix < si.num_splits; split_ix++){
        num_new_clusters = si.new_clusters[split_ix].size();
        k_star.clear();
        k_star.resize(num_new_clusters,-1);
        
        // first try to split all of the candidates but don't do any merges
        delete tmp_candidate;
        tmp_candidate = new Partition(init_particle);
        tmp_candidate->Split_Merge(si.split_k[split_ix], si.new_clusters[split_ix], k_star, ybar, T, A_block, rho, a1, a2, eta);
        sanity_flag = sanity_check(tmp_candidate);
        if(sanity_flag == false){
          Rcpp::Rcout << "[best_split]: Before even trying merges, something is wrong about this split!" << endl;
          Rcpp::Rcout << "   Attempted to split cluster " << si.split_k[split_ix] << " into " << si.new_clusters[split_ix].size() << " parts " << endl;
          for(int nc_ix = 0; nc_ix < si.new_clusters[split_ix].size(); nc_ix++){
            Rcpp::Rcout << "    new cluster " << nc_ix << " of size " << si.new_clusters[split_ix][nc_ix].size() << " and neighbor " << si.nearest_neighbor[split_ix][nc_ix] << " : " << endl;
            for(int ii = 0; ii < si.new_clusters[split_ix][nc_ix].size(); ii++){
              Rcpp::Rcout << si.new_clusters[split_ix][nc_ix][ii] << " " ;
            }
            Rcpp::Rcout << endl;
          }
          Rcpp::Rcout << "Before split candidate is" << endl;
          init_particle->Print_Partition(total_ss, T, nu_sigma, lambda_sigma);
          Rcpp::Rcout << "The candidate is now" << endl;
          tmp_candidate->Print_Partition(total_ss, T, nu_sigma, lambda_sigma);
          Rcpp::stop("Terminating");
        } // closes if checking whether there is something wrong with the partition
        tmp_objective = total_log_post(tmp_candidate, total_ss, T, nu_sigma, lambda_sigma);
        if(tmp_objective > max_objective){
          delete max_candidate;
          max_candidate = new Partition(tmp_candidate);
          max_objective = tmp_objective;
        }
        // now try to do all of the proposed merges
        for(int nc_ix = 0; nc_ix < num_new_clusters; nc_ix++){
          if(si.nearest_neighbor[split_ix][nc_ix] != -1){
            k_star[nc_ix] = si.nearest_neighbor[split_ix][nc_ix];
            delete tmp_candidate;
            tmp_candidate = new Partition(init_particle);
            tmp_candidate->Split_Merge(si.split_k[split_ix], si.new_clusters[split_ix], k_star, ybar, T, A_block, rho, a1, a2, eta);
            sanity_flag = sanity_check(tmp_candidate);
            if(sanity_flag == false){
              Rcpp::Rcout << "[best_split]: Before even trying merges, something is wrong about this split!" << endl;
              Rcpp::Rcout << "   Attempted to split cluster " << si.split_k[split_ix] << " into " << si.new_clusters[split_ix].size() << " parts " << endl;
              for(int nc_ix = 0; nc_ix < si.new_clusters[split_ix].size(); nc_ix++){
                Rcpp::Rcout << "    new cluster " << nc_ix << " of size " << si.new_clusters[split_ix][nc_ix].size() << " and neighbor " << si.nearest_neighbor[split_ix][nc_ix] << " : " << endl;
                for(int ii = 0; ii < si.new_clusters[split_ix][nc_ix].size(); ii++){
                  Rcpp::Rcout << si.new_clusters[split_ix][nc_ix][ii] << " " ;
                }
                Rcpp::Rcout << endl;
              }
              Rcpp::Rcout << "Before split candidate is" << endl;
              init_particle->Print_Partition(total_ss, T, nu_sigma, lambda_sigma);
              Rcpp::Rcout << "The candidate is now" << endl;
              tmp_candidate->Print_Partition(total_ss, T, nu_sigma, lambda_sigma);
              Rcpp::stop("Terminating");
            } // closes if checking whether there is something wrong with the partition
            tmp_objective = total_log_post(tmp_candidate, total_ss, T, nu_sigma, lambda_sigma);
            if(tmp_objective > max_objective){
              delete max_candidate;
              max_candidate = new Partition(tmp_candidate);
              max_objective = tmp_objective;
            }
          } // closes if checking to see if we need to update k_star
        } // closes loop over all of the new sub-clusters
      } // closes loop over each potential split
    } // close else part checking whether we do progressive merges or not
    candidate->Copy_Partition(max_candidate);
  } else{
    candidate->Copy_Partition(init_particle);
  } // closes if/else checking that we actually try splits at all
  
}

void best_merge(merge_info &mi, LPPartition candidate, const int current_l, const std::vector<LPPartition> particle_set, const std::vector<double> w, const arma::vec &ybar, const double total_ss, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double nu_sigma, const double lambda_sigma, const double eta, const double lambda)
{
  
  if(mi.num_merges > 0){
    LPPartition max_candidate = new Partition(particle_set[current_l]);
    //double max_objective = w[current_l]*total_log_post(max_candidate, a_sigma, nu_sigma) + lambda * Entropy(current_l, max_candidate, particle_set, w);
    double max_objective = w[current_l]*total_log_post(max_candidate, total_ss, T, nu_sigma, lambda_sigma) + lambda * Entropy(current_l, max_candidate, particle_set, w);
    LPPartition tmp_candidate = new Partition(particle_set[current_l]); // the running candidate
    double tmp_objective = 0.0;
    
    for(int m_ix = 0; m_ix < mi.num_merges; m_ix++){
      delete tmp_candidate;
      tmp_candidate = new Partition(particle_set[current_l]);
      tmp_candidate->Merge(mi.rec_k[m_ix], mi.donor_k[m_ix], ybar, T, A_block, rho, a1, a2, eta);
      tmp_objective = w[current_l]*total_log_post(tmp_candidate, total_ss, T, nu_sigma, lambda_sigma) + lambda*Entropy(current_l, tmp_candidate, particle_set, w);
      if(tmp_objective > max_objective){
        //Rcpp::Rcout << "[best merge]: max_candidate." << endl;
        delete max_candidate;
        max_candidate = new Partition(tmp_candidate);
        max_objective = tmp_objective;
      }
    }
    candidate->Copy_Partition(max_candidate);
  } else{
    candidate->Copy_Partition(particle_set[current_l]);
  }
  
}


void best_merge_map(merge_info &mi, LPPartition candidate, LPPartition init_particle, const arma::vec &ybar, const double total_ss, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double nu_sigma, const double lambda_sigma, const double eta)
{
  
  if(mi.num_merges > 0){
    LPPartition max_candidate = new Partition(init_particle);
    double max_objective = total_log_post(max_candidate, total_ss, T, nu_sigma, lambda_sigma);
    LPPartition tmp_candidate = new Partition(init_particle); // the running candidate
    double tmp_objective = 0.0;
    
    for(int m_ix = 0; m_ix < mi.num_merges; m_ix++){
      delete tmp_candidate;
      tmp_candidate = new Partition(init_particle);
      tmp_candidate->Merge(mi.rec_k[m_ix], mi.donor_k[m_ix], ybar, T, A_block, rho, a1, a2, eta);
      tmp_objective = total_log_post(tmp_candidate, total_ss, T, nu_sigma, lambda_sigma);
      if(tmp_objective > max_objective){
        //Rcpp::Rcout << "[best merge]: max_candidate." << endl;
        delete max_candidate;
        max_candidate = new Partition(tmp_candidate);
        max_objective = tmp_objective;
      }
    }
    candidate->Copy_Partition(max_candidate);
  } else{
    candidate->Copy_Partition(init_particle);
  }
  
}

bool sanity_check(LPPartition partition){
  bool flag = true;
  int K = partition->K;
  int n = partition->nObs;

  // check that the number of elements in all of the clusters adds up to n
  int running_count = 0;
  for(int k = 0; k < K; k++){
    running_count += partition->cluster_config[k];
  }
  if(running_count != n){
    flag = false;
  } else{
    //Rcpp::Rcout << "Sum of cluster sizes is not n!" << endl; this gets printed way too many times
    flag = true;
  }
  
  // check that we haven't duplicated an index
  std::vector<int> duplicate(partition->nObs,-1);
  for(int i = 0; i < partition->nObs; i++){
    duplicate[i] = -1;
  }
  for(int k = 0; k < K; k++){
    for(int ix = 0; ix < partition->cluster_config[k]; ix++){
      if(duplicate[partition->clusters[k][ix]] == -1) duplicate[partition->clusters[k][ix]] = 1;
      else{
        //Rcpp::Rcout << "Particle has duplicate indicies" << endl; this gets printed way too many times
        flag = false;
      }
    }
  }
  return flag;
}

// this function will find nearest neighbors of each new subcluster and also sort them based on how close they are to their neighbors
void get_subcluster_neighbor(std::vector<std::vector<int> > &init_new_clusters, std::vector<std::vector<int> > &new_clusters, std::vector<int> &k_star, const int split_k, LPPartition gamma_l, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2)
{
  int num_new_clusters = init_new_clusters.size();
  new_clusters.clear();
  new_clusters.resize(num_new_clusters);
  k_star.clear();
  k_star.resize(num_new_clusters,-1);
  std::vector<int> tmp_k_star(num_new_clusters,-1);
  arma::mat A_tmp = arma::zeros<arma::mat>(num_new_clusters, num_new_clusters);
  double tmp_alpha_bar = 0.0; // holds the new alphabar estiate for subcluster
  std::vector<int> tmp_nn; // holds potential nearest neighbors
  std::vector<double> tmp_dist; // holds distance to potential nearest neighbor
  arma::vec tmp_dist_vec = arma::zeros<vec>(1); // for sorting distances to nearest neighbors for single subcluster
  arma::uvec tmp_dist_indices(1); // for getting the index after sorting
  
  arma::vec dist_vec = arma::zeros<vec>(num_new_clusters); // holds distances from each subcluster to nearest neighbor
  arma::uvec dist_indices(num_new_clusters);
  
  // we need to handle the case where there are no other existing sub-clusters
  int orig_K = gamma_l->K;
  if(orig_K > 1){
    // now it makes sense to re-order the new sub-cluster by their distance to their neighbors
    for(int new_k = 0; new_k < num_new_clusters; new_k++){
      tmp_alpha_bar = alpha_bar_func(init_new_clusters[new_k], gamma_l, T, A_block, rho, a1, a2);
      tmp_nn.clear();
      tmp_dist.clear();
      for(int kk = 0; kk < orig_K; kk++){
        if(kk != split_k){
          A_tmp = Submatrix(A_block, init_new_clusters[new_k].size(), gamma_l->cluster_config[kk], init_new_clusters[new_k], gamma_l->clusters[kk]);
          if(any(vectorise(A_tmp) == 1)){
            tmp_nn.push_back(kk);
            tmp_dist.push_back(abs(tmp_alpha_bar - gamma_l->alpha_bar[kk]));
          }
        } // closes if checking that kk != split_k
      } // closes loop over the original cluster ids (kk)
      if(tmp_nn.size() > 0){
        // new subcluster new_k is adjacent to an existing cluster
        tmp_dist_vec.reset();
        tmp_dist_indices.reset();
        tmp_dist_vec.set_size(tmp_dist.size());
        tmp_dist_indices.set_size(tmp_dist.size());
        for(int kk = 0; kk < tmp_dist.size(); kk++){
          tmp_dist_vec(kk) = tmp_dist[kk];
        }
        tmp_dist_indices = arma::sort_index(tmp_dist_vec, "ascend");
        tmp_k_star[new_k] = tmp_nn[tmp_dist_indices(0)];
        dist_vec(new_k) = tmp_dist_vec(tmp_dist_indices(0)); // distance from subcluster new_k to its nearest neighbor
      } else{
        // new subcluster new_k is adjacent to no other cluster
        // we will set dist[new_k] = 0.0 so that when we sort them, these subclusters come first in the order
        tmp_k_star[new_k] = -1;
        dist_vec(new_k) = 0.0;
      }
    } // closes loop over the new sub-clusters
    // now the vector dist_vec contains the distance of each subcluster to nearest neighbor
    dist_indices = arma::sort_index(dist_vec, "ascend");
    for(int new_k = 0; new_k < num_new_clusters; new_k++){
      // we want get init_new_clusters[dist_indices(new_k)]
      for(int ii = 0; ii < init_new_clusters[dist_indices(new_k)].size(); ii++){
        new_clusters[new_k].push_back(init_new_clusters[dist_indices(new_k)][ii]);
      } // closes loop over the elements in the subcluster
      k_star[new_k] = tmp_k_star[dist_indices(new_k)];
    } // closes loop over the new subclusters
  } else {
    // gamma_l consisted of only one block. we may as well arrange the clusters based on their variance
    arma::vec alpha_hat_cluster = arma::zeros<arma::vec>(1); // will hold the alpha_hats from each clusters
    arma::vec cluster_var = arma::zeros<arma::vec>(num_new_clusters);
    arma::uvec var_indices(num_new_clusters);
    for(int new_k = 0; new_k < num_new_clusters; new_k++){
      alpha_hat_cluster.reset();
      alpha_hat_cluster.set_size(init_new_clusters[new_k].size());
      for(int ii = 0; ii < init_new_clusters[new_k].size(); ii++){
        alpha_hat_cluster(ii) = gamma_l->alpha_hat[init_new_clusters[new_k][ii]];
      }
      cluster_var[new_k] = arma::var(alpha_hat_cluster);
    }
    var_indices = arma::sort_index(cluster_var, "descend");
    for(int new_k = 0; new_k < num_new_clusters; new_k++){
      for(int ii = 0; ii < init_new_clusters[var_indices(new_k)].size(); ii++){
        new_clusters[new_k].push_back(init_new_clusters[var_indices(new_k)][ii]);
      }
      k_star[new_k] = -1;
    }
  } // closes else part checking whether original partition has one cluster or not
}


void format_particle_set(std::vector<LPPartition> particle_set, Rcpp::List &output_list)
{
  arma::vec tmp_vec = arma::zeros<vec>(1);
  Rcpp::NumericVector output_vec;
  Rcpp::List tmp_list;
  output_list = Rcpp::List(particle_set.size());
  
  
  for(int l = 0; l < particle_set.size(); l++){
    tmp_list = Rcpp::List(particle_set[l]->K);
    for(int k = 0; k < particle_set[l]->K; k++){
      tmp_vec = arma::zeros<vec>(particle_set[l]->cluster_config[k]);
      for(int i = 0; i  < particle_set[l]->cluster_config[k]; i++) tmp_vec(i) = particle_set[l]->clusters[k][i] + 1; // Remember that R is 1-indexed
      output_vec = Rcpp::wrap(arma::sort(tmp_vec, "ascend"));
      output_vec.attr("dim") = R_NilValue;
      tmp_list[k] = output_vec;
    }
    output_list[l] = tmp_list;
  }
  
}
