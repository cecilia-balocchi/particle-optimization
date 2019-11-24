//
//  initialize_particle_set.cpp
//  
//


#include "initialize_particle_set.h"


void initialize_particle_set(std::vector<LPPartition> &particle_set, const int L, const arma::vec ybar, const double total_ss, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double eta, const double nu_sigma, const double lambda_sigma, const int init_id, RNG &gen)
{
  
  int n = ybar.size();
  
  // create a partition with all elements in a single cluster
  std::vector<std::vector<int> > tmp_clusters(1);
  for(int i = 0; i < n; i++) tmp_clusters[0].push_back(i);
  
  LPPartition gamma_init = new Partition(n, tmp_clusters, ybar, T, A_block, rho, a1, a2, eta);
  if(init_id == 0){
    for(int l = 0; l < L; l++) particle_set[l]->Copy_Partition(gamma_init);
  } else{
    // run k-means
    int max_splits = floor(log( (double) n));
    int num_copies = floor( (double) L / ( (double) max_splits));
    int excess = L - max_splits * num_copies; // if excess > 1, we make an additional copy of the partitions with highest probability
    int n_k = n;
    arma::mat U = arma::zeros<arma::mat>(n_k, 1);
    arma::mat means = arma::zeros<arma::mat>(1, max_splits);
    double min_score = 0.0;
    bool status = true;
    
    arma::vec cluster_dist = arma::zeros<arma::vec>(max_splits); // holds distance from single observaiton to each cluster mean
    arma::uvec cluster_dist_indices(max_splits); // used to sort cluster_dist
    
    std::vector<std::vector<int> > init_new_clusters; // stores output of k-means. NOTE: these clusters may not be connected.
    arma::mat A_tmp = arma::zeros<arma::mat>(n_k,n_k); // store submatrices of A_block. used to determine connected components
    std::vector<std::vector<int> > connected_components; // connected components of the individual sub clusters found by k-means
    
    std::vector<std::vector<int> > tmp_new_clusters; // used for figuring out the connected components of the new clusters
    std::vector<std::vector<int> > new_clusters(1); // holds the final new clusters
    std::vector<int> k_star; // actual nearest neighbors of newly formed clusters
    
    arma::vec init_log_post(max_splits);
    arma::uvec init_log_post_indices(max_splits); // for sorting the log posterior values
    std::vector<double> log_post_prob(max_splits, 0.0);
    double tmp_log_post = 0.0;
    double max_log_post = 0.0;
    
    
    LPPartition gamma_new = new Partition();
    std::vector<LPPartition> unik_particles(max_splits); // holds the results from running K-mean
    
    for(int k = 0; k <= max_splits-1; k++) unik_particles[k] = new Partition();
    //Rcpp::Rcout << "initialized unik_particles." << unik_particles.size() << std::endl;
    
    // now we run k-means
    for(int num_splits = 1; num_splits <= max_splits; num_splits++){
      //Rcpp::Rcout << "num_splits = " << num_splits << std::endl;
      if(num_splits > 1){
        U.set_size(n_k,1);
        means.set_size(1, num_splits);
        min_score = 0.0;
        status = true;
        cluster_dist.set_size(num_splits);
        cluster_dist_indices.set_size(num_splits);
        for(int ii = 0; ii < n_k; ii++){
          //Rcpp::Rcout << " " << ii;
          U(ii,0) = ybar(ii);
        }
        kmeans_repeat(U, means, status, min_score, num_splits, 5000);
        if(status == true){
          init_new_clusters.clear();
          init_new_clusters.resize(num_splits);
          for(int ii = 0; ii < n_k; ii++){
            cluster_dist.zeros();
            for(int new_k = 0; new_k < num_splits; new_k++) cluster_dist(new_k) = arma::norm(U.row(ii).t() - means.col(new_k));
            cluster_dist_indices = arma::sort_index(cluster_dist, "ascend");
            init_new_clusters[cluster_dist_indices(0)].push_back(ii);
          }
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
          gamma_new = new Partition(n, new_clusters, ybar, T, A_block, rho, a1, a2, eta);
          //Rcpp::Rcout << "  updated gamma_new" << std::endl;
          unik_particles[num_splits-1]->Copy_Partition(gamma_new);
          tmp_log_post = total_log_post(unik_particles[num_splits-1], total_ss, T, nu_sigma, lambda_sigma);
          init_log_post(num_splits-1) = tmp_log_post;
          log_post_prob[num_splits - 1] = tmp_log_post;
          if(tmp_log_post > max_log_post) max_log_post = tmp_log_post;
        } // closes if/else checking that kmeans ran successfull
        
      } else{
        new_clusters.clear();
        new_clusters.resize(1);
        for(int i = 0; i < n; i++) new_clusters[0].push_back(i);
        delete gamma_new;
        gamma_new = new Partition(n, new_clusters, ybar, T, A_block, rho, a1, a2, eta);
        unik_particles[num_splits-1]->Copy_Partition(gamma_new);
        tmp_log_post = total_log_post(unik_particles[num_splits-1], total_ss, T, nu_sigma, lambda_sigma);
        init_log_post(num_splits-1) = tmp_log_post;
        log_post_prob[num_splits - 1] = tmp_log_post;
        max_log_post = tmp_log_post; // this is the first one encountered
        
      }
    } // closes loop running k-means
    //Rcpp::Rcout << "Finished k-means" << std::endl;
    // now we order by the log_post values
    init_log_post_indices = arma::sort_index(init_log_post, "descend");
    int l = 0;
    if(init_id == 1){
      // make copies of each partition found by k-means
      for(int k = 0; k < max_splits; k++){
        for(int ix = 0; ix < num_copies; ix++){
          particle_set[l]->Copy_Partition(unik_particles[init_log_post_indices(k)]);
          l++;
        }
        if(k < excess){
          // make an extra copy of the initializers with highest probability if needed
          particle_set[l]->Copy_Partition(unik_particles[init_log_post_indices(k)]);
          l++;
        }
      }
    } else{
      std::vector<double> post_prob(max_splits,0.0);
      std::vector<double> cdf(max_splits,0.0);
      double tmp_sum = 0.0;
      for(size_t ix = 0; ix < max_splits; ix++){
        log_post_prob[ix] -= max_log_post;
        post_prob[ix] = exp(log_post_prob[ix]);
        tmp_sum += post_prob[ix];
        cdf[ix] = tmp_sum;
      }
      for(size_t ix = 0; ix < max_splits; ix++){
        post_prob[ix] /= tmp_sum; // normalizes the CDF
        cdf[ix] /= tmp_sum ; // normalizes the CDF
      }
      
      std::vector<int> possible_ix;
      double unif = 0.0;
      int gamma_ix = 0;
      
      for(int l  = 0; l < L; l++){
        unif = gen.uniform(); // draw the uniform variable
        possible_ix.clear();
        for(int ix = 0; ix < max_splits; ix++){
          if(cdf[ix] <= unif) possible_ix.push_back(ix);
        }
        
        if(possible_ix.size() == 0) gamma_ix = 0; // if we somehow generated unif = 0...
        else if(possible_ix.size() == max_splits) gamma_ix = max_splits-1; // if we somehow generated unif == 1
        else{
          gamma_ix = possible_ix[possible_ix.size() - 1] + 1; // cdf[gamma_ix] < unif < cdf[gamma_ix + 1]
        }
        particle_set[l]->Copy_Partition(unik_particles[gamma_ix]);

        //Rcpp::Rcout << "l = " << l << "unif = " << unif << " gamma_ix = " << gamma_ix << std::endl;
        
        
        
      }
    } // closes if/else checking whether init_id == 1
  } // closes if/else checking whether init_id == 0
}

void initialize_particle(LPPartition &init_gamma, const arma::vec ybar, const double total_ss, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double eta, const double nu_sigma, const double lambda_sigma, RNG &gen)
{
  int n = ybar.size();
  
  std::vector<std::vector<int> > tmp_clusters(1);
  for(int i = 0; i < n; i++) tmp_clusters[0].push_back(i);
  
  LPPartition gamma_init = new Partition(n, tmp_clusters, ybar, T, A_block, rho, a1, a2, eta);
  
  int max_splits = floor(log( (double) n));
  int n_k = n;
  arma::mat U = arma::zeros<arma::mat>(n_k,1);
  arma::mat means = arma::zeros<arma::mat>(1, max_splits);
  double min_score = 0.0;
  bool status = true;
  
  arma::vec cluster_dist = arma::zeros<arma::vec>(max_splits); // holds distance from single observaiton to each cluster mean
  arma::uvec cluster_dist_indices(max_splits); // used to sort cluster_dist
  
  std::vector<std::vector<int> > init_new_clusters; // stores output of k-means. NOTE: these clusters may not be connected.
  arma::mat A_tmp = arma::zeros<arma::mat>(n_k,n_k); // store submatrices of A_block. used to determine connected components
  std::vector<std::vector<int> > connected_components; // connected components of the individual sub clusters found by k-means
  
  std::vector<std::vector<int> > tmp_new_clusters; // used for figuring out the connected components of the new clusters
  std::vector<std::vector<int> > new_clusters(1); // holds the final new clusters
  
  arma::vec init_log_post(max_splits);
  arma::uvec init_log_post_indices(max_splits); // for sorting the log posterior values
  std::vector<double> log_post_prob(max_splits, 0.0);
  double tmp_log_post = 0.0;
  double max_log_post = 0.0;
  
  LPPartition gamma_new = new Partition();
  std::vector<LPPartition> unik_particles(max_splits); // holds the results from running K-mean
  
  for(int k = 0; k <= max_splits-1; k++) unik_particles[k] = new Partition();
  
  for(int num_splits = 1; num_splits <= max_splits; num_splits++){
    //Rcpp::Rcout << "num_splits = " << num_splits << std::endl;
    if(num_splits > 1){
      U.set_size(n_k,1);
      means.set_size(1, num_splits);
      min_score = 0.0;
      status = true;
      cluster_dist.set_size(num_splits);
      cluster_dist_indices.set_size(num_splits);
      for(int ii = 0; ii < n_k; ii++){
        //Rcpp::Rcout << " " << ii;
        U(ii,0) = ybar(ii);
      }
      kmeans_repeat(U, means, status, min_score, num_splits, 5000);
      if(status == true){
        init_new_clusters.clear();
        init_new_clusters.resize(num_splits);
        for(int ii = 0; ii < n_k; ii++){
          cluster_dist.zeros();
          for(int new_k = 0; new_k < num_splits; new_k++) cluster_dist(new_k) = arma::norm(U.row(ii).t() - means.col(new_k));
          cluster_dist_indices = arma::sort_index(cluster_dist, "ascend");
          init_new_clusters[cluster_dist_indices(0)].push_back(ii);
        }
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
        gamma_new = new Partition(n, new_clusters, ybar, T, A_block, rho, a1, a2, eta);
        //Rcpp::Rcout << "  updated gamma_new" << std::endl;
        unik_particles[num_splits-1]->Copy_Partition(gamma_new);
        tmp_log_post = total_log_post(unik_particles[num_splits-1], total_ss, T, nu_sigma, lambda_sigma);
        init_log_post(num_splits-1) = tmp_log_post;
        log_post_prob[num_splits - 1] = tmp_log_post;
        if(tmp_log_post > max_log_post) max_log_post = tmp_log_post;
      } // closes if/else checking that kmeans ran successfull
      
    } else{
      new_clusters.clear();
      new_clusters.resize(1);
      for(int i = 0; i < n; i++) new_clusters[0].push_back(i);
      delete gamma_new;
      gamma_new = new Partition(n, new_clusters, ybar, T, A_block, rho, a1, a2, eta);
      unik_particles[num_splits-1]->Copy_Partition(gamma_new);
      tmp_log_post = total_log_post(unik_particles[num_splits-1], total_ss, T, nu_sigma, lambda_sigma);
      init_log_post(num_splits-1) = tmp_log_post;
      log_post_prob[num_splits - 1] = tmp_log_post;
      max_log_post = tmp_log_post; // this is the first one encountered
      
    }
  } // closes loop running k-means
    //Rcpp::Rcout << "Finished k-means" << std::endl;
    // now we order by the log_post values
  init_log_post_indices = arma::sort_index(init_log_post, "descend");
  delete init_gamma;
  init_gamma = new Partition(unik_particles[init_log_post_indices(0)]);

}
