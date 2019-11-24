
// partition_functions.h



#ifndef GUARD_PARTITION_FUNCTIONS_H
#define GUARD_PARTITION_FUNCTIONS_H
#include <RcppArmadillo.h>
#include <vector>


using namespace arma;
using namespace std;

class split_info{
public:
  int num_splits;
  std::vector<int> split_k; // holds the cluster id being split
  std::vector<std::vector<std::vector<int> > > new_clusters; // holds the new sub-clusters being created
  std::vector<std::vector<int> > nearest_neighbor; // hold the cluster id of the existing cluster closest to each new sub-cluster
  split_info() {num_splits = 0; split_k = std::vector<int>(1,0); new_clusters = std::vector<std::vector<std::vector<int> > > (1); nearest_neighbor = std::vector<std::vector<int> >(1, std::vector<int>(1,0));}
  ~split_info(){};
};

class merge_info{
public:
  int num_merges;
  std::vector<int> rec_k; // holds the cluster id of what is accepting the new elemets
  std::vector<int> donor_k; // holds cluster id of what is being merged
};

int Partition_Equal(Partition *partition1, Partition *partition2);

void get_unik_particles(std::vector<std::vector<int> > &particle_map, std::vector<double> &p_star, std::vector<int> &counts, std::vector<LPPartition> particle_set, std::vector<double> w);


double Entropy(unsigned current_l, Partition* candidate_particle, std::vector<LPPartition> particle_set, std::vector<double> w);

// need to account for the block group variances in the log likelihood


double total_log_post(LPPartition partition, const double total_ss, const int T, const double nu_sigma, const double lambda_sigma);
double total_log_like(LPPartition partition, const double total_ss, const int T, const double nu_sigma, const double lambda_sigma);
double total_log_prior(LPPartition partition);
double alpha_bar_func(std::vector<int> new_cluster, LPPartition gamma_l, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2);

void update_w(std::vector<LPPartition> particle_set, std::vector<double> &w, const int L, const double total_ss, const int T, const double nu_sigma, const double lambda_sigma, const double lambda);


void get_island(split_info &si, LPPartition gamma_l, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double island_frac);
void get_border(split_info &si, LPPartition gamma_l, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2);

void get_merge(merge_info &mi, LPPartition gamma_l, const arma::mat &A_block);
void get_spectral_split(split_info &si, LPPartition gamma_l, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const int reps);
void get_tail_split(split_info &si, LPPartition gamma_l, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double tail_frac);
void get_km_split(split_info &si, LPPartition gamma_l, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const int reps);


void get_local(split_info &si, LPPartition gamma_l, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2);



void best_split_orig(split_info &si, LPPartition candidate, const int current_l, const std::vector<LPPartition> particle_set, const std::vector<double> w, const arma::vec &ybar, const double total_ss, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double nu_sigma, const double lambda_sigma, const double eta, const double lambda);
void best_split(split_info &si, LPPartition candidate, const int current_l, const std::vector<LPPartition> particle_set, const std::vector<double> w, const arma::vec &ybar, const double total_ss, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double nu_sigma, const double lambda_sigma, const double eta, const double lambda, const bool merge_flag);

void best_split_map(split_info &si, LPPartition candidate, LPPartition init_particle, const arma::vec &ybar, const double total_ss, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double nu_sigma, const double lambda_sigma, const double eta, const bool merge_flag);

void best_merge(merge_info &mi, LPPartition candidate, const int current_l, const std::vector<LPPartition> particle_set, const std::vector<double> w, const arma::vec &ybar, const double total_ss, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double nu_sigma, const double lambda_sigma, const double eta, const double lambda);

void best_merge_map(merge_info &mi, LPPartition candidate, LPPartition init_particle, const arma::vec &ybar, const double total_ss, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double nu_sigma, const double lambda_sigma, const double eta);





bool sanity_check(LPPartition partition);

void get_subcluster_neighbor(std::vector<std::vector<int> > &init_new_clusters, std::vector<std::vector<int> > &new_clusters, std::vector<int> &kstar, const int split_k, LPPartition gamma_l, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2);

void format_particle_set(std::vector<LPPartition> particle_set, Rcpp::List &output_list);
#endif

