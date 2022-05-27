/*
 * particle.h
 *
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <dlib/matrix.h>
#include <armadillo>
using namespace std;
using namespace arma;

typedef class Particle* LPParticle;
class Particle
{
	//data members
public:
	int nObs;
	LPPartition partition_A;
	LPPartition partition_B;
	double* alpha;
	double* beta;
  arma::mat* I11;
  arma::mat* I12;
  arma::mat* I21;
  arma::mat* I22;
	// double sigma;
	// double tau_A;
	// double tau_B;

public:
  Particle(); // constructor
  Particle(LPParticle initial_particle); // constructor. create a new particle from an old one
  ~Particle(); // destructor

public:
  void get_alpha_beta();
  void get_alpha_beta_mle(double* alpha_mle, double* beta_mle);
  void get_alpha_beta_mle_normal(double* alpha_mle, double* beta_mle);
  void get_alpha_beta_mle_poisson(double* alpha_mle, double* beta_mle);
  void get_alpha_beta_poisson();
  void get_alpha_beta_iter();
  void get_parameter(bool A_or_B, int cluster_id);
  
  void Copy_Particle(LPParticle initial_particle);
  // void Initialize_Particle(int n);
  void Initialize_Particle(int n, Rcpp::List gamma_init_A = R_NilValue, Rcpp::List gamma_init_B = R_NilValue);
  void Initialize_Particle_nclusters(int n);
  
  void Partition_Split(bool A_or_B, int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2);
  void Partition_Split(bool A_or_B, int split_k, std::vector<int> new_cluster1, std::vector<int> new_cluster2, int size1, int size2);
  void Partition_KSplit(bool A_or_B, int split_k, std::vector<std::vector<int> > indices, std::vector<int> ns);
  void Partition_Merge(bool A_or_B, int k_1, int k_2);
  void Partition_Split_and_Merge(bool A_or_B, int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2, int k_star_1, int k_star_2);
  void Partition_Modify(bool A_or_B, int cl_ind);

  void Update_Partition(bool A_or_B, int current_l, std::vector<LPParticle> Particle_Set, std::vector<double> w, double lambda, double xi, bool restricted = false);
  void Update_Particle(int current_l, std::vector<LPParticle> Particle_Set, std::vector<double> w, double lambda, double xi, bool restricted = false); //done

  void Island_moves(bool A_or_B, double *pointer_to_maxobjective, int current_l, LPParticle& max_candidate, std::vector<LPParticle> Particle_Set, std::vector<double> w, double lambda, double xi, string *pt_maxstr, double percentage = 0.05);
  void Border_moves(bool A_or_B, double *pointer_to_maxobjective, int current_l, LPParticle& max_candidate, std::vector<LPParticle> Particle_Set, std::vector<double> w, double lambda, double xi, string *pt_maxstr);
  void Merge_moves(bool A_or_B, double *pointer_to_maxobjective, int current_l, LPParticle& max_candidate, std::vector<LPParticle> Particle_Set, std::vector<double> w, double lambda, double xi, string *pt_maxstr);
  
  void Spectral_SplitMerge_moves(bool A_or_B, double *pointer_to_maxobjective, int current_l, LPParticle& max_candidate, std::vector<LPParticle> Particle_Set, std::vector<double> w, double lambda, double xi, string *pt_maxstr, bool use_mle = false);
  void KM_SplitMerge_moves(bool A_or_B, double *pointer_to_maxobjective, int current_l, LPParticle& max_candidate, std::vector<LPParticle> Particle_Set, std::vector<double> w, double lambda, double xi, string *pt_maxstr, bool use_mle = false);
  void Tail_Split_moves(bool A_or_B, double *pointer_to_maxobjective, int current_l, LPParticle& max_candidate, std::vector<LPParticle> Particle_Set, std::vector<double> w, double lambda, double xi, string *pt_maxstr);
  
  std::vector<int> get_jstar(bool A_or_B, int split_k, int num_splits, std::vector< std::vector<int> > indices, std::vector<int> ns);
  std::vector<int> get_tail_jstar(bool A_or_B, int split_k, std::vector<std::vector<int> > tail_conncomp);
  
  double Total_Log_Post(); //done
  double LikelihoodYAB();
  double LikelihoodY();
  double LikelihoodY_MCMC_poisson();
  double LikelihoodY_MCMC_normal();
  arma::mat get_Omega(bool A_or_B);
  double fun_like_pois(dlib::matrix<double,0,1> alpha_beta, bool mle_bool = false);
  // double fun_like_negbin(dlib::matrix<double,0,1> alpha_beta);
  double fun_like_pois_singleunit(dlib::matrix<double,0,1> alpha_beta, int i);

  dlib::matrix<double,0,1> fun_deriv_like_pois(dlib::matrix<double,0,1> alpha_beta, bool mle_bool = false);
  dlib::matrix<double,0,1> fun_deriv_like_pois_singleunit(dlib::matrix<double,0,1> alpha_beta, int i);
  dlib::matrix<double> fun_hess_like_pois(dlib::matrix<double,0,1> alpha_beta, bool mle_bool = false);
  dlib::matrix<double> fun_hess_like_pois_singleunit(dlib::matrix<double,0,1> alpha_beta);

  arma::mat fun_hess_like_pois(arma::vec alpha, arma::vec beta);
  double optim(dlib::matrix<double,0,1> &alpha_beta);
  void optim_mle(dlib::matrix<double,0,1> &alpha_beta);
  double LikelihoodY_Laplace_Pois();
  double total_log_like();
  
  void Print_Particle(); //done
  void Print_Particle_Short(); //done
  void Print_Particle_ToFile(string file); //done
  
  void Find_Splits_nonconnected(bool A_or_B, int cluster_id, int n_cl, int* component_i, int* component_not_i, int **index1_ptr, int **index2_ptr, int &n1, int &n2);
  void Find_Splits_nonspatial(bool A_or_B, int cluster_id, int **index1_ptr, int **index2_ptr, int &n1, int &n2);
  bool Find_K_SpectralSplits(bool A_or_B, int split_k, int num_splits, std::vector<std::vector<int> >& indices, std::vector<int>& ns, bool use_mle = false); 
  bool Find_KM_Splits(bool A_or_B, int split_k, int num_splits, std::vector<std::vector<int> >& indices, std::vector<int>& ns, bool use_mle = false);
  void Find_TailSplit(bool A_or_B, int left0right2, int i, int split_k, std::vector<std::vector<int> > left_center_right, std::vector<std::vector<int> >& tail_conncomp, std::vector<std::vector<int> >& indices, std::vector<int>& ns);
  
  void Initial_K_Splits(int num_splits_A, int num_splits_B, bool use_mle = false);
  void Initial_KM_Splits(int num_splits_A, int num_splits_B, bool use_mle = false);
  double param_bar(bool A_or_B, int k);
  void get_leftcenterright(bool A_or_B, std::vector<std::vector<int> >& left_center_right, int split_k);

  void update_sigma2();
  // void Partition_K_Splits(bool A_or_B, int cluster_id, int num_splits);
  // void Find_SpectralSplits(bool A_or_B, int cluster_id, int **index1_ptr, int **index2_ptr, int &n1, int &n2);
  // int get_jstar(bool A_or_B, int k);
};



#endif /* PARTICLE_H_ */

