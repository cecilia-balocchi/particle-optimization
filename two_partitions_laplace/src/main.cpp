#include <stdio.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(dlib)]]
#include <dlib/optimization.h>
#include <dlib/global_optimization.h>

// // [[Rcpp::depends(RcppProgress)]]
// #include <progress.hpp>

#include <vector>

#include "partition.h"
#include "particle.h"
#include "partition_functions.h"
#include "various_functions.h"
#include "update_particle.h"
#include "particle_functions.h"


arma::mat Y;
arma::mat X;
arma::mat A_block;
arma::mat Lambda;
arma::mat sX;
double rho;
double a1;
double a2;
double b1;
double b2;
double a_cohes = 200;
double b_cohes = 10;
double nu_sigma = 1;
double lambda_sigma = 1;
double alpha_sigma = nu_sigma/2; 
double beta_sigma = nu_sigma*lambda_sigma/2; 
double eta;
double sigma_py;
double lambda;
double xi = 0.0;
int max_iter;
bool resampling;
int* sigma;
int priorA;
int priorB;
int opt_method;
int opt_Y;
double A_or_B_first;
int coordasc_iter = 2;
unsigned seed;
double split_frac;
Rcpp::List empty(0);
double* alpha_mle;
double* beta_mle;
bool sampleA;
bool sampleB;
int kmeans_rep = 500;
double sigma2_laplace = 0.5;
bool update_sigma2;
using namespace Rcpp;
using namespace std;
using namespace arma;


// [[Rcpp::export]]
Rcpp::List ensm_cluster_mean(arma::mat Y_input, arma::mat X_input, const arma::mat A_block_input, const int L, 
  const double lambda_input = 1.0, const double eta_input = 1.0, const double sigma_py_input = 0.0, const int max_iter_input = 10, 
  double a1_input = 0.01, double a2_input = 10, double b1_input = 0.001, double b2_input = 1,  
  double alpha_sigma_input = 0.5, double beta_sigma_input = 0.5, 
  int opt_method_input = 0, int opt_Y_input = 1, int priorA_input = 0, int priorB_input = 0, bool Kmeans_initialize = false, int k_A = 0, int k_B = 0, 
  Rcpp::Nullable<Rcpp::NumericVector> ks_A_ = R_NilValue, Rcpp::Nullable<Rcpp::NumericVector> ks_B_ = R_NilValue, //IntegerVector ks_A = 0, IntegerVector ks_B = 0, 
  const double rho_input = 0.8, const double eps = 1e-3, double split_frac_input = 0.1, double A_or_B_first_input = 1.0,
  bool resampling_input = false, Rcpp::List gamma_init_A = R_NilValue, Rcpp::List gamma_init_B = R_NilValue, 
  Rcpp::Nullable<Rcpp::NumericVector> w_init_ = R_NilValue, bool sampleA_input = true, bool sampleB_input = true, 
  bool last_islands = false, bool update_sigma2_ = true, double starting_sigma2_ = 0.5)
{
  priorA = priorA_input;
  priorB = priorB_input;
  opt_method = opt_method_input;
  opt_Y = opt_Y_input;
  A_or_B_first = A_or_B_first_input;
  // We should not be using opt_Y = 0 anymore! (It computes LikelihoodYAB)
  if(opt_Y == 0){
    cout << "WARNING: opt_Y = 0 so LikelihoodYAB is used. Choose instead opt_Y = 1 for LikelihoodY " <<
      "and opt_Y = 2 for using independence of alpha and beta." << endl;
  }

  Y = Y_input;
  X = X_input;
  // We should check that the data is actually normalized.
  bool X_orthogonal = is_orthogonal(X);
  if(X_orthogonal){
    if(opt_method != 0){
      cout << "WARNING: data is orthogonal and opt_method is different from 0; you are not using the fastest method!" << endl;
    }
    if((opt_Y != 2) & (opt_Y != 4)){
      cout << "WARNING: data is orthogonal and opt_Y is different from 2; you are not using the fastest method!" << endl;
    }
  } else {
    if(opt_method == 0){
      cout << "WARNING: data is NOT orthogonal and opt_method is equal to 0; the method you are using is not correct!" << endl;
    }
    if(opt_Y == 2){
      cout << "WARNING: data is NOT orthogonal and opt_Y is equal to 2; the method you are using is not correct!" << endl;
    }
  }
  A_block = A_block_input;
  if(L == 1){
    update_sigma2 = true;
    if(update_sigma2_ == false){
      update_sigma2 = false;
    }
    // basically:
    // update_sigma2 = update_sigma2_;
  } else {
    update_sigma2 = false;
  }
  sigma2_laplace = starting_sigma2_;

  a1 = a1_input;
  a2 = a2_input;
  b1 = b1_input;
  b2 = b2_input;
  std::vector<double> hyper_par(4);
  hyper_par[0] = a1;
  hyper_par[1] = a2;
  hyper_par[2] = b1;
  hyper_par[3] = b2;

  alpha_sigma = alpha_sigma_input;
  beta_sigma = beta_sigma_input;

  sampleA = sampleA_input;
  sampleB = sampleB_input;
  rho = rho_input;
  lambda = lambda_input;
  eta = eta_input;
  sigma_py = sigma_py_input;
  max_iter = max_iter_input;
  resampling = resampling_input;
  split_frac = split_frac_input;
  int n = Y.n_rows;

  vector<LPParticle> Particle_Set(L);
  vector<double> w(L);
  double max_log_post = 0.0;

  // unsigned seed = time(0);
  // std::default_random_engine eng(seed);

  // initialize with Kmeans:
  // using the MLE find pairs of partitions with (k,h) clusters for k,h = 1 .. Kmeans (=log(n))
  // (the Kmeans is made "spatial" by finding connected subclusters)
  // create a cdf with the log posterior values and sample L values
  // the particle elements are initialized with possibly different values

  LPParticle Gamma_0 = new Particle;

  // **** Here I am starting tests for optim and optim_mle (used in get_alpha_beta_mle) ****

  // Gamma_0->Initialize_Particle(n);
  // alpha_mle = new double[n];
  // beta_mle = new double[n];
  // cout << "before get_alpha_beta_mle" << endl;
  // Gamma_0->get_alpha_beta_mle(alpha_mle, beta_mle);
  // cout << "after get_alpha_beta_mle" << endl;
  // for(int i = 0; i < n; i++){
  //   cout << alpha_mle[i] << " ";
  // }
  // cout << endl << endl;
  // for(int i = 0; i < n; i++){
  //   cout << beta_mle[i] << " ";
  // }
  // cout << endl << endl;


  // dlib::matrix<double,0,1> alpha_beta(2*(Y.n_rows));
  // for(int i =0; i < Y.n_rows; i++){
  //   alpha_beta(i) = 0;
  //   alpha_beta(i+Y.n_rows) = 0;
  // }
  // cout << "before optim" << endl;
  // double res = Gamma_0->optim(alpha_beta); 
  // cout << "after optim" << endl;
  // for(int i = 0; i < n; i++){
  //   cout << alpha_beta(i) << " ";
  // }
  // cout << endl << endl;
  // for(int i = 0; i < n; i++){
  //   cout << alpha_beta(n+i) << " ";
  // }
  // cout << endl << endl;
  
  // // If we want to check the objective function in a given datapoint 
  // // The one below is the glm estimate computed in R for i = 2, .. , 11
  // dlib::matrix<double,0,1> alpha_beta(2*n);
  // alpha_beta = 1.426294,
  //             1.226410,
  //             1.295612,
  //             1.299088,
  //             3.576279,
  //             3.464746,
  //             3.458611,
  //             3.413097,
  //             3.318283,
  //             3.458568,
  //             -0.04237306,
  //             0.07017038,
  //             0.08952695,
  //             -0.02063098,
  //             -0.02433473,
  //             0.08854025,
  //             0.09932463,
  //             0.06429234,
  //             0.15467500,
  //             -0.97444112;
  // cout << "glm: " <<  Gamma_0->fun_like_pois(alpha_beta) << endl; 
  // cout << res << endl;

  // **** End of tests for optim and optim_mle (used in get_alpha_beta_mle) ****

  if(Kmeans_initialize){
    Gamma_0->Initialize_Particle(n);
    alpha_mle = new double[n];
    beta_mle = new double[n];
    Gamma_0->get_alpha_beta_mle(alpha_mle, beta_mle);
    
    int Kmeans = (int)floor(eta*log(n));
    cout << "Beginning Kmeans_initialize ";
    cout << "with k from 1 to Kmeans = " << Kmeans;
    int Km_comb = Kmeans*Kmeans;
    if(!((gamma_init_A.size() == 0) & (gamma_init_B.size() == 0))){
      // if we provide a starting point, let's add it to the Initial_Particle_Set 
      Km_comb++;
      cout << " and a starting point.";
    }
    cout << std::endl;
    // Initial_Particle_Set contains all the combinations of kmeans partitions:
    vector<LPParticle> Initial_Particle_Set(Km_comb);
    std::vector<double> log_post_prob(Km_comb, 0.0);
    double tmp_log_post = 0.0;
    int l = 0;
    for(int i = 0; i < Kmeans; i++){
      for(int j = 0; j < Kmeans; j++){
        Initial_Particle_Set[l] = new Particle(Gamma_0);
        Initial_Particle_Set[l]->Initial_KM_Splits(i+1, j+1);
        tmp_log_post = Initial_Particle_Set[l]->Total_Log_Post();
        log_post_prob[l] = tmp_log_post;
        if(l == 0){
          max_log_post = tmp_log_post;
        } else {
          if(tmp_log_post > max_log_post){
            max_log_post = tmp_log_post;
          }
        }
        l++;
      }
    }
    if(!((gamma_init_A.size() == 0) & (gamma_init_B.size() == 0))){
      Initial_Particle_Set[l] = new Particle(Gamma_0);
      Initial_Particle_Set[l]->Initialize_Particle(n, gamma_init_A, gamma_init_B);
      tmp_log_post = Initial_Particle_Set[l]->Total_Log_Post();
      log_post_prob[l] = tmp_log_post;
      if(tmp_log_post > max_log_post){
        max_log_post = tmp_log_post;
      }
    }
    // create now the cdf
    std::vector<double> post_prob(Km_comb,0.0);
    std::vector<double> cdf(Km_comb,0.0);
    double tmp_sum = 0.0;
    for(size_t ix = 0; ix < Km_comb; ix++){
      log_post_prob[ix] -= max_log_post;
      post_prob[ix] = exp(log_post_prob[ix]);
      tmp_sum += post_prob[ix];
      cdf[ix] = tmp_sum;
    }
    for(size_t ix = 0; ix < Km_comb; ix++){
      post_prob[ix] /= tmp_sum; // normalizes the CDF
      cdf[ix] /= tmp_sum ; // normalizes the CDF
    }

    // print
    std::vector<int> possible_ix;
    double unif = 0.0;
    int gamma_ix = 0;
    // cout << "max_log_post = " << max_log_post << std::endl;
    
    // cout << "Log-post is : ";
    // for(size_t ix = 0; ix < Km_comb; ix++) cout << " " << log_post_prob[ix];
    // cout << std::endl << std::endl;
    
    cout << "posterior is : ";
    for(size_t ix = 0; ix < Km_comb; ix++) cout << " " << post_prob[ix];
    cout << std::endl << std::endl;
    
    // cout << "cdf if : ";
    // for(size_t ix = 0; ix < Km_comb; ix++) cout << " " << cdf[ix];
    // cout << std::endl << std::endl;
    
    // sample the L starting points of the particle set
    for(int l  = 0; l < L; l++){
      unif = R::runif(0.0, 1.0);
      possible_ix.clear();
      for(int ix = 0; ix < Km_comb; ix++){
        if(cdf[ix] <= unif) possible_ix.push_back(ix);
      }
      
      if(possible_ix.size() == 0) gamma_ix = 0; // if we somehow generated unif = 0...
      else if(possible_ix.size() == Km_comb) gamma_ix = Km_comb-1; // if we somehow generated unif == 1
      else{
        gamma_ix = possible_ix[possible_ix.size() - 1] + 1; // cdf[gamma_ix] < unif < cdf[gamma_ix + 1]
      }
      Particle_Set[l] = new Particle(Initial_Particle_Set[gamma_ix]);
      
      cout << "l = " << l << " unif = " << unif << " gamma_ix = " << gamma_ix << std::endl;
    }
    for(int i = 0; i < Km_comb; i++){
      delete Initial_Particle_Set[i];
    }
    cout << "Initialized with Kmeans." << endl << endl;
  } else {
    // if no starting point are provided, just initialize with one cluster
    if((gamma_init_A.size() == 0) & (gamma_init_B.size() == 0)){
      Gamma_0->Initialize_Particle(n);
      alpha_mle = new double[n];
      beta_mle = new double[n];
      Gamma_0->get_alpha_beta_mle(alpha_mle, beta_mle);
    } else {
      Gamma_0->Initialize_Particle(n, gamma_init_A, gamma_init_B);
      cout << "Initialized with custom partition." << endl;
      // alpha_mle and beta_mle need a partition where everyone is in one cluster
      LPParticle Gamma_1 = new Particle;
      Gamma_1->Initialize_Particle(n);
      alpha_mle = new double[n];
      beta_mle = new double[n];
      Gamma_1->get_alpha_beta_mle(alpha_mle, beta_mle);
      delete Gamma_1;
    }
    // if we need to do some additional splits
    if((k_A > 1) || (k_B > 1)){
      Gamma_0->Initial_KM_Splits(k_A, k_B);
      cout << "Starting point: KM_Splits" << endl;
      Gamma_0->Print_Particle();
    }
    for(int l = 0; l < L; l++){
      Particle_Set[l] = new Particle(Gamma_0);
    }
    bool multiple_splits = false;
    std::vector<int> ksA(L,0);
    std::vector<int> ksB(L,0);
    if (ks_A_.isNotNull()) {
      NumericVector ks_A(ks_A_.get());
      multiple_splits = true;
      for(int l = 0; l < L; l++){
        ksA[l] = ks_A[l];
      }
    }
    if (ks_B_.isNotNull()) {
      NumericVector ks_B(ks_B_.get());
      multiple_splits = true;
      for(int l = 0; l < L; l++){
        ksB[l] = ks_B[l];
      }
    }
    if(multiple_splits){
      for(int l = 0; l < L; l++){
        Particle_Set[l]->Initial_KM_Splits(ksA[l], ksB[l]);
      }
      cout << "Starting point: Different KM_Splits" << endl;
    }
  }
  delete Gamma_0;
  if (w_init_.isNotNull()) {
    NumericVector w_init(w_init_.get());
    if(w_init.size() == L){
      for(int l = 0; l < L; l++){
        w[l] = w_init[l];
      }
      cout << "Initialized with custom weights." << endl;
    } else {
      for(int l = 0; l < L; l++){
        w[l] = (double) 1/L;
      } 
    }
  } else {
    for(int l = 0; l < L; l++){
      w[l] = (double) 1/L;
    }  
  }
  std::vector< std::vector<double> > w_history;
  w_history.push_back(w);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // this is the code we had in ensm_particle:
  int iter = 0;
  int flag = 0;
  int conv_counter = 0; // counter to see how many particles remain unchanged

  double old_objective = 0.0;
  double objective = 0.0;

  objective = lambda * Entropy(0, Particle_Set[0], Particle_Set, w) + xi * VI_Avg(0, Particle_Set[0], Particle_Set);
  for(int l = 0; l < L; l++){
    objective += w[l] * Particle_Set[l]->Total_Log_Post();
  }
  LPParticle old_particle = new Particle(Particle_Set[0]);
  
  while((iter < max_iter) & (flag == 0)){

    // if(Progress::check_abort()){
    //   return -1;
    // }
    // try to update the particle
    // compute the old objective value
    cout << "[particle_spatial]: Starting iter = " << iter << std::endl;
    old_objective = objective;
    
    
    if(resampling){
      Resampling(Particle_Set, w);
    }

    // sweep over the particle set
    conv_counter = 0;
    for(int l = 0; l < L; l++){
      cout << "-- [particle_spatial]: updating particle " << l << " --" << endl;
      // free up the old_particle
      delete old_particle; //not used anymore
      try{
        old_particle = new Particle(Particle_Set[l]);
      }
      catch(const std::bad_alloc& e){
        cout << "EXCEPTION IN PARTICLE SPATIAL"  << e.what() << endl;
      }
      Particle_Set[l]->Update_Particle(l, Particle_Set, w, lambda, xi);
      conv_counter += Particle_Equal(Particle_Set[l], old_particle); // check if the current particle is the best
    }

    // now let us update w
    update_w(Particle_Set, w, lambda, xi);
    w_history.push_back(w);

    objective = lambda * Entropy(0, Particle_Set[0], Particle_Set, w) + xi * VI_Avg(0, Particle_Set[0], Particle_Set); // compute the entropy
    for(int l = 0; l < L; l++){
      objective += w[l] * Particle_Set[l]->Total_Log_Post();
    }

    // cout << endl;
    // cout << "[ensm_partition]: obj = " << setprecision(8) << objective << "    old_obj = " << setprecision(8) << old_objective << endl;
    // cout << "[ensm_partition]: percent increase in objective = " << setprecision(6) << 100*fabs((objective - old_objective)/old_objective) << endl;
    // cout << "[ensm_partition]: number of stationary particles = " << conv_counter << endl;
    // cout << "[ensm_partition]: importance weights :  " << endl;
    cout << "[ensm_partition]: importance weights :  ";
    for(int l = 0; l < L; l++){
      cout << setprecision(6) << w[l] << "   " ;
    }
    cout << endl;
    // check for convergence
    if(objective < old_objective){
      cout << "WARNING THE OBJECTIVE DECREASED" << endl;
    }
    cout << endl;

    // Note: at iter (t+1), objective was computed with sigma2_laplace computed at time (t) and
    // old_objective was computed with sigma2_laplace computed at time (t-1), 
    // so it's normal for the likelihood to decrease

    // flag = 0;
    if(conv_counter == L){ //((conv_counter == L) || ( fabs((objective - old_objective)/old_objective) < 1e-6))
      flag = 1;
    }

    if(update_sigma2){
      // can only happen if L = 1
      cout << "sigma2_laplace before update: " << sigma2_laplace << " -- ";
      Particle_Set[0]->update_sigma2();
      cout << "sigma2_laplace after update: " << sigma2_laplace << endl;
    }

    iter++;
  }
  cout << "[ensm_partition]: Total number of iterations = " << iter << endl;
  delete old_particle;
  // this is to check all islands before ending the algorithm
  if(last_islands){
    cout << "Beginning last_islands - ";
    bool A_or_B;
    LPParticle max_candidate;
    double max_objective;
    string max_str;
    old_objective = objective;
    for(int l = 0; l < L; l++){
      max_str = "Not changed";
      A_or_B = true;
      max_candidate = new Particle(Particle_Set[l]);
      max_objective = w[l]*max_candidate->Total_Log_Post() + lambda * Entropy(l, max_candidate, Particle_Set, w) + xi * VI_Avg(l, max_candidate, Particle_Set);
      Particle_Set[l]->Island_moves(A_or_B, &max_objective, l, max_candidate, Particle_Set, w, lambda, xi, &max_str, 1.0);
      Particle_Set[l]->Copy_Particle(max_candidate);
      // delete max_candidate;
      A_or_B = false;
      // max_candidate = new Particle(Particle_Set[l]);
      // max_objective = w[l]*max_candidate->Total_Log_Post() + lambda * Entropy(l, max_candidate, Particle_Set, w) + xi * VI_Avg(l, max_candidate, Particle_Set);
      Particle_Set[l]->Island_moves(A_or_B, &max_objective, l, max_candidate, Particle_Set, w, lambda, xi, &max_str, 1.0);
      Particle_Set[l]->Copy_Particle(max_candidate);
      delete max_candidate;
    }
    objective = lambda * Entropy(0, Particle_Set[0], Particle_Set, w) + xi * VI_Avg(0, Particle_Set[0], Particle_Set); // compute the entropy
    for(int l = 0; l < L; l++){
      objective += w[l] * Particle_Set[l]->Total_Log_Post();
    }
    if(objective > old_objective){
      cout << "Objective increased with last_islands" << std::endl;
    }
    cout << " - Finished last_islands." << std::endl;

    if(update_sigma2){
      Particle_Set[0]->update_sigma2();
    }
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  

  
  // prepare the output containers
  Rcpp::List tmp_list;
  arma::vec tmp_vec = arma::zeros<vec>(1);
  Rcpp::List particle_set_outA(L);
  Rcpp::List particle_set_outB(L);
  Rcpp::NumericVector output_vec;
  
  for(int l = 0; l < L; l++){
    tmp_list = List(Particle_Set[l]->partition_A->K);
    for(int k = 0; k < Particle_Set[l]->partition_A->K; k++){
      tmp_vec = arma::zeros<vec>(Particle_Set[l]->partition_A->cluster_config[k]);
      for(int i = 0; i < Particle_Set[l]->partition_A->cluster_config[k]; i++){
        tmp_vec[i] = Particle_Set[l]->partition_A->clusters[k][i] + 1; // remember that in R we are 1-indexed
      }
      output_vec = Rcpp::wrap(tmp_vec);
      output_vec.attr("dim") = R_NilValue;
      tmp_list[k] = output_vec;
    }
    particle_set_outA[l] = tmp_list;
    tmp_list = List(Particle_Set[l]->partition_B->K);
    for(int k = 0; k < Particle_Set[l]->partition_B->K; k++){
      tmp_vec = arma::zeros<vec>(Particle_Set[l]->partition_B->cluster_config[k]);
      for(int i = 0; i < Particle_Set[l]->partition_B->cluster_config[k]; i++){
        tmp_vec[i] = Particle_Set[l]->partition_B->clusters[k][i] + 1; // remember that in R we are 1-indexed
      }
      output_vec = Rcpp::wrap(tmp_vec);
      output_vec.attr("dim") = R_NilValue;
      tmp_list[k] = output_vec;
    }
    particle_set_outB[l] = tmp_list;
  }
  
  arma::mat alpha_particle = arma::zeros<mat>(n, L); // estimates for each unique particle
  arma::mat beta_particle = arma::zeros<mat>(n, L); // estimates for each unique particle
  for(int l = 0; l < L; l++){
    for(int i = 0; i < n; i++){
      alpha_particle(i,l) = Particle_Set[l]->alpha[i];
      beta_particle(i,l) = Particle_Set[l]->beta[i];
    }
  }
  vector<double> log_posterior(L);
  vector<double> log_likelihood(L);
  for(int l = 0; l < L; l++){
    log_posterior[l] = Particle_Set[l]->Total_Log_Post();
    log_likelihood[l] = Particle_Set[l]->total_log_like();
  }

  // find unique values of the particles and of w (p_star)
  std::vector<int> unik_index;
  int num_unik_particles = 0;
  std::vector<int> particle_assignment(L);
  double tmp_log_post, tmp_p_star, tmp_norm;
  max_log_post = 0.0;
  for(int l = 0; l < L; l++){
    if(l == 0){
      unik_index.push_back(l);
      num_unik_particles++;
      particle_assignment[l] = 0;
      
      max_log_post = log_posterior[l];
    } else {
      bool is_duplicate = false;
      for(int ul = 0; ul < num_unik_particles; ul++){
        // if it's a duplicate we don't have to check anymore
        if(!is_duplicate){ 
          if(Particle_Equal(Particle_Set[l], Particle_Set[ unik_index[ul] ]) == 1){ 
            is_duplicate = true;
            particle_assignment[l] = ul;
          }
        }
      } // end for
      if(is_duplicate == false){
        unik_index.push_back(l);
        num_unik_particles++;
        particle_assignment[l] = num_unik_particles - 1;
      }

      if(log_posterior[l] > max_log_post){
        max_log_post = log_posterior[l];
      }
    }
  }
  std::vector<double> p_star;
  if(num_unik_particles == L){
    p_star = w;
  } else {
    tmp_norm = 0.0;
    for(int ul = 0; ul < num_unik_particles; ul++){
      tmp_log_post = log_posterior[ particle_assignment[ul] ] - max_log_post;
      tmp_p_star = exp(1/lambda * tmp_log_post); 
      tmp_norm += tmp_p_star;
      p_star.push_back(tmp_p_star);
    }
    for(int ul = 0; ul < num_unik_particles; ul++){
      p_star[ul] /= tmp_norm;
    }
  }
  



  for(int l = 0; l < L; l++){
    delete Particle_Set[l];
  }
  
  Rcpp::List results;
  results["particle_set_A"] = particle_set_outA;
  results["particle_set_B"] = particle_set_outB;
  results["w"] = w;
  results["w_history"] = w_history;
  results["p_star"] = p_star; // new
  results["unik_index"] = unik_index; // new
  results["particle_assignment"] = particle_assignment; // new
  results["alpha_particle"] = alpha_particle;
  results["beta_particle"] = beta_particle;
  results["log_posterior"] = log_posterior;
  results["log_likelihood"] = log_likelihood;
  results["hyper_parameters"] = hyper_par;
  results["sigma2_laplace"] = sigma2_laplace;
  return results;
  
  // return 0;
}