#include <stdio.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "partition.h"
#include "particle.h"
#include "partition_functions.h"
#include "particle_functions.h"
#include "various_functions.h"

#include <vector>

using namespace Rcpp;
using namespace std;
using namespace arma;

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
double sigma_py = 0.0;
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


// [[Rcpp::export]]
Rcpp::List particle_summary(arma::mat Y_input,
                            arma::mat X_input,
                            arma::mat A_block_input,
                            Rcpp::List gamma_init_A,
                            Rcpp::List gamma_init_B,
                            double a1_input = 0.01, double a2_input = 10,
                            double b1_input = 0.001, double b2_input = 1,
                            double alpha_sigma_input = 0.5, double beta_sigma_input = 0.5,
                            int priorA_input = 0, int priorB_input = 0,
                            double eta_input = 1.0,
                            double rho_input = 0.99){
  Y = Y_input;
  X = X_input;
  A_block = A_block_input;

  a1 = a1_input;
  a2 = a2_input;
  b1 = b1_input;
  b2 = b2_input;

  alpha_sigma = alpha_sigma_input;
  beta_sigma = beta_sigma_input;
  rho = rho_input;
  eta = eta_input;
  
  priorA = priorA_input;
  priorB = priorB_input;
  
  opt_method = 0;
  opt_Y = 2;

  int n = Y.n_rows;
  // int T = Y.n_cols;

  LPParticle Gamma_0 = new Particle;
  Gamma_0->Initialize_Particle(n, gamma_init_A, gamma_init_B);

  Rcpp::List particle_set_outA(1);
  Rcpp::List particle_set_outB(1);
  particle_set_outA[0] = gamma_init_A;
  particle_set_outB[0] = gamma_init_B;

  std::vector<double> alpha(n);
  std::vector<double> beta(n);
  for(int i = 0; i < n; i++){
    alpha[i] = Gamma_0->alpha[i];
    beta[i] = Gamma_0->beta[i];
  }

  Rcpp::List results;
  results["particle_set_A"] = particle_set_outA;
  results["particle_set_B"] = particle_set_outB;
  results["log_like"] = Gamma_0->total_log_like();
  results["log_prior_A"] = total_log_prior(Gamma_0->partition_A);
  results["log_prior_B"] = total_log_prior(Gamma_0->partition_B);
  results["log_post"] = Gamma_0->Total_Log_Post();
  results["alpha"] = alpha;
  results["beta"] = beta;
  return results;
}
