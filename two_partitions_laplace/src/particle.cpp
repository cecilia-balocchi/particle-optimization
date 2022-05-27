/*
 * particle.cpp
 *
 */
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <dlib/matrix.h>
#include <dlib/optimization.h>
#include <dlib/global_optimization.h>
using namespace Rcpp;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <limits.h>
#include <armadillo>
#include <time.h>
#include <random>
#include "partition.h"
#include "particle.h"
#include "various_functions.h"
#include "partition_functions.h"
#include "particle_functions.h"
using namespace std;
using namespace arma;

extern arma::mat Y;
extern arma::mat X;
extern arma::mat A_block;
// Set the hyper-parameters
extern double rho;
// extern double a;
// extern double b;
extern double a1;
extern double a2;
extern double b1;
extern double b2;
extern double alpha_sigma;
extern double beta_sigma;
extern double eta;
extern int prior;
extern int opt_method;
extern int opt_Y;
extern int coordasc_iter;
extern double A_or_B_first;
extern int* sigma;
extern double split_frac;
extern double* alpha_mle;
extern double* beta_mle;
extern bool sampleA;
extern bool sampleB;
extern int kmeans_rep;
extern double sigma2_laplace;


// extern double a_taualpha;
// extern double b_taualpha;
// extern double a_taubeta;
// extern double b_taubeta;

Particle::Particle(){
  nObs = 0;
  partition_A = NULL;
  partition_B = NULL;
  alpha = NULL;
  beta = NULL;
  I11 = NULL;
  I12 = NULL;
  I21 = NULL;
  I22 = NULL;
  // sigma = 1.0;
  // tau_A = 1.0;
  // tau_B = 1.0;
  return;
}

Particle::Particle(LPParticle initial_particle){
  nObs = initial_particle->nObs;
  partition_A = new Partition(initial_particle->partition_A, 1);
  partition_B = new Partition(initial_particle->partition_B, 0);
  
  alpha = new double[nObs];
  beta = new double[nObs];
  for(int i = 0; i < nObs; i++){
    alpha[i] = initial_particle->alpha[i];
    beta[i] = initial_particle->beta[i];
  }
  if(opt_method == 1){
    I11 = new arma::mat(nObs, nObs, fill::zeros);
    I12 = new arma::mat(nObs, nObs, fill::zeros);
    I21 = new arma::mat(nObs, nObs, fill::zeros);
    I22 = new arma::mat(nObs, nObs, fill::zeros);

    *I11 = 0 + *(initial_particle->I11);
    *I12 = 0 + *(initial_particle->I12);
    *I21 = 0 + *(initial_particle->I21);
    *I22 = 0 + *(initial_particle->I22);
  }
  return;
}

Particle::~Particle(){
  delete partition_A;
  delete partition_B;
  delete[] alpha;
  delete[] beta;
  partition_A = NULL;
  partition_B = NULL;
  alpha = NULL;
  beta = NULL;

  if(opt_method == 1){
    delete I11;
    delete I12;
    delete I21;
    delete I22;
    
    I11 = NULL;
    I12 = NULL;
    I21 = NULL;
    I22 = NULL;
  }

  return;
}

// copy particle is useful to use within a function of the class,
// when you want to delete this and the assign it a value
void Particle::Copy_Particle(LPParticle initial_particle){
  nObs = initial_particle->nObs;
  delete partition_A;
  delete partition_B;
  partition_A = new Partition(initial_particle->partition_A,1);
  partition_B = new Partition(initial_particle->partition_B,0);
  delete[] alpha;
  delete[] beta;
  alpha = new double[nObs];
  beta = new double[nObs];
  for(int i = 0; i < nObs; i++){
    alpha[i] = initial_particle->alpha[i];
    beta[i] = initial_particle->beta[i];
  }
  if(opt_method == 1){
    delete I11;
    delete I12;
    delete I21;
    delete I22;

    I11 = new arma::mat(nObs, nObs, fill::zeros);
    I12 = new arma::mat(nObs, nObs, fill::zeros);
    I21 = new arma::mat(nObs, nObs, fill::zeros);
    I22 = new arma::mat(nObs, nObs, fill::zeros);

    *I11 = 0 + *(initial_particle->I11);
    *I12 = 0 + *(initial_particle->I12);
    *I21 = 0 + *(initial_particle->I21);
    *I22 = 0 + *(initial_particle->I22);
  }
  return;
}

// void Particle::Initialize_Particle(int n){
//   nObs = n;
//   partition_A = new Partition;
//   partition_B = new Partition;
//   partition_A->Initialize_Partition(n, 1);
//   partition_B->Initialize_Partition(n, 0);
  
//   alpha = new double[nObs];
//   beta = new double[nObs];

//   if(opt_method == 1){
//     I11 = new arma::mat(nObs, nObs, fill::zeros);
//     I12 = new arma::mat(nObs, nObs, fill::zeros);
//     I21 = new arma::mat(nObs, nObs, fill::zeros);
//     I22 = new arma::mat(nObs, nObs, fill::zeros);
//   }
//   get_alpha_beta();
//   return;
// }

void Particle::Initialize_Particle(int n, Rcpp::List gamma_init_A, Rcpp::List gamma_init_B){
  nObs = n;
  partition_A = new Partition;
  partition_B = new Partition;
  
  if(gamma_init_A.size() != 0){
    partition_A->Initialize_Partition(n, gamma_init_A, 1);
  } else {
    partition_A->Initialize_Partition(n, 1);
  }

  if(gamma_init_B.size() != 0){
    partition_B->Initialize_Partition(n, gamma_init_B, 0);
  } else {
    partition_B->Initialize_Partition(n, 0);
  }
  
  alpha = new double[nObs];
  beta = new double[nObs];

  if(opt_method == 1){
    I11 = new arma::mat(nObs, nObs, fill::zeros);
    I12 = new arma::mat(nObs, nObs, fill::zeros);
    I21 = new arma::mat(nObs, nObs, fill::zeros);
    I22 = new arma::mat(nObs, nObs, fill::zeros);
  }
  if(opt_Y == 4){
    get_alpha_beta_poisson();
  } else {
    get_alpha_beta();
  }
  return;
}


void Particle::Initialize_Particle_nclusters(int n){
  nObs = n;
  partition_A = new Partition;
  partition_B = new Partition;
  partition_A->Initialize_Partition_nclusters(n, 1);
  partition_B->Initialize_Partition_nclusters(n, 0);
  
  alpha = new double[nObs];
  beta = new double[nObs];

  if(opt_method == 1){
    I11 = new arma::mat(nObs, nObs, fill::zeros);
    I12 = new arma::mat(nObs, nObs, fill::zeros);
    I21 = new arma::mat(nObs, nObs, fill::zeros);
    I22 = new arma::mat(nObs, nObs, fill::zeros);
  }
  if(opt_Y == 4){
    get_alpha_beta_poisson();
  } else {
    get_alpha_beta();
  }
  return;
}

void Particle::get_alpha_beta(){
  if(opt_Y == 4) 
    cout << "WARNING: Poisson likelihood selected but estimating alpha and beta using normal likelihood!" << endl;
  // double a1 = 100;
  // double a2 = 100;
  // double b1 = 100;
  // double b2 = 100;

  int T = Y.n_cols;
  rowvec x = X.row(0);

  mat E = eye<mat>(nObs,nObs);
  mat O = ones<mat>(nObs,nObs);
  mat Z = zeros<mat>(nObs,nObs);
  mat O_k, A_block_k, Omega_k, D, M;
  mat Omega_a = zeros<mat>(nObs,nObs);
  mat Omega_b = zeros<mat>(nObs,nObs);
  // I am updating Omega_a (the inverse of the covariance matrix of alpha) cluster-wise
  for(int k = 0; k < partition_A->K; k++){
    int cluster_size = partition_A->cluster_config[k]; 
    if(cluster_size == 1){
      Omega_a(partition_A->clusters[k][0], partition_A->clusters[k][0]) = 1/(a1/(1-rho)+a2);
    } else {
      O_k = ones<mat>(cluster_size,cluster_size);
      // Creating M: the precision matrix of the CAR model
      A_block_k = Submatrix(A_block, cluster_size, cluster_size, partition_A->clusters[k], partition_A->clusters[k]);
      vec row_sums = zeros<vec>(cluster_size);
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          row_sums(i) += A_block_k(i,j);
        }
      }
      D = diagmat(row_sums);
      M = rho * (D - A_block_k);
      for(int i = 0; i < cluster_size; i++){
        M(i,i) += 1 - rho;
      }
      // This is the Woodbury matrix identity for inversion
      // Omega_k = (M/a1) - (M/a1) * O_k * (M/a1) / (sum(sum(M))/a1 + 1/a2);
      Omega_k = (M/a1) - ((1-rho)/a1)*((1-rho)/a1) / ( cluster_size*(1-rho)/a1 + 1/a2) * O_k;
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          Omega_a(partition_A->clusters[k][i], partition_A->clusters[k][j]) = Omega_k(i,j);
        }
      }
    }
  }
  // I am updating Omega_b (the inverse of the covariance matrix of beta) cluster-wise
  for(int k = 0; k < partition_B->K; k++){
    int cluster_size = partition_B->cluster_config[k];   
    if(cluster_size == 1){
      Omega_b(partition_B->clusters[k][0], partition_B->clusters[k][0]) = 1/(b1/(1-rho)+b2);
    } else { 
      O_k = ones<mat>(cluster_size,cluster_size);
      // Creating M: the precision matrix of the CAR model
      A_block_k = Submatrix(A_block, cluster_size, cluster_size, partition_B->clusters[k], partition_B->clusters[k]);
      vec row_sums = zeros<vec>(cluster_size);
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          row_sums(i) += A_block_k(i,j);
        }
      }
      D = diagmat(row_sums);
      M = rho * (D - A_block_k);
      for(int i = 0; i < cluster_size; i++){
        M(i,i) += 1 - rho;
      }
      // This is the Woodbury matrix identity for inversion
      // Omega_k = (M/b1) - (M/b1) * O_k * (M/b1) / (sum(sum(M))/b1 + 1/b2);
      Omega_k = (M/b1) - ((1-rho)/b1)*((1-rho)/b1) / ( cluster_size*(1-rho)/b1 + 1/b2) * O_k;
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          Omega_b(partition_B->clusters[k][i], partition_B->clusters[k][j]) = Omega_k(i,j);
        }
      }
    }
  }

  rowvec x2 = pow(x,2);

  vec xty1 = sum(Y,1);
  vec xty2 = sum(X%Y,1);
  arma::vec XtY = arma::join_cols(xty1,xty2);
  
  // mat Sigma0inv = arma::join_cols(arma::join_rows(Omega_a, Z),arma::join_rows(Z, Omega_b));
  // mat XtX = arma::join_cols(arma::join_rows(E*T, E*sum(x)),arma::join_rows(E*sum(x), E*sum(x2)));
  // mat P = inv_sympd(XtX + Sigma0inv); // this should be replaced with the block inversion formula
  mat P;
  if(opt_method == 1){
    P = block_inverse_ret(E*T + Omega_a, E*sum(x), E*sum(x), E*sum(x2)+Omega_b, I11, I12, I21, I22);  
  } else {
    mat L11, L12, L21, L22;
    P = block_inverse_ret(E*T + Omega_a, E*sum(x), E*sum(x), E*sum(x2)+Omega_b, &L11, &L12, &L21, &L22);  
  }
  
  vec alpha_beta = P * XtY;
  
  for(int i = 0; i < nObs; i++){
    alpha[i] = alpha_beta(i);
    beta[i] = alpha_beta(i+nObs);
  }
  return;
}

void Particle::get_alpha_beta_mle_normal(double* alpha_mle, double* beta_mle){
  double a1 = 100;
  double a2 = 100;
  double b1 = 100;
  double b2 = 100;

  int T = Y.n_cols;
  rowvec x = X.row(0);

  mat E = eye<mat>(nObs,nObs);
  mat O = ones<mat>(nObs,nObs);
  mat Z = zeros<mat>(nObs,nObs);
  mat O_k, A_block_k, Omega_k, D, M;
  mat Omega_a = zeros<mat>(nObs,nObs);
  mat Omega_b = zeros<mat>(nObs,nObs);
  // I am updating Omega_a (the inverse of the covariance matrix of alpha) cluster-wise
  for(int k = 0; k < partition_A->K; k++){
    int cluster_size = partition_A->cluster_config[k]; 
    if(cluster_size == 1){
      Omega_a(partition_A->clusters[k][0], partition_A->clusters[k][0]) = 1/(a1/(1-rho)+a2);
    } else {
      O_k = ones<mat>(cluster_size,cluster_size);
      // Creating M: the precision matrix of the CAR model
      A_block_k = Submatrix(A_block, cluster_size, cluster_size, partition_A->clusters[k], partition_A->clusters[k]);
      vec row_sums = zeros<vec>(cluster_size);
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          row_sums(i) += A_block_k(i,j);
        }
      }
      D = diagmat(row_sums);
      M = rho * (D - A_block_k);
      for(int i = 0; i < cluster_size; i++){
        M(i,i) += 1 - rho;
      }
      // This is the Woodbury matrix identity for inversion
      // Omega_k = (M/a1) - (M/a1) * O_k * (M/a1) / (sum(sum(M))/a1 + 1/a2);
      Omega_k = (M/a1) - ((1-rho)/a1)*((1-rho)/a1) / ( cluster_size*(1-rho)/a1 + 1/a2) * O_k;
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          Omega_a(partition_A->clusters[k][i], partition_A->clusters[k][j]) = Omega_k(i,j);
        }
      }
    }
  }
  // I am updating Omega_b (the inverse of the covariance matrix of beta) cluster-wise
  for(int k = 0; k < partition_B->K; k++){
    int cluster_size = partition_B->cluster_config[k];   
    if(cluster_size == 1){
      Omega_b(partition_B->clusters[k][0], partition_B->clusters[k][0]) = 1/(b1/(1-rho)+b2);
    } else { 
      O_k = ones<mat>(cluster_size,cluster_size);
      // Creating M: the precision matrix of the CAR model
      A_block_k = Submatrix(A_block, cluster_size, cluster_size, partition_B->clusters[k], partition_B->clusters[k]);
      vec row_sums = zeros<vec>(cluster_size);
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          row_sums(i) += A_block_k(i,j);
        }
      }
      D = diagmat(row_sums);
      M = rho * (D - A_block_k);
      for(int i = 0; i < cluster_size; i++){
        M(i,i) += 1 - rho;
      }
      // This is the Woodbury matrix identity for inversion
      // Omega_k = (M/b1) - (M/b1) * O_k * (M/b1) / (sum(sum(M))/b1 + 1/b2);
      Omega_k = (M/b1) - ((1-rho)/b1)*((1-rho)/b1) / ( cluster_size*(1-rho)/b1 + 1/b2) * O_k;
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          Omega_b(partition_B->clusters[k][i], partition_B->clusters[k][j]) = Omega_k(i,j);
        }
      }
    }
  }

  rowvec x2 = pow(x,2);

  vec xty1 = sum(Y,1);
  vec xty2 = sum(X%Y,1);
  arma::vec XtY = arma::join_cols(xty1,xty2);
  
  // mat Sigma0inv = arma::join_cols(arma::join_rows(Omega_a, Z),arma::join_rows(Z, Omega_b));
  // mat XtX = arma::join_cols(arma::join_rows(E*T, E*sum(x)),arma::join_rows(E*sum(x), E*sum(x2)));
  // mat P = inv_sympd(XtX + Sigma0inv); // this should be replaced with the block inversion formula
  mat P;
  if(opt_method == 1){
    P = block_inverse_ret(E*T + Omega_a, E*sum(x), E*sum(x), E*sum(x2)+Omega_b, I11, I12, I21, I22);  
  } else {
    mat L11, L12, L21, L22;
    P = block_inverse_ret(E*T + Omega_a, E*sum(x), E*sum(x), E*sum(x2)+Omega_b, &L11, &L12, &L21, &L22);  
  }
  
  vec alpha_beta = P * XtY;
  
  for(int i = 0; i < nObs; i++){
    alpha_mle[i] = alpha_beta(i);
    beta_mle[i] = alpha_beta(i+nObs);
  }
  return;
}


void Particle::get_alpha_beta_mle_poisson(double* alpha_mle, double* beta_mle){
  dlib::matrix<double,0,1> alpha_beta(2*(Y.n_rows));
  for(int i =0; i < Y.n_rows; i++){
    alpha_beta(i) = 0;
    alpha_beta(i+Y.n_rows) = 0;
  }
  optim_mle(alpha_beta);

  for(int i = 0; i < nObs; i++){
    alpha_mle[i] = alpha_beta(i);
    beta_mle[i] = alpha_beta(i+nObs);
  }
  return;
}

void Particle::get_alpha_beta_mle(double* alpha_mle, double* beta_mle){
  if(opt_Y == 4){
    get_alpha_beta_mle_poisson(alpha_mle, beta_mle);
  } else {
    get_alpha_beta_mle_normal(alpha_mle, beta_mle);
  }
}

void Particle::get_alpha_beta_iter(){
  LPPartition partition;
  for(int iter = 0; iter < coordasc_iter; iter++){
    partition = partition_A;
    for(int k = 0; k < partition->K; k++){
      get_parameter(1,k);
    }
    partition = partition_B;
    for(int k = 0; k < partition->K; k++){
      get_parameter(0,k);
    }
  }
}

void Particle::get_alpha_beta_poisson(){
  dlib::matrix<double,0,1> alpha_beta(2*(Y.n_rows));
  for(int i =0; i < Y.n_rows; i++){
    alpha_beta(i) = 0;
    alpha_beta(i+Y.n_rows) = 0;
  }
  optim(alpha_beta);

  for(int i = 0; i < nObs; i++){
    alpha[i] = alpha_beta(i);
    beta[i] = alpha_beta(i+nObs);
  }
  return;
}

void Particle::get_parameter(bool A_or_B, int cluster_id){
  if(opt_Y == 4) 
    cout << "WARNING: Poisson likelihood selected but estimating alpha and beta using normal likelihood!" << endl;
  
  // double a1 = 100;
  // double a2 = 100;
  // double b1 = 100;
  // double b2 = 100;
  
  if(A_or_B){
    int T = Y.n_cols;
    int cluster_size = partition_A->cluster_config[cluster_id];
    rowvec x = X.row(0);
    double txx = T; // this is correct for A_or_B = 1
    
    if(cluster_size == 1){
      vec y_adj(T);
      for(int t = 0; t < T; t++){
        // y_adj(t) = Y(partition_A->clusters[cluster_id][0],t); // they are equivalent when X is orthogonal
        y_adj(t) = Y(partition_A->clusters[cluster_id][0],t) - x(t)*beta[partition_A->clusters[cluster_id][0]];
      }
      double Omega_beta_k = 1/(a1 + a2);
      
      double tXY = sum(y_adj);
      alpha[partition_A->clusters[cluster_id][0]] = 1/(txx + Omega_beta_k) * tXY;
    } else {
      mat O = ones<mat>(cluster_size, cluster_size);
      mat E = eye<mat>(cluster_size, cluster_size);
      mat A_block_k = Submatrix(A_block, cluster_size, cluster_size, partition_A->clusters[cluster_id], partition_A->clusters[cluster_id]);
      // get rowsums
      vec row_sums = zeros<vec>(cluster_size);
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          row_sums(i) += A_block_k(i,j);
        }
      }
      mat D = diagmat(row_sums);
      mat A_star_k = D - A_block_k;
      mat M = rho * A_star_k;
      for(int i = 0; i < cluster_size; i++){
        M(i,i) += 1 - rho;
      }
      // mat Sigma0inv = (M/a1) - (M/a1) * O * (M/a1) / (sum(sum(M))/a1 + 1/a2);
      mat Sigma0inv = (M/a1) - ((1-rho)/a1)*((1-rho)/a1) / ( cluster_size*(1-rho)/a1 + 1/a2) * O;
      mat Y_adj = zeros<mat>(cluster_size,T);
      for(int i = 0; i < cluster_size; i++){
        // Y_adj.row(i) = Y.row(partition_A->clusters[cluster_id][i]); // they are equivalent when X is orthogonal
        Y_adj.row(i) = Y.row(partition_A->clusters[cluster_id][i]) - x*beta[partition_A->clusters[cluster_id][i]];
      }
      vec tXY = sum(Y_adj, 1);
      
      mat tXX = E * txx;
      mat temp = inv_sympd(tXX + Sigma0inv);
      vec tmp_alpha = temp * tXY;
      for(int i = 0; i < cluster_size; i++){
        alpha[partition_A->clusters[cluster_id][i]] = tmp_alpha(i);
      }
    }
  } else {
    int T = Y.n_cols;
    int cluster_size = partition_B->cluster_config[cluster_id];
    rowvec x = X.row(0);
    double txx = as_scalar(x * x.t());

    if(cluster_size == 1){
      vec y_adj(T);
      for(int t = 0; t < T; t++){
        // y_adj(t) = Y(partition_B->clusters[cluster_id][0],t); // they are equivalent when X is orthogonal
        y_adj(t) = Y(partition_B->clusters[cluster_id][0],t) - alpha[partition_B->clusters[cluster_id][0]];
      }
      double Omega_beta_k = 1/(b1 + b2);
      
      double tXY = as_scalar(x * y_adj);
      beta[partition_B->clusters[cluster_id][0]] = 1/(txx + Omega_beta_k) * tXY;
    } else {
      mat O = ones<mat>(cluster_size, cluster_size);
      mat E = eye<mat>(cluster_size, cluster_size);
      mat A_block_k = Submatrix(A_block, cluster_size, cluster_size, partition_B->clusters[cluster_id], partition_B->clusters[cluster_id]);
      // get rowsums
      vec row_sums = zeros<vec>(cluster_size);
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          row_sums(i) += A_block_k(i,j);
        }
      }
      mat D = diagmat(row_sums);
      mat A_star_k = D - A_block_k;
      mat M = rho * A_star_k;
      for(int i = 0; i < cluster_size; i++){
        M(i,i) += 1 - rho;
      }
      // mat Sigma0inv = (M/b1) - (M/b1) * O * (M/b1) / (sum(sum(M))/b1 + 1/b2);
      mat Sigma0inv = (M/b1) - ((1-rho)/b1)*((1-rho)/b1) / ( cluster_size*(1-rho)/b1 + 1/b2) * O;
      mat Y_adj = zeros<mat>(cluster_size,T);
      mat XtY = zeros<mat>(cluster_size,T);
      // this could change: we could have only the vector tXY and inside the loop..
      for(int i = 0; i < cluster_size; i++){
        // Y_adj.row(i) = Y.row(partition_B->clusters[cluster_id][i]); // they are equivalent when X is orthogonal
        Y_adj.row(i) = Y.row(partition_B->clusters[cluster_id][i]) - alpha[partition_B->clusters[cluster_id][i]];
        XtY.row(i) = x % Y_adj.row(i);
        // here we could do tXY(i) = sum(x % Y_adj.row(i));
      }
      vec tXY = sum(XtY, 1);
      
      mat tXX = E * txx;
      mat temp = inv_sympd(tXX + Sigma0inv);
      vec tmp_beta = temp * tXY;
      for(int i = 0; i < cluster_size; i++){
        beta[partition_B->clusters[cluster_id][i]] = tmp_beta(i);
      }
    }
  }
  return;
}

// if A_or_B = True we do A
void Particle::Partition_Split(bool A_or_B, int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2){
  int orig_K;
  if(A_or_B == true){
    orig_K = partition_A->K;
    partition_A->Split(split_k, new_cluster1, new_cluster2, size1, size2, A_or_B);
  } else {
    orig_K = partition_B->K;
    partition_B->Split(split_k, new_cluster1, new_cluster2, size1, size2, A_or_B);
  }
  if(opt_Y == 4){
    get_alpha_beta_poisson();
  } else {
    if(opt_method == 0){
      // just do what we were doing before
      for(int k = 0; k < orig_K+1; k++){
        if(k == split_k){ // need to re-compute
          get_parameter(A_or_B, k);
        } else if(k == orig_K){
          get_parameter(A_or_B, k);
        } 
      }
    } else if(opt_method == 1){
      // optimize alpha_beta given the two partitions
      get_alpha_beta();
    } else if(opt_method == 2){
      // optimize alpha_beta given the two partitions, by doing coordinate ascend
      get_alpha_beta_iter();
    }
  }
  return;
}

void Particle::Partition_Split(bool A_or_B, int split_k, std::vector<int> new_cluster1, std::vector<int> new_cluster2, int size1, int size2){
  int orig_K;
  if(A_or_B == true){
    orig_K = partition_A->K;
    partition_A->Split(split_k, new_cluster1, new_cluster2, size1, size2, A_or_B);
  } else {
    orig_K = partition_B->K;
    partition_B->Split(split_k, new_cluster1, new_cluster2, size1, size2, A_or_B);
  }
  if(opt_Y == 4){
    get_alpha_beta_poisson();
  } else {
    if(opt_method == 0){
      // just do what we were doing before
      for(int k = 0; k < orig_K+1; k++){
        if(k == split_k){ // need to re-compute
          get_parameter(A_or_B, k);
        } else if(k == orig_K){
          get_parameter(A_or_B, k);
        } 
      }
    } else if(opt_method == 1){
      // optimize alpha_beta given the two partitions
      get_alpha_beta();
    } else if(opt_method == 2){
      // optimize alpha_beta given the two partitions, by doing coordinate ascend
      get_alpha_beta_iter();
    }
  }
  return;
}

void Particle::Partition_KSplit(bool A_or_B, int split_k, std::vector<std::vector<int> > indices, std::vector<int> ns){
  int num_splits = indices.size();
  if(num_splits != (int)ns.size()){
    cout << "Partition_KSplit: ERROR, ns.size should be equal to indices.size()" << endl;
  }
  int orig_K;
  if(A_or_B == true){
    orig_K = partition_A->K;
    partition_A->KSplit(split_k, num_splits, indices, ns, A_or_B);
  } else {
    orig_K = partition_B->K;
    partition_B->KSplit(split_k, num_splits, indices, ns, A_or_B);
  }
  if(opt_Y == 4){
    get_alpha_beta_poisson();
  } else {
    if(opt_method == 0){
      // just re-compute for the new clusters
      int new_K = orig_K + num_splits - 1;
      for(int k = 0; k < new_K; k++){
        if(k == split_k){ 
          get_parameter(A_or_B, k);
        } else if(k >= orig_K){
          get_parameter(A_or_B, k);
        } 
      }
    } else if(opt_method == 1){
      // optimize alpha_beta given the two partitions
      get_alpha_beta();
    } else if(opt_method == 2){
      // optimize alpha_beta given the two partitions, by doing coordinate ascend
      get_alpha_beta_iter();
    }
  }
  return;
}


void Particle::Partition_Merge(bool A_or_B, int k_1, int k_2){
  // int k_max = max(k_1, k_2);
  int k_min = min(k_1, k_2);
  int orig_K;
  if(A_or_B == true){
    orig_K = partition_A->K;
    partition_A->Merge(k_1, k_2, A_or_B);
  } else {
    orig_K = partition_B->K;
    partition_B->Merge(k_1, k_2, A_or_B);
  }
  if(opt_Y == 4){
    get_alpha_beta_poisson();
  } else {
    if(opt_method == 0){
      // just do what we were doing before
      for(int k = 0; k < orig_K; k++){
        if(k == k_min){
          get_parameter(A_or_B, k);
        }
      }
    } else if(opt_method == 1){
      // optimize alpha_beta given the two partitions
      get_alpha_beta();
    } else if(opt_method == 2){
      // optimize alpha_beta given the two partitions, by doing coordinate ascend
      get_alpha_beta_iter();
    }
  }
  return;
}


void Particle::Partition_Split_and_Merge(bool A_or_B, int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2, int k_star_1, int k_star_2){
  // this should be good, so they are pointing at the same thing.
  LPPartition partition;
  if(A_or_B){
    partition = partition_A;
  } else {
    partition = partition_B;
  }
  if(k_star_1 != k_star_2){
    Partition_Split(A_or_B, split_k, new_cluster1, new_cluster2, size1, size2);
    if((split_k == k_star_1) & (split_k != k_star_2)){ // leave new_cluster1 alone and just attempt to merge new_cluster2 with k_star_2
      Partition_Merge(A_or_B, partition->K-1, k_star_2); // remember new_cluster2's label is K-1
    } else if( (split_k != k_star_1) & (split_k == k_star_2)){
      Partition_Merge(A_or_B, split_k, k_star_1);
    } else if((split_k != k_star_1) & (split_k != k_star_2) & (k_star_2 > max(split_k, k_star_1))){
      Partition_Merge(A_or_B, split_k, k_star_1);
      Partition_Merge(A_or_B, partition->K, k_star_2 - 1);
    } else if((split_k != k_star_1) & (split_k != k_star_2) & (k_star_2 < max(split_k, k_star_1))){
      Partition_Merge(A_or_B, split_k, k_star_1);
      Partition_Merge(A_or_B, partition->K, k_star_2);
    }
  }
  return;
}

void Particle::Partition_Modify(bool A_or_B, int cl_ind){
  // this should be good, so they are pointing at the same thing.
  LPPartition partition;
  if(A_or_B){
    partition = partition_A;
  } else {
    partition = partition_B;
  }
  // cl_ind is the index of the cluster that needs to be modified
  if(partition->cluster_config[cl_ind] > 1)
  {
    int n = partition->cluster_config[cl_ind];
    int *index = partition->clusters[cl_ind];
    mat A_cluster = Submatrix(A_block, n, n, index, index);
    int *components = new int[n];
    int *count = new int;
    Connected_Components(A_cluster, n, components, count);
    if(*count > 1)
    {
      //I will use split iteratively, first split 0 from !=0
      //then among the remaining, split 1 from !=1, etc
      *count = *count - 1;
      int *new_components, *index1, *index2, *i_ind, n1, n2;
      for(int tmp = 0; tmp < *count; tmp++)
      {
        index1 = new int[n];
        index2 = new int[n];
        i_ind = new int[n];
        n1 = 0;
        n2 = 0;
        for(int i = 0; i < n; i++){
          if(components[i] == tmp){
            index1[n1] = index[i];
            n1++;
          } else {
            index2[n2] = index[i];
            i_ind[n2] = i;
            n2++;
          }
        }
        Partition_Split(A_or_B, cl_ind, index1, index2, n1, n2);
        if(tmp > 0)
          delete[] index;
        index = index2;
        n = n2;
        new_components = new int[n2];
        for(int j = 0; j < n2; j++)
          new_components[j] = components[i_ind[j]];
        delete[] components;
        components = new_components;
        cl_ind = partition->K-1;
        delete[] index1;
      }
      delete[] index2;
    }
    delete[] components;
    delete count;
  }
}

// find splits needs more work, because of all the computation
// same for ksplits

void Particle::Update_Partition(bool A_or_B, int current_l, std::vector<LPParticle> Particle_Set, std::vector<double> w, double lambda, double xi, bool restricted){
  string max_str = "Not changed";
  LPParticle max_candidate;

  max_candidate = new Particle(this); // holds the running "best local candidate"
  double max_objective = 0.0;
  max_objective = w[current_l]*max_candidate->Total_Log_Post() + lambda * Entropy(current_l, max_candidate, Particle_Set, w) + xi * VI_Avg(current_l, max_candidate, Particle_Set);
  if(!restricted){
    // Islands - 5% in each cluster
    // cout << "before island" << endl;
    Island_moves(A_or_B, &max_objective, current_l, max_candidate, Particle_Set, w, lambda, xi, &max_str);
    // cout << "before border" << endl;
    Border_moves(A_or_B, &max_objective, current_l, max_candidate, Particle_Set, w, lambda, xi, &max_str);
    // Merges: Sweep over all clusters and find the ones that are adjacent
    // cout << "before merge" << endl;
    Merge_moves(A_or_B, &max_objective, current_l, max_candidate, Particle_Set, w, lambda, xi, &max_str);
    // we will loop over the clusters, run spectral clustering on the beta_hat's within the cluster, and propose various split and merge candidates
    // cout << "before splitmerge" << endl;
    Spectral_SplitMerge_moves(A_or_B, &max_objective, current_l, max_candidate, Particle_Set, w, lambda, xi, &max_str);
    Spectral_SplitMerge_moves(A_or_B, &max_objective, current_l, max_candidate, Particle_Set, w, lambda, xi, &max_str, true);
    // cout << "before tailsplit" << endl;
    Tail_Split_moves(A_or_B, &max_objective, current_l, max_candidate, Particle_Set, w, lambda, xi, &max_str);
  }
  // cout << "before km" << endl;
  KM_SplitMerge_moves(A_or_B, &max_objective, current_l, max_candidate, Particle_Set, w, lambda, xi, &max_str);
  KM_SplitMerge_moves(A_or_B, &max_objective, current_l, max_candidate, Particle_Set, w, lambda, xi, &max_str, true);
  // cout << "after km" << endl;
  Copy_Particle(max_candidate);
  cout << max_str << " " << max_objective << " ";
  cout << endl;
  delete max_candidate;
  return;
}

std::vector<int> Particle::get_jstar(bool A_or_B, int split_k, int num_splits, std::vector< std::vector<int> > indices, std::vector<int> ns){
  LPPartition partition;
  double* parameter;
  if(A_or_B){
    partition = partition_A;
    parameter = alpha;
  } else {
    partition = partition_B;
    parameter = beta;
  }

  std::vector<std::vector<int> > adj_clusters;
  std::vector<int> jstar(num_splits);
  double beta_bar, tmp, tmp_min;
  adj_clusters.resize(num_splits);
  arma::mat tmp_A;
  for(int j = 0; j < num_splits; j++){
    arma::vec param_cluster_vec(ns[j]);
    for(int i = 0; i < ns[j]; i++){
      param_cluster_vec(i) = parameter[ indices[j][i] ];
    }
    beta_bar = arma::mean(param_cluster_vec);
    for(int kk = 0; kk < partition->K; kk++){
      if(kk != split_k){
        tmp_A = Submatrix(A_block, ns[j], partition->cluster_config[kk], indices[j], partition->clusters[kk]);
        if(any(vectorise(tmp_A == 1.0))){
          adj_clusters[j].push_back(kk);
        }
      }
    }
    if((int)adj_clusters[j].size() > 0){ // find jstar (the closest in mean adjacent cluster)
      tmp_min = abs(param_bar(A_or_B, adj_clusters[j][0])-beta_bar);
      jstar[j] = adj_clusters[j][0];
      for(int ii = 1; ii < (int)adj_clusters[j].size(); ii++){
        tmp = abs(param_bar(A_or_B, adj_clusters[j][ii])-beta_bar);
        if(tmp < tmp_min){
          tmp_min = tmp;
          jstar[j] = adj_clusters[j][ii];
        }
      }
    } else { // -1 means no adjacent cluster
      jstar[j] = -1;
    }
  }
  return jstar;
}

std::vector<int> Particle::get_tail_jstar(bool A_or_B, int split_k, std::vector<std::vector<int> > tail_conncomp){
  LPPartition partition;
  double* parameter;
  if(A_or_B){
    partition = partition_A;
    parameter = alpha;
  } else {
    partition = partition_B;
    parameter = beta;
  }
  std::vector<int> left_jstar;
  arma::vec param_cluster;
  double param_mean;
  std::vector<double> tmp_dist; // holds distance from the new sub-clusters to adjacent existing clusters
  std::vector<double> tmp_nn; // holds indices of potential nearest neighbors for new sub-clusters
  arma::vec tmp_dist_vec; // arma vector version of tmp_dist
  arma::uvec tmp_dist_indices; // for soring tmp_dist_vec
  arma::mat A_tmp;
  
  for(int new_k = 0; new_k < (int)tail_conncomp.size(); new_k++){
    // now figure out the connected components
    param_cluster.clear();
    param_cluster.set_size( tail_conncomp[new_k].size() );
    for(int j = 0; j < (int)tail_conncomp[new_k].size(); j++){
      param_cluster(j) = parameter[ tail_conncomp[new_k][j] ];
    }
    param_mean = arma::mean(param_cluster);

    tmp_nn.clear();
    tmp_dist.clear();
    for(int kk = 0; kk < partition->K; kk++){
      if(kk != split_k){
        A_tmp = Submatrix(A_block, tail_conncomp[new_k].size(), partition->cluster_config[kk], tail_conncomp[new_k], partition->clusters[kk]);
        if(any(vectorise(A_tmp == 1.0))){
          tmp_nn.push_back(kk);
          param_cluster.clear();
          param_cluster.set_size( partition->cluster_config[kk] );
          for(int j = 0; j < partition->cluster_config[kk]; j++){
            param_cluster(j) = parameter[ partition->clusters[kk][j] ];
          }
          tmp_dist.push_back(abs(param_mean - arma::mean(param_cluster)));
        }
      }
    }
    if(tmp_nn.size() > 0){
      tmp_dist_vec.reset();
      tmp_dist_indices.reset();
      tmp_dist_vec.set_size(tmp_dist.size());
      tmp_dist_indices.set_size(tmp_dist.size());
      for(int kk = 0; kk < (int)tmp_dist.size(); kk++){
        tmp_dist_vec(kk) = tmp_dist[kk];
      }
      tmp_dist_indices = sort_index(tmp_dist_vec,"ascend");
      left_jstar.push_back(tmp_nn[tmp_dist_indices(0)]);
    } else{
      left_jstar.push_back(-1);
    }
  }
  return left_jstar;
}

void Particle::Island_moves(bool A_or_B, double *pointer_to_maxobjective, int current_l, LPParticle& max_candidate, std::vector<LPParticle> Particle_Set, std::vector<double> w, double lambda, double xi, string *pt_maxstr, double percentage){
  // cout << "START ISLAND" << endl;
  LPPartition partition;
  double* parameter;
  if(A_or_B){
    partition = partition_A;
    parameter = alpha;
  } else {
    partition = partition_B;
    parameter = beta;
  }
  double tmp_objective = 0.0; 
  int split_k, size1, size2;
  // int* new_cluster1;
  // int* new_cluster2;
  // new_cluster1 = new int[1];
  // new_cluster2 = new int[1];
  LPParticle tmp_candidate = new Particle(this);
  
  for(int k = 0; k < partition->K; k++){
    // cout << "*Start cluster " << k << endl;
    // find which elements are in the top - bottom 5%
    if(partition->cluster_config[k] > 1){
      
      double* beta_hat_copy = new double[ partition->cluster_config[k] ];
      for(int i = 0; i < partition->cluster_config[k]; i++){
        beta_hat_copy[i] = parameter[ partition->clusters[k][i] ];
      }
      sort(beta_hat_copy, beta_hat_copy + partition->cluster_config[k]);
      int threshold = (int) ceil(partition->cluster_config[k] * percentage);
      // std::vector<int> index_topbottom(2*threshold);
      std::vector<int> index_topbottom;
      index_topbottom.reserve(2*threshold);
      // int *index_topbottom = new int[2*threshold];
      // int count = 0;
      for(int i = 0; i < threshold; i++){
        for(int j = 0; j < partition->cluster_config[k]; j++){
          if(parameter[ partition->clusters[k][j] ] == beta_hat_copy[i]){
            index_topbottom.push_back(partition->clusters[k][j]);
            // count ++;
          }
        }
        for(int j = 0; j < partition->cluster_config[k]; j++){
          if(parameter[ partition->clusters[k][j] ] == beta_hat_copy[(partition->cluster_config[k])-1-i]){
            index_topbottom.push_back(partition->clusters[k][j]);
            // count ++;
          }
        }
      }

      // index_topbottom contains the index in the whole city like {1, ..255}, not the index in the cluster
      delete[] beta_hat_copy;
      // cout << "*Done index - 2*threshold: " << 2*threshold << " and index_topbottom.size() " << index_topbottom.size() << endl;
      for(int c = 0; c < 2*threshold; c++){
        // cout << "**Start element " << c << endl;
        int i = index_topbottom[c];
        // cout << "after i" << endl;
        split_k = k;
        delete tmp_candidate; // delete the old value of tmp_candidate .. i can either use Copy_Partition or do a delete[] and new call
        // cout << "after delete tmp_candidate" << endl;
        tmp_candidate = new Particle(this);
        // cout << "after new Particle" << endl;
        size1 = partition->cluster_config[split_k] - 1; 
        size2 = 1;
        std::vector<int> new_cluster1;
        std::vector<int> new_cluster2(size2);
        new_cluster1.reserve(size1);
        new_cluster2[0] = i;
        for(int ii = 0; ii < partition->cluster_config[split_k]; ii++){
          if(partition->clusters[split_k][ii] != i){
            new_cluster1.push_back(partition->clusters[split_k][ii]);
          }
        }
        // now actually perform the split
        tmp_candidate->Partition_Split(A_or_B, split_k, new_cluster1, new_cluster2, size1, size2);
        tmp_candidate->Partition_Modify(A_or_B, split_k);
        tmp_objective = w[current_l]*tmp_candidate->Total_Log_Post() + lambda * Entropy(current_l, tmp_candidate, Particle_Set, w) + xi * VI_Avg(current_l, tmp_candidate, Particle_Set);
        if(tmp_objective > *pointer_to_maxobjective){
          delete max_candidate;
          try{
            max_candidate = new Particle(tmp_candidate);
          }
          catch(const std::bad_alloc& e){
            cout << "######## EXCEPTION 4 ########"  << e.what() << endl;
          }
          *pointer_to_maxobjective = tmp_objective;
          *pt_maxstr = "Island";
        }
      }
    }
  }
  delete tmp_candidate;
  // delete[] new_cluster1;
  // delete[] new_cluster2;
  // cout << "END ISLAND" << endl;
  return;
}

void Particle::Border_moves(bool A_or_B, double *pointer_to_maxobjective, int current_l, LPParticle& max_candidate, std::vector<LPParticle> Particle_Set, std::vector<double> w, double lambda, double xi, string *pt_maxstr){
  LPPartition partition;
  // double* parameter;
  if(A_or_B){
    partition = partition_A;
    // parameter = alpha;
  } else {
    partition = partition_B;
    // parameter = beta;
  }
  
  double tmp_objective = 0.0; 
  int split_k, k_star_1, k_star_2;
  std::vector<int> neighboring_clusters; // holds the labels of clusters adjacent to block-group i
  int neighbor_counter = 0;
  int size1, size2;
  int* new_cluster1;
  int* new_cluster2;
  new_cluster1 = new int[1];
  new_cluster2 = new int[1];
  LPParticle tmp_candidate = new Particle(this);

  for(int i = 0; i < partition->nObs; i++){
    split_k = partition->cluster_assignment[i];
    if(partition->cluster_config[split_k] > 1){ // only attempt if there are at least 2 elements in that cluster. Otherwise this is just a Merge move
      neighbor_counter = 0;
      neighboring_clusters.clear();
      for(int k = 0; k < partition->K;k++){
        if(k != split_k){ // loop over the elements of cluster k to see if any are adjacent to i
          for(int j = 0; j < partition->cluster_config[k]; j++){
            if(A_block(i, partition->clusters[k][j]) == 1){ // j^th blockgroup of cluster k is adjacent to block group i
              neighboring_clusters.push_back(k); // add k to the vector of adjacent
              neighbor_counter++;
              break; // once we know that i is adjacent to something in cluster k we can stop
            }
          }
        }
      }
      if(neighbor_counter > 0){
        // we are splitting i away from its cluster
        // let us define all of the necessary components for that split
        size1 = partition->cluster_config[split_k] - 1;
        size2 = 1;
        delete[] new_cluster1;
        delete[] new_cluster2;
        try{
          new_cluster1 = new int[size1];
          new_cluster2 = new int[size2];
        }
        catch(const std::bad_alloc& e){
          cout << "######## EXCEPTION 34b ########"  << e.what() << endl;
        }
        
        new_cluster2[0] = i;
        int counter = 0;
        for(int ii = 0; ii < partition->cluster_config[split_k]; ii++){
          if(partition->clusters[split_k][ii] != i){
            new_cluster1[counter] = partition->clusters[split_k][ii];
            counter++;
          }
        }
        // new_cluster1 contains all of the things in cluster split_k EXCEPT for i
        // it will be "merged" with itself (i.e. k_star_1 = split_k)
        k_star_1 = split_k;
        for(int nc = 0; nc < neighbor_counter; nc++){
          k_star_2 = neighboring_clusters[nc]; // this is the new cluster that we are adding i to!
          delete tmp_candidate;
          tmp_candidate = new Particle(this);
          
          tmp_candidate->Partition_Split_and_Merge(A_or_B, split_k, new_cluster1, new_cluster2, size1, size2, k_star_1, k_star_2);
          tmp_candidate->Partition_Modify(A_or_B, split_k);
          tmp_objective = w[current_l]*tmp_candidate->Total_Log_Post() + lambda * Entropy(current_l, tmp_candidate, Particle_Set, w) + xi * VI_Avg(current_l, tmp_candidate, Particle_Set);
          // cout << tmp_objective << endl;
          if(tmp_objective > *pointer_to_maxobjective){
            delete max_candidate;
            try{
              max_candidate = new Particle(tmp_candidate);
            }
            catch(const std::bad_alloc& e){
              cout << "######## EXCEPTION 6 ########"  << e.what() << endl;
            }
            *pointer_to_maxobjective = tmp_objective;
            *pt_maxstr = "Border";
          }
        } // closes loop over the neighboring clusters
      } // closes if statement that checks if there are neighboring clusters
    } // closes if cluster_config>1
  } // closes loop over the blockgroups
  delete tmp_candidate;
  delete[] new_cluster1;
  delete[] new_cluster2;
  return;
}

void Particle::Merge_moves(bool A_or_B, double *pointer_to_maxobjective, int current_l, LPParticle& max_candidate, std::vector<LPParticle> Particle_Set, std::vector<double> w, double lambda, double xi, string *pt_maxstr){
  // cout << "START MERGE" << endl;
  LPPartition partition;
  // double* parameter;
  if(A_or_B){
    partition = partition_A;
    // parameter = alpha;
  } else {
    partition = partition_B;
    // parameter = beta;
  }
  
  arma::mat tmp_A;
  bool adj_clusters;
  double tmp_objective = 0.0; 
  LPParticle tmp_candidate = new Particle(this);

  if(partition->K > 1){ // only makes sense to merge clusters if there are at least 2
    for(int k = 0; k < partition->K - 1 ; k++){
      for(int kk = k+1; kk < partition->K; kk++){
        // cout << "Start Clusters " << k << " and " << kk << endl;
        tmp_A = Submatrix(A_block, partition->cluster_config[k], partition->cluster_config[kk], partition->clusters[k], partition->clusters[kk]);
        // check if any element of tmp_A is equal to
        adj_clusters = any(vectorise(tmp_A == 1.0));
        if(adj_clusters){ // propose merging clusters k and kk!
          delete tmp_candidate;
          tmp_candidate = new Particle(this);
          
          // cout << "before merge" << endl;
          tmp_candidate->Partition_Merge(A_or_B,k,kk);
          // cout << "after merge" << endl;
          tmp_objective = w[current_l]*tmp_candidate->Total_Log_Post() + lambda * Entropy(current_l, tmp_candidate, Particle_Set, w) + xi * VI_Avg(current_l, tmp_candidate, Particle_Set);
          // cout << "after tmp_objective" << endl;
          // cout << tmp_objective << endl;
          if(tmp_objective > *pointer_to_maxobjective){
            delete max_candidate;
            try{
              max_candidate = new Particle(tmp_candidate);
            }
            catch(const std::bad_alloc& e){
              cout << "######## EXCEPTION 8 ########"  << e.what() << endl;
            }
            *pointer_to_maxobjective = tmp_objective;
            *pt_maxstr = "Merge";
          } // closes if statement that updates max_candidate
        } // closes if statement that proposes the merge
        // cout << "Start Clusters " << k << " and " << kk << endl;
      } // closes inner loop over remaining clusters
    } // closes outer loop over clusters
  }
  delete tmp_candidate;
  // cout << "END MERGE" << endl;
  return;
}

void Particle::Spectral_SplitMerge_moves(bool A_or_B, double *pointer_to_maxobjective, int current_l, LPParticle& max_candidate, std::vector<LPParticle> Particle_Set, std::vector<double> w, double lambda, double xi, string *pt_maxstr, bool use_mle){
  LPPartition partition;
  // double* parameter;
  if(A_or_B){
    partition = partition_A;
    // parameter = alpha;
  } else {
    partition = partition_B;
    // parameter = beta;
  }
  double tmp_objective = 0.0; 
  int max_splits, n_cl;
  std::vector<std::vector<int> > indices; // indices contains numbers between 0 and nObs-1
  std::vector<int> ns;
  LPParticle tmp_candidate = new Particle(this); 
  
  for(int split_k = 0; split_k < partition->K; split_k++){
    if(partition->cluster_config[split_k] > 2){
      n_cl = partition->cluster_config[split_k];
      
      if(ceil(split_frac * sqrt(n_cl)) <= 2) max_splits = 3;
      else max_splits = ceil(split_frac * sqrt(n_cl));

      for(int num_splits = 2; num_splits < max_splits; num_splits++){
        if(n_cl > num_splits){
          delete tmp_candidate;
          tmp_candidate = new Particle(this);
          
          tmp_objective = 0.0;
          // now just try splitting
          tmp_candidate->Find_K_SpectralSplits(A_or_B, split_k, num_splits, indices, ns, use_mle);
          tmp_objective = w[current_l]*tmp_candidate->Total_Log_Post() + lambda * Entropy(current_l, tmp_candidate, Particle_Set, w) + xi * VI_Avg(current_l, tmp_candidate, Particle_Set);

          if(tmp_objective > *pointer_to_maxobjective){
            delete max_candidate;
            max_candidate = new Particle(tmp_candidate);
            *pointer_to_maxobjective = tmp_objective;
            *pt_maxstr = "Split";
            if(use_mle){
              *pt_maxstr = "SplitMle";
            }
          }
          // now do the split and merge part
          if(partition->K>1){
            std::vector<int> jstar = get_jstar(A_or_B, split_k, num_splits, indices, ns);
            int to_be_merged;
            for(int j = 0; j < num_splits; j++){
              if(jstar[j] != -1){
                // the new cluster j is now in position 'to_be_merged':
                if(j == 0){
                  to_be_merged = split_k;
                } else {
                  to_be_merged = partition->K + j - 1;
                }
                delete tmp_candidate;
                tmp_candidate = new Particle(this);
                tmp_candidate->Partition_KSplit(A_or_B, split_k, indices, ns);
                tmp_candidate->Partition_Merge(A_or_B, to_be_merged, jstar[j]);
                tmp_objective = w[current_l]*tmp_candidate->Total_Log_Post() + lambda * Entropy(current_l, tmp_candidate, Particle_Set, w) + xi * VI_Avg(current_l, tmp_candidate, Particle_Set);
                if(tmp_objective > *pointer_to_maxobjective){
                  delete max_candidate;
                  max_candidate = new Particle(tmp_candidate);
                  *pointer_to_maxobjective = tmp_objective;
                  *pt_maxstr = "SplitMerge";
                  if(use_mle){
                    *pt_maxstr = "SplitMergeMle";
                  }
                }
              }
            } // close FOR over new clusters
          } // close IF for split-merge
          for(int c = 0; c < (int)indices.size(); c++){
            indices[c].clear();  
          }
          indices.clear();
        } // close IF n_cl > num_splits
      } // close FOR num_split
    } // close IF at least two elements
  } // close FOR for split_k
  delete tmp_candidate;
  return;
}

void Particle::KM_SplitMerge_moves(bool A_or_B, double *pointer_to_maxobjective, int current_l, LPParticle& max_candidate, std::vector<LPParticle> Particle_Set, std::vector<double> w, double lambda, double xi, string *pt_maxstr, bool use_mle){
  LPPartition partition;
  // double* parameter;
  if(A_or_B){
    partition = partition_A;
    // parameter = alpha;
  } else {
    partition = partition_B;
    // parameter = beta;
  }
  double tmp_objective = 0.0; 
  int max_splits, n_cl;
  std::vector<std::vector<int> > indices; // indices contains numbers between 0 and nObs-1
  std::vector<int> ns;
  LPParticle tmp_candidate = new Particle(this);
  
  for(int split_k = 0; split_k < partition->K; split_k++){
    if(partition->cluster_config[split_k] > 2){
      n_cl = partition->cluster_config[split_k];
      
      if(ceil(split_frac * sqrt(n_cl)) <= 2) max_splits = 3;
      else max_splits = ceil(split_frac * sqrt(n_cl));
      
      for(int num_splits = 2; num_splits < max_splits; num_splits++){
        delete tmp_candidate;
        tmp_candidate = new Particle(this);
        
        tmp_objective = 0.0;
        // now just try splitting
        tmp_candidate->Find_KM_Splits(A_or_B, split_k, num_splits, indices, ns, use_mle);
        tmp_objective = w[current_l]*tmp_candidate->Total_Log_Post() + lambda * Entropy(current_l, tmp_candidate, Particle_Set, w) + xi * VI_Avg(current_l, tmp_candidate, Particle_Set);

        if(tmp_objective > *pointer_to_maxobjective){
          delete max_candidate;
          max_candidate = new Particle(tmp_candidate);
          *pointer_to_maxobjective = tmp_objective;
          *pt_maxstr = "KM Split";
          if(use_mle){
            *pt_maxstr = "KM SplitMle";
          }
        }
        // now do the split and merge part 
        
        if(partition->K>1){
          std::vector<int> jstar = get_jstar(A_or_B, split_k, num_splits, indices, ns);
          int to_be_merged;
          for(int j = 0; j < num_splits; j++){
            if(jstar[j] != -1){
              // the new cluster j is now in position 'to_be_merged':
              if(j == 0){
                to_be_merged = split_k;
              } else {
                to_be_merged = partition->K + j - 1;
              }
              delete tmp_candidate;
              tmp_candidate = new Particle(this);
              tmp_candidate->Partition_KSplit(A_or_B, split_k, indices, ns);
              tmp_candidate->Partition_Merge(A_or_B, to_be_merged, jstar[j]);
              tmp_objective = w[current_l]*tmp_candidate->Total_Log_Post() + lambda * Entropy(current_l, tmp_candidate, Particle_Set, w) + xi * VI_Avg(current_l, tmp_candidate, Particle_Set);
              if(tmp_objective > *pointer_to_maxobjective){
                delete max_candidate;
                max_candidate = new Particle(tmp_candidate);
                *pointer_to_maxobjective = tmp_objective;
                *pt_maxstr = "KM SplitMerge";
                if(use_mle){
                  *pt_maxstr = "KM SplitMergeMle";
                }
              }
            }
          } // close FOR over new clusters
        } // close IF for split-merge
        
        for(int c = 0; c < (int)indices.size(); c++){
          indices[c].clear();  
        }
        indices.clear();
      } // close FOR num_split
    } // close IF at least two elements
  } // close FOR for split_k
  delete tmp_candidate;
  return;
}

void Particle::Tail_Split_moves(bool A_or_B, double *pointer_to_maxobjective, int current_l, LPParticle& max_candidate, std::vector<LPParticle> Particle_Set, std::vector<double> w, double lambda, double xi, string *pt_maxstr){
  LPPartition partition;
  // double* parameter;
  if(A_or_B){
    partition = partition_A;
    // parameter = alpha;
  } else {
    partition = partition_B;
    // parameter = beta;
  }
  double tmp_objective = 0.0; 
  
  std::vector<std::vector<int> > indices;
  std::vector<int> ns;

  std::vector<std::vector<int> > left_center_right;
  int to_be_merged;

  // to determine the connected components
  std::vector<std::vector<int> > left_conncomp; // holds the connected components of sub_left_tail
  std::vector<std::vector<int> > right_conncomp; // holds the connected components of sub_left_tail
  std::vector<int> left_k_star;
  std::vector<int> right_k_star;
  LPParticle tmp_candidate = new Particle(this);

  for(int split_k = 0; split_k < partition->K; split_k++){
    if(partition->cluster_config[split_k] > 2){ 
      
      get_leftcenterright(A_or_B, left_center_right, split_k);

      // now try to remove just the left-tail
      left_conncomp.clear();
      left_conncomp.push_back(std::vector<int>(1, left_center_right[0][0]));
      for(int i = 1; i < (int)left_center_right[0].size(); i++){
        delete tmp_candidate;
        tmp_candidate = new Particle(this);
        tmp_candidate->Find_TailSplit(A_or_B, 0, i, split_k, left_center_right, left_conncomp, indices, ns);
        tmp_objective = w[current_l]*tmp_candidate->Total_Log_Post() + lambda * Entropy(current_l, tmp_candidate, Particle_Set, w) + xi * VI_Avg(current_l, tmp_candidate, Particle_Set);
        if(tmp_objective > *pointer_to_maxobjective){
          delete max_candidate;
          max_candidate = new Particle(tmp_candidate);
          *pointer_to_maxobjective = tmp_objective;
          *pt_maxstr = "left TailSplit";
        }

        // The merge part only merges one of the splits to its best neighboring cluster.
        if(partition->K > 1){
          left_k_star = get_tail_jstar(A_or_B, split_k, left_conncomp);
          for(int j = 0; j < (int)left_conncomp.size(); j++){
            if(left_k_star[j] != -1){
              if(j == 0){
                to_be_merged = split_k;
              } else {
                to_be_merged = partition->K + j - 1; // because indices contains left_conncomp first
              }
              delete tmp_candidate;
              tmp_candidate = new Particle(this);
              tmp_candidate->Partition_KSplit(A_or_B, split_k, indices, ns);
              tmp_candidate->Partition_Merge(A_or_B, to_be_merged, left_k_star[j]);
              tmp_objective = w[current_l]*tmp_candidate->Total_Log_Post() + lambda * Entropy(current_l, tmp_candidate, Particle_Set, w) + xi * VI_Avg(current_l, tmp_candidate, Particle_Set);
              if(tmp_objective > *pointer_to_maxobjective){
                delete max_candidate;
                max_candidate = new Particle(tmp_candidate);
                *pointer_to_maxobjective = tmp_objective;
                *pt_maxstr = "left TailSplitMerge";
              }
            }
          }
        }

        for(int c = 0; c < (int)indices.size(); c++){
          indices[c].clear();  
        }
        indices.clear();
        ns.clear();
      } // end left-tail
      // now try to remove just the right tail
      right_conncomp.clear();
      right_conncomp.push_back(std::vector<int>(1, left_center_right[2][0])); 
      for(int i = 1; i < (int)left_center_right[2].size(); i++){
        delete tmp_candidate;
        tmp_candidate = new Particle(this);
        tmp_candidate->Find_TailSplit(A_or_B, 2, i, split_k, left_center_right, right_conncomp, indices, ns);
        tmp_objective = w[current_l]*tmp_candidate->Total_Log_Post() + lambda * Entropy(current_l, tmp_candidate, Particle_Set, w) + xi * VI_Avg(current_l, tmp_candidate, Particle_Set);
        if(tmp_objective > *pointer_to_maxobjective){
          delete max_candidate;
          max_candidate = new Particle(tmp_candidate);
          *pointer_to_maxobjective = tmp_objective;
          *pt_maxstr = "right TailSplit";
        }

        if(partition->K > 1){
          right_k_star = get_tail_jstar(A_or_B, split_k, right_conncomp);
          for(int j = 0; j < (int)right_conncomp.size(); j++){
            if(right_k_star[j] != -1){
              if(j == 0){
                to_be_merged = split_k;
              } else {
                to_be_merged = partition->K + j - 1; // because indices contains right_conncomp first
              }
              delete tmp_candidate;
              tmp_candidate = new Particle(this);
              tmp_candidate->Partition_KSplit(A_or_B, split_k, indices, ns);
              tmp_candidate->Partition_Merge(A_or_B, to_be_merged, right_k_star[j]);
              tmp_objective = w[current_l]*tmp_candidate->Total_Log_Post() + lambda * Entropy(current_l, tmp_candidate, Particle_Set, w) + xi * VI_Avg(current_l, tmp_candidate, Particle_Set);
              if(tmp_objective > *pointer_to_maxobjective){
                delete max_candidate;
                max_candidate = new Particle(tmp_candidate);
                *pointer_to_maxobjective = tmp_objective;
                *pt_maxstr = "right TailSplitMerge";
              }
            }
          }
        }
      }
    }
  }
  delete tmp_candidate;
  return;
}

void Particle::Update_Particle(int current_l, std::vector<LPParticle> Particle_Set, std::vector<double> w, double lambda, double xi, bool restricted){
  if(A_or_B_first == 1.0){
    if(sampleA){
      Update_Partition(1, current_l, Particle_Set, w, lambda, xi, restricted);  
    }
    if(sampleB){
      Update_Partition(0, current_l, Particle_Set, w, lambda, xi, restricted);  
    }
  } else if(A_or_B_first == 0.0){
    if(sampleB){
      Update_Partition(0, current_l, Particle_Set, w, lambda, xi, restricted);  
    }
    if(sampleA){
      Update_Partition(1, current_l, Particle_Set, w, lambda, xi, restricted);  
    }
  } else if(A_or_B_first == 0.5){
    if(R::runif(0,1) < 0.5){
      if(sampleA){
        Update_Partition(1, current_l, Particle_Set, w, lambda, xi, restricted);  
      }
      if(sampleB){
        Update_Partition(0, current_l, Particle_Set, w, lambda, xi, restricted);  
      }
    } else {
      if(sampleB){
        Update_Partition(0, current_l, Particle_Set, w, lambda, xi, restricted);  
      }
      if(sampleA){
        Update_Partition(1, current_l, Particle_Set, w, lambda, xi, restricted);  
      }
    }
  }
  // alpha and beta are updated when we update the partition
}

double Particle::Total_Log_Post(){
  double tot = 0;
  if(opt_Y == 1){
    tot += LikelihoodY();
  } else if(opt_Y == 0){
    tot += LikelihoodYAB();
  } else if(opt_Y == 2){
    tot += total_log_like();
  } else if(opt_Y == 3){
    tot += LikelihoodY_MCMC_poisson();
  } else if(opt_Y == 4){
    tot += LikelihoodY_Laplace_Pois();
  }

  // you need the prior for Palpha and Pbeta
  tot += total_log_prior(partition_A);
  tot += total_log_prior(partition_B);
  return tot;
}

double Particle::LikelihoodYAB(){
  int T = Y.n_cols;
  rowvec x = X.row(0);
  double dof = nObs*T/2 + nObs; // sarebbe: + nObs/2 + nObs/2
  
  vec Y_vec = zeros<vec>(nObs*T);
  for(int i = 0; i < nObs; i++){
    for(int t = 0; t < T; t++){
      Y_vec(i*T + t) = Y(i,t) - alpha[i] - x(t)*beta[i];
    }
  }
  double quad_form = as_scalar(Y_vec.t() * Y_vec);
  double log_determ = 0.0;
  for(int k = 0; k < partition_A->K; k++){
    int cluster_size = partition_A->cluster_config[k];
    if(cluster_size == 1){
      double param = alpha[ partition_A->clusters[k][0] ];
      quad_form += param*param/(a1/(1-rho)+a2);
      log_determ += -log(a1/(1-rho)+a2);
    } else {
      mat O = ones<mat>(cluster_size, cluster_size);
      vec param_vec = zeros<vec>(cluster_size);
      for(int i = 0; i < cluster_size; i++){
        param_vec(i) = alpha[ partition_A->clusters[k][i] ];
      }
      mat A_block_k = Submatrix(A_block, cluster_size, cluster_size, partition_A->clusters[k], partition_A->clusters[k]);
      vec row_sums = zeros<vec>(cluster_size);
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          row_sums(i) += A_block_k(i,j);
        }
      }
      mat D = diagmat(row_sums);
      mat A_star_k = D - A_block_k;
      mat M = rho * A_star_k;
      for(int i = 0; i < cluster_size; i++){
        M(i,i) += 1 - rho;
      }
      mat Sigma0inv = (M/a1) - (M/a1) * O * (M/a1) / (sum(sum(M))/a1 + 1/a2);
      double Sigma0inv_det = 0;
      double Sigma0inv_det_sgn = 0;
      log_det(Sigma0inv_det, Sigma0inv_det_sgn, Sigma0inv);
      log_determ += Sigma0inv_det;
      quad_form += as_scalar(param_vec.t() * Sigma0inv * param_vec);
    }
  }
  for(int k = 0; k < partition_B->K; k++){
    int cluster_size = partition_B->cluster_config[k];
    if(cluster_size == 1){
      double param = beta[ partition_B->clusters[k][0] ];
      quad_form += param*param/(b1/(1-rho)+b2);
      log_determ += -log(b1/(1-rho)+b2);
    } else {
      mat O = ones<mat>(cluster_size, cluster_size);
      vec param_vec = zeros<vec>(cluster_size);
      for(int i = 0; i < cluster_size; i++){
        param_vec(i) = beta[ partition_B->clusters[k][i] ];
      }
      mat A_block_k = Submatrix(A_block, cluster_size, cluster_size, partition_B->clusters[k], partition_B->clusters[k]);
      vec row_sums = zeros<vec>(cluster_size);
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          row_sums(i) += A_block_k(i,j);
        }
      }
      mat D = diagmat(row_sums);
      mat A_star_k = D - A_block_k;
      mat M = rho * A_star_k;
      for(int i = 0; i < cluster_size; i++){
        M(i,i) += 1 - rho;
      }
      mat Sigma0inv = (M/b1) - (M/b1) * O * (M/b1) / (sum(sum(M))/b1 + 1/b2);
      
      double Sigma0inv_det = 0;
      double Sigma0inv_det_sgn = 0;
      log_det(Sigma0inv_det, Sigma0inv_det_sgn, Sigma0inv);
      log_determ += Sigma0inv_det;
      quad_form += as_scalar(param_vec.t() * Sigma0inv * param_vec);
    }
  }
  double numerator = lgamma(alpha_sigma + dof) + alpha_sigma * log(beta_sigma);
  double denominator = lgamma(alpha_sigma) + (alpha_sigma + dof) * log(beta_sigma + quad_form/2);
  double log_likeYAB = 0.5 * log_determ + numerator - denominator;
  return log_likeYAB;
}

double Particle::LikelihoodY(){
  int T = Y.n_cols;
  rowvec x = X.row(0);
  rowvec x2 = pow(x,2);
  mat E = eye<mat>(nObs,nObs);
  mat O_k, A_block_k, Omega_k, D, M;
  mat Omega_a = zeros<mat>(nObs,nObs);
  mat Omega_b = zeros<mat>(nObs,nObs);
  // I am updating Omega_a (the inverse of the covariance matrix of alpha) cluster-wise
  for(int k = 0; k < partition_A->K; k++){
    int cluster_size = partition_A->cluster_config[k]; 
    if(cluster_size == 1){
      Omega_a(partition_A->clusters[k][0], partition_A->clusters[k][0]) = 1/(a1/(1-rho)+a2);
    } else {
      O_k = ones<mat>(cluster_size,cluster_size);
      // Creating M: the precision matrix of the CAR model
      A_block_k = Submatrix(A_block, cluster_size, cluster_size, partition_A->clusters[k], partition_A->clusters[k]);
      vec row_sums = zeros<vec>(cluster_size);
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          row_sums(i) += A_block_k(i,j);
        }
      }
      D = diagmat(row_sums);
      M = rho * (D - A_block_k);
      for(int i = 0; i < cluster_size; i++){
        M(i,i) += 1 - rho;
      }
      // This is the Woodbury matrix identity for inversion
      // Omega_k = (M/a1) - (M/a1) * O_k * (M/a1) / (sum(sum(M))/a1 + 1/a2);
      // Woodbury matrix identity + remember that M*1 = (1-rho)*1
      Omega_k = (M/a1) - ((1-rho)/a1)*((1-rho)/a1) / ( cluster_size*(1-rho)/a1 + 1/a2) * O_k;
      // this is formula (S1) [Sept 10, 2021]
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          Omega_a(partition_A->clusters[k][i], partition_A->clusters[k][j]) = Omega_k(i,j);
        }
      }
    }
  }
  // I am updating Omega_b (the inverse of the covariance matrix of beta) cluster-wise
  for(int k = 0; k < partition_B->K; k++){
    int cluster_size = partition_B->cluster_config[k];   
    if(cluster_size == 1){
      Omega_b(partition_B->clusters[k][0], partition_B->clusters[k][0]) = 1/(b1/(1-rho)+b2);
    } else { 
      O_k = ones<mat>(cluster_size,cluster_size);
      // Creating M: the precision matrix of the CAR model
      A_block_k = Submatrix(A_block, cluster_size, cluster_size, partition_B->clusters[k], partition_B->clusters[k]);
      vec row_sums = zeros<vec>(cluster_size);
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          row_sums(i) += A_block_k(i,j);
        }
      }
      D = diagmat(row_sums);
      M = rho * (D - A_block_k);
      for(int i = 0; i < cluster_size; i++){
        M(i,i) += 1 - rho;
      }
      // This is the Woodbury matrix identity for inversion
      // Omega_k = (M/b1) - (M/b1) * O_k * (M/b1) / (sum(sum(M))/b1 + 1/b2);
      // Woodbury matrix identity + remember that M*1 = (1-rho)*1
      Omega_k = (M/b1) - ((1-rho)/b1)*((1-rho)/b1) / ( cluster_size*(1-rho)/b1 + 1/b2) * O_k;
      // this is formula (S1) [Sept 10, 2021]
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          Omega_b(partition_B->clusters[k][i], partition_B->clusters[k][j]) = Omega_k(i,j);
        }
      }
    }
  }
  mat P, P11, P12, P21, P22;
  if(opt_method == 1){
    P11 = *I11;
    P12 = *I12;
    P21 = *I21;
    P22 = *I22;
    P = arma::join_cols(arma::join_rows(P11, P12),arma::join_rows(P21, P22));
  } else {
    P = block_inverse_ret(E*T + Omega_a, E*sum(x), E*sum(x), E*sum(x2)+Omega_b, &P11, &P12, &P21, &P22);
    // this is formula the parenthesis in (S4) [Sept 10, 2021]
  }
  arma::vec xty1 = sum(Y,1);
  arma::vec xty2 = sum(X%Y,1);
  arma::vec XtY = arma::join_cols(xty1,xty2);

  double quad_form = accu(Y%Y);
  quad_form += - as_scalar(XtY.t() * P * XtY);
  // this is formula (S4) [Sept 10, 2021]

  double Omega_log_det = 0;
  double Omega_log_det_sgn = 0;
  // Using Sylvester's + Woodbury determinant theorem: det(I + X Sigma Xt)^-1 = det(I - X P Xt) = det(I - P XtX) - so we need to multiply for 0.5
  arma::mat Sigma_det = arma::join_cols(arma::join_rows(E - (T*P11 + sum(x)*P12), -(sum(x)*P11+sum(x2)*P12)),arma::join_rows(-(T*P21 + sum(x)*P22), E-(sum(x)*P21+sum(x2)*P22)));
  // this is formula (S6) [Sept 10, 2021]
  arma::log_det(Omega_log_det, Omega_log_det_sgn, Sigma_det);
  // potentially you can try the formula for block matrices
  // Omega_log_det = block_log_det(E - (T*P11 + sum(x)*P12), -(sum(x)*P11+sum(x2)*P12), -(T*P21 + sum(x)*P22), E-(sum(x)*P21+sum(x2)*P22));


  double numerator = lgamma(alpha_sigma +  ( (double) nObs * T/2)) + alpha_sigma * log(beta_sigma);
  double denominator = lgamma(alpha_sigma) + (alpha_sigma + ( (double)nObs * T/2)) * log(beta_sigma + quad_form/2);
  return 0.5 * Omega_log_det - ( (double)nObs * T/2)*log(2*M_PI)+ numerator - denominator;
  // this is formula (S3) [Sept 10, 2021]
}

double Particle::LikelihoodY_MCMC_poisson(){
  int n_MCMC = 100;
  int n_MCMC_batch = 100; // faster if smaller n_MCMC_batch
  int T= Y.n_cols;
  rowvec x = X.row(0);
  mat O_k, A_block_k, Sigma_k, R_Sigma, D, Omega_CAR;
  vec row_sums;
  arma::vec log_post_MCMC(n_MCMC);
  arma::vec sums(n_MCMC_batch);

  // sample alpha and beta
  arma::mat alpha(nObs, n_MCMC);
  arma::mat beta(nObs, n_MCMC);
  arma::vec sigma2 = ones<vec>(n_MCMC);
  mat RSigma_a = zeros<mat>(nObs,nObs);
  mat RSigma_b = zeros<mat>(nObs,nObs);
  
  // we need the whole alpha/beta vector before computing tmp and log_post_MCMC
  // so we cannot break down the likelihood calculation (to speed up)
  for(int k = 0; k < partition_A->K; k++){
    int cluster_size = partition_A->cluster_config[k]; 
    if(cluster_size == 1){
      RSigma_a(partition_A->clusters[k][0], partition_A->clusters[k][0]) = sqrt(a1/(1-rho)+a2);
    } else {
      O_k = ones<mat>(cluster_size,cluster_size);
      // Creating Omega_CAR: the precision matrix of the CAR model
      A_block_k = Submatrix(A_block, cluster_size, cluster_size, partition_A->clusters[k], partition_A->clusters[k]);
      row_sums = zeros<vec>(cluster_size);
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          row_sums(i) += A_block_k(i,j);
        }
      }
      D = diagmat(row_sums);
      Omega_CAR = rho * (D - A_block_k);
      for(int i = 0; i < cluster_size; i++){
        Omega_CAR(i,i) += 1 - rho;
      }
      
      Sigma_k = inv_sympd(Omega_CAR) * a1 + O_k * a2;
      R_Sigma = chol( Sigma_k, "lower" ); // with layout = "lower": X = R * R.t()

      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          RSigma_a(partition_A->clusters[k][i], partition_A->clusters[k][j]) = R_Sigma(i,j);
        }
      }

    }
  }
  for(int k = 0; k < partition_B->K; k++){
    int cluster_size = partition_B->cluster_config[k]; 
    if(cluster_size == 1){
      RSigma_b(partition_B->clusters[k][0], partition_B->clusters[k][0]) = sqrt(b1/(1-rho)+b2);
    } else {
      O_k = ones<mat>(cluster_size,cluster_size);
      // Creating Omega_CAR: the precision matrix of the CAR model
      A_block_k = Submatrix(A_block, cluster_size, cluster_size, partition_B->clusters[k], partition_B->clusters[k]);
      row_sums = zeros<vec>(cluster_size);
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          row_sums(i) += A_block_k(i,j);
        }
      }
      D = diagmat(row_sums);
      Omega_CAR = rho * (D - A_block_k);
      for(int i = 0; i < cluster_size; i++){
        Omega_CAR(i,i) += 1 - rho;
      }
      
      Sigma_k = inv_sympd(Omega_CAR) * b1 + O_k * b2;
      R_Sigma = chol( Sigma_k, "lower" ); // with layout = "lower": X = R * R.t()

      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          RSigma_b(partition_B->clusters[k][i], partition_B->clusters[k][j]) = R_Sigma(i,j);
        }
      }

    }
  }

  for(int batch = 0; batch < n_MCMC_batch; batch++){
    sigma2 = 1/randg<vec>(n_MCMC, distr_param(alpha_sigma,1/beta_sigma)); 
    alpha = RSigma_a * randn<mat>(nObs, n_MCMC);
    beta = RSigma_b * randn<mat>(nObs, n_MCMC);
    alpha = alpha.each_row() % sqrt(sigma2).t();
    beta = beta.each_row() % sqrt(sigma2).t();
    arma::rowvec tmp = zeros<rowvec>(n_MCMC);
    for(int t = 0; t < T; t++){
      // remember that lambda = exp(a + b*x)    &     Y.col(t).t() * (...) is a matrix product
      // sum(M, 0) means summing elements in columns (colSums)
      tmp += -sum(exp(alpha + beta * x[t]), 0) + ( Y.col(t).t() * (alpha + beta * x[t]) );  
    }
    for(int m = 0;m < n_MCMC; m++){
      for(int t = 0; t < T; t++){
        // tmp(m) += - sum(exp(alpha.col(m) + beta.col(m) * x[t])) + sum(Y.col(t).t() * (alpha.col(m) + beta.col(m) * x[t]));
        tmp(m) += - sum(exp(alpha.col(m) + beta.col(m) * x[t])) + dot(Y.col(t), (alpha.col(m) + beta.col(m) * x[t]));
      }
    }
    log_post_MCMC = tmp.t(); // we don't normalize this
    sums(batch) = logsumexp(log_post_MCMC);
  }
  double log_like = logsumexp(sums) - log(n_MCMC_batch*n_MCMC);
  return log_like;
}


double Particle::LikelihoodY_MCMC_normal(){
  int n_MCMC = 10000;
  int n_MCMC_batch = 100; // faster if smaller n_MCMC_batch
  int T= Y.n_cols;
  rowvec x = X.row(0);
  rowvec x2 = pow(x,2);
  mat O_k, A_block_k, Sigma_k, R_Sigma, D, Omega_CAR;
  vec row_sums;
  arma::vec log_post_MCMC(n_MCMC);
  arma::vec sums(n_MCMC_batch);
  

  // sample alpha and beta
  arma::mat alpha(nObs, n_MCMC);
  arma::mat beta(nObs, n_MCMC);
  arma::vec sigma2 = ones<vec>(n_MCMC);
  mat RSigma_a = zeros<mat>(nObs,nObs);
  mat RSigma_b = zeros<mat>(nObs,nObs);
  
  // we need the whole alpha/beta vector before computing tmp and log_post_MCMC
  // so we cannot break down the likelihood calculation (to speed up)
  for(int k = 0; k < partition_A->K; k++){
    int cluster_size = partition_A->cluster_config[k]; 
    if(cluster_size == 1){
      RSigma_a(partition_A->clusters[k][0], partition_A->clusters[k][0]) = sqrt(a1/(1-rho)+a2);
    } else {
      O_k = ones<mat>(cluster_size,cluster_size);
      // Creating Omega_CAR: the precision matrix of the CAR model
      A_block_k = Submatrix(A_block, cluster_size, cluster_size, partition_A->clusters[k], partition_A->clusters[k]);
      row_sums = zeros<vec>(cluster_size);
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          row_sums(i) += A_block_k(i,j);
        }
      }
      D = diagmat(row_sums);
      Omega_CAR = rho * (D - A_block_k);
      for(int i = 0; i < cluster_size; i++){
        Omega_CAR(i,i) += 1 - rho;
      }
      
      Sigma_k = inv_sympd(Omega_CAR) * a1 + O_k * a2;
      R_Sigma = chol( Sigma_k, "lower" ); // with layout = "lower": X = R * R.t()

      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          RSigma_a(partition_A->clusters[k][i], partition_A->clusters[k][j]) = R_Sigma(i,j);
        }
      }

    }
  }
  for(int k = 0; k < partition_B->K; k++){
    int cluster_size = partition_B->cluster_config[k]; 
    if(cluster_size == 1){
      RSigma_b(partition_B->clusters[k][0], partition_B->clusters[k][0]) = sqrt(b1/(1-rho)+b2);
    } else {
      O_k = ones<mat>(cluster_size,cluster_size);
      // Creating Omega_CAR: the precision matrix of the CAR model
      A_block_k = Submatrix(A_block, cluster_size, cluster_size, partition_B->clusters[k], partition_B->clusters[k]);
      row_sums = zeros<vec>(cluster_size);
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          row_sums(i) += A_block_k(i,j);
        }
      }
      D = diagmat(row_sums);
      Omega_CAR = rho * (D - A_block_k);
      for(int i = 0; i < cluster_size; i++){
        Omega_CAR(i,i) += 1 - rho;
      }
      
      Sigma_k = inv_sympd(Omega_CAR) * b1 + O_k * b2;
      R_Sigma = chol( Sigma_k, "lower" ); // with layout = "lower": X = R * R.t()

      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          RSigma_b(partition_B->clusters[k][i], partition_B->clusters[k][j]) = R_Sigma(i,j);
        }
      }

    }
  }

  for(int batch = 0; batch < n_MCMC_batch; batch++){
    // cout << batch << " ";
    sigma2 = 1/randg<vec>(n_MCMC, distr_param(alpha_sigma,1/beta_sigma)); 
    alpha = RSigma_a * randn<mat>(nObs, n_MCMC);
    beta = RSigma_b * randn<mat>(nObs, n_MCMC);
    alpha = alpha.each_row() % sqrt(sigma2).t();
    beta = beta.each_row() % sqrt(sigma2).t();
    arma::rowvec tmp = zeros<rowvec>(n_MCMC);
    for(int m = 0;m < n_MCMC; m++){
      for(int t = 0; t < T; t++){
        tmp(m) += - sum(  pow(Y.col(t) - alpha.col(m) - beta.col(m) * x[t],2) )/(2*sigma2(m));
      }
    }
    log_post_MCMC = tmp.t() - log(2*M_PI*sigma2)*nObs*T/2;  // alternative to M_PI: datum::pi
    sums(batch) = logsumexp(log_post_MCMC);
  }
  double log_like = logsumexp(sums) - log(n_MCMC_batch*n_MCMC);
  return log_like;
}

arma::mat Particle::get_Omega(bool A_or_B){
  mat E = eye<mat>(nObs,nObs);
  mat O_k, A_block_k, Omega_k, D, M;
  mat Omega = zeros<mat>(nObs,nObs);

  LPPartition partition;
  double c1, c2;
  if(A_or_B){
    partition = partition_A;
    c1 = a1;
    c2 = a2;
  } else {
    partition = partition_B;
    c1 = b1;
    c2 = b2;
  }
  
  for(int k = 0; k < partition->K; k++){
    int cluster_size = partition->cluster_config[k]; 
    if(cluster_size == 1){
      Omega(partition->clusters[k][0], partition->clusters[k][0]) = 1.0/(c1/(1.0-rho)+c2);
    } else {
      O_k = ones<mat>(cluster_size,cluster_size);
      // Creating M: the precision matrix of the CAR model
      A_block_k = Submatrix(A_block, cluster_size, cluster_size, partition->clusters[k], partition->clusters[k]);
      vec row_sums = zeros<vec>(cluster_size);
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          row_sums(i) += A_block_k(i,j);
        }
      }
      D = diagmat(row_sums);
      M = rho * (D - A_block_k);
      for(int i = 0; i < cluster_size; i++){
        M(i,i) += 1 - rho;
      }
      // This is the Woodbury matrix identity for inversion
      // Omega_k = (M/c1) - (M/c1) * O_k * (M/c1) / (sum(sum(M))/c1 + 1/c2);
      // Woodbury matrix identity + remember that M*1 = (1-rho)*1
      Omega_k = (M/c1) - ((1-rho)/c1)*((1-rho)/c1) / ( cluster_size*(1-rho)/c1 + 1/c2) * O_k;
      // this is formula (S1) [Sept 10, 2021]
      for(int i = 0; i < cluster_size; i++){
        for(int j = 0; j < cluster_size; j++){
          Omega(partition->clusters[k][i], partition->clusters[k][j]) = Omega_k(i,j);
        }
      }
    }
  }
  return Omega;
}

// // check http://dlib.net/matrix_ex.cpp.html
// double Particle::fun_like_negbin(dlib::matrix<double,0,1> alpha_beta){
//   // NOT USED AT THE MOMENT.
//   double r = 1.2; // todo
//   int nObs = Y.n_rows;
//   arma::vec alpha(nObs);
//   arma::vec beta(nObs);
//   rowvec x = X.row(0);

//   for(int i = 0; i < nObs; i++){
//     alpha[i] = alpha_beta(i);
//     beta[i] = alpha_beta(i+nObs);
//   }
//   arma::mat Xgamma = alpha * (x*0 + 1) + beta * x; // this should be a matrix (column x row)
//   arma::mat tmp = Xgamma - log(r + exp(Xgamma));

//   double loglik = accu(Y % tmp); // sum(element-wise multiplication)
//   double quad_form = 0.0;

//   mat Omega_a = get_Omega(true);
//   mat Omega_b = get_Omega(false);
  
//   quad_form += as_scalar(alpha.t() * Omega_a * alpha);
//   quad_form += as_scalar(beta.t() * Omega_b * beta);
//   return loglik - quad_form/2;
// }

double Particle::fun_like_pois(dlib::matrix<double,0,1> alpha_beta, bool mle_bool){
  // cout << "fun ";
  // for(int i = 0; i < nObs; i++){
  //   cout << alpha_beta(i) << " ";
  // }
  // cout << endl;
  // for(int i = 0; i < nObs; i++){
  //   cout << alpha_beta(nObs+i) << " ";
  // }
  // cout << endl;
  // //////

  int nObs = Y.n_rows;
  arma::vec alpha(nObs);
  arma::vec beta(nObs);
  rowvec x = X.row(0);

  for(int i = 0; i < nObs; i++){
    alpha[i] = alpha_beta(i);
    beta[i] = alpha_beta(i+nObs);
  }
  arma::mat Xgamma = alpha * (x*0 + 1) + beta * x; // this should be a matrix (column x row)
  double loglik = accu(-exp(Xgamma) + Y % Xgamma); // sum(element-wise multiplication)
  // if(!isfinite(loglik)){
  //   return nan(""); // did not change anything
  // }
  double quad_form = 0.0;

  if(!mle_bool){
    arma::mat Omega_a = get_Omega(true);
    arma::mat Omega_b = get_Omega(false);
    
    // log determinants not needed because the partition is fixed here (and we optimize for alpha and beta)
    quad_form += as_scalar(alpha.t() * Omega_a * alpha); 
    quad_form += as_scalar(beta.t() * Omega_b * beta);
  }
  ///
  // cout << loglik - quad_form/2 << " - ";
  return loglik - quad_form/(2*sigma2_laplace);
}

double Particle::fun_like_pois_singleunit(dlib::matrix<double,0,1> alpha_beta, int i){
  double alpha, beta;
  rowvec x = X.row(0);
  rowvec y = Y.row(i);

  alpha = alpha_beta(0);
  beta = alpha_beta(1);
  arma::rowvec Xgamma = alpha * (x*0 + 1) + beta * x; 
  double loglik = accu(-exp(Xgamma) + y % Xgamma); // sum(element-wise multiplication)
  
  return loglik;
}

dlib::matrix<double,0,1> Particle::fun_deriv_like_pois(dlib::matrix<double,0,1> alpha_beta, bool mle_bool){
  // cout << "_der " << alpha_beta(0) << " ";
  int nObs = Y.n_rows;
  arma::vec alpha(nObs);
  arma::vec beta(nObs);
  rowvec x = X.row(0);

  for(int i = 0; i < nObs; i++){
    alpha[i] = alpha_beta(i);
    beta[i] = alpha_beta(i+nObs);
  }
  arma::mat Xgamma = alpha * (x*0 + 1) + beta * x; // this should be a matrix (column x row)
  
  arma::vec dalpha = arma::sum(-exp(Xgamma) + Y,1);
  arma::vec dbeta = arma::sum((-exp(Xgamma) + Y) % X,1);
  
  if(!mle_bool){
    arma::mat Omega_a = get_Omega(true);
    arma::mat Omega_b = get_Omega(false);
    
    dalpha = dalpha - Omega_a * alpha / sigma2_laplace;
    dbeta = dbeta - Omega_b * beta / sigma2_laplace;
  }
  
  dlib::matrix<double,0,1> dgamma;
  dgamma.set_size(2*nObs);
  for(int i = 0; i < nObs; i++){
    dgamma(i) = dalpha[i];
    dgamma(nObs+i) = dbeta[i];
  }
  // cout << dgamma(0) << " - ";
  return dgamma;
}

dlib::matrix<double,0,1> Particle::fun_deriv_like_pois_singleunit(dlib::matrix<double,0,1> alpha_beta, int i){
  // cout << "deriv " << alpha_beta(0) << " " << alpha_beta(1) << " - ";
  double alpha, beta;
  rowvec x = X.row(0);
  rowvec y = Y.row(i);
  
  alpha = alpha_beta(0);
  beta = alpha_beta(1);

  arma::rowvec Xgamma = alpha * (x*0 + 1) + beta * x; 

  dlib::matrix<double,0,1> dgamma(2);
  dgamma(0) = accu(-exp(Xgamma) + y);
  dgamma(1) = accu((-exp(Xgamma) + y) % x);

  return dgamma;
}

dlib::matrix<double> Particle::fun_hess_like_pois_singleunit(dlib::matrix<double,0,1> alpha_beta){
  // cout << "hess " << alpha_beta(0) << " ";
  rowvec x = X.row(0);
  
  dlib::matrix<double> Hess;
  Hess.set_size(2,2);

  arma::rowvec Xgamma = alpha_beta(0) * (x*0 + 1) + alpha_beta(1) * x; 

  Hess(0,0) = accu(-exp(Xgamma));
  Hess(1,0) = accu(-exp(Xgamma) % x);
  Hess(0,1) = Hess(1,0);
  Hess(1,1) = accu(-exp(Xgamma) % x % x);
  return Hess;
}

dlib::matrix<double> Particle::fun_hess_like_pois(dlib::matrix<double,0,1> alpha_beta, bool mle_bool){
  // cout << "!hess " << alpha_beta(0) << " ";
  int nObs = Y.n_rows;

  dlib::matrix<double,0,1> alpha; // column vector
  dlib::matrix<double,0,1> beta;
  alpha.set_size(nObs);
  beta.set_size(nObs);
  // to subset alpha_beta:
  alpha = dlib::rowm(alpha_beta,dlib::range(0,nObs-1));
  beta = dlib::rowm(alpha_beta,dlib::range(nObs,2*nObs-1));

  dlib::matrix<double,1,0> x; // row vector
  x.set_size(Y.n_cols);
  dlib::matrix<double,1,0> xsquared;
  xsquared.set_size(Y.n_cols);
  for(int j = 0; j < Y.n_cols; j++){
    x(j) = X(0,j);
    xsquared(j) = X(0,j)*X(0,j);
  }

  dlib::matrix<double> eXgamma = dlib::exp(alpha * (0*x + 1) + beta * x); // this should be a matrix (column x row)
  dlib::matrix<double> tmp2 = -eXgamma * trans(x*0+1); // DO NOT USED dlib::sum_rows(-eXgamma); -- returns the sum of the COLUMNS!!
  dlib::matrix<double> tmp3 = -eXgamma * trans(x);
  dlib::matrix<double> tmp4 = -eXgamma * trans(xsquared);

  dlib::matrix<double> Hess;
  Hess.set_size(2*nObs,2*nObs);
  set_all_elements(Hess,0);
  for(int i = 0; i < nObs; i++){
    Hess(i,i) = tmp2(i);
    Hess(i,nObs+i) = tmp3(i);
    Hess(nObs+i,i) = tmp3(i);
    Hess(nObs+i,nObs+i) = tmp4(i);
  }

  if(!mle_bool){
    arma::mat Omega_a = get_Omega(true);
    arma::mat Omega_b = get_Omega(false);
    dlib::matrix<double> Omega_a2,Omega_b2;
    Omega_a2.set_size(nObs,nObs);
    Omega_b2.set_size(nObs,nObs);
    for(int i = 0; i < nObs; i++){
      for(int j = i; j < nObs; j++){
        Omega_a2(i,j) = Omega_a(i,j);
        Omega_a2(j,i) = Omega_a(i,j);
        Omega_b2(i,j) = Omega_b(i,j);
        Omega_b2(j,i) = Omega_b(i,j);
      }
    }
    dlib::set_subm(Hess, dlib::range(0,nObs-1), dlib::range(0,nObs-1))           =   dlib::subm(Hess, dlib::range(0,nObs-1),dlib::range(0,nObs-1))           - Omega_a2 / sigma2_laplace;
    dlib::set_subm(Hess, dlib::range(nObs,2*nObs-1), dlib::range(nObs,2*nObs-1)) =   dlib::subm(Hess, dlib::range(nObs,2*nObs-1),dlib::range(nObs,2*nObs-1)) - Omega_b2 / sigma2_laplace;
  }
  // cout << endl;
  // for(int i = 0; i < 2*nObs; i++){
  //   for(int j = 0; j < 2*nObs; j++){
  //     cout << Hess(i,j) << " ";
  //   }
  //   cout << endl;
  // }
  // cout << endl;
  return Hess;
}

arma::mat Particle::fun_hess_like_pois(arma::vec alpha, arma::vec beta){
  // int nObs = Y.n_rows;
  rowvec x = X.row(0);
  rowvec xsquared = x % x;
  
  arma::mat Xgamma = alpha * (x*0 + 1) + beta * x; // this should be a matrix (column x row)
  arma::vec tmp2 = sum(-exp(Xgamma),1); // sum over the rows
  arma::vec tmp3 = -exp(Xgamma) * x.t();
  arma::vec tmp4 = -exp(Xgamma) * xsquared.t();
  
  arma::mat Hess(2*nObs,2*nObs,fill::zeros);
  for(int i = 0; i < nObs; i++){
    Hess(i,i) = tmp2[i];
    Hess(i,nObs+i) = tmp3[i];
    Hess(nObs+i,i) = tmp3[i];
    Hess(nObs+i,nObs+i) = tmp4[i];
  }

  arma::mat Omega_a = get_Omega(true);
  arma::mat Omega_b = get_Omega(false);
  

  // X.submat( first_row, first_col, last_row, last_col )
  Hess.submat( 0, 0, nObs-1, nObs-1 ) = Hess.submat( 0, 0, nObs-1, nObs-1 ) - Omega_a / sigma2_laplace;
  Hess.submat( nObs, nObs, 2*nObs-1, 2*nObs-1 ) = Hess.submat( nObs, nObs, 2*nObs-1, 2*nObs-1 ) - Omega_b / sigma2_laplace;

  return Hess;
}

// there are two steps and we might want to break them down. 
// the p(gamma, theta \vert Y) is a product of the prior, divided by the normal approximation
// so the laplace is used to get the normal approximation at the denominator


double Particle::LikelihoodY_Laplace_Pois(){
  // you are calling get_Omega within the optimization function all the times!!
  // you should just call it here and pass it as parameter **TODO**
  arma::vec alpha(nObs);
  arma::vec beta(nObs);
  rowvec x = X.row(0);
  double sigma2 = sigma2_laplace;

  // we want this to return the optimal alpha_beta
  dlib::matrix<double,0,1> alpha_beta(2*(Y.n_rows));
  for(int i =0; i < Y.n_rows; i++){
    alpha_beta(i) = 0;
    alpha_beta(i+Y.n_rows) = 0;
  }
  optim(alpha_beta);

  for(int i = 0; i < nObs; i++){
    alpha[i] = alpha_beta(i);
    beta[i] = alpha_beta(i+nObs);
  }
  arma::mat Xgamma = alpha * (x*0 + 1) + beta * x; // this should be a matrix (column x row)
  // this is \log p(Y \vert \theta, \sigma, \gamma)
  double loglik = accu(-exp(Xgamma) + Y % Xgamma); // sum(element-wise multiplication) 
  // now the prior for theta: \log p(\theta \vert \sigma, \gamma)
  arma::mat Omega_a = get_Omega(true);
  arma::mat Omega_b = get_Omega(false);
  double log_det_a, log_det_b, sign_det;
  arma::log_det(log_det_a, sign_det, Omega_a);
  arma::log_det(log_det_b, sign_det, Omega_b);
  double logprior = 0.5 * log_det_a - 0.5 * log(sigma2)  - 0.5 * as_scalar(alpha.t() * Omega_a * alpha)/sigma2; 
  logprior += 0.5 * log_det_b - 0.5 * log(sigma2) - 0.5 * as_scalar(beta.t() * Omega_b * beta)/sigma2; 

  // we can choose if we want to compute the likelihood of sigma2 here or together with the log prior of gamma **TODO**
  
  // now the gaussian approximation of the denominator
  // p(\theta \vert \gamma, sigma2, Y)
  // since it's computed on the MAP alpha_beta, we only have the normalizing constant.

  // Following two lines do not work, given that it's not log_det and gives over/underflow
  // dlib::matrix<double> hessian2 = fun_hess_like_pois(alpha_beta);
  // double denom_gauss_approx = 0.5 * log(det(hessian)); 
  arma::mat hessian = fun_hess_like_pois(alpha, beta);
  double log_det_hess;
  arma::log_det(log_det_hess, sign_det, hessian);

  return (loglik + logprior - 0.5 * log_det_hess);
  // we are missing sigma! we need a new branch where sigma is fixed together with the particle. **TODO**
}


double Particle::optim(dlib::matrix<double,0,1> &alpha_beta){
  // // https://stackoverflow.com/questions/56069273/passing-member-functions-as-arguments-to-optimizer-passing-member-functions-to
  // double max_objective = 1e+10; // careful: should be large enough. think about it more carefully (TODO)

  // // First attempt: using bfgs_search_strategy (approximate newton)
  // // Does not really work and gives "The objective function generated non-finite outputs" for some i (like i = 11)
  // double tmp = find_max(dlib::bfgs_search_strategy(),
  //                       dlib::objective_delta_stop_strategy(1e-7),
  //                       [this](const dlib::matrix<double,0,1>& a) { 
  //                           return this->fun_like_pois(a,false);
  //                       }, [this](const dlib::matrix<double,0,1>& a) { 
  //                           return this->fun_deriv_like_pois(a,false);
  //                       }, alpha_beta, max_objective);

  // // First attempt B: This did not improve
  // double tmp = find_max_using_approximate_derivatives(dlib::bfgs_search_strategy(),
  //                                        dlib::objective_delta_stop_strategy(1e-7),
  //                                        [this](const dlib::matrix<double,0,1>& a) { 
  //                                           return this->fun_like_pois(a);
  //                                         }, starting_point, 1000);

  // // Second attempt: using the newton search strategy
  // // For small n it runs but sometimes finds non-optimal values. For large n, "The objective function generated non-finite outputs"
  // double tmp = find_max(dlib::newton_search_strategy(
  //             [this](const dlib::matrix<double,0,1>& a) {
  //                 return this->fun_hess_like_pois(a,false);
  //             }),
  //            dlib::objective_delta_stop_strategy(1e-8),
  //            [this](const dlib::matrix<double,0,1>& a) { 
  //                 return this->fun_like_pois(a,false);
  //             }, [this](const dlib::matrix<double,0,1>& a) { 
  //                 return this->fun_deriv_like_pois(a,false);
  //             }, alpha_beta, max_objective); 

  // Third attempt: using the constrained search
  // First we specify what the upper and lower limit of the box search are:
  dlib::matrix<double,0,1> lower(2*(Y.n_rows));
  dlib::matrix<double,0,1> upper(2*(Y.n_rows));
  // TODO: choose better constrains
  for(int i =0; i < Y.n_rows; i++){
    lower(i) = -10;
    lower(i+Y.n_rows) = -5;
    upper(i) = 10;
    upper(i+Y.n_rows) = 5;
  }
  double tmp = find_max_box_constrained(dlib::newton_search_strategy(
              [this](const dlib::matrix<double,0,1>& a) {
                  return this->fun_hess_like_pois(a,false);
              }),
             dlib::objective_delta_stop_strategy(1e-6),
             // dlib::gradient_norm_stop_strategy(1e-3), // did not work
             [this](const dlib::matrix<double,0,1>& a) { 
                  return this->fun_like_pois(a,false);
              }, [this](const dlib::matrix<double,0,1>& a) { 
                  return this->fun_deriv_like_pois(a,false);
              }, 
              alpha_beta, 
              lower, upper); 

  
  return tmp;
}

void Particle::optim_mle(dlib::matrix<double,0,1> &alpha_beta){
  cout << "optim_mle" << endl;
  // double max_objective = 1e+10; // careful: should be large enough. think about it more carefully (TODO)

  // // First attempt: using the multivariate objective, bfgs_search_strategy (approximate newton)
  // find_max(dlib::bfgs_search_strategy(),
  //         dlib::objective_delta_stop_strategy(1e-7),
  //         [this](const dlib::matrix<double,0,1>& a) {
  //             return this->fun_like_pois(a,true);
  //         }, [this](const dlib::matrix<double,0,1>& a) {
  //             return this->fun_deriv_like_pois(a,true);
  //         }, alpha_beta, max_objective);
  // // gives wrong results. sometimes non-finite values. 
  // // I tried using approximate_derivatives to see if my derivative function was wrong. But same results (find_max_using_approximate_derivatives)
  
  // Second attempt: using the multivariate objective, newton search strategy
  // find_max(dlib::newton_search_strategy(
  //             [this](const dlib::matrix<double,0,1>& a) {
  //                 return this->fun_hess_like_pois(a,true);
  //             }),
  //            dlib::objective_delta_stop_strategy(1e-5),
  //            [this](const dlib::matrix<double,0,1>& a) { 
  //                 return this->fun_like_pois(a,true);
  //             }, [this](const dlib::matrix<double,0,1>& a) { 
  //                 return this->fun_deriv_like_pois(a,true);
  //             }, alpha_beta, max_objective); 
  // Works for n = 5. Otherwise:
  // The objective function generated non-finite outputs. 
  // Because it's searching values of alpha = 24 or 30 which become Inf when we do exp(Xgamma)
  // It might be because it's inverting the full hessian matrix and there might be numerical errors? (the inverson seems to be done using LU factorization)
  
  // // Third attempt: using the multivariate objective, constrained search
  // // First we specify what the upper and lower limit of the box search are:
  dlib::matrix<double,0,1> lower(2*(Y.n_rows));
  dlib::matrix<double,0,1> upper(2*(Y.n_rows));
  for(int i =0; i < Y.n_rows; i++){
    lower(i) = -10;
    lower(i+Y.n_rows) = -5;
    upper(i) = 10;
    upper(i+Y.n_rows) = 5;
  }
  find_max_box_constrained(dlib::newton_search_strategy(
              [this](const dlib::matrix<double,0,1>& a) {
                  return this->fun_hess_like_pois(a,true);
              }),
             dlib::objective_delta_stop_strategy(1e-6),
             // dlib::gradient_norm_stop_strategy(1e-3), // did not work
             [this](const dlib::matrix<double,0,1>& a) { 
                  return this->fun_like_pois(a,true);
              }, [this](const dlib::matrix<double,0,1>& a) { 
                  return this->fun_deriv_like_pois(a,true);
              }, 
              alpha_beta, 
              lower, upper); 
  // // find_max_box_constrained was kind of working with n = 10, but goes to the boundary when n = 100
  // // using gradient_norm_stop_strategy makes the procedure very long (infinite loop?)

  // Fourth attempt: using univariate optimization. 
  // dlib::matrix<double,0,1> alpha_beta_tmp(2);
  // for(int i = 0; i < nObs; i++){
  //   // cout << i << " ";
  //   alpha_beta_tmp(0) = 1;
  //   alpha_beta_tmp(1) = -0.1;
  //   // // Fourth (a): works up to 10, then "The objective function generated non-finite outputs"
  //   // find_max_using_approximate_derivatives(dlib::bfgs_search_strategy(),
  //   //                                        dlib::objective_delta_stop_strategy(1e-4),
  //   //                                        [this, i](const dlib::matrix<double,0,1>& a) { // dlib::matrix<double,0,1> alpha_beta
  //   //                                             return this->fun_like_pois_singleunit(a,i);
  //   //                                         }, alpha_beta_tmp,  max_objective);
  //   // // Fourth (b): works up to 10, then "The objective function generated non-finite outputs"
  //   // find_max(dlib::bfgs_search_strategy(),
  //   //          dlib::objective_delta_stop_strategy(1e-4),
  //   //          [this, i](const dlib::matrix<double,0,1>& a) { 
  //   //               return this->fun_like_pois_singleunit(a,i);
  //   //           }, [this,i](const dlib::matrix<double,0,1>& a) { 
  //   //               return this->fun_deriv_like_pois_singleunit(a,i);
  //   //           }, alpha_beta_tmp, max_objective);
  //   // // Fourth (c): works! even for n = 400.
  //   find_max(dlib::newton_search_strategy(
  //             [this](const dlib::matrix<double,0,1>& a) {
  //                 return this->fun_hess_like_pois_singleunit(a);
  //             }),
  //            dlib::objective_delta_stop_strategy(1e-4), // CAREFUL with tolerance -- should not be too small
  //            [this, i](const dlib::matrix<double,0,1>& a) { 
  //                 return this->fun_like_pois_singleunit(a,i);
  //             }, [this,i](const dlib::matrix<double,0,1>& a) { 
  //                 return this->fun_deriv_like_pois_singleunit(a,i);
  //             }, alpha_beta_tmp, max_objective); // CAREFUL WITH THE STARTING POINT -- should not be alpha_beta 
  //   alpha_beta(i) = alpha_beta_tmp(0);
  //   alpha_beta(i+nObs) = alpha_beta_tmp(1);
  // }
  
  // return tmp;
  return;
}

double Particle::total_log_like(){
  int T = Y.n_cols;
  double n = nObs;
  double log_like;
  double log_det = 0.0;
  double quad_form = accu(Y%Y);
  for(int k = 0; k < partition_A->K; k++){
    log_det += partition_A->log_det_Omegay[k];
    quad_form += partition_A->quad_forms[k];
  }
  for(int k = 0; k < partition_B->K; k++){
    log_det += partition_B->log_det_Omegay[k];
    quad_form += partition_B->quad_forms[k];
  }
  log_like = lgamma(alpha_sigma +  ( n * T/2)) + alpha_sigma * log(beta_sigma); // numerator
  log_like += - lgamma(alpha_sigma) - (alpha_sigma + ( n * T/2)) * log(beta_sigma + quad_form/2); // denominator
  log_like += 0.5 * log_det - (n * T/2)*log(2*M_PI);
  return log_like;
}

void Particle::Print_Particle_Short(){
  cout << "Partition A: " << endl;
  partition_A->Print_Partition_Short();
  cout << "Partition B: " << endl;
  partition_B->Print_Partition_Short();
  // cout << "alpha: "<< endl;
  // for(int i = 0; i < nObs; i++){
  //   cout << alpha[i] << " ";
  // }
  // cout << endl;
  // cout << "beta: "<< endl;
  // for(int i = 0; i < nObs; i++){
  //   cout << beta[i] << " ";
  // }
  // cout << endl;
  
  //cout << "hyper-parameters: " << sigma << " " << tau_A << " " << tau_B << endl;
}

void Particle::Print_Particle(){
  cout << "Partition A: " << endl;
  partition_A->Print_Partition();
  cout << "Partition B: " << endl;
  partition_B->Print_Partition();
  cout << "Log-posterior: ";
  if(opt_Y == 1){
    cout << LikelihoodY() << endl;
  } else if(opt_Y == 0){
    cout << LikelihoodYAB() << endl;
  } else if(opt_Y == 2){
    cout << total_log_like() << endl;
  } else if(opt_Y == 3){
    cout << LikelihoodY_MCMC_poisson() << endl;
  } else if(opt_Y == 4){
    cout << LikelihoodY_Laplace_Pois() << endl;
  }
  cout << "Log-posterior: " << Total_Log_Post() << endl;
  // cout << "alpha: "<< endl;
  // for(int i = 0; i < nObs; i++){
  //   cout << alpha[i] << " ";
  // }
  // cout << endl;
  // cout << "beta: "<< endl;
  // for(int i = 0; i < nObs; i++){
  //   cout << beta[i] << " ";
  // }
  // cout << endl;
}

void Particle::Print_Particle_ToFile(string file){
  partition_A->Print_Partition_ToFile(file);
  partition_B->Print_Partition_ToFile(file);
  return;
}

void Particle::Find_Splits_nonconnected(bool A_or_B, int cluster_id, int n_cl, int* component_i, int* component_not_i, int **index1_ptr, int **index2_ptr, int &n1, int &n2){
  LPPartition partition;
  double* parameter;
  if(A_or_B){
    partition = partition_A;
    parameter = alpha;
  } else {
    partition = partition_B;
    parameter = beta;
  }
  
  arma::mat A_block_cluster = Submatrix(A_block, n_cl, n_cl, component_i, component_i);
  arma::mat param_sim = zeros<mat>(n_cl, n_cl);
  arma::vec param_cluster(n_cl); // an armadillo vector holding the beta-hats
  for(int i = 0; i < n_cl; i++){
    param_cluster(i) = parameter[component_i[i]];
  }
  double beta_hat_var = arma::var(param_cluster); // variance of the beta_hats within the cluster
  // populate the param_similarity matrix
  double error = 0.0;
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0.0,beta_hat_var);
  for(int i = 0; i < n_cl - 1; i++){
    for(int j = i; j < n_cl; j++){
      // error = distribution(generator);
      param_sim(i,j) = exp(-1 * (param_cluster(i) - param_cluster(j) + error) * (param_cluster(i) - param_cluster(j) + error)/(2 * beta_hat_var));
      param_sim(j,i) = exp(-1 * (param_cluster(i) - param_cluster(j) + error) * (param_cluster(i) - param_cluster(j) + error)/(2 * beta_hat_var));
    }
  }
  arma::mat diag_ncl(n_cl,n_cl,fill::eye);
  arma::mat W_beta_cl =  diag_ncl + param_sim % A_block_cluster;
  arma::mat Dinv_sqrt = arma::diagmat(1/sqrt(arma::sum(W_beta_cl, 1)));
  arma::mat L = diag_ncl - Dinv_sqrt * W_beta_cl * Dinv_sqrt;
  arma::vec eigval; // the smallest eigenvalues are the first two
  arma::mat eigvec;
  eig_sym(eigval, eigvec, L);
  mat U = eigvec.cols(0,1);
  U = arma::diagmat(1/sqrt(arma::sum(arma::square(U), 1))) * U;
  arma::mat means;
  // kmeans(means, U.t(), 2, random_subset, 10, false);
  bool status = arma::kmeans(means, U.t(), 2, random_subset, 10, false);
  if(status == false)
    cout << "clustering failed" << endl;
  int * membership = which_is_nearest_k(means, U.t());
  (*index1_ptr) = new int[n_cl];
  (*index2_ptr) = new int[partition->cluster_config[cluster_id]];
  n1 = 0;
  n2 = 0;
  for(int i = 0; i < n_cl; i++){
    if(membership[i] == 0){
      (*index1_ptr)[n1] = component_i[i];
      (n1)++;
    }
    else {
      (*index2_ptr)[n2] = component_i[i];
      (n2)++;
    }
  }
  delete[] membership;
  for(int i = 0; i < (partition->cluster_config[cluster_id]-n_cl); i++){
    (*index2_ptr)[n2] = component_not_i[i];
    (n2)++;
  }
  Partition_Split(A_or_B, cluster_id, *index1_ptr, *index2_ptr, n1, n2);
}

void Particle::Find_Splits_nonspatial(bool A_or_B, int cluster_id, int **index1_ptr, int **index2_ptr, int &n1, int &n2){
  LPPartition partition;
  double* parameter;
  if(A_or_B){
    partition = partition_A;
    parameter = alpha;
  } else {
    partition = partition_B;
    parameter = beta;
  }
  int n_cl = partition->cluster_config[cluster_id];

  arma::vec param_cluster(n_cl); // an armadillo vector holding the beta-hats
  for(int i = 0; i < n_cl; i++){
    param_cluster(i) = parameter[partition->clusters[cluster_id][i]];
  }

  arma::mat means;
  bool status = arma::kmeans(means, param_cluster.t(), 2, random_subset, 10, false);
  if(status == false)
    cout << "clustering failed" << endl;
  int * membership = which_is_nearest_k(means, param_cluster.t());
  (*index1_ptr) = new int[n_cl];
  (*index2_ptr) = new int[n_cl];
  n1 = 0;
  n2 = 0;
  for(int i = 0; i < n_cl; i++){
    if(membership[i] == 0){
      (*index1_ptr)[n1] = partition->clusters[cluster_id][i];
      (n1)++;
    }
    else {
      (*index2_ptr)[n2] = partition->clusters[cluster_id][i];
      (n2)++;
     }
   }
   delete[] membership;
   Partition_Split(A_or_B, cluster_id, *index1_ptr, *index2_ptr, n1, n2);
}

bool Particle::Find_K_SpectralSplits(bool A_or_B, int split_k, int num_splits, std::vector<std::vector<int> >& indices, std::vector<int>& ns, bool use_mle){
  // indices now contains numbers between 0 and nObs-1
  if(num_splits <= 1){
    return false;
  } else {
    LPPartition partition;
    double* parameter;
    if(A_or_B){
      partition = partition_A;
      parameter = alpha;
      if(use_mle){
        parameter = alpha_mle;
      }
    } else {
      partition = partition_B;
      parameter = beta;
      if(use_mle){
        parameter = beta_mle;
      }
    }
    // int SCiter = 2;
    int n_cl = partition->cluster_config[split_k];

    /* populate the similarity_matrixilarity matrix */
    arma::vec param_cluster(n_cl);
    for(int i = 0; i < n_cl; i++){
      param_cluster(i) = parameter[partition->clusters[split_k][i]];
    }
    double param_var = arma::var(param_cluster);
    
    arma::mat similarity_matrix = zeros<mat>(n_cl, n_cl);
    /* // for randomness, not used for now 
    std::default_random_engine generator; 
    std::normal_distribution<double> distribution(0.0,param_var); 
    */
    double error = 0.0; // for randomness, not used for now
    for(int i = 0; i < n_cl - 1; i++){
      for(int j = i; j < n_cl; j++){
        // error = distribution(generator); // for randomness, not used for now
        similarity_matrix(i,j) = exp(-1 * (param_cluster(i) - param_cluster(j) + error) * (param_cluster(i) - param_cluster(j) + error)/(2 * param_var));
        similarity_matrix(j,i) = exp(-1 * (param_cluster(i) - param_cluster(j) + error) * (param_cluster(i) - param_cluster(j) + error)/(2 * param_var));
      }
    }
    
    /* Here goes spectral clustering */
    arma::mat A_block_cluster = Submatrix(A_block, n_cl, n_cl, partition->clusters[split_k], partition->clusters[split_k]);
    arma::mat diag_n(n_cl,n_cl,fill::eye);
    arma::mat W_beta_cl =  diag_n + similarity_matrix % A_block_cluster;
    arma::mat Dinv_sqrt = arma::diagmat(1/sqrt(arma::sum(W_beta_cl, 1)));
    arma::mat L = diag_n - Dinv_sqrt * W_beta_cl * Dinv_sqrt;
    arma::vec eigval; // the smallest eigenvalues are the first two
    arma::mat eigvec;
    eig_sym(eigval, eigvec, L);
    mat U = eigvec.cols(0,num_splits-1);
    U = arma::diagmat(1/sqrt(arma::sum(arma::square(U), 1))) * U;
    arma::mat means;

    indices.clear();
    indices.resize(num_splits);
    ns.clear();
    ns.resize(num_splits);

    double min_score; // actually I don't think I'll use it
    int * membership;
    bool status = false;
    int loop_check = 0;
    int kmeans_rep_tmp = kmeans_rep;
    if(num_splits > 3){
      kmeans_rep_tmp = (num_splits - 2)*kmeans_rep;
    }
    while((!status) && (loop_check < 10)){
      status = kmeans_repeat(means, U, num_splits, kmeans_rep_tmp, min_score);
      if(status == false){
        cout << "clustering failed" << endl;
        // go directly to next iteration: i.e. do kmeans again
      }
      loop_check++;
    }
    // membership is long n_cl and has elements between 0 and num_splits-1
    membership = which_is_nearest_k(means, U.t());
        
    for(int j = 0; j<num_splits; j++){
      ns[j] = 0;
    }

    // CHANGED: now indices has the actual indices from 0 to nObs-1
    // before indices had elements between 0 and n_cl-1
    for(int i = 0; i < n_cl; i++){
      indices[membership[i]].push_back( partition->clusters[split_k][i] );
      (ns[membership[i]])++;
    }
  
    int *components;
    int *count;
    for(int j = 0; j < num_splits; j++){
      mat A_cluster = Submatrix(A_block, ns[j], ns[j], indices[j], indices[j]);
      components = new int[ ns[j] ];
      count = new int;
      Connected_Components(A_cluster, ns[j],components,count);
      if((*count)!=1){
        // instead of iterating let's change indices:
        // save the connected components of indices[j] as a vector
        std::vector<std::vector<int> > new_indices;
        new_indices.resize(*count);
        for(int i = 0; i < ns[j]; i++){
          new_indices[ components[i] ].push_back(indices[j][i]);
        }
        // we can eliminate the current element of indices[j] and ns[j]
        indices[j].clear();
        // now we update the previous element indices[j] and add new ones
        for(int c = 0; c < (*count); c++){
          if(c == 0){
            indices[j] = new_indices[0];
            ns[j] = new_indices[0].size();
          } else {
            indices.push_back(new_indices[c]);
            ns.push_back(new_indices[c].size());
          }
        }
      }
      delete[] components;
      delete count;
    }
    delete[] membership;
  
    Partition_KSplit(A_or_B, split_k, indices, ns);
    return true;
  }
}

bool Particle::Find_KM_Splits(bool A_or_B, int split_k, int num_splits, std::vector<std::vector<int> >& indices, std::vector<int>& ns, bool use_mle){
  // indices now contains numbers between 0 and nObs-1
  if(num_splits <= 1){
    return false;
  } else {
    LPPartition partition;
    double* parameter;
    if(A_or_B){
      partition = partition_A;
      parameter = alpha;
      if(use_mle){
        parameter = alpha_mle;
      }
    } else {
      partition = partition_B;
      parameter = beta;
      if(use_mle){
        parameter = beta_mle;
      }
    }
    // int SCiter = 2;
    int n_cl = partition->cluster_config[split_k];
    arma::mat U = arma::zeros<mat>(n_cl, 1); // holds the data passed to k-means
    U.set_size(n_cl, 1);
    for(int i = 0; i < n_cl; i++){
      U(i,0)  = parameter[partition->clusters[split_k][i]];
    }
    arma::mat means;

    indices.clear();
    indices.resize(num_splits);
    ns.clear();
    ns.resize(num_splits);

    double min_score;
    int * membership;
    bool status = false;
    int loop_check = 0;
    while(!status){
      // status = arma::kmeans(means, U.t(), num_splits, random_subset, 10, false);
      status = kmeans_repeat(means, U, num_splits, kmeans_rep, min_score);
      if(status == false){
        cout << "clustering failed" << endl;
        // go directly to next iteration: i.e. do kmeans again
      }
      loop_check++;
      if(loop_check > 3){
        return -1;
      }
    }
    // membership is long n_cl and has elements between 0 and num_splits-1
    membership = which_is_nearest_k(means, U.t());
    
    for(int j = 0; j<num_splits; j++){
      ns[j] = 0;
    }

    // indices has the actual indices from 0 to nObs-1
    for(int i = 0; i < n_cl; i++){
      indices[membership[i]].push_back( partition->clusters[split_k][i] );
      (ns[membership[i]])++;
    }

    int *components;
    int *count;
    for(int j = 0; j < num_splits; j++){
      mat A_cluster = Submatrix(A_block, ns[j], ns[j], indices[j], indices[j]);
      components = new int[ ns[j] ];
      count = new int;
      Connected_Components(A_cluster, ns[j],components,count);
      if((*count)!=1){
        // instead of iterating let's change indices:
        // save the connected components of indices[j] as a vector
        std::vector<std::vector<int> > new_indices;
        new_indices.resize(*count);
        for(int i = 0; i < ns[j]; i++){
          new_indices[ components[i] ].push_back(indices[j][i]);
        }
        // we can eliminate the current element of indices[j] and ns[j]
        indices[j].clear();
        // now we update the previous element indices[j] and add new ones
        for(int c = 0; c < (*count); c++){
          if(c == 0){
            indices[j] = new_indices[0];
            ns[j] = new_indices[0].size();
          } else {
            indices.push_back(new_indices[c]);
            ns.push_back(new_indices[c].size());
          }
        }
      }
      delete[] components;
      delete count;
    }
    delete[] membership;

    Partition_KSplit(A_or_B, split_k, indices, ns);
    return true;
  }
}

void Particle::Find_TailSplit(bool A_or_B, int left0right2, int i, int split_k, std::vector<std::vector<int> > left_center_right, std::vector<std::vector<int> >& tail_conncomp, std::vector<std::vector<int> >& indices, std::vector<int>& ns){
  // adds element i to left_tail and finds its connected component (tail_conncomp), then splits
  std::vector<int> remain; // contains all of the indices that are not split
  std::vector<std::vector<int> > remain_conncomp;
  for(int ii = i+1; ii < (int)left_center_right[left0right2].size(); ii++){
    remain.push_back(left_center_right[left0right2][ii]);
  }
  for(int ii = 0; ii < (int)left_center_right[1].size(); ii++){
    remain.push_back(left_center_right[1][ii]);
  }
  for(int ii = 0; ii < (int)left_center_right[2-left0right2].size(); ii++){
    remain.push_back(left_center_right[2-left0right2][ii]);
  }
  remain_conncomp = Alternative_Connected_Components(remain);
  Alternative_Connected_Components(left_center_right[left0right2][i], tail_conncomp);
  
  // prepare indices for KSplit
  indices.clear();
  ns.clear();
  for(int cc = 0; cc < (int)tail_conncomp.size(); cc++){
    indices.push_back(tail_conncomp[cc]);
    ns.push_back(tail_conncomp[cc].size());
  }
  for(int cc = 0; cc < (int)remain_conncomp.size(); cc++){
    indices.push_back(remain_conncomp[cc]);
    ns.push_back(remain_conncomp[cc].size());
  }
  Partition_KSplit(A_or_B, split_k, indices, ns);
  return;
}

void Particle::Initial_K_Splits(int num_splits_A, int num_splits_B, bool use_mle){
  std::vector<std::vector<int> > indices;
  std::vector<int> ns;
  if(num_splits_A > 1)
    Find_K_SpectralSplits(1, 0, num_splits_A, indices, ns, use_mle);
  if(num_splits_B > 1)
    Find_K_SpectralSplits(0, 0, num_splits_B, indices, ns, use_mle);
  return;
}


void Particle::Initial_KM_Splits(int num_splits_A, int num_splits_B, bool use_mle){
  std::vector<std::vector<int> > indices;
  std::vector<int> ns;
  if(num_splits_A > 1)
    Find_KM_Splits(1, 0, num_splits_A, indices, ns, use_mle);
  if(num_splits_B > 1)
    Find_KM_Splits(0, 0, num_splits_B, indices, ns, use_mle);
  return;
}


double Particle::param_bar(bool A_or_B, int k){
  LPPartition partition;
  double* parameter;
  if(A_or_B){
    partition = partition_A;
    parameter = alpha;
  } else {
    partition = partition_B;
    parameter = beta;
  }
  arma::vec param_cl(partition->cluster_config[k]);
  for(int i = 0; i < partition->cluster_config[k]; i++){
    param_cl(i) = parameter[ partition->clusters[k][i] ];
  }
  return arma::mean(param_cl);
}

void Particle::get_leftcenterright(bool A_or_B, std::vector<std::vector<int> >& left_center_right, int split_k){
  LPPartition partition;
  double* parameter;
  if(A_or_B){
    partition = partition_A;
    parameter = alpha;
  } else {
    partition = partition_B;
    parameter = beta;
  }
  int n_k = partition->cluster_config[split_k];
  int size_tail = ceil(n_k*split_frac);

  arma::vec param_cluster(1);
  arma::uvec param_indices(1);
  param_cluster.set_size(n_k); // holds the alpha_hat's within the cluster
  for(int i = 0; i < n_k; i++){
    param_cluster(i) = parameter[partition->clusters[split_k][i]];
  }
  param_indices = arma::sort_index(param_cluster, "ascend"); // sort in ascending order
  
  left_center_right.clear();
  left_center_right.resize(3);
  for(int i = 0; i < size_tail; i++){
    left_center_right[0].push_back(partition->clusters[split_k][param_indices(i)]);
    left_center_right[2].push_back(partition->clusters[split_k][param_indices(n_k-1 - i)]);
  }
  for(int i = size_tail; i < n_k-size_tail; i++){
    left_center_right[1].push_back(partition->clusters[split_k][param_indices(i)]);
  }
  return;
}


void Particle::update_sigma2(){
  arma::vec alpha(nObs);
  arma::vec beta(nObs);
  dlib::matrix<double,0,1> alpha_beta(2*(Y.n_rows));
  for(int i =0; i < Y.n_rows; i++){
    alpha_beta(i) = 0;
    alpha_beta(i+Y.n_rows) = 0;
  }
  optim(alpha_beta);
  for(int i = 0; i < nObs; i++){
    alpha[i] = alpha_beta(i);
    beta[i] = alpha_beta(i+nObs);
  }

  arma::mat Omega_a = get_Omega(true);
  arma::mat Omega_b = get_Omega(false);

  double quad_form = as_scalar(alpha.t() * Omega_a * alpha + beta.t() * Omega_b * beta);
  double mode = (quad_form + beta_sigma)/(alpha_sigma + ((double) nObs * Y.n_cols/2) + 1);
  // double mean = (quad_form + beta_sigma)/(alpha_sigma + ((double) nObs * Y.n_cols/2) - 1);
  sigma2_laplace = mode;
  return;
}
