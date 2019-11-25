/*
 * partition.cpp
 *
 */

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

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
#include "various_functions.h"
using namespace std;
using namespace arma;

extern arma::mat Y;
extern arma::mat X;
extern arma::mat A_block;
extern arma::mat Lambda;
extern arma::mat sX;
// Set the hyper-parameters
extern double rho;
// extern double a;
// extern double b;
extern double a1;
extern double a2;
extern double b1;
extern double b2;
extern double a_cohes;
extern double b_cohes;
extern double alpha_sigma;
extern double beta_sigma;
extern double eta;
extern double sigma_py;
extern int priorA;
extern int priorB;
extern int* sigma;

//class Partition::Begins
Partition::Partition(){
  nObs = 0;
  K = 0;
  cluster_config = NULL;
  clusters = NULL;
  cluster_assignment = NULL;
  pairwise_assignment = NULL;
  log_prior = NULL;
  log_cohesion = 0.0;
  beta_hat = NULL;
  log_det_Omegay = NULL;
  quad_forms = NULL;
  return;
}
// initializes a partition
Partition::Partition(LPPartition initial_partition, bool A_or_B){
  int prior;
  if(A_or_B){
    prior = priorA;
  } else {
    prior = priorB;
  }
  nObs = initial_partition->nObs;
  K = initial_partition->K;
  cluster_config = new int[K];
  clusters = new int*[K];
  cluster_assignment = new int[nObs];
  pairwise_assignment = new int*[nObs];
  log_prior = new double[K];
  // beta_hat = new double[nObs];
  beta_hat = NULL;
  log_det_Omegay = new double[K];
  quad_forms = new double[K];

  for(int k = 0; k < K; k++){
    cluster_config[k] = initial_partition->cluster_config[k];
    clusters[k] = new int[cluster_config[k]];
    for(int i = 0; i < cluster_config[k]; i++){
      clusters[k][i] = initial_partition->clusters[k][i];
    }
    log_det_Omegay[k] = initial_partition->log_det_Omegay[k];
    quad_forms[k] = initial_partition->quad_forms[k];
  }

  for(int i = 0; i < nObs; i++){
    cluster_assignment[i] = initial_partition->cluster_assignment[i];
    pairwise_assignment[i] = new int[nObs];
    for(int j = 0; j < nObs; j++){
      pairwise_assignment[i][j] = initial_partition->pairwise_assignment[i][j];
    }
  }
  get_prior_all(prior, sigma);
  return;
}

Partition::~Partition(){
  int i;
  delete[] cluster_config; cluster_config = NULL;
  for(i = 0; i < K; i++){
    delete[] clusters[i]; clusters[i] = NULL;
  }
  delete[] clusters; clusters = NULL;
  delete[] cluster_assignment; cluster_assignment = NULL;
  for(i = 0; i < nObs; i++){
    delete[] pairwise_assignment[i]; pairwise_assignment[i] = NULL;
  }
  delete[] pairwise_assignment; pairwise_assignment = NULL;
  delete[] log_prior; log_prior = NULL;
  delete[] log_det_Omegay; log_det_Omegay = NULL;
  delete[] quad_forms; quad_forms = NULL;
  // delete[] beta_hat; 
  beta_hat = NULL;
  return;
}

void Partition::Copy_Partition(LPPartition initial_partition, bool A_or_B){
  delete[] cluster_config;
  delete[] log_prior;
  delete[] log_det_Omegay; 
  delete[] quad_forms; 

  int prior;
  if(A_or_B){
    prior = priorA;
  } else {
    prior = priorB;
  }

  for(int k = 0; k < K; k++){
    delete[] clusters[k];
  }
  delete[] clusters;
  K = initial_partition->K;

  cluster_config = new int[K];
  log_prior = new double[K];
  log_det_Omegay = new double[K];
  quad_forms = new double[K];
  clusters = new int*[K];

  for(int i = 0; i < nObs; i++){
    cluster_assignment[i] = initial_partition->cluster_assignment[i];
    for(int j = 0 ; j < nObs; j++){
      pairwise_assignment[i][j] = initial_partition->pairwise_assignment[i][j];
    }
  }
  for(int k = 0; k < K; k++){
    cluster_config[k] = initial_partition->cluster_config[k];
    clusters[k] = new int[cluster_config[k]];
    for(int i = 0; i < cluster_config[k]; i++){
      clusters[k][i] = initial_partition->clusters[k][i];
    }
    log_det_Omegay[k] = initial_partition->log_det_Omegay[k];
    quad_forms[k] = initial_partition->quad_forms[k];
  }
  get_prior_all(prior, sigma);
  return;
}

void Partition::Initialize_Partition(int n, bool A_or_B){
  int prior;
  if(A_or_B){
    prior = priorA;
  } else {
    prior = priorB;
  }
  nObs = n;
  K = 1;
  cluster_config = new int[K];
  cluster_config[0] = nObs; // initial partition has 1 cluster of size nObs
  clusters = new int*[K];
  cluster_assignment = new int[nObs];

  log_det_Omegay = new double[K];
  quad_forms = new double[K];
  log_prior = new double[K];

  for(int k = 0; k < K; k++){
    clusters[k] = new int[cluster_config[k]];
  }
  // now since K = 1, we only have one cluster:
  for(int i = 0; i < nObs; i++){
    clusters[0][i] = i;
  }
  for(int i = 0; i < nObs; i++){
    cluster_assignment[i] = 0; // initial partition has 1 cluster
  }
  get_pairwise();
  
  get_prior_all(prior, sigma);
  get_likelihood(A_or_B, 0);
  return;
}

void Partition::Initialize_Partition(int n, Rcpp::List gamma_init, bool A_or_B){
  int prior;
  if(A_or_B){
    prior = priorA;
  } else {
    prior = priorB;
  }
  nObs = n;
  K = gamma_init.size();
  cluster_config = new int[K];
  cluster_assignment = new int[nObs];
  clusters = new int*[K];
  log_det_Omegay = new double[K];
  quad_forms = new double[K];
  log_prior = new double[K];
  beta_hat = NULL;
  
  Rcpp::NumericVector tmp_vec;
  for(int k = 0; k < K; k++){
    tmp_vec = gamma_init[k];
    cluster_config[k] = tmp_vec.size();
    clusters[k] = new int[cluster_config[k]];
    for(int i = 0; i < cluster_config[k]; i++){
      clusters[k][i] = tmp_vec[i] - 1; // R is 1-indexed and C++ is 0-indexed
      cluster_assignment[clusters[k][i]] = k;
    }
  }
  get_prior_all(prior, sigma);
  get_pairwise();
  for(int k = 0; k < K; k++){
    get_likelihood(A_or_B, k);
  }
  return;
}

void Partition::Initialize_Partition_nclusters(int n, bool A_or_B){
  int prior;
  if(A_or_B){
    prior = priorA;
  } else {
    prior = priorB;
  }
  nObs = n;
  K = n;
  cluster_config = new int[K];
  clusters = new int*[K];
  for(int k = 0; k < K; k++){
    cluster_config[k] = 1;
    clusters[k] = new int;
    clusters[k][0] = k;
  }
  cluster_assignment = new int[nObs];
  for(int i = 0; i < nObs; i++){
    cluster_assignment[i] = i; // initial partition has 1 cluster
  }
  log_det_Omegay = new double[K];
  quad_forms = new double[K];
  log_prior = new double[K];
  // beta_hat = new double[nObs];
  get_pairwise();
  get_prior_all(prior, sigma);
  for(int k = 0; k < K; k++){
    get_likelihood(A_or_B, k);
  }
  return;
}

/* // code that I don't use anymore
void Partition::Initialize_Partition3(int n, int k, int* cl_sizes, int* cl_limits){
  nObs = n;
  K = k;
  cluster_config = new int[K];
  clusters = new int*[K];
  cluster_assignment = new int[nObs];
  
  log_prior = new double[K];
  // beta_hat = new double[nObs];

  for(int k = 0; k < K; k++){
    cluster_config[k] = cl_sizes[k];
    clusters[k] = new int[cluster_config[k]];
    int i_tmp = 0;
    if(k == 0){
      for(int i = 0; i < nObs; i++){
        if(i < cl_limits[k]){
          clusters[k][i_tmp] = i;
          cluster_assignment[i] = k;
          i_tmp++;
        }
      }
    } else if(k == K){
      for(int i = 0; i < nObs; i++){
        if(i >= cl_limits[k-1]){
          clusters[k][i_tmp] = i;
          cluster_assignment[i] = k;
          i_tmp++;
        }
      }
    } else{
      for(int i = 0; i < nObs; i++){
        if(i >= cl_limits[k-1] && i < cl_limits[k]){
          clusters[k][i_tmp] = i;
          cluster_assignment[i] = k;
          i_tmp++;
        }
      }
    }
  }

  get_pairwise();

  get_prior_all(prior, sigma);
  return;
}

// useful just to try a handful for partitions
// corresponding to the example with square grid
// here n = 100
void Partition::Initialize_Partition2(int id){

  nObs = 100;

  if(id == 0){
  K = 3;
  cluster_config = new int[K];
  cluster_config[0] = 25;
  cluster_config[1] = 25;
  cluster_config[2] = 50;
  clusters = new int*[K];
  clusters[0] = new int[25];
  clusters[1] = new int[25];
  clusters[2] = new int[50];

  log_prior = new double[K];
  // beta_hat = new double[nObs];
  cluster_assignment = new int[nObs];

  int counter0 = 0;
  int counter1 = 0;
  for(int i = 0; i < 5; i++){
    for(int j = 0; j < 5; j++){
    clusters[0][counter0] = 10*i + j;
    clusters[1][counter1] = 10*(i+5) + j;
    counter0++;
    counter1++;
    cluster_assignment[10*i+j] = 0;
    cluster_assignment[10*(i+5) + j] = 1;
    }
  }
  int counter2 = 0;
  for(int i = 0; i < 10; i++){
    for(int j = 5; j < 10; j++){
    clusters[2][counter2] = 10*i + j;
    cluster_assignment[10*i + j] = 2;
    counter2++;
    }
  }
  } else if(id == 1){
  // everything belongs to a single cluster
  //Initialize_Partition(100, Y,X, A_block, rho, a, b, alpha, nu, eta);
  Initialize_Partition(100);
  } else if(id == 2){
  // put original clusters 1 and 2 together
  nObs = 100;
    K = 2;
    cluster_config = new int[K];
    cluster_config[0] = 50;
    cluster_config[1] = 50;
    clusters = new int*[K];
    clusters[0] = new int[50];
    clusters[1] = new int[50];
    log_prior = new double[K];
    // beta_hat = new double[nObs];

    cluster_assignment = new int[nObs];

        int counter0 = 0;

    for(int i = 0; i < 5; i++){
      for(int j = 0; j < 5; j++){
        clusters[0][counter0] = 10*i + j;
        counter0++;
        clusters[0][counter0] = 10*(i+5) + j;
        counter0++;
        cluster_assignment[10*i+j] = 0;
        cluster_assignment[10*(i+5) + j] = 0;
      }
    }
    int counter1 = 0;
    for(int i = 0; i < 10; i++){
      for(int j = 5; j < 10; j++){
        clusters[1][counter1] = 10*i + j;
        cluster_assignment[10*i + j] = 1;
        counter1++;
      }
    }
  } else if(id == 3){
    // everything is in its own cluster
    nObs = 100;
    K = 100;
    cluster_config = new int[K];
    clusters = new int*[K];
    cluster_assignment = new int[nObs];
    log_prior = new double[K];
    // beta_hat = new double[nObs];

    for(int i = 0; i < nObs; i++){
      cluster_config[i] = 1;
      clusters[i] = new int[1];
      clusters[i][0] = i;
      cluster_assignment[i] = i;
    }
  } else if(id == 4){
    nObs = 100;
    K = 2;
    cluster_config = new int[K];
    cluster_config[0] = 60;
    cluster_config[1] = 40;
    clusters = new int*[K];
    clusters[0] = new int[60];
    clusters[1] = new int[40];
    log_prior = new double[K];
    // beta_hat = new double[nObs];
    cluster_assignment = new int[nObs];

      int counter0 = 0;
    for(int i = 0; i < 10; i++){
      for(int j = 0; j < 2; j++){
      clusters[0][counter0] = 10*i + j;
      counter0++;
      clusters[0][counter0] = 10*i + (j+4);
      counter0++;
      clusters[0][counter0] = 10*i + (j+8);
      counter0++;
      cluster_assignment[10*i+j] = 0;
      cluster_assignment[10*i + (j+4)] = 0;
      cluster_assignment[10*i + (j+8)] = 0;
      }
    }
    int counter1 = 0;
    for(int i = 0; i < 10; i++){
    for(int j = 2; j < 4; j++){
      clusters[1][counter1] = 10*i + j;
      cluster_assignment[10*i + j] = 1;
      counter1++;
    }
      for(int j = 6; j < 8; j++){
      clusters[1][counter1] = 10*i + j;
      cluster_assignment[10*i + j] = 1;
      counter1++;
    }
    }
  }
  if(id != 1){
    get_pairwise();
    get_prior_all(prior, sigma);
  }
  return;
}

void Partition::Initialize_Partition_FromFile(int n, string file){
  nObs = n;
  cluster_assignment = new int[nObs];
  Read_Partition_FromFile(file, n);
  cout << "read" << endl;
  K = 0;
  vector<int> checked;
  vector<int> cl_sizes;
  int cl_size;
  int checked_size;
  bool i_in_check;
  for(int i = 0; i < nObs; i++){
    i_in_check = false;
    // check if i has already been checked
    checked_size = checked.size();
    if(checked_size>0){
      for(int jj = 0; jj < checked_size; jj++){
        if(checked[jj] == i){
          i_in_check = true;
        }
      }
    }
    if(!i_in_check){
      cl_size = 0;
      for(int ii = 0; ii < nObs; ii++){
        if(pairwise_assignment[i][ii] == 1){
          // ii is included in this loop
          checked.push_back(ii);
          cl_size++;
          cluster_assignment[ii] = K;
        }
      }
      // checked.push_back(i);
      cl_sizes.push_back(cl_size);
      // cluster_assignment[i] = K;
      K++;
    }
  }
  cout << checked_size << endl;
  cluster_config = new int[K];
  clusters = new int*[K];
  for(int k = 0; k < K; k++){
    cluster_config[k] = cl_sizes[k];
    clusters[k] = new int[cluster_config[k]];
  }
  int count;
  for(int k = 0; k < K; k++){
    count = 0;
    for(int i = 0; i < nObs; i++){
      if(cluster_assignment[i] == k){
        clusters[k][count] = i;
        count++;
      }
    }
  }
  log_prior = new double[K];
  get_prior_all(prior, sigma);
  return;
}
*/
// Function to get the pairwise assignments
void Partition::get_pairwise(){
  pairwise_assignment = new int*[nObs];
  for(int i = 0; i < nObs; i++){
    pairwise_assignment[i] = new int[nObs];
  }
  for(int i = 0; i < nObs; i++){
    for(int j = i; j < nObs; j++){
      if(cluster_assignment[i] == cluster_assignment[j]){
        pairwise_assignment[i][j] = 1;
        pairwise_assignment[j][i] = 1;
      } else{
        pairwise_assignment[i][j] = 0;
        pairwise_assignment[j][i] = 0;
      }
    }
  }
  return;

}
// Function to print out the partition
void Partition::Print_Partition(){
  std::cout << "Number of clusters K: " << K << std::endl;
  std::cout << "Size of clusters:";
  for(int k = 0; k < K; k++){
    std::cout << cluster_config[k] << " ";
  }
  std::cout << std::endl;
  std::cout << "Clusters:" << std::endl;
  for(int k = 0; k < K; k++){
    std::cout << "Cluster " << k  << " : ";
    for(int j = 0; j < cluster_config[k]; j++){
      std::cout << clusters[k][j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Log-prior:" ;
  for(int k = 0; k < K; k++){
    std::cout <<  log_prior[k] << " ";
  }
  std::cout << std::endl;
  return;
}
void Partition::Print_Partition_Short(){
  std::cout << "Number of clusters K: " << K << std::endl;
  std::cout << "Size of clusters:";
  for(int k = 0; k < K; k++){
  std::cout << cluster_config[k] << " ";
  }
  std::cout << std::endl;
  return;
}
void Partition::Print_Partition_ToFile(string file){
  ofstream myfile;
  myfile.open(file.c_str(), std::ofstream::app);
  // myfile.open(file, std::ios_base::app);
  for(int i = 0; i < nObs-1; i++){
    for(int j = i+1; j < nObs; j++){
      myfile << pairwise_assignment[i][j] << ",";
    }
  }
  myfile << endl;
  myfile.close();
  return;
}
void Partition::Read_Partition_FromFile(string file, int n){
  ifstream filein ( file.c_str() );
  string value;
  int value_int;

  // mat mat_temp(n,n,fill::randu);
  int aa = 0;
  pairwise_assignment = new int*[nObs];
  for(int i = 0; i < nObs; i++){
    pairwise_assignment[i] = new int[nObs];
  }
  cout << "before loop" << endl;
  for(int i = 0; i < nObs-1; i++){
    pairwise_assignment[i][i] = 1;
    // cout << "i" << endl;
    // mat_temp(i,i) = 1;
    for(int j = i+1; j < nObs; j++){
      getline ( filein, value, ',' );
      value_int = stoi(value);
      pairwise_assignment[j][i] = value_int;
      pairwise_assignment[i][j] = value_int;
      // mat_temp(i,j) = value_int;
      // mat_temp(j,i) = value_int;
      // cout << value_int << " ";
      aa++;
    }
    // cout << endl;
  }
  cout << aa << endl;
  // cout << mat_temp << endl;
  filein.close();
  return;
}

void Partition::Print_Means(){
  std::cout << "Posterior means of beta:" << std::endl;
  for(int i = 0; i < nObs; i++){
    std::cout << setprecision(6) << beta_hat[i] << " ";
  }
  std::cout << std::endl;
}


void Partition::log_pi_ep(int cluster_id){
  log_prior[cluster_id] = lgamma(cluster_config[cluster_id]);
}
void Partition::log_pi_epy(int cluster_id){
  log_prior[cluster_id] = lgamma(cluster_config[cluster_id] - sigma_py) - lgamma(1-sigma_py);
}

void Partition::log_cohes(double a_cohes, double b_cohes){
  int nkk = 0;
  for(int k = 0; k < K; k++){
    mat A_block_k = Submatrix(A_block, cluster_config[k], cluster_config[k], clusters[k], clusters[k]);
    nkk += arma::accu(A_block_k)/2;
  }
  int n_tot = arma::accu(A_block)/2;
  double logc = lbinomial(n_tot, nkk) + lbeta(nkk + a_cohes,n_tot - nkk + b_cohes) - lbeta(a_cohes,b_cohes);
  log_cohesion = logc;
  return;
}

void Partition::get_prior_onlylogprior(int cluster_id, int prior){
  if(prior == 0){
    // DP/PY prior
    log_pi_epy(cluster_id);
  } else if(prior == 1){
    // Areal PPMx
    log_pi_ep(cluster_id);
  } else if(prior == 2){
    // Dahl
    log_prior[cluster_id] = 0;
  } else if(prior == 3){
    // uniform
    log_prior[cluster_id] = 0.0;
  } else if(prior == 4){
    // Orbanz?
    mat A_block_k = Submatrix(A_block, cluster_config[cluster_id], cluster_config[cluster_id], clusters[cluster_id], clusters[cluster_id]);
    log_prior[cluster_id] = log(eta) + lgamma(cluster_config[cluster_id]) + arma::accu(A_block_k);
  }
  return;
}
void Partition::get_prior_onlycohesion(int prior, int* sigma){
  if(prior == 0){
    // DP/PY prior
    log_cohesion = -lgamma(eta + nObs) + lgamma(eta);
    for(int k = 0; k < K; k++){
      log_cohesion += log(eta+k*sigma_py);
    }
  } else if(prior == 1){
    // Areal PPMx
    log_cohes(a_cohes, b_cohes);
  } else if(prior == 2){
    // Dahl
    log_Dahl_cohes(sigma);
  } else if(prior == 3){
    // uniform
    log_cohesion = 0.0;
  } else if(prior == 4){
    log_cohesion = 0.0;
  }
  return;
}

void Partition::get_prior_all(int prior, int* sigma){
  if(prior == 0){
    // DP/PY prior
    for(int k = 0; k < K; k++){
      log_pi_ep(k);
    }
    log_cohesion = -lgamma(eta + nObs) + lgamma(eta);
    for(int k = 0; k < K; k++){
      log_cohesion += log(eta+k*sigma_py);
    }
  } else if(prior == 1){
    // Areal PPMx
    for(int k = 0; k < K; k++){
      log_pi_ep(k);
      // log_prior[k] = log(eta) + lgamma(cluster_config[k]);
    }
    log_cohes(a_cohes, b_cohes);
    log_cohesion += K*log(eta);
  } else if(prior == 2){
    // Dahl
    for(int k = 0; k < K; k++){
      log_prior[k] = 0.0;
    }
    log_Dahl_cohes(sigma);
    log_cohesion += K*log(eta);
  } else if(prior == 3){
    // uniform
    for(int k = 0; k < K; k++){
      log_prior[k] = 0.0;
    }
    log_cohesion = 0.0;
  } else if(prior == 4){
    for(int k = 0; k < K; k++){
      mat A_block_k = Submatrix(A_block, cluster_config[k], cluster_config[k], clusters[k], clusters[k]);
      log_prior[k] = log(eta) + lgamma(cluster_config[k]) + arma::accu(A_block_k);
    }
    log_cohesion = 0.0;
  }
  return;
}

void Partition::log_Dahl_cohes(int* sigma){
  // eta is public, Lambda is public
  // sigma should be long nObs
  int* z = new int[nObs];
  double lp = 0.0;
  int q = 0; // number of current clusters
  for(int t = 0; t < nObs; t++){
    int s = sigma[t];
    if(t == 0){
      z[t] = cluster_assignment[s];
      q ++;
    } else {
      int cl = cluster_assignment[s];
      vector<int> same_cl;
      for(int r = 0; r < t; r++){
        if(z[r] == cl){
          same_cl.push_back(sigma[r]);
        }
      }
      if(same_cl.size() == 0){// new cluster
        lp += log(eta) -log(eta + t); // not t-1 because 0-index
        q++;
      } else { 
        double n_const = 0.0;
        for(int j = 0; j < t; j++){
          n_const += Lambda(s, sigma[j]);
        }
        double l_const = 0.0;
        for(int i = 0; i < (int)same_cl.size(); i++){
          l_const += Lambda(s, same_cl[i]);
        }
        lp += log(t) -log(eta + t);
        lp += log(l_const) - log(n_const);
      }
      z[t] = cl;
    }
  }
  delete[] z;
  log_cohesion = lp;
  return;
}

void Partition::get_likelihood(bool A_or_B, int cluster_id){
  int T = Y.n_cols;
  rowvec x = X.row(0);
  double c1,c2;
  double txx;
  if(A_or_B){
    c1 = a1;
    c2 = a2;
    txx = T;
  } else {
    c1 = b1;
    c2 = b2;
    rowvec x2 = pow(x,2);
    txx = sum(x2);
  }
  int cluster_size = cluster_config[cluster_id];
  mat E = eye<mat>(cluster_size,cluster_size);
  mat Omega_k = zeros<mat>(cluster_size,cluster_size);
  mat O_k, A_block_k, D, M;
  mat Y_k = zeros<mat>(cluster_size,T);
  // Omega is the inverse of the covariance matrix of the parameter (for one cluster)
  if(cluster_size == 1){
    Omega_k(0,0) = 1/(c1/(1-rho)+c2);
    Y_k = Y.row(clusters[cluster_id][0]);
  } else {
    O_k = ones<mat>(cluster_size,cluster_size);
    // Creating M: the precision matrix of the CAR model
    A_block_k = Submatrix(A_block, cluster_size, cluster_size, clusters[cluster_id], clusters[cluster_id]);
    vec row_sums = zeros<vec>(cluster_size);
    for(int i = 0; i < cluster_size; i++){
      row_sums(i) = accu(A_block_k.row(i));
    }
    D = diagmat(row_sums);
    M = rho * (D - A_block_k);
    for(int i = 0; i < cluster_size; i++){
      M(i,i) += 1 - rho;
      Y_k.row(i) = Y.row(clusters[cluster_id][i]);
    }
    // This is the Woodbury matrix identity for inversion
    // Omega_k = (M/c1) - (M/c1) * O_k * (M/c1) / (sum(sum(M))/c1 + 1/c2);
    // Woodbury matrix identity + remember that M*1 = (1-rho)*1
    Omega_k = (M/c1) - ((1-rho)/c1) * ((1-rho)/c1) * O_k / ( cluster_size*(1-rho)/c1 + 1/c2);
  }
  
  arma::vec XtY;
  if(A_or_B){
    XtY = sum(Y_k,1);
  } else {
    XtY = Y_k * x.t();
  }
  mat P_k = inv_sympd(E*txx + Omega_k);
  quad_forms[cluster_id] = - as_scalar(XtY.t() * P_k * XtY);
  
  double Omega_log_det = 0;
  double Omega_log_det_sgn = 0;
  // Using Sylvester's + Woodbury determinant theorem: det(I + X Sigma Xt)^-1 = det(I - X P Xt) = det(I - P XtX) - so we need to multiply for 0.5
  arma::mat Sigma_det = E - txx * P_k;
  arma::log_det(Omega_log_det, Omega_log_det_sgn, Sigma_det);
  log_det_Omegay[cluster_id] = Omega_log_det;
  return;
}


// Function will split the cluster into two parts
// new cluster1 will still be called cluster split_k
// new cluster2 will be called cluster K+1
//void Partition::Split(int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2, mat Y, mat X, mat A_block, double rho, double a, double b, double alpha, double nu, double eta){
void Partition::Split(int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2, bool A_or_B){
  // create a temporary copy of the main attributes
  int prior;
  if(A_or_B){
    prior = priorA;
  } else {
    prior = priorB;
  }
  int orig_K = K;
  int* orig_cluster_config = new int[orig_K];
  int** orig_clusters = new int*[orig_K];
  int* orig_cluster_assignment = new int[nObs];
  double* orig_log_prior = new double[orig_K];
  double* orig_log_det_Omegay = new double[orig_K];
  double* orig_quad_forms = new double[orig_K];
  
  for(int k = 0; k < orig_K; k++){
    orig_cluster_config[k] = cluster_config[k];
    orig_clusters[k] = new int[cluster_config[k]];
    for(int i = 0; i < cluster_config[k]; i++){
      orig_clusters[k][i] = clusters[k][i];
    }
    orig_log_prior[k] = log_prior[k];
    orig_log_det_Omegay[k] = log_det_Omegay[k];
    orig_quad_forms[k] = quad_forms[k];
  }
  for(int i = 0; i < nObs; i++){
    orig_cluster_assignment[i] = cluster_assignment[i];
  }
  // clear up the memory from the original values and re-initialize
  // no need to delete pairwsie_assignment or cluster_assignment; the sizes are fixed
  K = orig_K+1; // we now have K+1 clusters
  delete[] cluster_config; cluster_config = NULL;
  cluster_config = new int[K];

  for(int k = 0; k < K-1; k++){
    delete[] clusters[k]; clusters[k] = NULL;
  }
  delete[] clusters; clusters = NULL;
  clusters = new int*[K];
  delete[] log_prior; log_prior = NULL;
  delete[] quad_forms; quad_forms = NULL;
  delete[] log_det_Omegay; log_det_Omegay = NULL;
  log_prior = new double[K];
  quad_forms = new double[K];
  log_det_Omegay = new double[K];
  for(int i = 0; i < nObs; i++){
    delete[] pairwise_assignment[i]; pairwise_assignment[i] = NULL;
  }
  delete[] pairwise_assignment; pairwise_assignment = NULL;
  // how big are the new clusters?

  // update cluster_config
  for(int k = 0; k < K; k++){
    if(k == split_k){
      cluster_config[k] = size1;
    } else if(k == K-1){
      cluster_config[k] = size2;
    } else{
      cluster_config[k] = orig_cluster_config[k];
    }
  }
  for(int k = 0; k < K; k++){
    clusters[k] = new int[cluster_config[k]];
    if(k == split_k){
      for(int i = 0; i < size1; i++){
        clusters[k][i] = new_cluster1[i];
      }

    } else if(k == K - 1){ // remember the 0-indexing...
      for(int i = 0; i < size2;i++){
        clusters[k][i] = new_cluster2[i];
      }
    } else{
      for(int i = 0; i < cluster_config[k]; i++){
        clusters[k][i] = orig_clusters[k][i];
      }
    }
  }

  // now update new_cluster_assignments.
  for(int i = 0; i < nObs; i++){
    if(orig_cluster_assignment[i] != split_k){
      cluster_assignment[i] = orig_cluster_assignment[i];
    }
  }
  for(int ii = 0; ii < size1; ii++){
    cluster_assignment[new_cluster1[ii]] = split_k;
  }
  for(int ii = 0; ii < size2; ii++){
    cluster_assignment[new_cluster2[ii]] = K - 1; // remember, cluster labels go from 0 to K-1.
  }
  // update the pairwise allocations
  get_pairwise();
  //update the log-likelihood and log-prior now
  for(int k = 0; k < K; k++){
    if(k == split_k){ // need to re-compute
      get_prior_onlylogprior(k, prior);
      get_likelihood(A_or_B, k);
    } else if(k == K-1){
      get_prior_onlylogprior(k, prior);
      get_likelihood(A_or_B, k);
    } else{
      log_prior[k] = orig_log_prior[k];
      quad_forms[k] = orig_quad_forms[k];
      log_det_Omegay[k] = orig_log_det_Omegay[k];
    }
  }
  get_prior_onlycohesion(prior, sigma);
  // free up memory by deleting the local copies
  delete[] orig_cluster_config;
  for(int kk = 0; kk < orig_K; kk++){ //original_clusters has length orig_K (K and not K+1.)
    delete[] orig_clusters[kk];
  }
  delete[] orig_clusters;
  delete[] orig_cluster_assignment;
  delete[] orig_log_prior;
  delete[] orig_quad_forms;
  delete[] orig_log_det_Omegay;
  return;
}

void Partition::Split(int split_k, std::vector<int> new_cluster1, std::vector<int> new_cluster2, int size1, int size2, bool A_or_B){
  // create a temporary copy of the main attributes
  int prior;
  if(A_or_B){
    prior = priorA;
  } else {
    prior = priorB;
  }
  int orig_K = K;
  int* orig_cluster_config = new int[orig_K];
  int** orig_clusters = new int*[orig_K];
  int* orig_cluster_assignment = new int[nObs];
  double* orig_log_prior = new double[orig_K];
  double* orig_log_det_Omegay = new double[orig_K];
  double* orig_quad_forms = new double[orig_K];
  
  for(int k = 0; k < orig_K; k++){
    orig_cluster_config[k] = cluster_config[k];
    orig_clusters[k] = new int[cluster_config[k]];
    for(int i = 0; i < cluster_config[k]; i++){
      orig_clusters[k][i] = clusters[k][i];
    }
    orig_log_prior[k] = log_prior[k];
    orig_log_det_Omegay[k] = log_det_Omegay[k];
    orig_quad_forms[k] = quad_forms[k];
  }
  for(int i = 0; i < nObs; i++){
    orig_cluster_assignment[i] = cluster_assignment[i];
  }
  // clear up the memory from the original values and re-initialize
  // no need to delete pairwsie_assignment or cluster_assignment; the sizes are fixed
  K = orig_K+1; // we now have K+1 clusters
  delete[] cluster_config; cluster_config = NULL;
  cluster_config = new int[K];

  for(int k = 0; k < K-1; k++){
    delete[] clusters[k]; clusters[k] = NULL;
  }
  delete[] clusters; clusters = NULL;
  clusters = new int*[K];
  delete[] log_prior; log_prior = NULL;
  delete[] quad_forms; quad_forms = NULL;
  delete[] log_det_Omegay; log_det_Omegay = NULL;
  log_prior = new double[K];
  quad_forms = new double[K];
  log_det_Omegay = new double[K];
  for(int i = 0; i < nObs; i++){
    delete[] pairwise_assignment[i]; pairwise_assignment[i] = NULL;
  }
  delete[] pairwise_assignment; pairwise_assignment = NULL;
  // how big are the new clusters?

  // update cluster_config
  for(int k = 0; k < K; k++){
    if(k == split_k){
      cluster_config[k] = size1;
    } else if(k == K-1){
      cluster_config[k] = size2;
    } else{
      cluster_config[k] = orig_cluster_config[k];
    }
  }
  for(int k = 0; k < K; k++){
    clusters[k] = new int[cluster_config[k]];
    if(k == split_k){
      for(int i = 0; i < size1; i++){
        clusters[k][i] = new_cluster1[i];
      }

    } else if(k == K - 1){ // remember the 0-indexing...
      for(int i = 0; i < size2;i++){
        clusters[k][i] = new_cluster2[i];
      }
    } else{
      for(int i = 0; i < cluster_config[k]; i++){
        clusters[k][i] = orig_clusters[k][i];
      }
    }
  }

  // now update new_cluster_assignments.
  for(int i = 0; i < nObs; i++){
    if(orig_cluster_assignment[i] != split_k){
      cluster_assignment[i] = orig_cluster_assignment[i];
    }
  }
  for(int ii = 0; ii < size1; ii++){
    cluster_assignment[new_cluster1[ii]] = split_k;
  }
  for(int ii = 0; ii < size2; ii++){
    cluster_assignment[new_cluster2[ii]] = K - 1; // remember, cluster labels go from 0 to K-1.
  }
  // update the pairwise allocations
  get_pairwise();
  //update the log-likelihood and log-prior now
  for(int k = 0; k < K; k++){
    if(k == split_k){ // need to re-compute
      get_prior_onlylogprior(k, prior);
      get_likelihood(A_or_B, k);
    } else if(k == K-1){
      get_prior_onlylogprior(k, prior);
      get_likelihood(A_or_B, k);
    } else{
      log_prior[k] = orig_log_prior[k];
      quad_forms[k] = orig_quad_forms[k];
      log_det_Omegay[k] = orig_log_det_Omegay[k];
    }
  }
  get_prior_onlycohesion(prior, sigma);
  // free up memory by deleting the local copies
  delete[] orig_cluster_config;
  for(int kk = 0; kk < orig_K; kk++){ //original_clusters has length orig_K (K and not K+1.)
    delete[] orig_clusters[kk];
  }
  delete[] orig_clusters;
  delete[] orig_cluster_assignment;
  delete[] orig_log_prior;
  delete[] orig_quad_forms;
  delete[] orig_log_det_Omegay;
  return;
}

void Partition::KSplit(int split_k, int num_splits, std::vector<std::vector<int> > indices, std::vector<int> ns, bool A_or_B){
  int prior;
  if(A_or_B){
    prior = priorA;
  } else {
    prior = priorB;
  }
  int orig_K = K;
  int* orig_cluster_config = new int[orig_K];
  int** orig_clusters = new int*[orig_K];
  int* orig_cluster_assignment = new int[nObs];
  double* orig_log_prior = new double[orig_K];
  double* orig_log_det_Omegay = new double[orig_K];
  double* orig_quad_forms = new double[orig_K];

  for(int k = 0; k < orig_K; k++){
    orig_cluster_config[k] = cluster_config[k];
    orig_clusters[k] = new int[cluster_config[k]];
    for(int i = 0; i < cluster_config[k]; i++){
      orig_clusters[k][i] = clusters[k][i];
    }
    orig_log_prior[k] = log_prior[k];
    orig_log_det_Omegay[k] = log_det_Omegay[k];
    orig_quad_forms[k] = quad_forms[k];
  }
  for(int i = 0; i < nObs; i++){
    orig_cluster_assignment[i] = cluster_assignment[i];
  }

  /* Now we free the space */
  delete[] cluster_config; cluster_config = NULL;
  for(int i = 0; i < orig_K; i++){
    delete[] clusters[i]; clusters[i] = NULL;
  }
  delete[] clusters; clusters = NULL;
  delete[] log_prior; log_prior = NULL;
  delete[] quad_forms; quad_forms = NULL;
  delete[] log_det_Omegay; log_det_Omegay = NULL;
  for(int i = 0; i < nObs; i++){
    delete[] pairwise_assignment[i]; pairwise_assignment[i] = NULL;
  }
  delete[] pairwise_assignment; pairwise_assignment = NULL;
  // no need to delete cluster_assignment; the sizes are fixed

  /* Here we update each of them */
  K = orig_K - 1 + num_splits;
  cluster_config = new int[K];
  clusters = new int*[K];
  for(int k = 0; k < K; k++){
    if(k == split_k){
      cluster_config[k] = ns[0];
    } else if(k < orig_K){
      cluster_config[k] = orig_cluster_config[k];
    } else{
      cluster_config[k] = ns[k - orig_K + 1];
    }
  }
  for(int k = 0; k < K; k++){
    clusters[k] = new int[cluster_config[k]];
    if(k == split_k){
      for(int i = 0; i < cluster_config[k]; i++){
        clusters[k][i] = indices[0][i];
      }
    } else if(k < orig_K){
      for(int i = 0; i < cluster_config[k]; i++){
        clusters[k][i] = orig_clusters[k][i];
      }
    } else{
      for(int i = 0; i < cluster_config[k]; i++){
        clusters[k][i] = indices[k - orig_K + 1][i];
      }
    }
  }
  for(int i = 0; i < nObs; i++){
    if(orig_cluster_assignment[i] != split_k){
      cluster_assignment[i] = orig_cluster_assignment[i];
    }
  }
  for(int kk = 0; kk < num_splits; kk++){
    if(kk == 0){
      for(int i = 0; i < ns[kk]; i++){
        cluster_assignment[ indices[kk][i] ] = split_k;
      }
    } else {
      for(int i = 0; i < ns[kk]; i++){
        cluster_assignment[ indices[kk][i] ] = kk + orig_K - 1;
      }  
    }
  }
  get_pairwise();
  log_prior = new double[K];
  quad_forms = new double[K];
  log_det_Omegay = new double[K];
  for(int k = 0; k < K; k++){
    if(k == split_k){ // need to re-compute
      get_prior_onlylogprior(k, prior);
      get_likelihood(A_or_B, k);
    } else if(k < orig_K){
      log_prior[k] = orig_log_prior[k];
      quad_forms[k] = orig_quad_forms[k];
      log_det_Omegay[k] = orig_log_det_Omegay[k];
    } else{
      get_prior_onlylogprior(k, prior);
      get_likelihood(A_or_B, k);
    }
  }
  get_prior_onlycohesion(prior, sigma);
  
  // if(opt_prior == 0){
  //   // just do what we were doing before
  //   for(int k = 0; k < K; k++){
  //     if(A_or_B){
  //       get_parameter(A_or_B, k);
  //     } else {
  //       get_parameter(A_or_B, k);
  //     }
  //   }
  // } else if(opt_prior == 1){
  //   // optimize alpha_beta given the two partitions
  //   get_alpha_beta();
  // } else if(opt_prior == 2){
  //   // optimize alpha_beta given the two partitions, by doing coordinate ascend
  //   get_alpha_beta_iter();
  // }
  delete[] orig_cluster_config;
  for(int kk = 0; kk < orig_K; kk++){ //original_clusters has length orig_K (K and not K+1.)
    delete[] orig_clusters[kk];
  }
  delete[] orig_clusters;
  delete[] orig_cluster_assignment;
  delete[] orig_log_prior;
  delete[] orig_quad_forms;
  delete[] orig_log_det_Omegay;
  return;
}

//void Partition::Merge(int k_1, int k_2, mat Y, mat X, mat A_block, double rho, double a, double b, double alpha, double nu, double eta){
void Partition::Merge(int k_1, int k_2, bool A_or_B){
  int prior;
  if(A_or_B){
    prior = priorA;
  } else {
    prior = priorB;
  }
  int k_max = max(k_1, k_2);
  int k_min = min(k_1, k_2);
  int new_cluster_size = cluster_config[k_min] + cluster_config[k_max];

    // make a pointer to the new merged cluster
  int* new_merged_cluster = new int[new_cluster_size];
  for(int i = 0; i < cluster_config[k_min]; i++){
    new_merged_cluster[i] = clusters[k_min][i];
  }
  for(int i = 0; i < cluster_config[k_max]; i++){
    new_merged_cluster[cluster_config[k_min] + i] = clusters[k_max][i];
  }
  // Update cluster_assignment
  // for original cluster k_max: this now becomes k_min
  // for clusters with original label greater than k_max, we need to decrement by 1
  int tmp_assignment = 0;
  for(int i = 0; i < nObs; i++){
    if(cluster_assignment[i] > k_max){ // we decrement cluster label by 1
      tmp_assignment = cluster_assignment[i];
      cluster_assignment[i] = tmp_assignment - 1;
    } else if(cluster_assignment[i] == k_max){
      cluster_assignment[i] = k_min;
    }
  }
  // make a temporary copy of clusters and cluster_config
  int orig_K = K;
  int* orig_cluster_config = new int[orig_K];
  int** orig_clusters = new int*[orig_K];
  int* orig_cluster_assignment = new int[nObs];
  double* orig_log_prior = new double[orig_K];
  double* orig_log_det_Omegay = new double[orig_K];
  double* orig_quad_forms = new double[orig_K];
  
  for(int kk = 0; kk < orig_K; kk++){
    orig_cluster_config[kk] = cluster_config[kk];
    orig_clusters[kk] = new int[cluster_config[kk]];
    for(int i = 0; i < cluster_config[kk]; i++){
      orig_clusters[kk][i] = clusters[kk][i];
    }
    orig_log_prior[kk] = log_prior[kk];
    orig_log_det_Omegay[kk] = log_det_Omegay[kk];
    orig_quad_forms[kk] = quad_forms[kk];
  }
  for(int i = 0; i < nObs; i++){
    orig_cluster_assignment[i] = cluster_assignment[i];
  }
  delete[] cluster_config;
  cluster_config = new int[K-1];

  for(int kk = 0; kk < orig_K; kk++){
    delete[] clusters[kk];
  }
  delete[] clusters;
  clusters = new int*[K-1];
  delete[] log_prior;
  log_prior = new double[K-1];
  delete[] quad_forms; quad_forms = NULL;
  delete[] log_det_Omegay; log_det_Omegay = NULL;
  quad_forms = new double[K-1];
  log_det_Omegay = new double[K-1];
  // looping over the OLD labels
  // remember the labels don't change until we get to k_max
  // this loop visits every cluster EXCEPT k_max
  for(int kk = 0; kk < orig_K; kk++){
    if(kk == k_min){
      cluster_config[kk] = new_cluster_size;
      clusters[kk] = new int[cluster_config[kk]];
      for(int i = 0; i < cluster_config[kk]; i++){
        clusters[kk][i] = new_merged_cluster[i];
      }
      get_prior_onlylogprior(kk, prior);
      get_likelihood(A_or_B, kk);
    } else if(kk < k_max){
      cluster_config[kk] = orig_cluster_config[kk];
      clusters[kk] = new int[cluster_config[kk]];
      for(int i = 0; i < cluster_config[kk]; i++){
        clusters[kk][i] = orig_clusters[kk][i];
      }
      log_prior[kk] = orig_log_prior[kk];
      quad_forms[kk] = orig_quad_forms[kk];
      log_det_Omegay[kk] = orig_log_det_Omegay[kk];
    } else if(kk > k_max){
      cluster_config[kk-1] = orig_cluster_config[kk];
      clusters[kk-1] = new int[cluster_config[kk-1]];// changed cluster_config[kk] into cluster_config[kk-1]
      for(int i = 0; i < cluster_config[kk-1]; i++){ // changed cluster_config[kk] into cluster_config[kk-1]
        clusters[kk-1][i] = orig_clusters[kk][i];
      }
      log_prior[kk-1] = orig_log_prior[kk];
      quad_forms[kk-1] = orig_quad_forms[kk];
      log_det_Omegay[kk-1] = orig_log_det_Omegay[kk];
    }
  }
  for(int i = 0; i < nObs; i++){
    delete[] pairwise_assignment[i]; pairwise_assignment[i] = NULL;
  }
  delete[] pairwise_assignment; pairwise_assignment = NULL;
  
  // update K
  K = orig_K - 1;

  get_pairwise();
  get_prior_onlycohesion(prior, sigma);
  
  // clean-up the memory
  delete[] orig_cluster_config;
  for(int kk = 0; kk < orig_K; kk++){ //original_clusters has length K and not K+1.
    delete[] orig_clusters[kk];
  }
  delete[] orig_clusters;
  delete[] orig_cluster_assignment;
  delete[] orig_log_prior;
  delete[] orig_quad_forms;
  delete[] orig_log_det_Omegay;
  delete[] new_merged_cluster;
}

//void Partition::Split_and_Merge(int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2, int k_star_1, int k_star_2, mat Y, mat X, mat A_block, double rho, double a, double b, double alpha, double nu, double eta){
void Partition::Split_and_Merge(int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2, int k_star_1, int k_star_2, bool A_or_B){

  if(k_star_1 != k_star_2){
    Split(split_k, new_cluster1, new_cluster2, size1, size2, A_or_B);
    if((split_k == k_star_1) & (split_k != k_star_2)){ // leave new_cluster1 alone and just attempt to merge new_cluster2 with k_star_2
      Merge(K-1, k_star_2, A_or_B); // remember new_cluster2's label is K-1
    } else if( (split_k != k_star_1) & (split_k == k_star_2)){
      // leave new_clsuter2 by itself and merge k_star_1 with new_cluster1
      // remember new_cluster1's label is still split_k;
      Merge(split_k, k_star_1, A_or_B);
    } else if((split_k != k_star_1) & (split_k != k_star_2) & (k_star_2 > max(split_k, k_star_1))){
      // We need to perform two merges
      Merge(split_k, k_star_1, A_or_B);
      // k_star_2's label is now decremented by 1
      // new_cluster2's label is now K
      Merge(K, k_star_2 - 1, A_or_B);
    } else if((split_k != k_star_1) & (split_k != k_star_2) & (k_star_2 < max(split_k, k_star_1))){
      Merge(split_k, k_star_1, A_or_B);
      Merge(K, k_star_2, A_or_B);
    }
  }
  return;
}

// Function to modify a partition if any of the clusters are disconnected
void Partition::Modify(int cl_ind, bool A_or_B){
  // cl_ind is the index of the cluster that needs to be modified
  if(cluster_config[cl_ind] > 1)
  {
    int n = cluster_config[cl_ind];
    int *index = clusters[cl_ind];
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
        Split(cl_ind, index1, index2, n1, n2, A_or_B);
        if(tmp > 0)
          delete[] index;
        index = index2;
        n = n2;
        new_components = new int[n2];
        for(int j = 0; j < n2; j++)
          new_components[j] = components[i_ind[j]];
        delete[] components;
        components = new_components;
        cl_ind = K-1;
        delete[] index1;
      }
      delete[] index2;
    }
    delete[] components;
    delete count;
  }

}

// Eventually we can just put this inside of update_particle
// void Partition::Find_Splits(int cluster_id, int **index1_ptr, int **index2_ptr, int &n1, int &n2, bool A_or_B){
//   int n_cl = cluster_config[cluster_id];
//   //double* beta_hat_cluster = new double[n_cl]; // we may get rid of this
//   arma::mat A_block_cluster = Submatrix(A_block, n_cl, n_cl, clusters[cluster_id], clusters[cluster_id]);
//   arma::mat beta_sim = zeros<mat>(n_cl, n_cl);
//   arma::vec beta_hat_cluster(n_cl); // an armadillo vector holding the beta-hats
//   for(int i = 0; i < n_cl; i++){
//     beta_hat_cluster(i) = beta_hat[clusters[cluster_id][i]];
//   }
//   double beta_hat_var = arma::var(beta_hat_cluster); // variance of the beta_hats within the cluster
//   // populate the beta_similarity matrix
//   double error = 0.0;
//   std::default_random_engine generator;
//   std::normal_distribution<double> distribution(0.0,beta_hat_var);
//   for(int i = 0; i < n_cl - 1; i++){
//     for(int j = i; j < n_cl; j++){
//       // error = distribution(generator);
//       beta_sim(i,j) = exp(-1 * (beta_hat_cluster(i) - beta_hat_cluster(j) + error) * (beta_hat_cluster(i) - beta_hat_cluster(j) + error)/(2 * beta_hat_var));
//       beta_sim(j,i) = exp(-1 * (beta_hat_cluster(i) - beta_hat_cluster(j) + error) * (beta_hat_cluster(i) - beta_hat_cluster(j) + error)/(2 * beta_hat_var));
//     }
//   }
//   arma::mat diag_ncl(n_cl,n_cl,fill::eye);
//   arma::mat W_beta_cl =  diag_ncl + beta_sim % A_block_cluster;
//   arma::mat Dinv_sqrt = arma::diagmat(1/sqrt(arma::sum(W_beta_cl, 1)));
//   arma::mat L = diag_ncl - Dinv_sqrt * W_beta_cl * Dinv_sqrt;
//   arma::vec eigval; // the smallest eigenvalues are the first two
//   arma::mat eigvec;
//   eig_sym(eigval, eigvec, L);
//   mat U = eigvec.cols(0,1);
//   U = arma::diagmat(1/sqrt(arma::sum(arma::square(U), 1))) * U;
//   arma::mat means;
//   // kmeans(means, U.t(), 2, random_subset, 10, false);
//   bool status = arma::kmeans(means, U.t(), 2, random_subset, 10, false);
//   if(status == false)
//     cout << "clustering failed" << endl;
//   int * membership = which_is_nearest_k(means, U.t());
//   (*index1_ptr) = new int[n_cl];
//   (*index2_ptr) = new int[n_cl];
//   n1 = 0;
//   n2 = 0;
//   for(int i = 0; i < n_cl; i++){
//     if(membership[i] == 0){
//       (*index1_ptr)[n1] = clusters[cluster_id][i];
//       (n1)++;
//     }
//     else {
//       (*index2_ptr)[n2] = clusters[cluster_id][i];
//       (n2)++;
//      }
//    }
//    delete[] membership;
//    Split(cluster_id, *index1_ptr, *index2_ptr, n1, n2, A_or_B);
//    // delete[] index1;
//    // delete[] index2;
// }

