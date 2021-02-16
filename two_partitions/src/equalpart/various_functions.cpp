/*
 * various_functions.cpp
 *
 *  
 */
#include <RcppArmadillo.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <vector>
#include <random>
#include <armadillo>
#include "various_functions.h"
using namespace std;
using namespace arma;

extern arma::mat Y;
extern arma::mat X;
extern arma::mat A_block;


mat Submatrix(mat M, int n_row, int n_col, int* row_index, int* col_index){
  // creates a copy.
  mat N(n_row, n_col);
  for(int i = 0; i < n_row; i++){
    for(int j = 0; j < n_col; j++){
      N(i,j) = M(row_index[i],col_index[j]);
    }
  }
  return N;
}

mat Submatrix(mat M, int n_row, int n_col, std::vector<int> row_index, int* col_index){
  // creates a copy.
  mat N(n_row, n_col);
  for(int i = 0; i < n_row; i++){
    for(int j = 0; j < n_col; j++){
      N(i,j) = M(row_index[i],col_index[j]);
    }
  }
  return N;
}

mat Submatrix(mat M, int n_row, int n_col, std::vector<int> row_index, std::vector<int> col_index){
  // creates a copy.
  mat N(n_row, n_col);
  for(int i = 0; i < n_row; i++){
    for(int j = 0; j < n_col; j++){
      N(i,j) = M(row_index[i],col_index[j]);
    }
  }
  return N;
}

void Connected_Components(mat M, int n, int* components, int* count)
{
  // components is a vector where element i corresponds to the component of 
  // the ith element in the cluster considered. (from 0 to *count-1)
  // Mark all the vertices as not visited
  bool *visited = new bool[n];
  *count = 0;
  for(int v = 0; v < n; v++)
    visited[v] = false;
  for(int v = 0; v < n; v++)
  {
    if(visited[v] == false)
    {
      DFSUtil(M, n, v, visited, components, count);
      (*count)++;
    }
  }
  delete[] visited;
}
void DFSUtil(mat M, int n, int v, bool* visited, int* components, int* count)
{
  visited[v] = true;
  components[v] = *count;
  for(int i = 0; i < n; i++)
    if(M(v,i) == 1.0)
      if(visited[i] == false)
        DFSUtil(M, n, i, visited, components, count);
}

std::vector<std::vector<int> > Alternative_Connected_Components(std::vector<int> remain){
  std::vector<std::vector<int> > current_conncomp; 
  std::vector<std::vector<int> > tmp_new_conncomp; // temporarily used to build the connected components.
  current_conncomp.push_back(std::vector<int>(1,remain[0]));
  arma::mat A_tmp; 
  for(int ii = 1; ii < (int)remain.size(); ii++){
    // consider the next element in remain and check if it's connected to any of the current_conncomp
    tmp_new_conncomp.push_back(std::vector<int>(1, remain[ii]));
    // loop over the existing connected components
    for(int conncomp = 0; conncomp < (int)current_conncomp.size(); conncomp++){
      A_tmp = Submatrix(A_block, tmp_new_conncomp[0].size(), current_conncomp[conncomp].size(), tmp_new_conncomp[0], current_conncomp[conncomp]);
      if(any(vectorise(A_tmp) == 1.0)){ // something in current_conncomp[conncomp] is adjacent to remain[ii] so we should combine these clusters
        for(int ix = 0; ix < (int)current_conncomp[conncomp].size(); ix++){
          tmp_new_conncomp[0].push_back(current_conncomp[conncomp][ix]);
        }
      } else{ // current_conncomp[conncomp] remains its own distinct component
        tmp_new_conncomp.push_back(current_conncomp[conncomp]);
      }
      current_conncomp[conncomp].clear();
    } // closes loop over elements of current_conncomp
    // update current_conncomp: copy tmp_new_conncomp
    current_conncomp.clear();
    for(int cc = 0; cc < (int)tmp_new_conncomp.size(); cc++){
      current_conncomp.push_back(tmp_new_conncomp[cc]);
      tmp_new_conncomp[cc].clear();
    }
    tmp_new_conncomp.clear();
  } // closes loop over elements of remain used to determine connected components of remain
  return current_conncomp;
}

void Alternative_Connected_Components(int element, std::vector<std::vector<int> >& current_conncomp){
  // Adds element to the current_conncomp and updates it.
  std::vector<std::vector<int> > tmp_new_conncomp; // temporarily used to build the connected components.
  tmp_new_conncomp.push_back(std::vector<int>(1, element));
  arma::mat A_tmp; 
  // loop over the existing connected components
  for(int cc = 0; cc < (int)current_conncomp.size(); cc++){
    A_tmp = Submatrix(A_block, tmp_new_conncomp[0].size(), current_conncomp[cc].size(), tmp_new_conncomp[0], current_conncomp[cc]);
    if(any(vectorise(A_tmp) == 1.0)){ // something in current_conncomp[cc] is adjacent to element so we should combine these clusters
      for(int ix = 0; ix < (int)current_conncomp[cc].size(); ix++){
        tmp_new_conncomp[0].push_back(current_conncomp[cc][ix]);
      }
    } else{ // current_conncomp[cc] remains its own distinct component
      tmp_new_conncomp.push_back(current_conncomp[cc]);
    }
    current_conncomp[cc].clear();
  } // closes loop over elements of current_conncomp
  current_conncomp.clear();
  // update current_conncomp: copy tmp_new_conncomp
  for(int cc = 0; cc < (int)tmp_new_conncomp.size(); cc++){
    current_conncomp.push_back(tmp_new_conncomp[cc]);
    tmp_new_conncomp[cc].clear();
  }
  tmp_new_conncomp.clear();
  
  return ;
}

/*
int* which_is_nearest(arma::mat centroids, arma::mat data){
  // this is written for only two centroids
  int *membership;
  membership = new int[data.n_cols];
  int n0 = 0;
  int n1 = 0;
  for(int point_ind = 0; point_ind < data.n_cols; point_ind++){
    double dist0 = arma::norm(data.col(point_ind) - centroids.col(0));
    double dist1 = arma::norm(data.col(point_ind) - centroids.col(1));

    if(dist0 < dist1){
      membership[point_ind] = 0;
      n0++;
    } else if(dist0 > dist1) {
      membership[point_ind] = 1;
      n1++;
    } else { // if the distance is the same
      if(n0 == 0){
        membership[point_ind] = 0;
        n0++;
      } else if (n1 ==0){
        membership[point_ind] = 1;
        n1++;
      } else {
        membership[point_ind] = 0;
        n0++;
      }
    }
  }
  return membership;
}
*/

int* which_is_nearest_k(arma::mat centroids, arma::mat data){
  // centroids should be means, data should be U.t()
  // this is written for a general number of centroids
  int num_splits = centroids.n_cols;
  int n_k = data.n_cols;
  int *membership;
  membership = new int[n_k];
  arma::vec distance = zeros<vec>(num_splits);
  arma::uvec indices(num_splits);
  // cout << "--membership:--" << endl;
  for(int i = 0; i < n_k; i++){
    distance.zeros();
    for(int cent_index = 0; cent_index < num_splits; cent_index++){
      distance(cent_index) = arma::norm(data.col(i) - centroids.col(cent_index));
    }
    indices = arma::sort_index(distance, "ascend");
    membership[i] = indices(0);
    // cout << membership[i] << " ";
  }
  // cout << endl;
  return membership;
}
/*
int* which_is_nearest_k(arma::mat centroids, arma::mat data){
  // this is written for only two centroids
  int *membership;
  membership = new int[data.n_cols];
  double dists [centroids.n_cols];
  // cout << "--membership:--" << endl;
  for(int point_ind = 0; point_ind < data.n_cols; point_ind++){
    double min_val;
    int min_cent;
    for(int cent_ind = 0; cent_ind < centroids.n_cols; cent_ind++){
      dists[cent_ind] = arma::norm(data.col(point_ind) - centroids.col(cent_ind));
      if(cent_ind == 0){
        min_cent = cent_ind;
        min_val = dists[cent_ind];
      } else {
        if(dists[cent_ind] < min_val){
          min_cent = cent_ind;
          min_val = dists[cent_ind];
        }
      }
    }
    for(int cent_ind = 0; cent_ind < centroids.n_cols; cent_ind++){
      if(min_cent == cent_ind){
        membership[point_ind] = cent_ind;
      }
    }
    // cout << membership[point_ind] << " ";
  }
  // cout << endl;
  return membership;
}
*/

mat Distance_matrix(double* beta_hat, int nObs){
  mat dist(nObs,nObs, fill::zeros);
  for(int i = 0; i < nObs; i++)
  {
    for(int j = i+1; j < nObs; j++)
    {
      dist(i,j) = abs(beta_hat[i] - beta_hat[j]);
      dist(j,i) = dist(i,j);
    }
  }
  return dist;
}

double lbeta(double a, double b){
  return lgamma(a) + lgamma(b) - lgamma(a+b);
}

double lbinomial(int n, int k){
  return lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
}

double lmgamma(int p, double x){
  if(p == 1){
    return lgamma(x);
  } else {
    double c = 0.5*(p-1) * log(datum::pi) + lgamma(x);
    return c + lmgamma(p-1,x-0.5);
  }
}

arma::mat block_inverse(arma::mat A11, arma::mat A12, arma::mat A21, arma::mat A22){
  arma::mat A11inv = inv_sympd(A11);
  arma::mat temp = A22 - A21 * (A11inv * A12);
  arma::mat Ainv22 = inv_sympd(temp);
  arma::mat Ainv11 = A11inv + A11inv * A12 * Ainv22 * A21 * A11inv;
  arma::mat Ainv12 = - A11inv * A12 * Ainv22;
  arma::mat Ainv21 = - Ainv22 * A21 * A11inv;
  return arma::join_cols(arma::join_rows(Ainv11, Ainv12),arma::join_rows(Ainv21, Ainv22));
}

arma::mat block_inverse_ret(arma::mat A11, arma::mat A12, arma::mat A21, arma::mat A22, arma::mat* B11, arma::mat* B12, arma::mat* B21, arma::mat* B22){
  // the blocks of the inverse are passed by reference and modified inside the function; moreover the full matrix is returned
  arma::mat A11inv = inv_sympd(A11);
  arma::mat temp = A22 - A21 * (A11inv * A12);
  (*B22) = inv_sympd(temp);
  *B11 = A11inv + A11inv * A12 * (*B22) * A21 * A11inv;
  *B12 = - A11inv * A12 * (*B22);
  *B21 = - (*B22) * A21 * A11inv;
  return arma::join_cols(arma::join_rows(*B11, *B12),arma::join_rows(*B21, *B22));
}

double block_log_det(arma::mat A11, arma::mat A12, arma::mat A21, arma::mat A22){
  double logdet = 0;
  double tmp_log_det, tmp_log_det_sgn;
  arma::log_det(tmp_log_det, tmp_log_det_sgn, A11);
  logdet += tmp_log_det;
  arma::mat temp = A22 - A21 * inv_sympd(A11) * A12;
  arma::log_det(tmp_log_det, tmp_log_det_sgn, temp);
  logdet += tmp_log_det;
  return logdet;
}

bool is_orthogonal(arma::mat X){
  if( accu( abs( sum(X, 1) ) ) < 1e-10 ){
    return TRUE;
  } else {
    return FALSE;
  }
}

bool kmeans_repeat(arma::mat &means, arma::mat &U, int num_splits, int reps, double &min_score){
  // pick the kmeans that minimizes the sum of the distances squared of points to their centroids
  int n = U.n_rows; // each row of U is an observation. note that arma::kmeans wants the observations stored as columns.
  int d = U.n_cols; // dimension of the data
  arma::mat tmp_means(d, num_splits); // arma::kmeans returns the centroids stored as column vectors
  bool status;
  bool final_status = false; // if it returns false,
  int * membership;
  arma::mat cluster_means = arma::zeros<arma::mat>(num_splits, d); // each row is the mean of elements in cluster
  arma::vec cluster_counts = arma::zeros<arma::vec>(num_splits);
  
  double score = 0.0;
  min_score = 0.0;
  for(int r = 0; r < reps; r++){
    score = 0.0;
    status = arma::kmeans(tmp_means, U.t(), num_splits, arma::random_subset, 10, false);
    
    cluster_means.zeros();
    cluster_counts.zeros();
    
    if(status == false){
      cout << "kmeans failed";
    } else{
      final_status = true;
      membership = which_is_nearest_k(tmp_means, U.t());
      for(int i = 0; i < n; i++){
        cluster_means.row(membership[i]) += U.row(i); // in cluster membership[i] we save the average of all the U rows for that cluster
        ++(cluster_counts(membership[i]));
      }
      
      // at this point, rows of cluster_means contains SUM of the elements in the cluster and not their mean
      for(int k = 0; k < num_splits; k++){
        cluster_means.row(k) /= cluster_counts(k);
      }
      score = 0.0;
      for(int i = 0; i < n; i++){
        score += pow(arma::norm( U.row(i) - cluster_means.row(membership[i]) ),2);
      }
      if(r == 0){
        min_score = score;
        means = tmp_means;
      } else if(score < min_score){
        min_score = score;
        means = tmp_means;
        // Rcpp::Rcout << "r = " << r << " new min score = " << min_score << endl;
      }
    }
  }
  //Rcpp::Rcout << "[kmeans_repeat]: min_score = " << min_score << std::endl;
  return final_status;
}
