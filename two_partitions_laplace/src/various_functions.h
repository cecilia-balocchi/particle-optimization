/*
 * various_functions.h
 *
 */


#ifndef VARIOUS_FUNCTIONS_H_
#define VARIOUS_FUNCTIONS_H_
#include <RcppArmadillo.h>
#include <dlib/matrix.h>

#include <stdio.h>

arma::mat Submatrix(arma::mat M, int n_rows, int n_cols, int* row_index, int* col_index);
arma::mat Submatrix(arma::mat M, int n_rows, int n_cols, std::vector<int> row_index, int* col_index);
arma::mat Submatrix(arma::mat M, int n_rows, int n_cols, std::vector<int> row_index, std::vector<int> col_index);
void Connected_Components(arma::mat M, int n, int* components, int* count);
std::vector<std::vector<int> > Alternative_Connected_Components(std::vector<int> remain);
void Alternative_Connected_Components(int element, std::vector<std::vector<int> >& current_conncomp);
void DFSUtil(arma::mat M, int n, int v, bool* visited, int* components, int* count);
// int* which_is_nearest(arma::mat centroids, arma::mat data);
int* which_is_nearest_k(arma::mat centroids, arma::mat data);
arma::mat Distance_matrix(double* beta_hat, int nObs);
double lbeta(double a, double b);
double lbinomial(int n, int k);
double lmgamma(int p, double x);
arma::mat block_inverse(arma::mat A11, arma::mat A12, arma::mat A21, arma::mat A22);
arma::mat block_inverse_ret(arma::mat A11, arma::mat A12, arma::mat A21, arma::mat A22, arma::mat* B11, arma::mat* B12, arma::mat* B21, arma::mat* B22);
double block_log_det(arma::mat A11, arma::mat A12, arma::mat A21, arma::mat A22);
double logsumexp(arma::vec l);
bool is_orthogonal(arma::mat X);
bool kmeans_repeat(arma::mat &means, arma::mat &U, int num_splits, int reps, double &min_score);
#endif /* VARIOUS_FUNCTIONS_H_ */
