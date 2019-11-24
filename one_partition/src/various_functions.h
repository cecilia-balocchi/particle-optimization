//
//  various_functions.h


#ifndef GUARD_VARIOUS_FUNCTIONS_H
#define GUARD_VARIOUS_FUNCTIONS_H
#include <RcppArmadillo.h>

#include <stdio.h>
arma::mat Submatrix(arma::mat M, int n_rows, int n_cols, std::vector<int> row_index, std::vector<int> col_index);
void Connected_Components(const arma::mat M, const int n, int* components, int* count);
void DFSUtil(const arma::mat M, const int n, int v, bool* visited, int* components, int* count);
arma::mat Distance_matrix(arma::vec alpha_hat, int nObs);

void new_Connected_Components(const arma::mat &M, const int n, std::vector<int> &init_components, std::vector<std::vector<int> >&components);

std::vector<std::vector<int> > Alternative_Connected_Components(std::vector<int> remain, const arma::mat &A_block);
void Alternative_Connected_Components(int element, std::vector<std::vector<int> >& left_new_clusters, const arma::mat &A_block);


void kmeans_repeat(arma::mat &U, arma::mat &means, bool &final_status, double &min_score, int num_splits, int reps);
void kmeans_plus_plus(arma::mat &U, arma::mat &means, bool &final_status, double &min_score, int num_splits,int reps);

#endif
