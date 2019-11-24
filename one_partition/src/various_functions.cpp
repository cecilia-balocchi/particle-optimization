//
//  various_functions.cpp


#include "various_functions.h"
using std::endl;

arma::mat Submatrix(arma::mat M, int n_row, int n_col, std::vector<int> row_index, std::vector<int> col_index)
{
  arma::mat N(n_row, n_col);
  for(int i = 0; i < n_row; i++){
    for(int j = 0; j < n_col; j++){
      N(i,j) = M(row_index[i], col_index[j]);
    }
  }
  return N;
}
// will need to update these potentially

void Connected_Components(const arma::mat M, const int n, int* components, int* count)
{
  // Mark all the vertices as not visited
  bool *visited = new bool[n];
  *count = 0;
  for(int v = 0; v < n; v++)
    visited[v] = false;
  //Rcpp::Rcout << "[Connected_Components]: Starting main loop" << arma::endl;
  for(int v = 0; v < n; v++)
  {
    //Rcpp::Rcout << "visiting " << v << arma::endl;
    if(visited[v] == false)
    {
      DFSUtil(M, n, v, visited, components, count);
      (*count)++;
    }
  }
  delete[] visited;
}

// to be used in the split functions.
// reads in the sub-cluster found by the split (init_components) and forms the connected components (in the vector of vectors components)
void new_Connected_Components(const arma::mat &M, const int n, std::vector<int> &init_components, std::vector<std::vector<int> >&components)
{
  //Rcpp::Rcout << "[new_Connected_Components]: starting" << std::endl;
  int* count = new int;
  int* tmp_components = new int[n];
  bool* visited = new bool[n];
  
  *count = 0;
  for(int v = 0; v < n; v++){
    visited[v] = false;
  }
  //Rcpp::Rcout << "[new_Connected_Components] : about to start DFS" << std::endl;
  for(int v = 0; v < n; v++){
    if(visited[v] == false){
      DFSUtil(M,n,v,visited, tmp_components,count);
      (*count)++;
    }
  }
  //Rcpp::Rcout << "[new_Connected_Components]: Finished DFS" << std::endl;
  
  // use tmp_components to indx init_components
  components.clear();
  components.resize(*count);
  //int cluster_id = 0;
  for(int i = 0; i < n; i++){
    components[tmp_components[i]].push_back(init_components[i]);
  }
  delete count;
  delete[] tmp_components;
  delete[] visited;
}




void DFSUtil(const arma::mat M, const int n, int v, bool* visited, int* components, int* count)
{
  visited[v] = true;
  components[v] = *count;
  for(int i = 0; i < n; i++)
    if(M(v,i) == 1.0)
      if(visited[i] == false)
        DFSUtil(M, n, i, visited, components, count);
}

arma::mat Distance_matrix(arma::vec alpha_hat, int nObs){
  arma::mat dist = arma::zeros<arma::mat>(nObs, nObs);
  for(int i = 0; i < nObs; i++){
    for(int j = i+1; j < nObs; j++){
      dist(i,j) = abs(alpha_hat[i] - alpha_hat[j]);
      dist(j,i) = dist(i,j);
    }
  }
  return dist;
}







std::vector<std::vector<int> > Alternative_Connected_Components(std::vector<int> remain, const arma::mat &A_block){
  std::vector<std::vector<int> > remain_clusters; 
  std::vector<std::vector<int> > tmp_new_clusters; // temporarily used to build the connected components.
  remain_clusters.push_back(std::vector<int>(1,remain[0]));
  arma::mat A_tmp; 
  for(int ii = 1; ii < remain.size(); ii++){
    // consider the next element in remain and check if it's connected to any of the remain_clusters
    tmp_new_clusters.push_back(std::vector<int>(1, remain[ii]));
    // loop over the existing connected components
    for(int conncomp = 0; conncomp < remain_clusters.size(); conncomp++){
      A_tmp = Submatrix(A_block, tmp_new_clusters[0].size(), remain_clusters[conncomp].size(), tmp_new_clusters[0], remain_clusters[conncomp]);
      if(any(vectorise(A_tmp) == 1.0)){ // something in remain_clusters[conncomp] is adjacent to remain[ii] so we should combine these clusters
        for(int ix = 0; ix < remain_clusters[conncomp].size(); ix++){
          tmp_new_clusters[0].push_back(remain_clusters[conncomp][ix]);
        }
      } else{ // remain_clusters[conncomp] remains its own distinct component
        tmp_new_clusters.push_back(remain_clusters[conncomp]);
      }
      remain_clusters[conncomp].clear();
    } // closes loop over elements of remain_clusters
    // update remain_clusters: copy tmp_new_clusters
    remain_clusters.clear();
    for(int cc = 0; cc < tmp_new_clusters.size(); cc++){
      remain_clusters.push_back(tmp_new_clusters[cc]);
      tmp_new_clusters[cc].clear();
    }
    tmp_new_clusters.clear();
  } // closes loop over elements of remain used to determine connected components of remain
  return remain_clusters;
}



void Alternative_Connected_Components(int element, std::vector<std::vector<int> >& left_new_clusters, const arma::mat &A_block){
  std::vector<std::vector<int> > tmp_new_clusters; // temporarily used to build the connected components.
  tmp_new_clusters.push_back(std::vector<int>(1, element));
  arma::mat A_tmp; 
  // loop over the existing connected components
  for(int cc = 0; cc < left_new_clusters.size(); cc++){
    A_tmp = Submatrix(A_block, tmp_new_clusters[0].size(), left_new_clusters[cc].size(), tmp_new_clusters[0], left_new_clusters[cc]);
    if(any(vectorise(A_tmp) == 1.0)){ // something in left_new_clusters[cc] is adjacent to element so we should combine these clusters
      for(int ix = 0; ix < left_new_clusters[cc].size(); ix++){
        tmp_new_clusters[0].push_back(left_new_clusters[cc][ix]);
      }
    } else{ // left_new_clusters[cc] remains its own distinct component
      tmp_new_clusters.push_back(left_new_clusters[cc]);
    }
    left_new_clusters[cc].clear();
  } // closes loop over elements of left_new_clusters
  left_new_clusters.clear();
  // update left_new_clusters: copy tmp_new_clusters
  for(int cc = 0; cc < tmp_new_clusters.size(); cc++){
    left_new_clusters.push_back(tmp_new_clusters[cc]);
    tmp_new_clusters[cc].clear();
  }
  tmp_new_clusters.clear();
  
  return ;
}

// U is an n x d matrix, where rows represent observations
// means is a d x num_splits matrix, where columns represent cluster centroids
// final_status == false means that all of our runs of kmeans failed. final_status = true means we can proceed

// outside of kmeans_repeat, we can actually return the clusters.
// really we just need it to return means.

void kmeans_repeat(arma::mat &U, arma::mat &means, bool &final_status, double &min_score, int num_splits,int reps){
  int n = U.n_rows; // each row of U is an observation. note that arma::kmeans wants the observations stored as columns.
  int d = U.n_cols; // dimension of the data
  arma::mat tmp_means(d, num_splits); // arma::kmeans returns the centroids stored as column vectors
  bool status = true;
  final_status = false; // if it returns false,
  arma::vec cluster_dist = arma::zeros<arma::vec>(num_splits); // stores distance of each observation to each cluster centroid
  arma::uvec cluster_dist_indices(num_splits); // used to sort the distance from each point to each cluster centroid
  arma::vec cluster_assignment(n); // holds the assignment of each observation
  cluster_assignment.fill(-1); // initialize with really
  arma::mat cluster_means = arma::zeros<arma::mat>(num_splits, d); // each row is the mean of elements in cluster
  arma::vec cluster_counts = arma::zeros<arma::vec>(num_splits);
  
  double score = 0.0;
  min_score = 0.0;
  for(int r = 0; r < reps; r++){
    //Rcpp::Rcout << "r = " << r << endl;
    score = 0.0;
    status = arma::kmeans(tmp_means, U.t(), num_splits, arma::random_subset, 10, false);
    // reset the vectors for determining what's in each sub-cluster
    cluster_assignment.fill(-1);
    cluster_counts.zeros();
    cluster_means.zeros();
    cluster_dist.zeros();
    cluster_dist_indices.zeros();
    
    if(status == false){
      Rcpp::Rcout << "kmeans failed";
    } else{
      final_status = true;
      for(int i = 0; i < n; i++){
        //Rcpp::Rcout << "i = " << i << " : " ;
        cluster_dist.zeros();
        for(int k = 0; k < num_splits; k++){
          cluster_dist(k) = arma::norm(U.row(i).t() - tmp_means.col(k));
        }
        cluster_dist_indices = arma::sort_index(cluster_dist, "ascend");
        cluster_assignment(i) = cluster_dist_indices(0);
        cluster_means.row(cluster_dist_indices(0)) += U.row(i); //
        ++(cluster_counts(cluster_dist_indices(0)));
      }
      
      // at this point, rows of cluster_means contains SUM of the elements in the cluster and not their mean
      for(int k = 0; k < num_splits; k++){
        cluster_means.row(k) /= cluster_counts(k);
      }
      score = 0.0;
      for(int i = 0; i < n; i++){
        score += arma::norm(U.row(i) - cluster_means.row(cluster_assignment(i))) * arma::norm(U.row(i) - cluster_means.row(cluster_assignment(i)));
      }
      if(r == 0){
        min_score = score;
        //Rcpp::Rcout << " r = 0. re-setting min_score = " << min_score << std::endl;
        means = tmp_means;
      }
      else if(score < min_score){
        min_score = score;
        means = tmp_means;
        //Rcpp::Rcout << "[kmeans_repeat]:    r = " << r << " new min score = " << min_score << endl;
      }
    }
  }
  //Rcpp::Rcout << "[kmeans_repeat]: min_score = " << min_score << std::endl;
}

void kmeans_plus_plus(arma::mat &U, arma::mat &means, bool &final_status, double &min_score, int num_splits,int reps){
  
  int n = U.n_rows; // each row of U is an observation. note that arma::kmeans wants the observations stored as columns.
  int d = U.n_cols; // dimension of the data
  arma::mat tmp_means(d, num_splits); // arma::kmeans returns the centroids stored as column vectors
  bool status = true;
  final_status = false; // if it returns false,
  arma::vec cluster_dist = arma::zeros<arma::vec>(num_splits); // stores distance of each observation to each cluster centroid
  arma::uvec cluster_dist_indices(num_splits); // used to sort the distance from each point to each cluster centroid
  arma::vec cluster_assignment(n); // holds the assignment of each observation
  cluster_assignment.fill(-1); // initialize with really
  arma::mat cluster_means = arma::zeros<arma::mat>(num_splits, d); // each row is the mean of elements in cluster
  arma::vec cluster_counts = arma::zeros<arma::vec>(num_splits);
  
  // Form the probability matrix for selecting the starting centroids
  // cdf matrix is what we will compare to
  arma::mat cdf_mat = arma::mat(n,n); // matrix whose rows store the probability of selecting obs j given obs i is the initial see
  arma::vec tmp_dist_vec = arma::zeros<arma::vec>(n);
  double tmp_sum = 0.0;
  arma::vec tmp_cum_sum = arma::zeros<arma::vec>(n);
  for(int i = 0; i < n; i++){
    tmp_dist_vec.zeros();
    tmp_sum = 0.0;
    for(int j = 0; j < n; j++){
      tmp_dist_vec(j) = arma::norm(U.row(i) - U.row(j)) * arma::norm(U.row(i) - U.row(j));
      tmp_sum += tmp_dist_vec(j);
    }
    tmp_dist_vec /= tmp_sum; // normalize tmp_dist_vec to have sum 1
    tmp_cum_sum = arma::cumsum(tmp_dist_vec); // cumulative probabilities.
    cdf_mat.row(i) = tmp_cum_sum.t();
  }
  //Rcpp::Rcout << "[kmpp]: created cdf_mat" << endl;
  
  arma::vec tmp_unif(1);
  arma::vec index_selected(n); // keeps track of which points have been selected as initial points
  index_selected.fill(0);
  int init_index = -1;
  arma::uvec tmp_index(1);
  int index = 0; // used to find the index of the new point
  bool flag = false; // when flag == true, it means we have found a new index for a starting centroid
  double score = 0.0;
  min_score = 0.0;
  
  for(int r = 0; r < reps; r++){
    tmp_means.zeros();
    index_selected.zeros();
    // do the initialization
    // first pick the initial seed point
    //Rcpp::Rcout << "  r = " << r << endl;
    tmp_unif.randu(1);
    init_index = floor(n*tmp_unif(0));
    //Rcpp::Rcout << "    initial index = " << init_index << endl;
    index_selected(index) = 1;
    //Rcpp::Rcout << "    remaining indices: " << endl;
    tmp_means.col(0) = U.row(init_index).t();
    for(int i = 1; i < num_splits; i++){
      flag = false;
      //Rcpp::Rcout << "i = " << i << endl;
      while(flag == false){
        tmp_unif.randu(1);
        tmp_index = arma::find(cdf_mat.row(init_index) > tmp_unif(0), 1, "first"); // find first element in cdf_mat.row(init_index) exceeding tmp_unif(0)
        index = tmp_index(0);
        //Rcpp::Rcout << "tmp_unif = " << tmp_unif << "  cdf_mat(init_index, index) = " << cdf_mat(init_index, index) << endl;
        if(index_selected(index) == 0){
          flag = true;
          index_selected(index) = 1;
        } else {
          flag = false;
        }
        //Rcpp::Rcout << " " << index ;
      } // closes while loop
      tmp_means.col(i) = U.row(index).t();
      //Rcpp::Rcout << " " << index;
    } // closes loop that initializes tmp_mean
    score = 0.0;
    status = arma::kmeans(tmp_means, U.t(), num_splits, arma::keep_existing, 10, false);
    // reset the vectors for determining what's in each sub-cluster
    cluster_assignment.fill(-1);
    cluster_counts.zeros();
    cluster_means.zeros();
    cluster_dist.zeros();
    cluster_dist_indices.zeros();
    
    if(status == false){
      Rcpp::Rcout << "kmeans failed";
    } else{
      final_status = true;
      for(int i = 0; i < n; i++){
        //Rcpp::Rcout << "i = " << i << " : " ;
        cluster_dist.zeros();
        for(int k = 0; k < num_splits; k++){
          cluster_dist(k) = arma::norm(U.row(i).t() - tmp_means.col(k));
        }
        cluster_dist_indices = arma::sort_index(cluster_dist, "ascend");
        cluster_assignment(i) = cluster_dist_indices(0);
        cluster_means.row(cluster_dist_indices(0)) += U.row(i); //
        ++(cluster_counts(cluster_dist_indices(0)));
      }
      // at this point, rows of cluster_means contains SUM of the elements in the cluster and not their mean
      for(int k = 0; k < num_splits; k++){
        cluster_means.row(k) /= cluster_counts(k);
      }
      score = 0.0;
      for(int i = 0; i < n; i++){
        score += arma::norm(U.row(i) - cluster_means.row(cluster_assignment(i))) * arma::norm(U.row(i) - cluster_means.row(cluster_assignment(i)));
      }
      if(r == 0){
        min_score = score;
        //Rcpp::Rcout << " r = 0. re-setting min_score = " << min_score << std::endl;
        means = tmp_means;
      }
      else if(score < min_score){
        min_score = score;
        means = tmp_means;
        //Rcpp::Rcout << "r = " << r << " new min score = " << min_score << endl;
      }
    } // closes if/else checking whether kmeans failed
    
    
    //Rcpp::Rcout << endl;
  } // closes loop over r
  
  
  
  
}








