/*
 * partition.h
 *
 */


#ifndef PARTITION_H_
#define PARTITION_H_

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <string.h>
using namespace std;

typedef class Partition* LPPartition;
class Partition
{
	//data members
public:
	int nObs; // number of indices
	int K; // number of clusters
	int* cluster_config; // sizes of each cluster
  int** clusters; // the actual clusters  // might be easier to use vectors.
	int* cluster_assignment; // size nObs, tells us which cluster each index belongs to
	int** pairwise_assignment; // array of pairwise assignments
	double* log_prior; // size K. holds the log-prior evaluated on each cluster
	double* beta_hat; // will have size nObs. holds the point estimates of beta in each blockgroup
	double log_cohesion;

  double* log_det_Omegay;
  double* quad_forms;
  // methods
public:
	Partition(); // constructor
	Partition(LPPartition initial_partition, bool A_or_B); // constructor. create a new partition from an old one
	~Partition(); // destructor
public:
	void Copy_Partition(LPPartition initial_partition, bool A_or_B);
	void Initialize_Partition(int n, bool A_or_B);
  void Initialize_Partition(int n, Rcpp::List gamma_init, bool A_or_B);
  void Initialize_Partition_nclusters(int n, bool A_or_B);
	// void Initialize_Partition2(int id);
 //  void Initialize_Partition3(int n, int k, int* cl_sizes, int* cl_limits);
  void Initialize_Partition_FromFile(int n, string file);
	void get_pairwise();
	void log_pi_ep(int cluster_id);
  void log_pi_epy(int cluster_id);
  void log_cohes(double a_cohes, double b_cohes);
  void log_Dahl_cohes(int* sigma);
  // void get_prior_k(int cluster_id, int method,int* sigma);
  void get_prior_all(int method, int* sigma);
  void get_prior_onlylogprior(int cluster_id, int method);
  void get_prior_onlycohesion(int method, int* sigma);
  void Get_Likelihood(bool A_or_B);
  void get_likelihood(bool A_or_B, int cluster_id);
	void beta_postmean(int cluster_id);
	void Print_Partition();
  void Print_Partition_Short();
  void Print_Partition_ToFile(string file);
  void Read_Partition_FromFile(string file, int n);
  void Print_Means();
	void Split(int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2, bool A_or_B); // split cluster split_k into two parts: new_cluster1, new_cluster2
  void Split(int split_k, std::vector<int> new_cluster1, std::vector<int> new_cluster2, int size1, int size2, bool A_or_B);
  void KSplit(int split_k, int num_splits, std::vector<std::vector<int> > indices, std::vector<int> ns, bool A_or_B); // split cluster split_k into num_split parts
	void Merge(int k_1, int k_2, bool A_or_B); // merge cluster max(k_1, k_2) into min(k_1, k_2)

	void Split_and_Merge(int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2, int k_star_1, int k_star_2, bool A_or_B);

	void Modify(int cl_ind, bool A_or_B);

	// void Find_Splits(int cluster_id, int **index1, int **index2, int &n1, int &n2, bool A_or_B);
};



#endif /* PARTITION_H_ */
