/*
 * partition_functions.cpp
 *
 */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <vector>
#include <algorithm>
#include "partition.h"
#include "various_functions.h"
using namespace std;
using namespace arma;
extern double lambda;
extern double xi;

// None of the methods here actually need any of this data
// instead, they call methods which use the data and we already have
// the appropriate extern calls in the corresponding headers

extern arma::mat Y;
//extern arma::mat X;
extern arma::mat A_block; 


// Compare two partitions and see if they are equal
int Partition_Equal(Partition *partition1, Partition *partition2){
  int flag = 1;
    // simpler to just compare pairwise allocation
  for(int i = 0; i < partition1->nObs; i++){
  	for(int j = 0; j < partition1->nObs; j++){
  	  if(partition1->pairwise_assignment[i][j] != partition2->pairwise_assignment[i][j]){
  		flag = 0;
  		break;
  	  }
  	}
  	if(flag == 0) break;
  }
  return flag;
}

double beta_bar(Partition *partition, int k){
  arma::vec beta_cl(partition->cluster_config[k]);
  for(int i = 0; i < partition->cluster_config[k]; i++){
    beta_cl(i) = partition->beta_hat[ partition->clusters[k][i] ];
  }
  return arma::mean(beta_cl);
}

// function to compute entropy
// when we replace the (current_l)^th particle with the candidate_clusters
// NOT USED IN THE TWO PARTITION CODE! - CHECK THE PARTICLE VERSION!
/*
double Entropy(unsigned current_l, Partition* candidate_particle, std::vector<LPPartition> particle_set, std::vector<double> w){
  unsigned L = particle_set.size();
  // need to loop over to extract the unique partitions
  std::vector<LPPartition> unik_particles;
  std::vector<double> p_star;

  unik_particles.push_back(candidate_particle);
  p_star.push_back(w[current_l]);

  // in a sense, we are replacing particle_set[current_l] with candidate_particle
  // by adding it to the unik_particles vector


  int num_unik_particles = 1;
  int counter = 0;
  for(unsigned l = 0; l < L; l++){ // loop over current particle set
  	counter = 0;
  	if(l != current_l){
    	//std::cout << "[getUnik]: l = " << l << std::endl;
  	  for(unsigned ul = 0; ul < num_unik_particles; ul++){
  	    if(Partition_Equal(particle_set[l], unik_particles[ul]) == 1){
    	    // l^th partition is equal to the ul^th unique partition
    	    //std::cout << "particle " << l << " is equal to unik particle " << ul << std::endl;
  	      p_star[ul] += w[l]; // update p_star
  	      break;
        } else {
          counter++;
        }
  	  }
  	  //std::cout << "[getUnik]: counter = " << counter << std::endl;
  	  if(counter == num_unik_particles){
  	    //std::cout << "we found a new unique particle!" << std::endl;
  	    //particle_set[l]->Print_Partition();
  	    // we have found a new unique particle
  	    unik_particles.push_back(particle_set[l]);
  	    p_star.push_back(w[l]);
  	    num_unik_particles++;
      }
    }
  }
  double entropy = 0.0;
  double test;
  for(unsigned ul = 0; ul < num_unik_particles; ul++){
	  test = p_star[ul] * log(p_star[ul]);
	  if(!std::isnan(test)){
      entropy += test;
    }
    // entropy += p_star[ul] * log(p_star[ul]);
  }
  return -1.0 * entropy;

}
*/

double total_log_prior(LPPartition partition){
  double log_prior = partition->log_cohesion;
  for(int k = 0; k < partition->K; k++){
    log_prior += partition->log_prior[k];
  }
  return log_prior;
}

double Binder_Loss(LPPartition partition1, LPPartition partition2){
  int n = partition1->nObs;
  int K1 = partition1->K;
  int K2 = partition2->K;
  // need to get the matrix of counts n_ij which counts the number of indices that belong to cluster i in partition1 and cluster j in partition2
  arma::mat counts(K1, K2, fill::zeros);

  for(int i = 0; i < n; i++){
  counts(partition1->cluster_assignment[i], partition2->cluster_assignment[i])++;
  }
  //counts.print();

  double loss = 0.0;
  // Binder loss is quadratic in the elements of the matrix of counts
  // Need sum of squares of row sums + sum of squared columns sums and then subtract sum of squares of elements of counts
  // row sums of counts: just the configuration of partition1
  // column sums of counts: just the
  for(int k = 0; k < K1; k++){
  loss += 0.5 * ( (double) partition1->cluster_config[k])*( (double) partition1->cluster_config[k]);
  }
  for(int k = 0; k < K2; k++){
  loss += 0.5 * ( (double) partition2->cluster_config[k]) * ( (double) partition2->cluster_config[k]);
  }
  loss -= accu(pow(counts, 2));

  // Wade and Ghahramani consider a scaled Binder loss and we will as well
  // Property 4 in their paper asserts that this binder loss will always be between 0 and 1 - 1/n
  loss *= 2.0/( (double) n* (double) n);
  return loss;

}

double  VI_Loss(LPPartition partition1, LPPartition partition2){
  int n = partition1->nObs;
  int K1 = partition1->K;
  int K2 = partition2->K;
  // need to get the matrix of counts n_ij which counts the number of indices that belong to cluster i in partition1 and cluster j in partition2
  mat counts(K1, K2, fill::zeros);

  for(int i = 0; i < n; i++){
    counts(partition1->cluster_assignment[i], partition2->cluster_assignment[i])++;
  }
  double loss = 0.0;
  for(int k = 0; k < K1; k++){
    loss += ( (double) partition1->cluster_config[k])/( (double) n) * log( ((double) partition1->cluster_config[k])/( (double) n));
  }
  for(int k = 0; k < K2; k++){
    loss += ( (double) partition2->cluster_config[k])/( (double) n) * log(((double) partition2->cluster_config[k])/((double) n));
  }
  for(int k1 = 0; k1 < K1; k1 ++){
    for(int k2 = 0; k2 < K2; k2++){
      if(counts(k1, k2) != 0){ // 0 * log(0) = 0 so if any counts
        loss -= 2.0 * ( (double) counts(k1, k2)) / ((double) n) * log( ((double) counts(k1, k2))/ ( (double) n));
      }
    }
  }
  return loss;
}

// double VI_Avg(unsigned current_l, Partition* candidate_particle, std::vector<LPPartition> particle_set){
//   unsigned L = particle_set.size();
//   double tot = 0.0;
//   double dist;
//   for(int l = 0; l < L; l++){
//     for(int ll = l+1; ll < L; ll++){
//       if(l == current_l){
//         dist = VI_Loss(candidate_particle, particle_set[ll]);
//       } else if(ll == current_l){
//         dist = VI_Loss(particle_set[l], candidate_particle);
//       } else{
//         dist = VI_Loss(particle_set[l], particle_set[ll]);
//       }
//       tot = tot + dist;
//     }
//   }
//   return tot/(L*(L-1)/2);
// }

double Binder_DPP(int current_l, Partition* candidate_particle, std::vector<LPPartition> particle_set){

  // The first step is to get all of the unique particles when we replace the current_l^th particle with candidate particle
  int L = particle_set.size();
  // need to loop over to extract the unique partitions
  std::vector<LPPartition> unik_particles;
  std::vector<double> p_star;

  unik_particles.push_back(candidate_particle);
  int num_unik_particles = 1;
  int counter = 0;
  for(int l = 0; l < L; l++){ // loop over current particle set
  counter = 0;
  if(l != current_l){
    for(int ul = 0; ul < num_unik_particles; ul++){
      if(Partition_Equal(particle_set[l], unik_particles[ul]) == 1){
      // l^th partition is equal to the ul^th unique partition
      //std::cout << "particle " << l << " is equal to unik particle " << ul << std::endl;
      //p_star[ul] += w[l]; // update p_star
        break;
        } else {
        counter++;
        }
    }
    if(counter == num_unik_particles){
      unik_particles.push_back(particle_set[l]);
      num_unik_particles++;
      }
    }
  }

  mat kernel(num_unik_particles, num_unik_particles, fill::zeros);
  double dist = 0.0;
  double kernel_log_det = 0.0;
  double kernel_log_det_sgn = 0.0;
  for(int l = 0; l < num_unik_particles; l++){
    for(int ll = 0; ll < num_unik_particles; ll++){
      dist = Binder_Loss(unik_particles[l], unik_particles[ll]);
      kernel(l,ll) = exp(-0.5 * dist * dist);
    }
  }
  log_det(kernel_log_det, kernel_log_det_sgn, kernel);

/*
  double dist = 0.0;
  if(num_unik_particles == 1){
  dist = 0.0; // everything is equal to everything else
  } else {
  // now loop over the unique particles and compute pairwise distances
    for(unsigned l = 0; l < num_unik_particles - 1; l++){
    for(unsigned ll = l+1; ll < num_unik_particles; ll++){
      dist += Binder_Loss(unik_particles[l], unik_particles[ll]);
    }
    }
  }
*/
  return kernel_log_det;
}
double VI_DPP(int current_l, Partition* candidate_particle, std::vector<LPPartition> particle_set){
  // The first step is to get all of the unique particles when we replace the current_l^th particle with candidate particle
  int L = particle_set.size();
  // need to loop over to extract the unique partitions
  std::vector<LPPartition> unik_particles;
  std::vector<int> particle_counts;


  unik_particles.push_back(candidate_particle);
  particle_counts.push_back(1);
  int num_unik_particles = 1;
  int counter = 0;
  for(int l = 0; l < L; l++){ // loop over current particle set
  counter = 0;
  if(l != current_l){
    for(int ul = 0; ul < num_unik_particles; ul++){
      if(Partition_Equal(particle_set[l], unik_particles[ul]) == 1){
      // l^th partition is equal to the ul^th unique partition
      //std::cout << "particle " << l << " is equal to unik particle " << ul << std::endl;
      //p_star[ul] += w[l]; // update p_star
        particle_counts[ul]++;
        break;
        } else {
        counter++;
        }
    }
    if(counter == num_unik_particles){
      unik_particles.push_back(particle_set[l]);
      particle_counts.push_back(1);
      num_unik_particles++;
      }
    }
  }
  if(num_unik_particles == 1){
    return(-pow(10.0, 3.0));
  } else{

    mat kernel(num_unik_particles, num_unik_particles, fill::zeros);
    double dist = 0.0;
    double kernel_log_det = 0.0;
    double kernel_log_det_sgn = 0.0;
    for(int l = 0; l < num_unik_particles; l++){
      for(int ll = 0; ll < num_unik_particles; ll++){
        dist = VI_Loss(unik_particles[l], unik_particles[ll]);
        //kernel(l,ll) = exp(-1.0 * dist); // large distance means the particles are not similar, hence entry in kernel matrix should be small
        kernel(l,ll) = exp(-1.0 * dist /pow( particle_counts[l] * particle_counts[ll], 2.0));
      }
    }
    log_det(kernel_log_det, kernel_log_det_sgn, kernel);
    return(kernel_log_det);
  }
/*
  double dist = 0.0;
  if(num_unik_particles == 1){
  dist = 0.0; // everything is equal to everything else
  } else {
  // now loop over the unique particles and compute pairwise distances
    for(unsigned l = 0; l < num_unik_particles - 1; l++){
    for(unsigned ll = l+1; ll < num_unik_particles; ll++){
      dist += VI_Loss(unik_particles[l], unik_particles[ll]);
    }
    }
  }
  return dist;
*/
}

