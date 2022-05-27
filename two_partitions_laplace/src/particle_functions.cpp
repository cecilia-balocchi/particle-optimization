/*
 * particle_functions.cpp
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
#include <random>
#include <algorithm>
#include "partition.h"
#include "particle.h"
#include "partition_functions.h"
#include "particle_functions.h"
#include "various_functions.h"
using namespace std;
extern double lambda;
extern double xi;

// None of the methods here actually need any of this data
// instead, they call methods which use the data and we already have
// the appropriate extern calls in the corresponding headers

extern arma::mat Y;
//extern arma::mat X;
//extern arma::mat A_block; 


int Particle_Equal(Particle *particle1, Particle *particle2){
  int flag1, flag2;
  flag1 = Partition_Equal(particle1->partition_A,particle2->partition_A);
  flag2 = Partition_Equal(particle1->partition_B,particle2->partition_B);
  return (flag1*flag2);
}

double Entropy(int current_l, Particle* candidate_particle, std::vector<LPParticle> Particle_Set, std::vector<double> w){
  int L = Particle_Set.size();
  // need to loop over to extract the unique particles
  std::vector<LPParticle> unik_particles;
  std::vector<double> p_star;

  unik_particles.push_back(candidate_particle);
  p_star.push_back(w[current_l]);

  // in a sense, we are replacing Particle_Set[current_l] with candidate_particle
  // by adding it to the unik_particles vector

  int num_unik_particles = 1;
  int counter = 0;
  for(int l = 0; l < L; l++){ // loop over current particle set
    counter = 0;
    if(l != current_l){
      for(int ul = 0; ul < num_unik_particles; ul++){
        if(Particle_Equal(Particle_Set[l], unik_particles[ul]) == 1){
          // l^th particle is equal to the ul^th unique particle
          //std::cout << "particle " << l << " is equal to unik particle " << ul << std::endl;
          p_star[ul] += w[l]; // update p_star
          break;
        } else {
          counter++;
        }
      }
      if(counter == num_unik_particles){
        // we have found a new unique particle
        unik_particles.push_back(Particle_Set[l]);
        p_star.push_back(w[l]);
        num_unik_particles++;
      }
    }
  }
  double entropy = 0.0;
  double test;
  for(int ul = 0; ul < num_unik_particles; ul++){
    test = p_star[ul] * log(p_star[ul]);
    if(!std::isnan(test)){
      entropy += test; 
    }
  }
  return -1.0 * entropy;
}


double VI_Avg(int current_l, Particle* candidate_particle, std::vector<LPParticle> Particle_Set){
  int L = Particle_Set.size();
  if(L == 1){
    return 0.0;
  }
  double tot1 = 0.0;
  double tot2 = 0.0;
  double dist;
  for(int l = 0; l < L; l++){
    for(int ll = l+1; ll < L; ll++){
      if(l == current_l){
        dist = VI_Loss(candidate_particle->partition_A, Particle_Set[ll]->partition_A);
      } else if(ll == current_l){
        dist = VI_Loss(Particle_Set[l]->partition_A, candidate_particle->partition_A);
      } else{
        dist = VI_Loss(Particle_Set[l]->partition_A, Particle_Set[ll]->partition_A);
      }
      tot1 = tot1 + dist;
    }
  }
  for(int l = 0; l < L; l++){
    for(int ll = l+1; ll < L; ll++){
      if(l == current_l){
        dist = VI_Loss(candidate_particle->partition_B, Particle_Set[ll]->partition_B);
      } else if(ll == current_l){
        dist = VI_Loss(Particle_Set[l]->partition_B, candidate_particle->partition_B);
      } else{
        dist = VI_Loss(Particle_Set[l]->partition_B, Particle_Set[ll]->partition_B);
      }
      tot2 = tot2 + dist;
    }
  }
  return (tot1 + tot2)/(L*(L-1)/2);
}


void Resampling(std::vector<LPParticle>& Particle_Set, std::vector<double>& w){
  int L = Particle_Set.size();
  std::vector<LPParticle> old_Particle_Set(L);
  std::vector<double> old_w(L); 
  double w_threshold = 1e-50;
  int replacement;
  int number;
  double w_tot;
  std::default_random_engine generator;

  // first we need to save the old Particle_Set
  for(int l = 0; l < L; l++){
    old_Particle_Set[l] = new Particle(Particle_Set[l]);
    old_w[l] = w[l];
  }
  // resample the particles according to the weights
  // only if the current weights are less than a threshold
  w_tot = 0.0;
  replacement = 0;
  std::discrete_distribution<int> distribution (old_w.begin(), old_w.end());
  for(int l = 0; l < L; l++){
    if(old_w[l] < w_threshold){
      number = distribution(generator);
      Rcpp::Rcout << "Number sampled: " << number << endl;
      delete Particle_Set[l];
      Particle_Set[l] = new Particle(old_Particle_Set[(unsigned) number]);
      w[l] = old_w[(unsigned) number];
      w_tot += w[l];
      replacement++;
    } else {
      number = l; // otherwise I just use the previous one
      delete Particle_Set[l];
      Particle_Set[l] = new Particle(old_Particle_Set[(unsigned) number]);
      w[l] = old_w[(unsigned) number];
      w_tot += w[l];
    }
  }
  for(int l = 0; l < L; l++){
    delete old_Particle_Set[l];
  }
  for(int l = 0; l < L; l++){
    w[l] = w[l]/w_tot;
  }
  Rcpp::Rcout << "[particle_spatial]: New weights" << endl;
  for(int l = 0; l < L; l++){
    Rcpp::Rcout << w[l] << " ";
  }
  Rcpp::Rcout << endl << "[particle_spatial]: Number of replacements: " << replacement << endl << endl;
}


