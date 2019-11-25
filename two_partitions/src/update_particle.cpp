/*
 * update_particle.cpp
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
#include "particle.h"
#include "various_functions.h"
#include "partition_functions.h"
#include "particle_functions.h"
#include "update_particle.h"
extern arma::mat A_block;

using namespace std;

void update_w(std::vector<LPParticle> Particle_Set, std::vector<double>& w, double lambda, double xi){
  // First we identify the unik partitions
  double max_log_post = 0.0;
  double tmp_log_post = 0.0;
  double tmp_norm = 0.0;
  double tmp_p_star = 0.0;
  int L = Particle_Set.size();

  // double toll = pow (10.0, -10.0);

  // need to loop over to extract the unique partitions
  std::vector<LPParticle> unik_particles;
  unik_particles.push_back(Particle_Set[0]);

  // will tell us which unik particle each particle in the particle set is
  std::vector<int> particle_assignment(L); 
  particle_assignment[0] = 0;

  std::vector<double> p_star;
  std::vector<double> log_post;

  log_post.push_back(Particle_Set[0]->Total_Log_Post());
  p_star.push_back(0);

  max_log_post = Particle_Set[0]->Total_Log_Post();

  int num_unik_particles = 1;
  int counter;
  for(int l = 1; l < L; l++){ // loop over current particle set
    counter = 0;
    for(int ul = 0; ul < num_unik_particles; ul++){
      if(Particle_Equal(Particle_Set[l], unik_particles[ul]) == 1){
        // l^th partition is equal to the ul^th unique partition
        particle_assignment[l] = ul;
        break;
      } else {
        counter++;
      }
    }
    if(counter == num_unik_particles){
      // we have found a new unique particle
      unik_particles.push_back(Particle_Set[l]);
      particle_assignment[l] = num_unik_particles; // the labels are off-set by 1
      p_star.push_back(0.0); // for now, we will populate p_star with 0's
      tmp_log_post = Particle_Set[l]->Total_Log_Post();
      log_post.push_back(tmp_log_post);
      if(tmp_log_post > max_log_post){
        max_log_post = tmp_log_post;
      }
      num_unik_particles++;
    }
  }

  // how many particles are equal to this unique particle
  int particle_counts[num_unik_particles];
  int tmp_count;
  for(int ul = 0; ul < num_unik_particles; ul++){
    tmp_count = 0;
    for(int l = 0; l < L; l++){
      if(particle_assignment[l] == ul){
        tmp_count++;
      }
    }
    particle_counts[ul] = tmp_count;
  }

  for(int ul = 0; ul < num_unik_particles; ul++){
    tmp_log_post = log_post[ul] - max_log_post;
    tmp_p_star = exp(1/lambda * tmp_log_post); // introduce the 1/lambda
    tmp_norm += tmp_p_star;
    p_star[ul] = tmp_p_star;
  }
  for(int ul = 0; ul < num_unik_particles; ul++){
    p_star[ul] /= tmp_norm;
  }

  // threshold if some weights are too small
  // tmp_norm = 0.0;
  // for(int ul = 0; ul < num_unik_particles; ul++){
  //   if(p_star[ul] < toll){
  //     p_star[ul] = toll;
  //   }
  //   tmp_norm += p_star[ul];
  // }
  // for(int ul = 0; ul < num_unik_particles; ul++){
  //   p_star[ul] /= tmp_norm;
  // }

  for(int l = 0; l < L; l++){
    w[l] = (double) p_star[particle_assignment[l]]/particle_counts[particle_assignment[l]];
  }
  return;
}



