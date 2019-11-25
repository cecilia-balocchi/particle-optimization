/*
 * partition_functions.h
 *
 */
#include <RcppArmadillo.h>

#include <stdio.h>


int Partition_Equal(Partition *partition1, Partition *partition2);
double beta_bar(Partition *partition, int k);
// double Entropy(unsigned current_l, Partition* candidate_particle, std::vector<LPPartition> particle_set, std::vector<double> w);
double total_log_prior(LPPartition partition);
double Binder_Loss(LPPartition partition1, LPPartition partition2);
double VI_Loss(LPPartition partition1, LPPartition partition2);
// double VI_Avg(unsigned current_l, Partition* candidate_particle, std::vector<LPPartition> particle_set);

double Binder_DPP(int current_l, Partition* candidate_particle, std::vector<LPPartition> particle_set);
double VI_DPP(int current_l, Partition* candidate_particle, std::vector<LPPartition> particle_set);
