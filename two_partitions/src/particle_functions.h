/*
 * particle_functions.h
 *
 */

int Particle_Equal(Particle *particle1, Particle *particle2);
double Entropy(int current_l, Particle* candidate_particle, std::vector<LPParticle> particle_set, std::vector<double> w);
double VI_Avg(int current_l, Particle* candidate_particle, std::vector<LPParticle> Particle_Set);
void Resampling(std::vector<LPParticle>& Particle_Set, std::vector<double>& w);
