#ifndef MICROCANONICAL_SAMPLING_MAIN_H
#define MICROCANONICAL_SAMPLING_MAIN_H

#include "smash/particles.h"
#include "smash/forwarddeclarations.h"

int type_count(const smash::ParticleList &particles,
               const smash::ParticleTypePtr t);
void initialize_random_number_generator();
void load_smash_particles();

#endif  // MICROCANONICAL_SAMPLING_MAIN_H

