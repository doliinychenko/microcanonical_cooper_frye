#ifndef MICROCANONICAL_SAMPLING_MAIN_H
#define MICROCANONICAL_SAMPLING_MAIN_H

#include "smash/particletype.h"

#include "microcanonical_sampler.h"

int type_count(const MicrocanonicalSampler::SamplerParticleList &particles,
               const smash::ParticleTypePtr t, size_t cell_number);
bool is_sampled_type(const smash::ParticleTypePtr t);
void initialize_random_number_generator();
void load_smash_particles();

#endif  // MICROCANONICAL_SAMPLING_MAIN_H

