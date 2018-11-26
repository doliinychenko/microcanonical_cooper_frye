#ifndef MICROCANONICAL_SAMPLING_MAIN_H
#define MICROCANONICAL_SAMPLING_MAIN_H

#include "smash/particles.h"
#include "smash/forwarddeclarations.h"
#include "smash/quantumnumbers.h"

#include "threebody_integrals.h"

double compute_R2(double srts, double m1, double m2);
double compute_sum_R3(const smash::ParticleTypePtrList& sampled_types,
                      ThreeBodyIntegrals& three_body_int,
                      const smash::QuantumNumbers &cons);
double compute_sum_R2(const smash::ParticleTypePtrList& sampled_types,
                      const smash::QuantumNumbers &cons);
int type_count(const smash::ParticleList &particles,
               const smash::ParticleTypePtr t);
void sample_3body_phase_space(double srts,
                              smash::ParticleData &a,
                              smash::ParticleData &b,
                              smash::ParticleData &c);
bool quantum_numbers_match(const smash::ParticleTypePtrList& a,
                           const smash::QuantumNumbers& qn);
void random_two_to_three(smash::ParticleList &particles,
                         const smash::ParticleTypePtrList& sampled_types,
                         ThreeBodyIntegrals& three_body_int,
                         double V);
void random_three_to_two(smash::ParticleList &particles,
                         const smash::ParticleTypePtrList& sampled_types,
                         ThreeBodyIntegrals& three_body_int,
                         double V);

void renormalize_momenta(smash::ParticleList& plist,
                         smash::FourVector required_4mom);
void initialize_random_number_generator();
void load_smash_particles();

#endif  // MICROCANONICAL_SAMPLING_MAIN_H

