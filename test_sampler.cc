#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <streambuf>

#include <vector>
#include <set>

#include "gsl/gsl_sf_gamma.h"

#include "threebody_integrals.h"
#include "microcanonical_sampler.h"
#include "hydro_cells.h"
#include "main.h"

#include "smash/particles.h"
#include "smash/random.h"
#include "smash/decaymodes.h"
#include "smash/angles.h"
#include "smash/constants.h"
#include "smash/quantumnumbers.h"
#include "smash/pow.h"
#include "smash/isoparticletype.h"

using namespace smash;

void initialize_random_number_generator() {
  // Seed with a truly random 63-bit value, if possible
  std::random_device rd;
  static_assert(std::is_same<decltype(rd()), uint32_t>::value,
               "random_device is assumed to generate uint32_t");
  uint64_t unsigned_seed = (static_cast<uint64_t>(rd()) << 32) |
                            static_cast<uint64_t>(rd());
  // Discard the highest bit to make sure it fits into a positive int64_t
  int64_t seed = static_cast<int64_t>(unsigned_seed >> 1);
  random::set_seed(seed);
}

void load_smash_particles() {
  std::string smash_dir(std::getenv("SMASH_DIR"));
  if (smash_dir == "") {
    throw std::runtime_error("Failed to load SMASH particle types."
                             " SMASH_DIR is not set.");
  }
  std::cout << "Loading SMASH particle types and decay modes" << std::endl;
  std::ifstream particles_input_file(smash_dir + "/input/particles.txt");
  std::stringstream buffer;
  if (particles_input_file) {
    buffer << particles_input_file.rdbuf();
    ParticleType::create_type_list(buffer.str());
  } else {
    std::cout << "File with SMASH particle list not found." << std::endl;
  }
  std::ifstream decaymodes_input_file(smash_dir + "/input/decaymodes.txt");
  if (decaymodes_input_file) {
    buffer.clear();
    buffer.str(std::string());
    buffer << decaymodes_input_file.rdbuf();
    DecayModes::load_decaymodes(buffer.str());
    ParticleType::check_consistency();
  } else {
    std::cout << "File with SMASH decaymodes not found." << std::endl;
  }
}

int type_count(const MicrocanonicalSampler::SamplerParticleList &particles,
               const ParticleTypePtr t) {
  int cnt = 0;
  for (const auto &particle : particles) {
    if (particle.type == t) {
      cnt++;
    }
  }
  return cnt;
}

bool is_sampled_type(const smash::ParticleTypePtr t) {
  return t->is_hadron() && t->mass() < 1.0;
}

int main() {
  initialize_random_number_generator();
  load_smash_particles();

  HyperSurfacePatch hyper("../hydro_cells.dat", is_sampled_type);
  std::cout << hyper << std::endl;

  MicrocanonicalSampler sampler(is_sampled_type, 1);
  sampler.initialize(hyper);

  std::cout << "Warming up." << std::endl;
  constexpr int N_warmup = 10000;
  for (int i = 0; i < N_warmup; ++i) {
    // std::cout << i << std::endl;
    sampler.one_markov_chain_step(hyper);
  }
  std::cout << "Finished warming up." << std::endl;

  constexpr int N_decorrelate = 100;
  constexpr int N_printout = 10000;
  for (int j = 0; j < N_printout; j++) {
    for (int i = 0; i < N_decorrelate; ++i) {
      sampler.one_markov_chain_step(hyper);
    }
/*
    for (const auto &part : sampler.particles()) {
      std::cout << part.momentum.abs() << " " << part.momentum.x0() << std::endl;
    }
*/
  }

  MicrocanonicalSampler::QuantumNumbers cons(sampler.particles());
  std::cout << "Final momentum: " << cons.momentum << std::endl;
}
