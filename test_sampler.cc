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
               const ParticleTypePtr t, size_t cell_number) {
  int cnt = 0;
  for (const auto &particle : particles) {
    if (particle.type == t && cell_number == particle.cell_index) {
      cnt++;
    }
  }
  return cnt;
}

bool is_sampled_type(const smash::ParticleTypePtr t) {
  return t->is_hadron() && t->mass() < 0.2 && t->charge() == 0;
}

int main() {
  initialize_random_number_generator();
  load_smash_particles();
  constexpr bool quantum_statistics = false;

  HyperSurfacePatch hyper("../hydro_cells.dat", is_sampled_type,
                          quantum_statistics);
  std::cout << hyper << std::endl;

  MicrocanonicalSampler sampler(is_sampled_type, 0, quantum_statistics);
  sampler.initialize(hyper);

  std::cout << "Warming up." << std::endl;
  constexpr int N_warmup = 100000;
  for (int i = 0; i < N_warmup; ++i) {
    // std::cout << i << std::endl;
    sampler.one_markov_chain_step(hyper);
  }
  std::cout << "Finished warming up." << std::endl;

  ParticleTypePtrList sampled_types;
  for (const ParticleType &ptype : ParticleType::list_all()) {
    if (is_sampled_type(&ptype)) {
      sampled_types.push_back(&ptype);
    }
  }
  for (const ParticleTypePtr t : sampled_types) {
    std::cout << t->name() << "  ";
  }
  std::cout << std::endl;

  constexpr int N_decorrelate = 200;
  constexpr int N_printout = 100000;
  for (int j = 0; j < N_printout; j++) {
    for (int i = 0; i < N_decorrelate; ++i) {
      sampler.one_markov_chain_step(hyper);
    }
    for (const ParticleTypePtr t : sampled_types) {
      std::cout << " " << type_count(sampler.particles(), t, 0);
    }
    std::cout << std::endl;
  }

  MicrocanonicalSampler::QuantumNumbers cons(sampler.particles());
  assert((cons.momentum - hyper.pmu()).abs() < 1.e-6);
  assert(cons.B == hyper.B());
  assert(cons.S == hyper.S());
  assert(cons.Q == hyper.Q());
  // std::cout << "Final momentum diff: " << cons.momentum - hyper.pmu() << std::endl;
}
