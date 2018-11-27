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

int type_count(const ParticleList &particles, const ParticleTypePtr t) {
  int cnt = 0;
  for (const auto &particle : particles) {
    if (particle.type() == *t) {
      cnt++;
    }
  }
  return cnt;
}

int main() {
  initialize_random_number_generator();
  load_smash_particles();

  // Select Particle types to sample
  std::cout << "Selecting types to sample" << std::endl;
  ParticleTypePtrList sampled_types;
  for (const ParticleType &ptype : ParticleType::list_all()) {
    if (ptype.is_hadron() && ptype.mass() < 1.0) {
      sampled_types.push_back(&ptype);
    }
  }

  HyperSurfacePatch hyper("../hydro_cells.dat",
      [&](const ParticleTypePtr t) {
        return t->is_hadron() && t->mass() < 1.0;
      });
  std::cout << hyper << std::endl;

  MicrocanonicalSampler sampler([&](const ParticleTypePtr t) {
    return t->is_hadron() && t->mass() < 1.0;
  }, 1);

/*
  std::cout << "Constructing initial configuration." << std::endl;
  const double E_tot = 500.0;  // GeV
  const double V = 1000.0;    // fm^3

  ParticleList particles;
  ParticleData pi0{ParticleType::find(pdg::pi_z)};
  constexpr int Npart_init = 20;
  for (int i = 0; i < Npart_init; i++) {
    particles.push_back(pi0);
  }

  Angles phitheta;

  for (ParticleData &data : particles) {
    const double momentum_radial = 0.1;
    phitheta.distribute_isotropically();
    data.set_4momentum(data.type().mass(), phitheta.threevec() * momentum_radial);
    data.set_4position(FourVector());
  }

  FourVector required_total_4mom(E_tot, 0.0, 0.0, 0.0);
  renormalize_momenta(particles, required_total_4mom);

  std::cout << "Warming up." << std::endl;
  constexpr int N_warmup = 10000;
  for (int i = 0; i < N_warmup; ++i) {
    // std::cout << i << std::endl;
    if (random::uniform_int(0, 1) == 1) {
      random_two_to_three(particles, sampled_types, three_body_int, V);
    } else {
      random_three_to_two(particles, sampled_types, three_body_int, V);
    }
  }
  std::cout << "Finished warming up." << std::endl;
  for (const ParticleTypePtr t : sampled_types) {
    std::cout << t->name() << "  ";
  }
  std::cout << std::endl;
  constexpr int N_decorrelate = 100;
  constexpr int N_printout = 10000;
  for (int j = 0; j < N_printout; j++) {
    for (int i = 0; i < N_decorrelate; ++i) {
      if (random::uniform_int(0, 1) == 1) {
        random_two_to_three(particles, sampled_types, three_body_int, V);
      } else {
        random_three_to_two(particles, sampled_types, three_body_int, V);
      }
    }
*/
/*
    for (const auto &part : particles) {
      std::cout << part.momentum().abs() << " " << part.momentum().x0() << std::endl;
    }
*/
/*
    for (const ParticleTypePtr t : sampled_types) {
      std::cout << type_count(particles, t) << " ";
    }
    std::cout << std::endl;
  }

  QuantumNumbers cons(particles);
  std::cout << "Final momentum: " << cons.momentum() << std::endl;
*/
}
