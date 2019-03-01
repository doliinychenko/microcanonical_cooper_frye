#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <streambuf>

// #include "gsl/gsl_sf_gamma.h"

#include "kmeans_clustering.h"
#include "main.h"

#include "smash/random.h"
#include "smash/angles.h"
#include "smash/constants.h"
#include "smash/pow.h"

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
  std::string smash_dir("/home/dima/Work/SMASH/smash-devel/");
  // smash_dir = std::getenv("SMASH_DIR");
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

// Decides, which species are going to be sampled
bool is_sampled_type(const smash::ParticleTypePtr t) {
  return t->is_hadron() && t->mass() < 2.5;
}

void sample(std::string hypersurface_input_file,
            HyperSurfacePatch::InputFormat hypersurface_file_format,
            std::vector<ParticleTypePtr>printout_types,
            int N_warmup, int N_decorrelate, int N_printout) {
  constexpr bool quantum_statistics = false;

  HyperSurfacePatch hyper(hypersurface_input_file,
                          hypersurface_file_format,
                          is_sampled_type,
                          quantum_statistics);
  std::cout << hyper << std::endl;

  MicrocanonicalSampler sampler(is_sampled_type, 0, quantum_statistics);
  sampler.initialize(hyper);

  std::cout << "Warming up." << std::endl;
  for (int i = 0; i < N_warmup; ++i) {
    sampler.one_markov_chain_step(hyper);
  }
  std::cout << "Finished warming up." << std::endl;
  std::cout << sampler.particles().size() << " particles" << std::endl;
  sampler.print_rejection_stats();

  ParticleTypePtrList sampled_types;
  for (const ParticleType &ptype : ParticleType::list_all()) {
    if (is_sampled_type(&ptype)) {
      sampled_types.push_back(&ptype);
    }
  }
  for (const ParticleTypePtr t : sampled_types) {
    if (t->is_stable()) {
      std::cout << t->name() << "  ";
    }
  }
  std::cout << std::endl;
  std::cout << "Total sampled species: " << sampled_types.size() << std::endl;

  for (int j = 0; j < N_printout; j++) {
    for (int i = 0; i < N_decorrelate; ++i) {
      sampler.one_markov_chain_step(hyper);
    }
    std::map<std::pair<ParticleTypePtr, size_t>, size_t> counter;
    for (const auto &particle : sampler.particles()) {
      const std::pair<ParticleTypePtr, size_t> key(
                particle.type, particle.cell_index);
      if (counter.find(key) == counter.end()) {
        counter[key] = 0;
      }
      counter[key]++;
    }
    for (size_t icell = 0; icell < hyper.Ncells(); icell++) {
      for (const ParticleTypePtr t : printout_types) {
        const std::pair<ParticleTypePtr, size_t> key(t, icell);
        std::cout << " " << counter[key];
      }
    }
    std::cout << std::endl;
  }

  MicrocanonicalSampler::QuantumNumbers cons(sampler.particles());
  assert((cons.momentum - hyper.pmu()).abs() < 1.e-6);
  assert(cons.B == hyper.B());
  assert(cons.S == hyper.S());
  assert(cons.Q == hyper.Q());
  // sampler.print_rejection_stats();
}

int main() {
  initialize_random_number_generator();
  load_smash_particles();

  // test_clustering();
  // test_3body_integrals_precision();

  ParticleTypePtrList printout_types;
  for (const ParticleType &ptype : ParticleType::list_all()) {
    if (is_sampled_type(&ptype) && ptype.is_stable() && ptype.mass() < 2.0 &&
        ptype.baryon_number() != -1 && ptype.pdgcode().charmness() == 0) {
      printout_types.push_back(&ptype);
    }
  }
  for (const ParticleTypePtr t : printout_types) {
    std::cout << t->name() << "  ";
  }
  std::cout << std::endl;

  const int N_warmup = 1E6,
            N_decorrelate = 2E2,
            N_printout = 1E5;
  sample("../hydro_cells.dat",
         HyperSurfacePatch::InputFormat::DimaNaiveFormat,
         printout_types, N_warmup, N_decorrelate, N_printout);
}
