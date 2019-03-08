#include <fstream>
#include <iomanip>
#include <iostream>
#include <streambuf>
#include <string>

#include "kmeans_clustering.h"
#include "main.h"

#include "smash/angles.h"
#include "smash/constants.h"
#include "smash/pow.h"
#include "smash/random.h"
#include "smash/setup_particles_decaymodes.h"

using namespace smash;

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
            std::vector<ParticleTypePtr> printout_types, int N_warmup,
            int N_decorrelate, int N_printout) {
  constexpr bool quantum_statistics = false;

  HyperSurfacePatch hyper(hypersurface_input_file, hypersurface_file_format,
                          is_sampled_type, quantum_statistics);
  // auto patches = hyper.split(3);
  std::cout << hyper << std::endl;

  MicrocanonicalSampler sampler(is_sampled_type, 0, quantum_statistics);
  MicrocanonicalSampler::SamplerParticleList particles;
  sampler.initialize(hyper, particles);

  std::cout << "Warming up." << std::endl;
  for (int i = 0; i < N_warmup; ++i) {
    sampler.one_markov_chain_step(hyper, particles);
  }
  std::cout << "Finished warming up." << std::endl;
  std::cout << particles.size() << " particles" << std::endl;
  sampler.print_rejection_stats();

  for (const ParticleTypePtr t : printout_types) {
    std::cout << t->name() << "  ";
  }
  std::cout << std::endl;

  for (int j = 0; j < N_printout; j++) {
    for (int i = 0; i < N_decorrelate; ++i) {
      sampler.one_markov_chain_step(hyper, particles);
    }
    std::map<std::pair<ParticleTypePtr, size_t>, size_t> counter;
    for (const auto &particle : particles) {
      const std::pair<ParticleTypePtr, size_t> key(particle.type,
                                                   particle.cell_index);
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

  MicrocanonicalSampler::QuantumNumbers cons(particles);
  assert((cons.momentum - hyper.pmu()).abs() < 1.e-6);
  assert(cons.B == hyper.B());
  assert(cons.S == hyper.S());
  assert(cons.Q == hyper.Q());
  // sampler.print_rejection_stats();
}

int main() {
  smash::random::set_seed(smash::random::generate_63bit_seed());
  smash::load_default_particles_and_decaymodes();

  // test_clustering();
  // test_3body_integrals_precision();

  ParticleTypePtrList printout_types;
  for (const ParticleType &ptype : ParticleType::list_all()) {
    if (is_sampled_type(&ptype) && ptype.is_stable() && ptype.mass() < 2.0 &&
        ptype.baryon_number() != -1 && ptype.pdgcode().charmness() == 0) {
      printout_types.push_back(&ptype);
    }
  }

  const int N_warmup = 1E6, N_decorrelate = 2E2, N_printout = 1E5;
  sample("../hydro_cells.dat", HyperSurfacePatch::InputFormat::DimaNaiveFormat,
         printout_types, N_warmup, N_decorrelate, N_printout);
}
