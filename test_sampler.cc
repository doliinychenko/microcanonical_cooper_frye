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
  return t->is_hadron() && t->mass() < 1.5;
}

void sample(std::string hypersurface_input_file,
            HyperSurfacePatch::InputFormat hypersurface_file_format,
            std::vector<ParticleTypePtr> printout_types, int N_warmup,
            int N_decorrelate, int N_printout) {
  constexpr bool quantum_statistics = false;

  HyperSurfacePatch hyper(hypersurface_input_file, hypersurface_file_format,
                          is_sampled_type, quantum_statistics);
  std::cout << "Full hypersurface: " << hyper << std::endl;
  MicrocanonicalSampler sampler(is_sampled_type, 0, quantum_statistics);

  // constexpr size_t number_of_patches = 3;
  constexpr double E_patch = 5.0;  // GeV
  auto patches = hyper.split2(E_patch);
  size_t number_of_patches = patches.size();

  std::vector<MicrocanonicalSampler::SamplerParticleList> particles;
  particles.resize(number_of_patches);

  for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
    std::cout << "Initializing patch " << i_patch << std::endl;
    sampler.initialize(patches[i_patch], particles[i_patch]);
  }

  std::cout << "Warming up." << std::endl;
  #pragma omp parallel for
  for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
    std::cout << "Patch " << i_patch << std::endl;
    for (int i = 0; i < N_warmup; ++i) {
      sampler.one_markov_chain_step(patches[i_patch], particles[i_patch]);
    }
  }
  std::cout << "Finished warming up." << std::endl;
  size_t total_particles = 0;
  for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
    total_particles += particles[i_patch].size();
  }
  std::cout << total_particles << " particles" << std::endl;
  sampler.print_rejection_stats();

  for (const ParticleTypePtr t : printout_types) {
    std::cout << t->name() << "  ";
  }
  std::cout << std::endl;

  for (int j = 0; j < N_printout; j++) {
    #pragma omp parallel for
    for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
      for (int i = 0; i < N_decorrelate; ++i) {
        sampler.one_markov_chain_step(patches[i_patch], particles[i_patch]);
      }
    }
    std::map<std::pair<ParticleTypePtr, size_t>, size_t> counter;
    for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
      for (const auto &particle : particles[i_patch]) {
        const std::pair<ParticleTypePtr, size_t> key(particle.type, i_patch);
        if (counter.find(key) == counter.end()) {
          counter[key] = 0;
        }
        counter[key]++;
      }
    }
    for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
      for (const ParticleTypePtr t : printout_types) {
        const std::pair<ParticleTypePtr, size_t> key(t, i_patch);
        std::cout << " " << counter[key];
      }
    }
    std::cout << std::endl;
  }

  for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
    MicrocanonicalSampler::QuantumNumbers cons(particles[i_patch]);
    assert((cons.momentum - patches[i_patch].pmu()).abs() < 1.e-6);
    assert(cons.B == patches[i_patch].B());
    assert(cons.S == patches[i_patch].S());
    assert(cons.Q == patches[i_patch].Q());
  }
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
  //sample("../../hyper_from_Jan/hypersurface/spinodal_hyper_pbpb_elb3.5_39-2.f16",
  //       HyperSurfacePatch::InputFormat::Steinheimer,
  //       printout_types, N_warmup, N_decorrelate, N_printout);
}
