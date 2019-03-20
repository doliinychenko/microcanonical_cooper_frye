#include <getopt.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <streambuf>
#include <string>

#include "main.h"

#include "smash/angles.h"
#include "smash/constants.h"
#include "smash/pow.h"
#include "smash/random.h"
#include "smash/setup_particles_decaymodes.h"

using namespace smash;

void reproduce_arxiv_1902_09775() {
  auto species_to_sample = [](const smash::ParticleTypePtr t) {
    // Hadrons above this mass are not sampled. It also matters
    // for computing total energy and charges on the hypersurface.
    static constexpr double max_mass = 2.5;  // GeV
    return t->is_hadron() && t->mass() < max_mass;
  };

  constexpr bool quantum_statistics = false;
  const int N_warmup = 1E6, N_decorrelate = 2E2, N_printout = 1E5;
  const std::string hypersurface_input_file =
      "../hydro_cells_reproduce_arxiv1902.09775";
  const HyperSurfacePatch::InputFormat format =
      HyperSurfacePatch::InputFormat::DimaNaiveFormat;
  ParticleTypePtrList printout_types;
  for (const ParticleType &ptype : ParticleType::list_all()) {
    if (ptype.is_hadron() && ptype.is_stable() && ptype.mass() < 2.0 &&
        ptype.baryon_number() != -1 && ptype.pdgcode().charmness() == 0) {
      printout_types.push_back(&ptype);
    }
  }

  HyperSurfacePatch hyper(hypersurface_input_file, format,
                          species_to_sample, quantum_statistics);
  std::cout << "# Hypersurface: " << hyper << std::endl;
  MicrocanonicalSampler sampler(species_to_sample, 0, quantum_statistics);
  MicrocanonicalSampler::SamplerParticleList particles;
  std::cout << "# Initializing the sampler " << std::endl;
  sampler.initialize(hyper, particles);
  std::cout << "# Warming up." << std::endl;
  for (int i = 0; i < N_warmup; ++i) {
    sampler.one_markov_chain_step(hyper, particles);
  }
  std::cout << "# Finished warming up." << std::endl;
  std::cout << "# Total particles: " << particles.size() << std::endl;
  sampler.print_rejection_stats();
  for (const ParticleTypePtr t : printout_types) {
    std::cout << t->name() << "  ";
  }
  std::cout << std::endl;

  for (int j = 0; j < N_printout; j++) {
    for (int i = 0; i < N_decorrelate; ++i) {
      sampler.one_markov_chain_step(hyper, particles);
    }
    // Count every particle type in every cell
    std::map<std::pair<ParticleTypePtr, size_t>, size_t> counter;
    for (const auto &particle : particles) {
      const std::pair<ParticleTypePtr, size_t> key(particle.type,
                                                   particle.cell_index);
      if (counter.find(key) == counter.end()) {
        counter[key] = 0;
      }
      counter[key]++;
    }
    // Printout
    for (size_t i_cell = 0; i_cell < hyper.Ncells(); i_cell++) {
      for (const ParticleTypePtr t : printout_types) {
        const std::pair<ParticleTypePtr, size_t> key(t, i_cell);
        std::cout << " " << counter[key];
      }
    }
    std::cout << std::endl;
  }
  MicrocanonicalSampler::QuantumNumbers cons(particles);
  assert((cons.momentum - hyper.pmu()).abs() < 1.e-5);
  assert(cons.B == hyper.B());
  assert(cons.S == hyper.S());
  assert(cons.Q == hyper.Q());
}

void sample(std::string hypersurface_input_file,
            HyperSurfacePatch::InputFormat hypersurface_file_format,
            std::vector<ParticleTypePtr> printout_types, int N_warmup,
            int N_decorrelate, int N_printout,
            double max_mass, double E_patch) {
  constexpr bool quantum_statistics = false;
  /**
   * A function, which defines, which species will be sampled. For
   * a given species, if it returns true, the species will be sampled.
   */
  auto is_sampled_type = [&](const smash::ParticleTypePtr t) {
    return t->is_hadron() && t->mass() < max_mass;
  };

  HyperSurfacePatch hyper(hypersurface_input_file, hypersurface_file_format,
                          is_sampled_type, quantum_statistics);
  std::cout << "Full hypersurface: " << hyper << std::endl;
  MicrocanonicalSampler sampler(is_sampled_type, 0, quantum_statistics);

  auto patches = hyper.split(E_patch);
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

int main(int argc, char **argv) {
  smash::random::set_seed(smash::random::generate_63bit_seed());
  smash::load_default_particles_and_decaymodes();

  double Epatch = 10.0;  // GeV

  int opt = 0;
  while ((opt = getopt (argc, argv, "rte:")) != -1) {
    switch (opt) {
      case 'r':
        reproduce_arxiv_1902_09775();
        std::exit(EXIT_SUCCESS);
      case 't':
        test_3body_integrals_precision();
        std::exit(EXIT_SUCCESS);
      case 'e':
        Epatch = std::stod(optarg);
        std::cout << "Using patch energy " << Epatch << " GeV" << std::endl;
        break;
      default:
        break;
    }
  }

  ParticleTypePtrList printout_types;
  for (const ParticleType &ptype : ParticleType::list_all()) {
    if (ptype.is_hadron() && ptype.is_stable() && ptype.mass() < 2.0 &&
        ptype.baryon_number() != -1 && ptype.pdgcode().charmness() == 0) {
      printout_types.push_back(&ptype);
    }
  }

  const int N_warmup = 1E6, N_decorrelate = 2E2, N_printout = 1E4;
  constexpr double max_mass = 2.0;  // GeV
  sample("../../hyper_from_Jan/hypersurface/spinodal_hyper_pbpb_elb3.5_39-2.f16",
         HyperSurfacePatch::InputFormat::Steinheimer,
         printout_types, N_warmup, N_decorrelate, N_printout, max_mass, Epatch);
}
