#include <getopt.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <streambuf>
#include <string>
#include <cstring>

#include "main.h"
#include "microcanonical_sampler.h"

#include "smash/angles.h"
#include "smash/constants.h"
#include "smash/pow.h"
#include "smash/random.h"
#include "smash/setup_particles_decaymodes.h"
#include "smash/stringfunctions.h"

using namespace smash;

void reproduce_arxiv_1902_09775() {
  auto species_to_sample = [](const smash::ParticleTypePtr t) {
    // Hadrons above this mass are not sampled. It also matters
    // for computing total energy and charges on the hypersurface.
    static constexpr double max_mass = 2.5;  // GeV
    return t->is_hadron() && t->mass() < max_mass;
  };

  constexpr bool quantum_statistics = false;
  const size_t N_warmup = 1E6, N_decorrelate = 2E2, N_printout = 1E5;
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

  // The eta range {0.0, 0.0, 0.0} is completely dummy here
  HyperSurfacePatch hyper(hypersurface_input_file, format, {0.0, 0.0, 0.0},
                          species_to_sample, quantum_statistics);
  std::cout << "# Hypersurface: " << hyper << std::endl;
  MicrocanonicalSampler sampler(species_to_sample, 0, quantum_statistics);
  MicrocanonicalSampler::SamplerParticleList particles;
  std::cout << "# Initializing the sampler " << std::endl;
  sampler.initialize(hyper, particles);
  std::cout << "# Warming up." << std::endl;
  for (size_t i = 0; i < N_warmup; ++i) {
    sampler.one_markov_chain_step(hyper, particles);
  }
  std::cout << "# Finished warming up." << std::endl;
  std::cout << "# Total particles: " << particles.size() << std::endl;
  sampler.print_rejection_stats();
  for (const ParticleTypePtr t : printout_types) {
    std::cout << t->name() << "  ";
  }
  std::cout << std::endl;

  for (size_t j = 0; j < N_printout; j++) {
    for (size_t i = 0; i < N_decorrelate; ++i) {
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
  assert((cons.momentum - hyper.pmu()).abs() < 1.e-4);
  assert(cons.B == hyper.B());
  assert(cons.S == hyper.S());
  assert(cons.Q == hyper.Q());
}

void step_until_sufficient_decorrelation(
    MicrocanonicalSampler& sampler,
    const std::vector<HyperSurfacePatch>& patches,
    std::vector<MicrocanonicalSampler::SamplerParticleList>& particles,
    size_t min_steps_number, double required_decorrelation_degree) {
  size_t number_of_patches = patches.size();
  // Fixed amount of 2<->3. Forcing decorrelation in 2<->3 leads to a bias
  // in particle number distribution.
  #pragma omp parallel for
  for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
    for (size_t i = 0; i < min_steps_number; ++i) {
      sampler.one_markov_chain_step(patches[i_patch], particles[i_patch]);
    }
  }
  // Do 2<->2 until decorrelation. This may bias momentum correlations,
  // but does not bias multiplicity distribution.
  #pragma omp parallel for
  for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
    for (auto &particle : particles[i_patch]) {
      particle.decorrelated = false;
    }
    size_t non_decorr_counter;
    do {
      for (size_t i = 0; i < min_steps_number; ++i) {
        sampler.random_two_to_two(patches[i_patch], particles[i_patch]);
      }
      non_decorr_counter = std::count_if(
          particles[i_patch].begin(), particles[i_patch].end(),
          [](MicrocanonicalSampler::SamplerParticle p) { return !p.decorrelated; });
      // printf("Patch %lu: not decorrelated %lu/%lu\n", i_patch,
      //       non_decorr_counter, particles[i_patch].size());
    } while (non_decorr_counter >
             particles[i_patch].size() * required_decorrelation_degree);
  }
}

void sample(const std::string hypersurface_input_file,
            HyperSurfacePatch::InputFormat hypersurface_file_format,
            const std::array<double, 3> &eta_for_2Dhydro,
            const std::string output_file_name,
            const std::string patches_output_filename, size_t N_warmup,
            size_t N_decorrelate, size_t N_printout,
            double max_mass, double E_patch) {
  constexpr bool quantum_statistics = false;
  /**
   * A function, which defines, which species will be sampled. For
   * a given species, if it returns true, the species will be sampled.
   */
  auto is_sampled_type = [&](const smash::ParticleTypePtr t) {
    return t->is_hadron() && t->mass() < max_mass;
  };
  /**
   * Percentage of particles that are allowed to be untouched after a
   * decorrelation session: 0 is maximally strict, 1 is not strict
   * at all.
   */
  constexpr double sufficient_decorrelation = 0.01;

  HyperSurfacePatch hyper(hypersurface_input_file, hypersurface_file_format,
                          eta_for_2Dhydro,
                          is_sampled_type, quantum_statistics);
  std::cout << "Full hypersurface: " << hyper << std::endl;
  MicrocanonicalSampler sampler(is_sampled_type, 0, quantum_statistics);

  auto patches = hyper.split(E_patch);
  size_t number_of_patches = patches.size();

  // Save labelled patches
  if (patches_output_filename != "") {
     std::cout << "Dumping labeled hypersurface cells to "
               << patches_output_filename << std::endl;
    std::ofstream patches_output_file(patches_output_filename, std::ios::out);
    for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
      for (const auto a_cell : patches[i_patch].cells()) {
        const ThreeVector v = a_cell.u.velocity();
        patches_output_file << a_cell.r.x1() << " "
                            << a_cell.r.x2() << " "
                            << a_cell.r.x3() << " "
                            << a_cell.T << " "
                            << a_cell.muB << " "
                            << a_cell.muS << " "
                            << v.x1() << " " << v.x2() << " " << v.x3() << " "
                            << a_cell.dsigma.x0() << " "
                            << a_cell.dsigma.x1() << " "
                            << a_cell.dsigma.x2() << " "
                            << a_cell.dsigma.x3() << " "
                            << i_patch << std::endl;
      }
    }
  }

  std::vector<MicrocanonicalSampler::SamplerParticleList> particles;
  particles.resize(number_of_patches);

  for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
    std::cout << "Initializing patch " << i_patch << std::endl;
    sampler.initialize(patches[i_patch], particles[i_patch]);
  }

  std::cout << "Warming up." << std::endl;
  step_until_sufficient_decorrelation(sampler, patches, particles,
                                      N_warmup, sufficient_decorrelation);
  std::cout << "Finished warming up." << std::endl;
  size_t total_particles = 0;
  for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
    total_particles += particles[i_patch].size();
  }
  std::cout << total_particles << " particles" << std::endl;
  sampler.print_rejection_stats();

  std::ofstream output_file(output_file_name, std::ios::out);
  if (!output_file.is_open()) {
    throw std::runtime_error("Can't open the output file.");
  }

  for (size_t j = 0; j < N_printout; j++) {
    step_until_sufficient_decorrelation(sampler, patches, particles,
                                        N_decorrelate, sufficient_decorrelation);
    // print out
    if (j % 10000 == 0 && j != 0) {
      std::cout << "sample " << j << std::endl;
    }
    output_file << "# sample " << j << std::endl;
    for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
      for (const auto &particle : particles[i_patch]) {
        const auto &cell = patches[i_patch].cell(particle.cell_index);
        const FourVector p = particle.momentum;
        output_file << cell.r.x0() << " " << cell.r.x1() << " "
                    << cell.r.x2() << " " << cell.r.x3() << " "
                    << p.x0() << " " << p.x1() << " "
                    << p.x2() << " " << p.x3() << " "
                    << particle.type->pdgcode().get_decimal() << std::endl;
      }
    }
  }

  for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
    MicrocanonicalSampler::QuantumNumbers cons(particles[i_patch]);
    assert((cons.momentum - patches[i_patch].pmu()).abs() < 1.e-4);
    assert(cons.B == patches[i_patch].B());
    assert(cons.S == patches[i_patch].S());
    assert(cons.Q == patches[i_patch].Q());
  }
  sampler.print_rejection_stats();
}

void usage(const int rc, const std::string &progname) {
  std::printf("\nUsage: %s [option]\n\n", progname.c_str());
  std::printf(
      "  -h, --help              usage information\n\n"
      "  -p, --particles         <particles_file>,<decaymodes_file>\n"
      "                          list of particle species to be sampled in\n"
      "                          <particles_file> and their decays in the\n"
      "                          <decaymodes_file>, both in SMASH format\n"
      "        (see http://theory.gsi.de/~smash/doc/1.6/inputparticles.html)\n"
      "  -s, --surface           <surface_file>,<surface_format>\n"
      "                          File with the list of hypersurface elements\n"
      "                          and its format: Dima_format,\n"
      "                          Steinheimer_format, MUSIC_format,\n"
      "                          or VISH_format\n"
      "  -y, --eta_range         <eta_min>,<eta_max>,<d_eta>\n"
      "                          range of pseudorapidity for the case of\n"
      "                          2+1D hydro. Default is\n"
      "                          eta_min = -2, eta_max = 2, d_eta = 0.4\n"
      "  -n, --nevents           number of sampled instances to output\n"
      "  -r, --reproduce_1902_09775 reproduce results of arxiv::1902.09775\n"
      "  -t, --test                 run testing functions\n"
      "  -e, --energy_patch         set maximal energy [GeV] in the patch\n"
      "  -o, --outputfile           output file name, where sampled particles\n"
      "                             coordinates, momenta, and pdg ids\n"
      "                             will be printed out\n"
      "                             (default: ./sampled_particles.dat)\n"
      "  -l, --labeled_hyper_output name of the file to print out the same\n"
      "                             hypersurface cells that were read in, but\n"
      "                             with patch labels.\n\n");
  std::exit(rc);
}

int main(int argc, char **argv) {
  smash::random::set_seed(smash::random::generate_63bit_seed());
  char *particles_file = nullptr,
       *decaymodes_file = nullptr;
  double Epatch = 10.0;  // GeV
  size_t N_printout = 1E4;
  std::string output_file = "sampled_particles.dat",
              patches_output_filename = "";
  std::string hypersurface_input_file(
      "../../surface_from_Jan/spinodal_hyper_pbpb_elb3.5_28-5.f16");
  HyperSurfacePatch::InputFormat hypersurface_file_format =
      HyperSurfacePatch::InputFormat::Steinheimer;
  std::array<double, 3> eta_for_2Dhydro = {-2.0, 2.0, 0.4};


  constexpr option longopts[] = {
      {"help", no_argument, 0, 'h'},
      {"particles", required_argument, 0, 'p'},
      {"surface", required_argument, 0, 's'},
      {"eta_range", required_argument, 0, 'y'},
      {"nevents", required_argument, 0, 'n'},
      {"reproduce_1902_09775", no_argument, 0, 'r'},
      {"test", no_argument, 0, 't'},
      {"energy_patch", required_argument, 0, 'e'},
      {"outputfile", required_argument, 0, 'o'},
      {"labeled_hyper_output", required_argument, 0, 'l'},
      {nullptr, 0, 0, 0}};

  const std::string full_progname = std::string(argv[0]);
  const int i1 = full_progname.find_last_of("\\/") + 1,
            i2 = full_progname.size();
  const std::string progname = full_progname.substr(i1, i2);
  int opt = 0;
  while ((opt = getopt_long(argc, argv, "hp:s:y:n:rte:o:l:",
          longopts, nullptr)) != -1) {
    switch (opt) {
      case 'h':
        usage(EXIT_SUCCESS, progname);
        break;
      case 'p':
        {
          std::string arg_string(optarg);
          std::vector<std::string> pd_strings = split(arg_string, ',');
          if (pd_strings.size() != 2) {
            throw std::invalid_argument(
                "-p usage: particles_file,decaymodes_file");
          }
          particles_file = strdup(pd_strings[0].c_str());
          decaymodes_file = strdup(pd_strings[1].c_str());
          break;
        }
      case 's':
        {
          std::string arg_string(optarg);
          std::vector<std::string> args = split(arg_string, ',');
          if (args.size() != 2) {
            usage(EXIT_FAILURE, progname);
          }
          hypersurface_input_file = args[0];
          if (args[1] == "Dima_format") {
            hypersurface_file_format =
              HyperSurfacePatch::InputFormat::DimaNaiveFormat;
          } else if (args[1] == "Steinheimer_format") {
            hypersurface_file_format =
              HyperSurfacePatch::InputFormat::Steinheimer;
          } else if (args[1] == "MUSIC_format") {
            hypersurface_file_format =
              HyperSurfacePatch::InputFormat::MUSIC_ASCII_3plus1D;
          } else if (args[1] == "VISH_format") {
            hypersurface_file_format =
              HyperSurfacePatch::InputFormat::VISH_2files;
          } else {
            usage(EXIT_FAILURE, progname);
          }
          break;
        }
      case 'y':
        {
          std::string arg_string(optarg);
          std::vector<std::string> args = split(arg_string, ',');
          if (hypersurface_file_format !=
              HyperSurfacePatch::InputFormat::VISH_2files) {
            std::cout << "-y option only makes sense with 2+1D hydro."
                      << std::endl;
            usage(EXIT_FAILURE, progname);
          }
         if (args.size() != 3) {
            usage(EXIT_FAILURE, progname);
          }
          eta_for_2Dhydro[0] = std::stod(args[0]);
          eta_for_2Dhydro[1] = std::stod(args[1]);
          eta_for_2Dhydro[2] = std::stod(args[2]);
          break;
        }
      case 'n':
        N_printout = std::stoi(optarg);
        break;
      case 'r':
        smash::load_default_particles_and_decaymodes();
        reproduce_arxiv_1902_09775();
        std::exit(EXIT_SUCCESS);
      case 't':
        MicrocanonicalSampler::test_3body_integrals();
        MicrocanonicalSampler::test_3body_phase_space_sampling();
        std::exit(EXIT_SUCCESS);
      case 'e':
        Epatch = std::stod(optarg);
        std::cout << "Using patch energy " << Epatch << " GeV" << std::endl;
        break;
      case 'o':
        output_file = optarg;
        break;
      case 'l':
        patches_output_filename = optarg;
        break;
     default:
        usage(EXIT_FAILURE, progname);
    }
  }

  // Abort if there are unhandled arguments left.
  if (optind < argc) {
    std::cout << argv[0] << ": invalid argument -- '" << argv[optind] << "'\n";
    usage(EXIT_FAILURE, progname);
  }

  auto pd = smash::load_particles_and_decaymodes(particles_file, decaymodes_file);
  ParticleType::create_type_list(pd.first);
  DecayModes::load_decaymodes(pd.second);
  ParticleType::check_consistency();

  const size_t N_warmup = 1E6, N_decorrelate = 500;
  constexpr double max_mass = 2.5;  // GeV
  sample(hypersurface_input_file, hypersurface_file_format, eta_for_2Dhydro,
         output_file, patches_output_filename,
         N_warmup, N_decorrelate, N_printout, max_mass, Epatch);
}
