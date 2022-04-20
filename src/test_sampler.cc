#include <getopt.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <streambuf>
#include <string>
#include <cstring>
#include <random>

#include "microcanonical_sampler/main.h"
#include "microcanonical_sampler/microcanonical_sampler.h"
#include "microcanonical_sampler/sampler_particletype_list.h"
#include "microcanonical_sampler/statistics_summary.h"

#include "smash/angles.h"
#include "smash/constants.h"
#include "smash/pow.h"
#include "smash/random.h"
#include "smash/stringfunctions.h"

using namespace sampler;

namespace smash {
  random::Engine random::engine;
}

int64_t generate_63bit_seed() {
  std::random_device rd;
  static_assert(std::is_same<decltype(rd()), uint32_t>::value,
                "random_device is assumed to generate uint32_t");
  uint64_t unsigned_seed =
      (static_cast<uint64_t>(rd()) << 32) | static_cast<uint64_t>(rd());
  // Discard the highest bit to make sure it fits into a positive int64_t
  const int64_t seed = static_cast<int64_t>(unsigned_seed >> 1);
  return seed;
}

void reproduce_arxiv_1902_09775() {
  auto species_to_sample = [](ParticleTypePtr t) {
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
  std::vector<ParticleTypePtr> printout_types;
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
  for (ParticleTypePtr t : printout_types) {
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
      for (ParticleTypePtr t : printout_types) {
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
  auto is_sampled_type = [&](ParticleTypePtr t) {
    return t->is_hadron() && t->mass() < max_mass;
  };

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
        const smash::ThreeVector v = a_cell.u.velocity();
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
  step_until_sufficient_decorrelation(sampler, patches, particles, N_warmup);
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

  output_file <<
      "#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge\n"
      "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none none\n";
  Statistics stats = Statistics();
  for (size_t j = 0; j < N_printout; j++) {
    step_until_sufficient_decorrelation(sampler, patches, particles,
                                        N_decorrelate);
    // print out
    if (j % 10000 == 0 && j != 0) {
      std::cout << "sample " << j << std::endl;
    }
    output_file << "# event " << j << std::endl;
    int particle_counter = 0;
    for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
      for (const auto &particle : particles[i_patch]) {
        const auto &cell = patches[i_patch].cell(particle.cell_index);
        const smash::FourVector p = particle.momentum;
        output_file << cell.r.x0() << " " << cell.r.x1() << " "
                    << cell.r.x2() << " " << cell.r.x3() << " "
                    << particle.type->mass() << " "
                    << p.x0() << " " << p.x1() << " "
                    << p.x2() << " " << p.x3() << " "
                    << particle.type->pdgcode().get_decimal() << " "
                    << particle_counter << " "
                    << particle.type->pdgcode().charge() << std::endl;
        particle_counter++;
      }
    }
    output_file << "# event " << j << " end" << std::endl;
    stats.add_event(particles);
  }
  stats.printout("stat_test_E" + std::to_string(E_patch) + ".txt");

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
      "  -p, --particles         <particles_file>,<format>\n"
      "                          list of particle species to be sampled in\n"
      "                          <particles_file> and the file format,\n"
      "                          either \"SMASH\" or \"iSS\"\n"
      "                          (see README.md for more details)\n"
      "  -s, --surface           <surface_file>,<surface_format>\n"
      "                          File with the list of hypersurface elements\n"
      "                          and its format: Dima_format,\n"
      "                          Steinheimer_format, MUSIC_format,\n"
      "                          , VISH_format or vHLLE_format\n"
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
  smash::random::set_seed(generate_63bit_seed());
  std::string particles_file = "../smash/input/particles.txt";
  ParticleListFormat particles_file_format =
       ParticleListFormat::SMASH;
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
          std::vector<std::string> pd_strings = smash::split(arg_string, ',');
          if (pd_strings.size() != 2) {
            std::cout << "Need particle list file and format separated "
                  << "by comma: particle_list,format" << std::endl;
            usage(EXIT_FAILURE, progname);
          }
          particles_file = pd_strings[0];
          if (pd_strings[1] == "SMASH") {
            particles_file_format = ParticleListFormat::SMASH;
          } else if (pd_strings[1] == "iSS") {
            particles_file_format = ParticleListFormat::iSS;
          } else {
            usage(EXIT_FAILURE, progname);
          }
          break;
        }
      case 's':
        {
          std::string arg_string(optarg);
          std::vector<std::string> args = smash::split(arg_string, ',');
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
          } else if (args[1] == "vHLLE_format") {
            hypersurface_file_format =
              HyperSurfacePatch::InputFormat::vHLLE_ASCII;
          } else {
            usage(EXIT_FAILURE, progname);
          }
          break;
        }
      case 'y':
        {
          std::string arg_string(optarg);
          std::vector<std::string> args = smash::split(arg_string, ',');
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
        std::cout << "Number of events: " << N_printout << std::endl;
        break;
      case 'r':
        read_particle_list("../smash/input/particles.txt",
                           ParticleListFormat::SMASH);
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

  read_particle_list(particles_file, particles_file_format);

  const size_t N_warmup = 1E6, N_decorrelate = 500;
  constexpr double max_mass = 2.5;  // GeV
  sample(hypersurface_input_file, hypersurface_file_format, eta_for_2Dhydro,
         output_file, patches_output_filename,
         N_warmup, N_decorrelate, N_printout, max_mass, Epatch);
}
