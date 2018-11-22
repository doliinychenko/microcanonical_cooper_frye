#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <streambuf>

#include <vector>
#include <set>

#include "gsl/gsl_sf_gamma.h"

#include "3body_int.h"

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

void renormalize_momenta(ParticleList &particles,
                         const FourVector required_total_momentum) {

  // Centralize momenta
  QuantumNumbers conserved = QuantumNumbers(particles);
  std::cout << particles.size() << " particles. " << std::endl;
  std::cout << "Required 4-momentum: " << required_total_momentum << std::endl;
  std::cout << "Sampled 4-momentum: "  << conserved.momentum() << std::endl;
  const ThreeVector mom_to_add =
      (required_total_momentum.threevec() - conserved.momentum().threevec()) /
      particles.size();
  std::cout << "Adjusting momenta by " << mom_to_add << std::endl;
  for (auto &particle : particles) {
    particle.set_4momentum(particle.type().mass(),
                           particle.momentum().threevec() + mom_to_add);
  }

  // Boost every particle to the common center of mass frame
  conserved = QuantumNumbers(particles);
  const ThreeVector beta_CM_generated = conserved.momentum().velocity();
  const ThreeVector beta_CM_required = required_total_momentum.velocity();

  double E = 0.0;
  double E_expected = required_total_momentum.abs();
  for (auto &particle : particles) {
    particle.boost_momentum(beta_CM_generated);
    E += particle.momentum().x0();
  }
  // Renorm. momenta by factor (1+a) to get the right energy, binary search
  const double tolerance = really_small;
  double a, a_min, a_max, er;
  const int max_iter = 250;
  int iter = 0;
  if (E_expected >= E) {
    a_min = 0.0;
    a_max = 10000.0 * E_expected / E;
  } else {
    a_min = -1.0;
    a_max = 0.0;
  }
  do {
    a = 0.5 * (a_min + a_max);
    E = 0.0;
    for (const auto &particle : particles) {
      const double p2 = particle.momentum().threevec().sqr();
      const double E2 = particle.momentum().x0() * particle.momentum().x0();
      E += std::sqrt(E2 + a * (a + 2.0) * p2);
    }
    er = E - E_expected;
    if (er >= 0.0) {
      a_max = a;
    } else {
      a_min = a;
    }
    // std::cout << "Iteration " << iter << ": a = "
    //          << a << ", Δ = " << er << std::endl;
    iter++;
  } while (std::abs(er) > tolerance && iter < max_iter);

  std::cout << "Renormalizing momenta by factor 1+a, a = " << a << std::endl;
  for (auto &particle : particles) {
    particle.set_4momentum(particle.type().mass(),
                           (1 + a) * particle.momentum().threevec());
    particle.boost_momentum(-beta_CM_required);
  }
  QuantumNumbers conserved_final = QuantumNumbers(particles);
  std::cout << "Obtained total momentum: "
            << conserved_final.momentum() << std::endl;
}

void sample_3body_phase_space(double srts,
                              ParticleData &a,
                              ParticleData &b,
                              ParticleData &c) {
  const double m_a = a.type().mass(),
               m_b = b.type().mass(),
               m_c = c.type().mass();
  // sample mab from pCM(sqrt, mab, mc) pCM (mab, ma, mb) <= sqrts^2/4
  double mab, r, probability, pcm_ab, pcm;
  do {
    mab = random::uniform(m_a + m_b, srts - m_c);
    r = random::canonical();
    pcm = pCM(srts, mab, m_c);
    pcm_ab = pCM(mab, m_a, m_b);
    probability = pcm * pcm_ab * 4 / (srts * srts);
  } while (r > probability);
  Angles phitheta;
  phitheta.distribute_isotropically();
  c.set_4momentum(m_c, pcm * phitheta.threevec());
  const ThreeVector beta_cm =
      pcm * phitheta.threevec() / std::sqrt(pcm * pcm + mab * mab);

  phitheta.distribute_isotropically();
  a.set_4momentum(m_a,  pcm_ab * phitheta.threevec());
  b.set_4momentum(m_b, -pcm_ab * phitheta.threevec());
  a.boost_momentum(beta_cm);
  b.boost_momentum(beta_cm);
  // std::cout << a.momentum() + b.momentum() + c.momentum() << std::endl;
}

bool quantum_numbers_match(const ParticleTypePtrList& a,
                           const QuantumNumbers& qn) {
  int B = 0, S = 0, Q = 0;
  for (const auto ptype : a) {
    B += ptype->baryon_number();
    S += ptype->strangeness();
    Q += ptype->charge();
  }
  return (B == qn.baryon_number()) &&
         (S == qn.strangeness()) &&
         (Q == qn.charge());
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

double compute_R2(double srts, double m1, double m2) {
  return (srts > m1 + m2) ? pCM(srts, m1, m2) / (4.0 * M_PI * srts) : 0.0;
}

double compute_sum_R3(const ParticleTypePtrList& sampled_types,
              ThreeBodyIntegrals& three_body_int,
              const QuantumNumbers &cons) {
  // Compute sum of 3-body integrals for given quantum numbers
  const int Ntypes = sampled_types.size();
  FourVector momentum_total = cons.momentum();
  const double srts = momentum_total.abs();
  double sum_R3 = 0.0;
  // std::cout << "sqrts = " << srts << std::endl;
  for (int i1 = 0; i1 < Ntypes; ++i1) {
    for (int i2 = 0; i2 <= i1; ++i2) {
      for (int i3 = 0; i3 <= i2; ++i3) {
        if (quantum_numbers_match({sampled_types[i1],
                                   sampled_types[i2],
                                   sampled_types[i3]}, cons)) {
          const double tmp = three_body_int.value(srts,
                                sampled_types[i1]->mass(),
                                sampled_types[i2]->mass(),
                                sampled_types[i3]->mass());
/*          std::cout << tmp * 1E12 << std::endl;
          if (tmp < 0) {
            std::cout << srts << " " << sampled_types[i1]->mass() << " " <<
                         sampled_types[i2]->mass() << " " <<
                         sampled_types[i3]->mass() << ", R3_ch = " << tmp << std::endl;
          } 
*/
          assert(tmp >= 0);
          sum_R3 += tmp;
/*          
          std::cout << sampled_types[i1]->name()
                    << sampled_types[i2]->name()
                    << sampled_types[i3]->name() << " " << tmp * 1E6
                    << " at srts = " << srts << std::endl;
*/
        }
      }
    }
  }
  // std::cout << std::endl;
  return sum_R3;
}

double compute_sum_R2(const ParticleTypePtrList& sampled_types,
                      const QuantumNumbers &cons) {
  // Compute sum of 2-body integrals for given quantum numbers
  const int Ntypes = sampled_types.size();
  FourVector momentum_total = cons.momentum();
  const double srts = momentum_total.abs();
  double sum_R2 = 0.0;
  for (int i1 = 0; i1 < Ntypes; ++i1) {
    for (int i2 = 0; i2 <= i1; ++i2) {
      if (quantum_numbers_match({sampled_types[i1],
                                 sampled_types[i2]}, cons)) {
        const double m_i1 = sampled_types[i1]->mass(),
                     m_i2 = sampled_types[i2]->mass();
        const double R2 = compute_R2(srts, m_i1, m_i2);
        sum_R2 += R2;
      }
    }
  }
  return sum_R2;
}

void random_two_to_three(ParticleList &particles,
                         const ParticleTypePtrList& sampled_types,
                         ThreeBodyIntegrals& three_body_int,
                         double V) {
  const int Ntypes = sampled_types.size();
  const int Npart = particles.size();
  // Choose 2 random particles
  int in_index1, in_index2, i1, i2, i3;
  in_index1 = random::uniform_int(0, Npart - 1);
  do {
    in_index2 = random::uniform_int(0, Npart - 1);
  } while (in_index1 == in_index2);
  ParticleData in1 = particles[in_index1], in2 = particles[in_index2];
  QuantumNumbers cons({in1, in2});
  FourVector momentum_total = cons.momentum();
  const double srts = momentum_total.abs();
  const ThreeVector beta_cm = momentum_total.velocity();

  const double sum_R3 = compute_sum_R3(sampled_types, three_body_int, cons);
  const double sum_R2 = compute_sum_R2(sampled_types, cons);

  if (sum_R3 < 1.e-100) {
    // std::cout << "Not enough energy (" << srts << " GeV) for 2->3" << std::endl;
    return;
  }

  // Choose channel
  const double r = random::uniform(0.0, sum_R3);
  double partial_sum_R3 = 0.0;
  bool channel_found23 = false;
  for (i1 = 0; i1 < Ntypes && !channel_found23; i1++) {
    for (i2 = 0; i2 <= i1 && !channel_found23; i2++) {
      for (i3 = 0; i3 <= i2 && !channel_found23; i3++) {
        if (quantum_numbers_match({sampled_types[i1],
                                  sampled_types[i2],
                                  sampled_types[i3]}, cons)) {
          const double tmp = three_body_int.value(srts,
                                sampled_types[i1]->mass(),
                                sampled_types[i2]->mass(),
                                sampled_types[i3]->mass());
          assert(tmp >= 0.0);
          partial_sum_R3 += tmp;
          if (partial_sum_R3 >= r) {
            channel_found23 = true;
          }
          if (partial_sum_R3 > sum_R3) {
            throw std::runtime_error("2->3 partial sum exceeded total!");
          }
        }
      }
    }
  }
  i1--;
  i2--;
  i3--;
  if (!channel_found23) {
    return;
  }
  ParticleData out1{*(sampled_types[i1])},
               out2{*(sampled_types[i2])},
               out3{*(sampled_types[i3])};

  sample_3body_phase_space(srts, out1, out2, out3);
  out1.boost_momentum(-beta_cm);
  out2.boost_momentum(-beta_cm);
  out3.boost_momentum(-beta_cm);
  // std::cout << "In: " << in1.momentum() << " " << in2.momentum() << std::endl;
  // std::cout << "Out: " << out1.momentum() << " "
  //                      << out2.momentum() << " "
  //                      << out3.momentum() << std::endl;

  const double E1 = in1.momentum().x0(),
               E2 = in2.momentum().x0();
  const double E1_pr = out1.momentum().x0(),
               E2_pr = out2.momentum().x0(),
               E3_pr = out3.momentum().x0();

  const double phase_space_factor = sum_R3 / sum_R2;
  const double energy_factor = 2.0 * E1_pr * E2_pr * E3_pr / (E1 * E2);
  ParticleTypePtrList in_types {&(in1.type()), &(in2.type())},
                      out_types {&(out1.type()), &(out2.type()), &(out3.type())};
  double spin_factor = 1.0;
  for (const auto t : out_types) {
    spin_factor *= static_cast<double>(1 + t->spin());
  }
  for (const auto t : in_types) {
    spin_factor /= static_cast<double>(1 + t->spin());
  }

  double w = V * pow_int(1.0 / smash::hbarc, 3) * spin_factor *
             energy_factor * phase_space_factor;
  double identity_factor = 1.0;
  if (out1.type() == out2.type() && out2.type() == out3.type()) {
    identity_factor /= 6;
  } else if (out1.type() == out2.type() || out2.type() == out3.type()) {
    identity_factor /= 2;
  }
  if (in1.type() == in2.type()) {
    identity_factor *= 2;
  }
  w *= 3.0 / (particles.size() + 1) * identity_factor;
/*
  std::cout << "Trying 2->3 (" << in1.type().name() << in2.type().name() << "->"
            << out1.type().name() << out2.type().name() << out3.type().name()
            << ") with w = " << w << ", sym_factor: " << sym_factor
            << ", spin_factor = " << spin_factor << ", ph.sp. factor * 1e6= "
            << phase_space_factor * 1.e6 << ", srts = " << srts << std::endl;
*/
  w = std::min(1.0, w);
  if (random::canonical() < w) {
    particles[in_index1] = out1;
    particles[in_index2] = out2;
    particles.push_back(out3);

/*
    std::cout << in1.type().name() << "," << in2.type().name() << "->"
              << out1.type().name() << "," << out2.type().name()
              << "," << out3.type().name() << " ; " << srts << std::endl;
*/
  }
/*  std::cout << type_count(particles, &ParticleType::find(pdg::pi_z)) << " pi0 /"
            << type_count(particles, &ParticleType::find(pdg::pi_p)) << " pi+ /"
            << type_count(particles, &ParticleType::find(pdg::pi_m)) << " pi-" << std::endl;
*/

}

void random_three_to_two(ParticleList &particles,
                         const ParticleTypePtrList& sampled_types,
                         ThreeBodyIntegrals& three_body_int,
                         double V) {
  const int Ntypes = sampled_types.size();
  const int Npart = particles.size();
  if (Npart < 3) {
    return;
  }
  // Choose 3 random particles
  int in_index1, in_index2, in_index3, i1, i2;
  in_index1 = random::uniform_int(0, Npart - 1);
  do {
    in_index2 = random::uniform_int(0, Npart - 1);
  } while (in_index1 == in_index2);
  do {
    in_index3 = random::uniform_int(0, Npart - 1);
  } while (in_index3 == in_index1 || in_index3 == in_index2);
  ParticleData in1 = particles[in_index1],
               in2 = particles[in_index2],
               in3 = particles[in_index3];
  QuantumNumbers cons({in1, in2, in3});
  FourVector momentum_total = cons.momentum();
  const double srts = momentum_total.abs();
  const ThreeVector beta_cm = momentum_total.velocity();

  const double sum_R3 = compute_sum_R3(sampled_types, three_body_int, cons);
  const double sum_R2 = compute_sum_R2(sampled_types, cons);
  if (sum_R2 < 1.e-100) {
    // std::cout << "Not enough energy (" << srts << " GeV) or wrong quantum"
    //             " numbers for 3->2" << std::endl;
    return;
  }


  // Choose channel
  const double r = random::uniform(0.0, sum_R2);
  double partial_sum_R2 = 0.0;
  bool channel_found32 = false;
  for (i1 = 0; (i1 < Ntypes) && (!channel_found32); i1++) {
    for (i2 = 0; (i2 <= i1) && (!channel_found32); i2++) {
     if (quantum_numbers_match({sampled_types[i1],
                                 sampled_types[i2]}, cons)) {
        const double tmp = compute_R2(srts,
                            sampled_types[i1]->mass(),
                            sampled_types[i2]->mass());
        partial_sum_R2 += tmp;
//        std::cout << "Choosing: " << sampled_types[i1]->name() << sampled_types[i2]->name() << " "
//                  << tmp << "/" << sum_R2 << std::endl;
        if (partial_sum_R2 >= r) {
          channel_found32 = true;
        }
        if (partial_sum_R2 > sum_R2) {
          std::cout << partial_sum_R2 << " " << sum_R2 << std::endl;
          throw std::runtime_error("3->2, partial sum exceeded total!");
        }
 
      }
    }
  }
  i1--;
  i2--;
//  std::cout << "Chose: "  << sampled_types[i1]->name() << sampled_types[i2]->name() << std::endl;
  if (!channel_found32) {
    // std::cout << "No channel found, R2_sum = " << sum_R2 << std::endl;
    return;
  }
  ParticleData out1{*(sampled_types[i1])},
               out2{*(sampled_types[i2])};

  const double m1 = out1.type().mass(),
               m2 = out2.type().mass();
  const double p_cm = smash::pCM(srts, m1, m2);
  Angles phitheta;
  phitheta.distribute_isotropically();
  const ThreeVector mom = p_cm * phitheta.threevec();
  const double mom2 = mom.sqr();
  FourVector p1(std::sqrt(mom2 + m1*m1),  mom),
             p2(std::sqrt(mom2 + m2*m2), -mom);
  out1.set_4momentum(p1);
  out2.set_4momentum(p2);
  out1.boost_momentum(-beta_cm);
  out2.boost_momentum(-beta_cm);

  const double E1_pr = in1.momentum().x0(),
               E2_pr = in2.momentum().x0(),
               E3_pr = in3.momentum().x0();
  const double E1 = out1.momentum().x0(),
               E2 = out2.momentum().x0();
  const double phase_space_factor = sum_R3 / sum_R2;
  const double energy_factor = 2.0 * E1_pr * E2_pr * E3_pr / (E1 * E2);
  ParticleTypePtrList in_types {&(in1.type()), &(in2.type()), &(in3.type())},
                      out_types {&(out1.type()), &(out2.type())};
  double spin_factor = 1.0;
  for (const auto t : in_types) {
    spin_factor *= static_cast<double>(1 + t->spin());
  }
  for (const auto t : out_types) {
    spin_factor /= static_cast<double>(1 + t->spin());
  }

  double w = V * pow_int(1.0 / smash::hbarc, 3) * spin_factor *
             energy_factor * phase_space_factor;
  w = 1.0 / w;
  double identity_factor = 1.0;
  if (out1.type() == out2.type()) {
    identity_factor /= 2;
  }
  if (in1.type() == in2.type() && in2.type() == in3.type()) {
    identity_factor *= 6;
  } else if (in1.type() == in2.type() || in2.type() == in3.type()) {
    identity_factor *= 2;
  }
  w *= particles.size() / 3 * identity_factor;
/*
  std::cout << "Trying 3->2 (" << in1.type().name() << in2.type().name()
            << in3.type().name() << "->"
            << out1.type().name() << out2.type().name()
            << ") with w = " << w << ", sym_factor: " << sym_factor
            << ", spin_factor = " << spin_factor << ", ph.sp. factor = "
            << phase_space_factor << ", en. factor = " << energy_factor
            <<", srts = " << srts << std::endl;
*/
  w = std::min(1.0, w);

  if (random::canonical() < w) {
    particles[in_index1] = out1;
    particles[in_index2] = out2;
    particles.erase(particles.begin() + in_index3);
/*
    std::cout << in1.type().name() << "," << in2.type().name() 
              << "," << in3.type().name() << "->"
              << out1.type().name() << "," << out2.type().name()
              << " ; " << srts << std::endl;
*/
  }
/*
  std::cout << type_count(particles, &ParticleType::find(pdg::pi_z)) << " pi0 /"
            << type_count(particles, &ParticleType::find(pdg::pi_p)) << " pi+ /"
            << type_count(particles, &ParticleType::find(pdg::pi_m)) << " pi-" << std::endl;
*/
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

  // Prepare 3-body integrals
  std::cout << "Preparing 3-body integrals" << std::endl;
  ThreeBodyIntegrals three_body_int;
  std::vector<double> unique_masses;
  for (const ParticleTypePtr t : sampled_types) {
    unique_masses.push_back(t->mass());
  }
  std::sort(unique_masses.begin(), unique_masses.end());
  unique_masses.erase(std::unique(unique_masses.begin(), unique_masses.end()),
                      unique_masses.end() );
  const double N_masses = unique_masses.size();
  std::cout << N_masses << " different masses in total" << std::endl;
  for (unsigned int i = 0; i < N_masses; i++) {
    for (unsigned int j = 0; j <= i; j++) {
      for (unsigned int k = 0; k <= j; k++) {
        three_body_int.add(unique_masses[i], unique_masses[j], unique_masses[k]);
      }
    }
  }
  
  std::cout << "Constructing initial configuration." << std::endl;
  const double E_tot = 500.0;  // GeV
  const double V = 1000.0;    // fm^3

  ParticleList particles;
  ParticleData pi0{ParticleType::find(pdg::pi_z)};
  constexpr int Npart_init = 200;
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
/*
    for (const auto &part : particles) {
      std::cout << part.momentum().abs() << " " << part.momentum().x0() << std::endl;
    }
*/

    for (const ParticleTypePtr t : sampled_types) {
      std::cout << type_count(particles, t) << " ";
    }
    std::cout << std::endl;
  }

  QuantumNumbers cons(particles);
  std::cout << "Final momentum: " << cons.momentum() << std::endl;
}
