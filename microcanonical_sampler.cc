#include <vector>

#include "microcanonical_sampler.h"
#include "threebody_integrals.h"
#include "hydro_cells.h"

#include "smash/particletype.h"
#include "smash/random.h"
#include "smash/decaymodes.h"
#include "smash/angles.h"
#include "smash/constants.h"
#include "smash/quantumnumbers.h"
#include "smash/pow.h"
#include "smash/isoparticletype.h"
#include "smash/kinematics.h"

MicrocanonicalSampler::MicrocanonicalSampler(
    const std::function<bool(const smash::ParticleTypePtr)> &is_sampled,
    int debug_printout,
    bool quantum_statistics)
  : debug_printout_(debug_printout),
    quantum_statistics_(quantum_statistics) {
  // Initialize sampled types
  if (debug_printout_) {
    std::cout << "Initializing sampled types..." << std::endl;
  }
  for (const smash::ParticleType& ptype : smash::ParticleType::list_all()) {
    if (is_sampled(&ptype)) {
      sampled_types_.push_back(&ptype);
    }
  }

  if (debug_printout_) {
    std::cout << "Filling the map of channels..." << std::endl;
  }
  // Create lists of channels by qunatum numbers
  const int Ntypes = sampled_types_.size();
  for (int i1 = 0; i1 < Ntypes; ++i1) {
    for (int i2 = 0; i2 <= i1; ++i2) {
      std::array<smash::ParticleTypePtr, 2> y{sampled_types_[i1],
                                              sampled_types_[i2]};
      std::array<int,3> BSQ2{y[0]->baryon_number() + y[1]->baryon_number(),
                             y[0]->strangeness() +   y[1]->strangeness(),
                             y[0]->charge() +        y[1]->charge()};
      channels2_[BSQ2].push_back(y);
      for (int i3 = 0; i3 <= i2; ++i3) {
        std::array<smash::ParticleTypePtr, 3> x{sampled_types_[i1],
                                                sampled_types_[i2],
                                                sampled_types_[i3]};
        const int B = x[0]->baryon_number() + x[1]->baryon_number() +
                      x[2]->baryon_number();
        const int S = x[0]->strangeness() + x[1]->strangeness() +
                      x[2]->strangeness();
        const int Q = x[0]->charge() + x[1]->charge() + x[2]->charge();
        std::array<int,3> BSQ3{B,S,Q};
        channels3_[BSQ3].push_back(x);
      }
    }
  }

  // Sort lists of channels by sum of masses = threshold energy
  for (auto& channels_BSQ : channels3_) {
    std::sort(channels_BSQ.second.begin(), channels_BSQ.second.end(),
              [](const std::array<smash::ParticleTypePtr, 3> &a,
                 const std::array<smash::ParticleTypePtr, 3> &b) {
      return a[0]->mass() + a[1]->mass() + a[2]->mass() <
             b[0]->mass() + b[1]->mass() + b[2]->mass();
    });
  }
  for (auto& channels_BSQ : channels2_) {
    std::sort(channels_BSQ.second.begin(), channels_BSQ.second.end(),
              [](const std::array<smash::ParticleTypePtr, 2> &a,
                 const std::array<smash::ParticleTypePtr, 2> &b) {
      return a[0]->mass() + a[1]->mass() < b[0]->mass() + b[1]->mass();
    });
  }

  if (debug_printout_) {
    std::cout << "Total numbers of 2->3 channels with given quantum numbers. "
              << std::endl;
    int total_channels = 0;
    for (auto channels_BSQ : channels3_) {
      std::array<int, 3> BSQ = channels_BSQ.first;
      const int Nch = channels_BSQ.second.size();
      std::cout << "B = " << BSQ[0] << ", S = " << BSQ[1] << ", Q = " << BSQ[2]
                << ": Nch = " << Nch << std::endl;
      /*
      for (const auto t : channels_BSQ.second) {
        std::cout << t[0]->name() << t[1]->name() << t[2]->name() << " " <<
                     t[0]->mass() + t[1]->mass() + t[2]->mass() << ", ";
      }
      std::cout << std::endl;
      */
      total_channels += Nch;
    }
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "Total 2->3 channels: " << total_channels << std::endl;
    std::cout << "Total numbers of 3->2 channels with given quantum numbers. "
              << std::endl;
    total_channels = 0;
    for (auto channels_BSQ : channels2_) {
      std::array<int, 3> BSQ = channels_BSQ.first;
      const int Nch = channels_BSQ.second.size();
      std::cout << "B = " << BSQ[0] << ", S = " << BSQ[1] << ", Q = " << BSQ[2]
                << ": Nch = " << Nch << std::endl;
      /*
      for (const auto t : channels_BSQ.second) {
        std::cout << t[0]->name() << t[1]->name() << ", ";
      }
      std::cout << std::endl;
      */
      total_channels += Nch;
    }
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "Total 3->2 channels: " << total_channels << std::endl;
  }

  // Prepare 3-body integrals
  if (debug_printout_) {
    std::cout << "Preparing 3-body integrals..." << std::endl;
  }
  std::vector<double> unique_masses;
  for (const smash::ParticleTypePtr t : sampled_types_) {
    unique_masses.push_back(t->mass());
  }
  std::sort(unique_masses.begin(), unique_masses.end());
  unique_masses.erase(std::unique(unique_masses.begin(), unique_masses.end()),
                      unique_masses.end() );
  const size_t N_masses = unique_masses.size();
  if (debug_printout_) {
    std::cout << N_masses << " different masses in total" << std::endl;
  }
  for (size_t i = 0; i < N_masses; i++) {
    for (size_t j = 0; j <= i; j++) {
      for (size_t k = 0; k <= j; k++) {
        three_body_int_.add(unique_masses[i],
                            unique_masses[j],
                            unique_masses[k]);
      }
    }
  }
}

void MicrocanonicalSampler::initialize(const HyperSurfacePatch& hypersurface) {
  particles_.clear();
  const int B = hypersurface.B(),
            S = hypersurface.S(),
            Q = hypersurface.Q();
  smash::ParticleTypePtr baryon, antibaryon, strange_meson, antistrange_meson,
                         plus_nonstrange_meson, minus_nonstrange_meson,
                         neutral_meson;
  for (smash::ParticleTypePtr t : sampled_types_) {
    if (!baryon && t->baryon_number() > 0) {
      baryon = t;
    }
    if (!antibaryon && t->baryon_number() < 0) {
      antibaryon = t;
    }
    if (!strange_meson && t->baryon_number() == 0
                       && t->strangeness() > 0) {
      strange_meson = t;
    }
    if (!antistrange_meson && t->baryon_number() == 0
                           && t->strangeness() < 0) {
      antistrange_meson = t;
    }
    if (!plus_nonstrange_meson && t->baryon_number() == 0
        && t->strangeness() == 0 && t->charge() > 0) {
      plus_nonstrange_meson = t;
    }
    if (!minus_nonstrange_meson && t->baryon_number() == 0
        && t->strangeness() == 0 && t->charge() < 0) {
      minus_nonstrange_meson = t;
    }
    if (!neutral_meson && t->baryon_number() == 0 &&
        t->strangeness() == 0 && t->charge() == 0) {
      neutral_meson = t;
    }
  }
  int remaining_S = S, remaining_Q = Q;
  if (B > 0) {
    if (!baryon) {
      std::runtime_error("B > 0 requested, but hadron sorts"
                         " to sample do not include baryons");
    }
    for (int i = 0; i < B; i++) {
      particles_.push_back({smash::FourVector(), baryon, 0});
    }
    remaining_S -= baryon->strangeness() * B;
    remaining_Q -= baryon->charge() * B;
  } else if (B < 0) {
    if (!antibaryon) {
      std::runtime_error("B < 0 requested, but hadron sorts"
                         " to sample do not include antibaryons");
    }
    for (int i = 0; i < std::abs(B); i++) {
      particles_.push_back({smash::FourVector(), antibaryon, 0});
    }
    remaining_S -= antibaryon->strangeness() * std::abs(B);
    remaining_Q -= antibaryon->charge() * std::abs(B);
  }
  if (remaining_S > 0) {
    if (!strange_meson) {
      std::runtime_error("S > 0 needed, but hadron sorts"
                         " to sample do not include mesons with S > 0");
    }
    for (int i = 0; i < remaining_S; i++) {
      particles_.push_back({smash::FourVector(), strange_meson, 0});
    }
    remaining_Q -= strange_meson->charge() * remaining_S;
  } else if (remaining_S < 0) {
    if (!antistrange_meson) {
      std::runtime_error("S < 0 needed, but hadron sorts"
                         " to sample do not include mesons with S < 0");
    }
    for (int i = 0; i < std::abs(remaining_S); i++) {
      particles_.push_back({smash::FourVector(), antistrange_meson, 0});
    }
    remaining_Q -= antistrange_meson->charge() * std::abs(remaining_S);
  }
  if (remaining_Q > 0) {
    if (!plus_nonstrange_meson) {
      std::runtime_error("Q > 0 needed, but hadron sorts to sample"
                         " do not include nonstrange mesons with Q > 0");
    }
    for (int i = 0; i < remaining_Q; i++) {
      particles_.push_back({smash::FourVector(), plus_nonstrange_meson, 0});
    }
  } else if (remaining_Q < 0) {
    if (!minus_nonstrange_meson) {
      std::runtime_error("Q < 0 needed, but hadron sorts to sample"
                         " do not include nonstrange mesons with Q < 0");
    }
    for (int i = 0; i < std::abs(remaining_Q); i++) {
      particles_.push_back({smash::FourVector(), minus_nonstrange_meson, 0});
    }
  }
  if (neutral_meson) {
    for (int i = 0; i < 10; i++) {
      particles_.push_back({smash::FourVector(), neutral_meson, 0});
    }
  }

  // Assign some random momenta
  for (SamplerParticle& part : particles_) {
    const double m = part.type->mass();
    const double pabs = 0.1;
    smash::Angles phitheta;
    phitheta.distribute_isotropically();
    const smash::FourVector mom(std::sqrt(m*m + pabs*pabs),
                                pabs * phitheta.threevec());
    part.momentum = mom;
  }
  renormalize_momenta(hypersurface.pmu());

  QuantumNumbers cons(particles_);
  assert(cons.B == B);
  assert(cons.S == S);
  assert(cons.Q == Q);
}

void MicrocanonicalSampler::one_markov_chain_step(
                            const HyperSurfacePatch& hypersurface) {
  if (smash::random::uniform_int(0, 1)) {
    random_two_to_three(hypersurface);
  } else {
    random_three_to_two(hypersurface);
  }
}

double MicrocanonicalSampler::compute_R2(double srts, double m1, double m2) {
  return (srts > m1 + m2) ? smash::pCM(srts, m1, m2) / (4.0 * M_PI * srts)
                          : 0.0;
}

void MicrocanonicalSampler::sample_3body_phase_space(double srts,
                              SamplerParticle &a,
                              SamplerParticle &b,
                              SamplerParticle &c) {
  const double m_a = a.type->mass(),
               m_b = b.type->mass(),
               m_c = c.type->mass();
  // sample mab from pCM(sqrt, mab, mc) pCM (mab, ma, mb) <= sqrts^2/4
  double mab, r, probability, pcm_ab, pcm;
  do {
    mab = smash::random::uniform(m_a + m_b, srts - m_c);
    r = smash::random::canonical();
    pcm = smash::pCM(srts, mab, m_c);
    pcm_ab = smash::pCM(mab, m_a, m_b);
    probability = pcm * pcm_ab * 4 / (srts * srts);
  } while (r > probability);
  smash::Angles phitheta;
  phitheta.distribute_isotropically();
  c.momentum = smash::FourVector(std::sqrt(m_c * m_c + pcm * pcm),
                                 pcm * phitheta.threevec());
  const smash::ThreeVector beta_cm =
      pcm * phitheta.threevec() / std::sqrt(pcm * pcm + mab * mab);

  phitheta.distribute_isotropically();
  a.momentum = smash::FourVector(std::sqrt(m_a * m_a + pcm_ab * pcm_ab),
                                 pcm_ab * phitheta.threevec());
  b.momentum = smash::FourVector(std::sqrt(m_b * m_b + pcm_ab * pcm_ab),
                                -pcm_ab * phitheta.threevec());
  a.momentum = a.momentum.LorentzBoost(beta_cm);
  b.momentum = b.momentum.LorentzBoost(beta_cm);
  if (debug_printout_ == 3) {
    std::cout << "Sample 3-body phase space: total momentum = "
              << a.momentum + b.momentum + c.momentum << std::endl;
  }
  // assert(std::abs((a.momentum + b.momentum + c.momentum
  //                 - smash::FourVector(srts, 0, 0, 0)).abs()) < 1.e-12);
}

double MicrocanonicalSampler::mu_minus_E_over_T(const SamplerParticle& p,
    const HyperSurfacePatch& hypersurface) {
    const HyperSurfacePatch::hydro_cell c = hypersurface.cell(p.cell_index);
    const double mu = c.muB * p.type->baryon_number() +
                      c.muS * p.type->strangeness() +
                      c.muQ * p.type->charge();
    return (mu - p.momentum.Dot(c.u)) / c.T;
}

double MicrocanonicalSampler::compute_cells_factor(
    const SamplerParticleList& in, const SamplerParticleList& out,
    const HyperSurfacePatch& hypersurface) {
  double cells_factor = 0.0;
  if (!quantum_statistics_) {
    for (const SamplerParticle& p : out) {
      cells_factor += mu_minus_E_over_T(p, hypersurface);
    }
    for (const SamplerParticle& p : in) {
      cells_factor -= mu_minus_E_over_T(p, hypersurface);
    }
    cells_factor = std::exp(cells_factor);
  } else {
    cells_factor = 1.0;
    for (const SamplerParticle& p : in) {
      const double x = std::exp(- mu_minus_E_over_T(p, hypersurface));
      cells_factor *= p.type->is_meson() ? (x - 1) : (x + 1);
    }
    for (const SamplerParticle& p : out) {
      const double x = std::exp(- mu_minus_E_over_T(p, hypersurface));
      cells_factor /= p.type->is_meson() ? (x - 1) : (x + 1);
    }
  }
  return cells_factor;
}

double MicrocanonicalSampler::compute_spin_factor(
    const SamplerParticleList& in, const SamplerParticleList& out) {
  double spin_factor = 1.0;
  for (const SamplerParticle& part : out) {
    spin_factor *= part.type->pdgcode().spin_degeneracy();
  }
  for (const SamplerParticle& part : in) {
    spin_factor /= part.type->pdgcode().spin_degeneracy();
  }
  return spin_factor;
}

size_t MicrocanonicalSampler::N_available_channels3(std::array<int, 3>& BSQ,
                                                    double srts) {
  size_t N_available_channels = 0;
  if (channels3_.find(BSQ) != channels3_.end()) {
    for (const std::array<smash::ParticleTypePtr,3> &t : channels3_[BSQ]) {
      if (t[0]->mass() + t[1]->mass() + t[2]->mass() > srts) {
        break;
      }
      N_available_channels++;
    }
  }
  return N_available_channels;
}

size_t MicrocanonicalSampler::N_available_channels2(std::array<int, 3>& BSQ,
                                                    double srts) {
  size_t N_available_channels = 0;
  if (channels2_.find(BSQ) != channels2_.end()) {
    for (const std::array<smash::ParticleTypePtr,2> &t : channels2_[BSQ]) {
      if (t[0]->mass() + t[1]->mass() > srts) {
        break;
      }
      N_available_channels++;
    }
  }
  return N_available_channels;
}


void MicrocanonicalSampler::random_two_to_three(
                            const HyperSurfacePatch& hypersurface) {
  const size_t Npart = particles_.size();
  const size_t Ncells = hypersurface.Ncells();
  // Choose 2 random particles
  size_t in_index1, in_index2;
  in_index1 = smash::random::uniform_int<size_t>(0, Npart - 1);
  do {
    in_index2 = smash::random::uniform_int<size_t>(0, Npart - 1);
  } while (in_index1 == in_index2);
  SamplerParticleList in{particles_[in_index1], particles_[in_index2]};
  QuantumNumbers cons(in);
  const smash::FourVector momentum_total = cons.momentum;
  const double srts = momentum_total.abs();
  const smash::ThreeVector beta_cm = momentum_total.velocity();

  // Choose channel
  std::array<int, 3> BSQ{cons.B, cons.S, cons.Q};
  size_t N3 = N_available_channels3(BSQ, srts),
         N2 = N_available_channels2(BSQ, srts);

  if (N3 == 0 || N2 == 0) {
    if (debug_printout_) {
      std::cout << "Not enough energy ("
                << srts << " GeV) for 2->3" << std::endl;
    }
    rejected23_count_++;
    return;
  }

  SamplerParticleList out;
  out.clear();
  out.resize(3);
  std::array<size_t,3> cell =
      {smash::random::uniform_int<size_t>(0, Ncells - 1),
       smash::random::uniform_int<size_t>(0, Ncells - 1),
       smash::random::uniform_int<size_t>(0, Ncells - 1)};
  size_t i_channel = smash::random::uniform_int<size_t>(0, N3 - 1);
  std::array<smash::ParticleTypePtr, 3> out_types = channels3_[BSQ][i_channel];
  const double R3 = three_body_int_.value(srts,
       out_types[0]->mass(), out_types[1]->mass(), out_types[2]->mass());
  const double R2 = compute_R2(srts, in[0].type->mass(), in[1].type->mass());
  out[0] = {smash::FourVector(), out_types[0], cell[0]};
  out[1] = {smash::FourVector(), out_types[1], cell[1]};
  out[2] = {smash::FourVector(), out_types[2], cell[2]};


  sample_3body_phase_space(srts, out[0], out[1], out[2]);
  for (SamplerParticle& part : out) {
    part.momentum = part.momentum.LorentzBoost(-beta_cm);
  }
  if (debug_printout_ == 3) {
    std::cout << "In: " << in[0].momentum << " "
                        << in[1].momentum << std::endl;
    std::cout << "Out: " << out[0].momentum << " "
                         << out[1].momentum << " "
                         << out[2].momentum << std::endl;
  }
  // At this point the proposal function is set
  // Further is the calculation of probability to accept it
  const double phase_space_factor = R3 / R2 * N3 / N2;
  const double spin_factor = compute_spin_factor(in, out);

  double combinatorial_factor = 3.0 / (Npart + 1);
  if (out[0].type == out[1].type && out[1].type == out[2].type) {
    combinatorial_factor /= 6;
  } else if (out[0].type == out[1].type || out[1].type == out[2].type) {
    combinatorial_factor /= 2;
  }
  if (in[0].type == in[1].type) {
    combinatorial_factor *= 2;
  }

  double energy_factor = 2.0 * Ncells / smash::pow_int(smash::hbarc, 3);
  for (const SamplerParticle& p : out) {
    energy_factor *= p.momentum.Dot(hypersurface.cell(p.cell_index).dsigma);
  }
  for (const SamplerParticle& p : in) {
    energy_factor /= p.momentum.Dot(hypersurface.cell(p.cell_index).dsigma);
  }

  double cells_factor = compute_cells_factor(in, out, hypersurface);


  const double w = phase_space_factor * spin_factor * combinatorial_factor *
                   energy_factor * cells_factor;

  if (debug_printout_ == 1) {
    std::cout << "Trying 2->3 (" << in << "->" << out
              << ") with w = " << w
              << "; ph. sp. factor * 1e6 = " << phase_space_factor * 1.e6
              << "; spin_factor = " << spin_factor
              << "; comb. factor = " << combinatorial_factor
              << "; energy factor = " << energy_factor
              << "; cell factor = " << cells_factor
              << "; srts = " << srts << std::endl;
  }

  // No need to do it: w = std::min(1.0, w);
  if (smash::random::canonical() < w) {
    particles_[in_index1] = out[0];
    particles_[in_index2] = out[1];
    particles_.push_back(out[2]);
    accepted23_count_++;
    if (debug_printout_ == 1) {
      std::cout << in << "->" << out << ";" << srts << std::endl;
    }
  } else {
    rejected23_count_++;
  }
}


void MicrocanonicalSampler::random_three_to_two(
                            const HyperSurfacePatch& hypersurface) {
  const size_t Npart = particles_.size();
  const size_t Ncells = hypersurface.Ncells();
  if (Npart < 3) {
    return;
  }
  // Choose 3 random particles
  size_t in_index1, in_index2, in_index3;
  in_index1 = smash::random::uniform_int<size_t>(0, Npart - 1);
  do {
    in_index2 = smash::random::uniform_int<size_t>(0, Npart - 1);
  } while (in_index1 == in_index2);
  do {
    in_index3 = smash::random::uniform_int<size_t>(0, Npart - 1);
  } while (in_index3 == in_index1 || in_index3 == in_index2);
  SamplerParticleList in{particles_[in_index1],
                         particles_[in_index2],
                         particles_[in_index3]};
  QuantumNumbers cons(in);
  smash::FourVector momentum_total = cons.momentum;
  const double srts = momentum_total.abs();
  const smash::ThreeVector beta_cm = momentum_total.velocity();

  // Choose channel
  std::array<int, 3> BSQ{cons.B, cons.S, cons.Q};
  size_t N3 = N_available_channels3(BSQ, srts),
         N2 = N_available_channels2(BSQ, srts);
  if (N3 == 0 || N2 == 0) {
    if (debug_printout_) {
      std::cout << "Not enough energy (" << srts << " GeV) or wrong quantum"
                   " numbers for 3->2" << std::endl;
    }
    rejected32_count_++;
    return;
  }

  SamplerParticleList out;
  out.clear();
  out.resize(2);
  std::array<size_t,2> cell =
    {smash::random::uniform_int<size_t>(0, Ncells - 1),
     smash::random::uniform_int<size_t>(0, Ncells - 1)};
  size_t i_channel = smash::random::uniform_int<size_t>(0, N2 - 1);
  std::array<smash::ParticleTypePtr, 2> out_types = channels2_[BSQ][i_channel];
  const double R3 = three_body_int_.value(srts,
      in[0].type->mass(), in[1].type->mass(), in[2].type->mass());
  const double R2 =
      compute_R2(srts, out_types[0]->mass(), out_types[1]->mass());
  out[0] = {smash::FourVector(), out_types[0], cell[0]};
  out[1] = {smash::FourVector(), out_types[1], cell[1]};

  const double m1 = out[0].type->mass(),
               m2 = out[1].type->mass();
  const double p_cm = smash::pCM(srts, m1, m2);
  smash::Angles phitheta;
  phitheta.distribute_isotropically();
  const smash::ThreeVector mom = p_cm * phitheta.threevec();
  const double mom2 = p_cm*p_cm;
  smash::FourVector p1(std::sqrt(mom2 + m1*m1),  mom),
                    p2(std::sqrt(mom2 + m2*m2), -mom);
  out[0].momentum = p1.LorentzBoost(-beta_cm);
  out[1].momentum = p2.LorentzBoost(-beta_cm);

  // At this point the proposal function is set
  // Further is the calculation of probability to accept it
  const double phase_space_factor = R2 / R3 * N2 / N3;
  const double spin_factor = compute_spin_factor(in, out);

  double combinatorial_factor = static_cast<double>(Npart) / 3.0;
  if (out[0].type == out[1].type) {
    combinatorial_factor /= 2;
  }
  if (in[0].type == in[1].type && in[1].type == in[2].type) {
    combinatorial_factor *= 6;
  } else if (in[0].type == in[1].type || in[1].type == in[2].type) {
    combinatorial_factor *= 2;
  }

  double energy_factor = smash::pow_int(smash::hbarc, 3) / (2.0 * Ncells);
  for (const SamplerParticle& p : out) {
    energy_factor *= p.momentum.Dot(hypersurface.cell(p.cell_index).dsigma);
  }
  for (const SamplerParticle& p : in) {
    energy_factor /= p.momentum.Dot(hypersurface.cell(p.cell_index).dsigma);
  }

  double cells_factor = compute_cells_factor(in, out, hypersurface);


  const double w = phase_space_factor * spin_factor * combinatorial_factor *
                   energy_factor * cells_factor;

  if (debug_printout_ == 1) {
    std::cout << "Trying 3->2 (" << in << "->" << out
              << ") with w = " << w
              << "; ph. sp. factor * 1e6 = " << phase_space_factor * 1.e6
              << "; spin_factor = " << spin_factor
              << "; comb. factor = " << combinatorial_factor
              << "; energy factor = " << energy_factor
              << "; cell factor = " << cells_factor
              << "; srts = " << srts << std::endl;
  }

  // No need to do it: w = std::min(1.0, w);
  if (smash::random::canonical() < w) {
    particles_[in_index1] = out[0];
    particles_[in_index2] = out[1];
    particles_.erase(particles_.begin() + in_index3);
    accepted32_count_++;
    if (debug_printout_) {
      std::cout << in << "->" << out << ";" << srts << std::endl;
    }
  } else {
    rejected32_count_++;
  }
}


void MicrocanonicalSampler::renormalize_momenta(
                            const smash::FourVector& required_total_momentum) {

  // Centralize momenta
  QuantumNumbers conserved = QuantumNumbers(particles_);
  std::cout << particles_.size() << " particles. " << std::endl;
  std::cout << "Required 4-momentum: " << required_total_momentum << std::endl;
  std::cout << "Sampled 4-momentum: "  << conserved.momentum << std::endl;
  const smash::ThreeVector mom_to_add =
      (required_total_momentum.threevec() - conserved.momentum.threevec()) /
      particles_.size();
  std::cout << "Adjusting momenta by " << mom_to_add << std::endl;
  for (auto &particle : particles_) {
    const double m = particle.type->mass();
    const smash::ThreeVector mom3 = particle.momentum.threevec() + mom_to_add;
    const double Energy = std::sqrt(m*m + mom3.sqr());
    particle.momentum = smash::FourVector(Energy, mom3);
  }

  // Boost every particle to the common center of mass frame
  conserved = QuantumNumbers(particles_);
  const smash::ThreeVector beta_CM_generated = conserved.momentum.velocity();
  const smash::ThreeVector beta_CM_required = required_total_momentum.velocity();

  double E = 0.0;
  double E_expected = required_total_momentum.abs();
  for (auto &particle : particles_) {
    particle.momentum = particle.momentum.LorentzBoost(beta_CM_generated);
    E += particle.momentum.x0();
  }
  // Renorm. momenta by factor (1+a) to get the right energy, binary search
  const double tolerance = smash::really_small;
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
    for (const auto &particle : particles_) {
      const double p2 = particle.momentum.threevec().sqr();
      const double E2 = particle.momentum.x0() * particle.momentum.x0();
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
  for (auto &particle : particles_) {
    const double m = particle.type->mass();
    smash::ThreeVector mom3 = particle.momentum.threevec();
    mom3 *= (1 + a);
    const double Energy = std::sqrt(m*m + mom3.sqr());
    particle.momentum = smash::FourVector(Energy, mom3);
    particle.momentum = particle.momentum.LorentzBoost(-beta_CM_required);
  }
  QuantumNumbers conserved_final = QuantumNumbers(particles_);
  std::cout << "Obtained total momentum: "
            << conserved_final.momentum << std::endl;
}

void MicrocanonicalSampler::print_rejection_stats() {
  std::cout << "2->3 acceptance: " << accepted23_count_ << "/"
            << accepted23_count_ + rejected23_count_ << std::endl;
  std::cout << "3->2 acceptance: " << accepted32_count_ << "/"
            << accepted32_count_ + rejected32_count_ << std::endl;
}

std::ostream &operator<<(std::ostream &out,
    const MicrocanonicalSampler::SamplerParticleList &list) {
  const size_t N = list.size();
  out << list[0].type->name();
  for (size_t i = 1; i < N; i++) {
    out << "," << list[i].type->name();
  }
  return out;
}
