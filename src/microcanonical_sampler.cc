#include <fstream>
#include <vector>
#include <iostream>

#include "microcanonical_sampler/hydro_cells.h"
#include "microcanonical_sampler/microcanonical_sampler.h"

#include <gsl/gsl_sf_ellint.h>

#include "smash/angles.h"
#include "smash/constants.h"
#include "smash/kinematics.h"
#include "smash/pow.h"
#include "smash/random.h"

MicrocanonicalSampler::MicrocanonicalSampler(
    const std::function<bool(const ParticleTypePtr)> &is_sampled,
    int debug_printout, bool quantum_statistics)
    : debug_printout_(debug_printout), quantum_statistics_(quantum_statistics) {
  std::cout << "Initializing the sampler" << std::endl;
  // Initialize sampled types
  if (debug_printout_) {
    std::cout << "Initializing sampled types..." << std::endl;
  }
  for (const ParticleType &ptype : ParticleType::list_all()) {
    if (is_sampled(&ptype)) {
      sampled_types_.push_back(&ptype);
    }
  }

  if (debug_printout_) {
    std::cout << "Filling the map of channels..." << std::endl;
  }
  // Create lists of channels by quantum numbers
  const int Ntypes = sampled_types_.size();
  for (int i1 = 0; i1 < Ntypes; ++i1) {
    for (int i2 = 0; i2 <= i1; ++i2) {
      std::array<ParticleTypePtr, 2> y{sampled_types_[i1],
                                              sampled_types_[i2]};
      std::array<int, 3> BSQ2{y[0]->baryon_number() + y[1]->baryon_number(),
                              y[0]->strangeness() + y[1]->strangeness(),
                              y[0]->charge() + y[1]->charge()};
      channels2_[BSQ2].push_back(y);
      for (int i3 = 0; i3 <= i2; ++i3) {
        std::array<ParticleTypePtr, 3> x{
            sampled_types_[i1], sampled_types_[i2], sampled_types_[i3]};
        const int B = x[0]->baryon_number() + x[1]->baryon_number() +
                      x[2]->baryon_number();
        const int S =
            x[0]->strangeness() + x[1]->strangeness() + x[2]->strangeness();
        const int Q = x[0]->charge() + x[1]->charge() + x[2]->charge();
        std::array<int, 3> BSQ3{B, S, Q};
        channels3_[BSQ3].push_back(x);
      }
    }
  }

  // Sort lists of channels by sum of masses = threshold energy
  for (auto &channels_BSQ : channels3_) {
    std::sort(channels_BSQ.second.begin(), channels_BSQ.second.end(),
              [](const std::array<ParticleTypePtr, 3> &a,
                 const std::array<ParticleTypePtr, 3> &b) {
                return a[0]->mass() + a[1]->mass() + a[2]->mass() <
                       b[0]->mass() + b[1]->mass() + b[2]->mass();
              });
  }
  for (auto &channels_BSQ : channels2_) {
    std::sort(channels_BSQ.second.begin(), channels_BSQ.second.end(),
              [](const std::array<ParticleTypePtr, 2> &a,
                 const std::array<ParticleTypePtr, 2> &b) {
                return a[0]->mass() + a[1]->mass() <
                       b[0]->mass() + b[1]->mass();
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

  // Prepare threshold arrays
  for (auto channels_BSQ : channels3_) {
    std::array<int, 3> BSQ = channels_BSQ.first;
    thresholds3_[BSQ].reserve(channels_BSQ.second.size());
    for (const auto t : channels_BSQ.second) {
      thresholds3_[BSQ].push_back(t[0]->mass() + t[1]->mass() + t[2]->mass());
    }
  }
  for (auto channels_BSQ : channels2_) {
    std::array<int, 3> BSQ = channels_BSQ.first;
    thresholds2_[BSQ].reserve(channels_BSQ.second.size());
    for (const auto t : channels_BSQ.second) {
      thresholds2_[BSQ].push_back(t[0]->mass() + t[1]->mass());
    }
  }
}

void MicrocanonicalSampler::initialize(const HyperSurfacePatch &hypersurface,
                                       SamplerParticleList &particles) {
  // Find lightest types with every combination of quantum numbers
  std::map<std::array<int,3>, ParticleTypePtr> lightest_species_BSQ;
  for (ParticleTypePtr t : sampled_types_) {
    std::array<int,3> BSQ {t->baryon_number(), t->strangeness(), t->charge()};
    if (lightest_species_BSQ.find(BSQ) == lightest_species_BSQ.end() ||
        lightest_species_BSQ[BSQ]->mass() > t->mass()) {
      lightest_species_BSQ[BSQ] = t;
    }
  }

  int B = hypersurface.B(), S = hypersurface.S(), Q = hypersurface.Q();
  std::array<int,3> neutral_meson_BSQ{0, 0, 0};
  if (lightest_species_BSQ.find(neutral_meson_BSQ) !=
      lightest_species_BSQ.end()) {
    std::runtime_error("The sampling method requires"
        " presence of a neutral non-strange meson in the species list.");
  }
  std::map<ParticleTypePtr, size_t> specie_quantity;
  /* Apply heuristics to fill the required quantum numbers with possibly lower
   * total energy. This heuristics is not guaranteed to work in general case,
   * but it seems to be fine for a real list of hadronic species.
   */
  while (B || S || Q) {
    int signB = (B > 0) - (B < 0),
        signS = (S > 0) - (S < 0),
        signQ = (Q > 0) - (Q < 0);
    std::vector<std::array<int,3>> BSQ_to_try = {{signB, signS, signQ},
     {signB, signS, 0}, {signB, 0, signQ}, {0, signS, signQ},
     {signB, 0, 0}, {0, signS, 0}, {0, 0, signQ}};
    ParticleTypePtr t_minimizing =
        lightest_species_BSQ[neutral_meson_BSQ];
    for (const auto BSQ : BSQ_to_try) {
      if (lightest_species_BSQ.find(BSQ) != lightest_species_BSQ.end()) {
        t_minimizing = lightest_species_BSQ[BSQ];
        break;
      }
    }
    if (t_minimizing == lightest_species_BSQ[neutral_meson_BSQ]) {
      std::runtime_error("Initialization heuristics failed."
          " Can't combine given species to obtain required quantum numbers.");
    }
    int Bi = t_minimizing->baryon_number(),
        Si = t_minimizing->strangeness(),
        Qi = t_minimizing->charge();
    const std::array<int, 3> possible_mult =
      {(Bi != 0) ? B / Bi : 0, (Si != 0) ? S / Si : 0, (Qi != 0) ? Q / Qi : 0};
    int mult = std::numeric_limits<int>::max();
    for (const int x : possible_mult) {
      if (x > 0 && x < mult) {
        mult = x;
      }
    }
    if (mult == std::numeric_limits<int>::max()) {
      std::runtime_error("Initialization heuristics failed."
          " No positive multiplicity suggested.");
    }
    B -= Bi * mult;
    S -= Si * mult;
    Q -= Qi * mult;
    specie_quantity[t_minimizing] = mult;
    std::cout << t_minimizing->name() << ": " << mult << ", updated BSQ: "
              << B << " " << S << " " << Q << std::endl;
  }
  double m_sum = 0.0;
  for (const auto x : specie_quantity) {
    m_sum += x.first->mass() * x.second;
  }
  std::cout << "Sum of masses: " << m_sum << std::endl;
  for (const auto x : specie_quantity) {
    for (size_t i = 0; i < x.second; i++) {
      const size_t cell = smash::random::uniform_int(size_t(),
                              hypersurface.Ncells() - 1);
      particles.push_back({smash::FourVector(), x.first, cell, false});
    }
  }
  if (particles.size() < 2 && hypersurface.pmu().abs() >
      2 * lightest_species_BSQ[neutral_meson_BSQ]->mass()) {
    for (size_t i = 0; i < 2; i++) {
      const size_t cell = smash::random::uniform_int(size_t(),
                              hypersurface.Ncells() - 1);
      particles.push_back({smash::FourVector(),
                          lightest_species_BSQ[neutral_meson_BSQ], cell,
                          false});
    }
  }
  // Assign some random momenta
  for (SamplerParticle &part : particles) {
    const double m = part.type->mass();
    const double pabs = 0.01;
    smash::Angles phitheta;
    phitheta.distribute_isotropically();
    const smash::FourVector mom(std::sqrt(m * m + pabs * pabs),
                                pabs * phitheta.threevec());
    part.momentum = mom;
  }

  if (particles.size() > 0) {
    renormalize_momenta(hypersurface.pmu(), particles);
    QuantumNumbers cons(particles);
    assert(cons.B == hypersurface.B());
    assert(cons.S == hypersurface.S());
    assert(cons.Q == hypersurface.Q());
  }
}

void MicrocanonicalSampler::one_markov_chain_step(
    const HyperSurfacePatch &hypersurface,
    SamplerParticleList &particles) {
  if (smash::random::uniform_int(0, 1)) {
    random_two_to_three(hypersurface, particles);
  } else {
    random_three_to_two(hypersurface, particles);
  }
}

double
MicrocanonicalSampler::compute_R2(double srts, double m1, double m2) {
  return (srts > m1 + m2) ? smash::pCM(srts, m1, m2) / (4.0 * M_PI * srts)
                          : 0.0;
}

double MicrocanonicalSampler::compute_R3(double srts, double m1, double m2,
                                          double m3) {
  if (srts < m1 + m2 + m3) {
    return 0.0;
  }
  const double x1 = (m1 - m2) * (m1 - m2), x2 = (m1 + m2) * (m1 + m2),
               x3 = (srts - m3) * (srts - m3), x4 = (srts + m3) * (srts + m3);
  const double qmm = x3 - x1, qmp = x3 - x2, qpm = x4 - x1, qpp = x4 - x2;
  const double kappa = std::sqrt(qpm * qmp / (qpp * qmm));
  const double tmp = std::sqrt(qmm * qpp);
  const double c1 =
      4.0 * m1 * m2 * std::sqrt(qmm / qpp) * (x4 - m3 * srts + m1 * m2);
  const double c2 = 0.5 * (m1 * m1 + m2 * m2 + m3 * m3 + srts * srts) * tmp;
  const double c3 = 8 * m1 * m2 / tmp *
                    ((m1 * m1 + m2 * m2) * (m3 * m3 + srts * srts) -
                     2 * m1 * m1 * m2 * m2 - 2 * m3 * m3 * srts * srts);
  const double c4 =
      -8 * m1 * m2 / tmp * smash::pow_int(srts * srts - m3 * m3, 2);
  const double precision = 1.e-6;
  const double res =
      c1 * gsl_sf_ellint_Kcomp(kappa, precision) +
      c2 * gsl_sf_ellint_Ecomp(kappa, precision) +
      c3 * gsl_sf_ellint_Pcomp(kappa, -qmp / qmm, precision) +
      c4 * gsl_sf_ellint_Pcomp(kappa, -x1 * qmp / (x2 * qmm), precision);
  return res / (128.0 * M_PI * M_PI * M_PI * srts * srts);
}

void MicrocanonicalSampler::test_3body_integrals() {
  const double my_value1 = compute_R3(5.0, 1.2, 1.3, 1.4) * 1E8;
  const double mathematica_value1 = 0.000326309 * 1E8;
  const double my_value2 = compute_R3(1.0, 0.27, 0.3, 0.4) * 1E8;
  const double mathematica_value2 = 26.4663;
  std::cout << "Test 1: "
            << my_value1 << " == " << mathematica_value1 << std::endl;
  std::cout << "Test 2: "
            << my_value2 << " == " << mathematica_value2 << std::endl;
}


void MicrocanonicalSampler::sample_3body_phase_space(double srts,
                                                    SamplerParticle &a,
                                                    SamplerParticle &b,
                                                    SamplerParticle &c) {
  const double m_a = a.type->mass(), m_b = b.type->mass(), m_c = c.type->mass();
  // sample mab from pCM(sqrt, mab, mc) pCM (mab, ma, mb)
  double mab, probability, pcm_ab_sqr, pcm_sqr;
  do {
    mab = smash::random::uniform(m_a + m_b, srts - m_c);
    pcm_sqr = smash::pCM_sqr(srts, mab, m_c);
    // pCM_sqr decreases as a function of the second argument ->
    //         has maximum at minimal mab = m_a + m_b.
    const double pcm_sqr_max = smash::pCM_sqr(srts, m_a + m_b, m_c);
    pcm_ab_sqr = smash::pCM_sqr(mab, m_a, m_b);
    // pCM_sqr increases as a function of the first argument ->
    //         has maximum at maximal mab = srts - m_c.
    const double pcm_ab_sqr_max = smash::pCM_sqr(srts - m_c, m_a, m_b);

    probability = std::sqrt(pcm_sqr * pcm_ab_sqr /
                  (pcm_sqr_max * pcm_ab_sqr_max));
  } while (smash::random::canonical() > probability);
  const double pcm = std::sqrt(pcm_sqr), pcm_ab = std::sqrt(pcm_ab_sqr);
  smash::Angles phitheta;
  phitheta.distribute_isotropically();
  c.momentum = smash::FourVector(std::sqrt(m_c * m_c + pcm_sqr),
                                 pcm * phitheta.threevec());
  const smash::ThreeVector beta_cm =
      pcm * phitheta.threevec() / std::sqrt(pcm_sqr + mab * mab);

  phitheta.distribute_isotropically();
  a.momentum = smash::FourVector(std::sqrt(m_a * m_a + pcm_ab_sqr),
                                 pcm_ab * phitheta.threevec());
  b.momentum = smash::FourVector(std::sqrt(m_b * m_b + pcm_ab_sqr),
                                 -pcm_ab * phitheta.threevec());
  a.momentum = a.momentum.lorentz_boost(beta_cm);
  b.momentum = b.momentum.lorentz_boost(beta_cm);
  smash::FourVector mom_tot = a.momentum + b.momentum + c.momentum;
  mom_tot.set_x0(mom_tot.x0() - srts);
  assert(std::abs(mom_tot.sqr()) < 1.e-12);
}

double MicrocanonicalSampler::mu_minus_E_over_T(
    const SamplerParticle &p, const HyperSurfacePatch &hypersurface) {
  const HyperSurfacePatch::hydro_cell c = hypersurface.cell(p.cell_index);
  const double mu = c.muB * p.type->baryon_number() +
                    c.muS * p.type->strangeness() + c.muQ * p.type->charge();
  return (mu - p.momentum.Dot(c.u)) / c.T;
}

double MicrocanonicalSampler::compute_cells_factor(
    const SamplerParticleList &in, const SamplerParticleList &out,
    const HyperSurfacePatch &hypersurface) const {
  double cells_factor = 0.0;
  if (!quantum_statistics_) {
    for (const SamplerParticle &p : out) {
      cells_factor += mu_minus_E_over_T(p, hypersurface);
    }
    for (const SamplerParticle &p : in) {
      cells_factor -= mu_minus_E_over_T(p, hypersurface);
    }
    cells_factor = std::exp(cells_factor);
  } else {
    cells_factor = 1.0;
    for (const SamplerParticle &p : in) {
      const double x = std::exp(-mu_minus_E_over_T(p, hypersurface));
      cells_factor *= p.type->is_meson() ? (x - 1) : (x + 1);
    }
    for (const SamplerParticle &p : out) {
      const double x = std::exp(-mu_minus_E_over_T(p, hypersurface));
      cells_factor /= p.type->is_meson() ? (x - 1) : (x + 1);
    }
  }
  return cells_factor;
}

double
MicrocanonicalSampler::compute_spin_factor(const SamplerParticleList &in,
                                           const SamplerParticleList &out) {
  double spin_factor = 1.0;
  for (const SamplerParticle &part : out) {
    spin_factor *= part.type->pdgcode().spin_degeneracy();
  }
  for (const SamplerParticle &part : in) {
    spin_factor /= part.type->pdgcode().spin_degeneracy();
  }
  return spin_factor;
}

size_t MicrocanonicalSampler::N_available_channels3(std::array<int, 3> &BSQ,
                                                    double srts) {
  return (thresholds3_.find(BSQ) == thresholds3_.end())
             ? 0
             : std::upper_bound(thresholds3_[BSQ].begin(),
                                thresholds3_[BSQ].end(), srts) -
                   thresholds3_[BSQ].begin();
}

size_t MicrocanonicalSampler::N_available_channels2(std::array<int, 3> &BSQ,
                                                    double srts) {
  return (thresholds2_.find(BSQ) == thresholds2_.end())
             ? 0
             : std::upper_bound(thresholds2_[BSQ].begin(),
                                thresholds2_[BSQ].end(), srts) -
                   thresholds2_[BSQ].begin();
}

void MicrocanonicalSampler::random_two_to_three(
    const HyperSurfacePatch &hypersurface,
    SamplerParticleList &particles) {
  const size_t Npart = particles.size();
  const size_t Ncells = hypersurface.Ncells();
  // Choose 2 random particles
  size_t in_index1, in_index2;
  in_index1 = smash::random::uniform_int<size_t>(0, Npart - 1);
  do {
    in_index2 = smash::random::uniform_int<size_t>(0, Npart - 1);
  } while (in_index1 == in_index2);
  SamplerParticleList in{particles[in_index1], particles[in_index2]};
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
      std::cout << "Not enough energy (" << srts << " GeV) for 2->3 ("
                << in[0].type->name() << in[1].type->name()
                << "-> ?)" << std::endl;
    }
    #pragma omp atomic
    rejected23_count_++;
    return;
  }

  SamplerParticleList out;
  out.clear();
  out.resize(3);
  std::array<size_t, 3> cell = {
      smash::random::uniform_int<size_t>(0, Ncells - 1),
      smash::random::uniform_int<size_t>(0, Ncells - 1),
      smash::random::uniform_int<size_t>(0, Ncells - 1)};
  size_t i_channel = smash::random::uniform_int<size_t>(0, N3 - 1);
  std::array<ParticleTypePtr, 3> out_types = channels3_[BSQ][i_channel];
  const double R3 = compute_R3(
      srts, out_types[0]->mass(), out_types[1]->mass(), out_types[2]->mass());
  const double R2 = compute_R2(srts, in[0].type->mass(), in[1].type->mass());
  out[0] = {smash::FourVector(), out_types[0], cell[0], true};
  out[1] = {smash::FourVector(), out_types[1], cell[1], true};
  out[2] = {smash::FourVector(), out_types[2], cell[2], true};

  sample_3body_phase_space(srts, out[0], out[1], out[2]);
  for (SamplerParticle &part : out) {
    part.momentum = part.momentum.lorentz_boost(-beta_cm);
  }
  if (debug_printout_ == 3) {
    std::cout << "In: " << in[0].momentum << " " << in[1].momentum << std::endl;
    std::cout << "Out: " << out[0].momentum << " " << out[1].momentum << " "
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
  for (const SamplerParticle &p : out) {
    energy_factor *= p.momentum.Dot(hypersurface.cell(p.cell_index).dsigma);
  }
  for (const SamplerParticle &p : in) {
    energy_factor /= p.momentum.Dot(hypersurface.cell(p.cell_index).dsigma);
  }

  double cells_factor = compute_cells_factor(in, out, hypersurface);

  const double w = phase_space_factor * spin_factor * combinatorial_factor *
                   energy_factor * cells_factor;

  if (debug_printout_ == 1) {
    std::cout << "Trying 2->3 (" << in << "->" << out << ") with w = " << w
              << "; ph. sp. factor * 1e6 = " << phase_space_factor * 1.e6
              << "; spin_factor = " << spin_factor
              << "; comb. factor = " << combinatorial_factor
              << "; energy factor = " << energy_factor
              << "; cell factor = " << cells_factor << "; srts = " << srts
              << std::endl;
  }

  // No need to do it: w = std::min(1.0, w);
  if (smash::random::canonical() < w) {
    particles[in_index1] = out[0];
    particles[in_index2] = out[1];
    particles.push_back(out[2]);
    #pragma omp atomic
    accepted23_count_++;
    if (debug_printout_ == 1) {
      std::cout << in << "->" << out << ";" << srts << std::endl;
    }
  } else {
    #pragma omp atomic
    rejected23_count_++;
  }
}

void MicrocanonicalSampler::random_three_to_two(
    const HyperSurfacePatch &hypersurface,
    SamplerParticleList &particles) {
  const size_t Npart = particles.size();
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
  SamplerParticleList in{particles[in_index1], particles[in_index2],
                         particles[in_index3]};
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
      std::cout << "Not enough energy (" << srts
                << " GeV) or wrong quantum"
                   " numbers (B,S,Q = " << BSQ[0] << ", " << BSQ[1] << ", "
                << BSQ[2] << ") for 3->2"
                << std::endl;
    }
    #pragma omp atomic
    rejected32_count_++;
    return;
  }

  SamplerParticleList out;
  out.clear();
  out.resize(2);
  std::array<size_t, 2> cell = {
      smash::random::uniform_int<size_t>(0, Ncells - 1),
      smash::random::uniform_int<size_t>(0, Ncells - 1)};
  size_t i_channel = smash::random::uniform_int<size_t>(0, N2 - 1);
  std::array<ParticleTypePtr, 2> out_types = channels2_[BSQ][i_channel];
  const double R3 = compute_R3(
      srts, in[0].type->mass(), in[1].type->mass(), in[2].type->mass());
  const double R2 =
      compute_R2(srts, out_types[0]->mass(), out_types[1]->mass());
  out[0] = {smash::FourVector(), out_types[0], cell[0], true};
  out[1] = {smash::FourVector(), out_types[1], cell[1], true};

  const double m1 = out[0].type->mass(), m2 = out[1].type->mass();
  const double p_cm = smash::pCM(srts, m1, m2);
  smash::Angles phitheta;
  phitheta.distribute_isotropically();
  const smash::ThreeVector mom = p_cm * phitheta.threevec();
  const double mom2 = p_cm * p_cm;
  smash::FourVector p1(std::sqrt(mom2 + m1 * m1), mom),
      p2(std::sqrt(mom2 + m2 * m2), -mom);
  out[0].momentum = p1.lorentz_boost(-beta_cm);
  out[1].momentum = p2.lorentz_boost(-beta_cm);

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
  for (const SamplerParticle &p : out) {
    energy_factor *= p.momentum.Dot(hypersurface.cell(p.cell_index).dsigma);
  }
  for (const SamplerParticle &p : in) {
    energy_factor /= p.momentum.Dot(hypersurface.cell(p.cell_index).dsigma);
  }

  double cells_factor = compute_cells_factor(in, out, hypersurface);

  const double w = phase_space_factor * spin_factor * combinatorial_factor *
                   energy_factor * cells_factor;

  if (debug_printout_ == 1) {
    std::cout << "Trying 3->2 (" << in << "->" << out << ") with w = " << w
              << "; ph. sp. factor * 1e6 = " << phase_space_factor * 1.e6
              << "; spin_factor = " << spin_factor
              << "; comb. factor = " << combinatorial_factor
              << "; energy factor = " << energy_factor
              << "; cell factor = " << cells_factor << "; srts = " << srts
              << std::endl;
  }

  // No need to do it: w = std::min(1.0, w);
  if (smash::random::canonical() < w) {
    particles[in_index1] = out[0];
    particles[in_index2] = out[1];
    particles.erase(particles.begin() + in_index3);
    #pragma omp atomic
    accepted32_count_++;
    if (debug_printout_) {
      std::cout << in << "->" << out << ";" << srts << std::endl;
    }
  } else {
    #pragma omp atomic
    rejected32_count_++;
  }
}

void MicrocanonicalSampler::random_two_to_two(
    const HyperSurfacePatch &hypersurface,
    SamplerParticleList &particles,
    bool elastic) {
  const size_t Npart = particles.size();
  const size_t Ncells = hypersurface.Ncells();
  // Choose 2 random particles
  size_t in_index1, in_index2;
  in_index1 = smash::random::uniform_int<size_t>(0, Npart - 1);
  do {
    in_index2 = smash::random::uniform_int<size_t>(0, Npart - 1);
  } while (in_index1 == in_index2);
  SamplerParticleList in{particles[in_index1], particles[in_index2]};
  QuantumNumbers cons(in);
  const smash::FourVector momentum_total = cons.momentum;
  const double srts = momentum_total.abs();
  const smash::ThreeVector beta_cm = momentum_total.velocity();

  // Choose channel
  std::array<int, 3> BSQ{cons.B, cons.S, cons.Q};
  size_t N2 = N_available_channels2(BSQ, srts);

  SamplerParticleList out;
  out.clear();
  out.resize(2);
  const std::array<size_t, 2> cell = {
      smash::random::uniform_int<size_t>(0, Ncells - 1),
      smash::random::uniform_int<size_t>(0, Ncells - 1)};
  size_t i_channel = smash::random::uniform_int<size_t>(0, N2 - 1);
  std::array<ParticleTypePtr, 2> out_types{in[0].type, in[1].type};
  if (!elastic) {
    out_types = channels2_[BSQ][i_channel];
  }
  const double R2in = compute_R2(srts, in[0].type->mass(), in[1].type->mass());
  const double R2out = compute_R2(srts, out_types[0]->mass(), out_types[1]->mass());
  out[0] = {smash::FourVector(), out_types[0], cell[0], true};
  out[1] = {smash::FourVector(), out_types[1], cell[1], true};

  const double m1 = out[0].type->mass(), m2 = out[1].type->mass();
  const double p_cm = smash::pCM(srts, m1, m2);
  smash::Angles phitheta;
  phitheta.distribute_isotropically();
  const smash::ThreeVector mom = p_cm * phitheta.threevec();
  const double mom2 = p_cm * p_cm;
  out[0].momentum = smash::FourVector(std::sqrt(mom2 + m1 * m1),  mom);
  out[1].momentum = smash::FourVector(std::sqrt(mom2 + m2 * m2), -mom);
  for (SamplerParticle &part : out) {
    part.momentum = part.momentum.lorentz_boost(-beta_cm);
  }
  if (debug_printout_ == 3) {
    std::cout << "In: " << in[0].momentum << " " << in[1].momentum << std::endl;
    std::cout << "Out: " << out[0].momentum << " " << out[1].momentum << " "
              << std::endl;
  }
  // At this point the proposal function is set
  // Further is the calculation of probability to accept it
  const double phase_space_factor = R2out / R2in;
  const double spin_factor = compute_spin_factor(in, out);

  double combinatorial_factor = 1.0;
  if (out[0].type == out[1].type) {
    combinatorial_factor /= 2;
  }
  if (in[0].type == in[1].type) {
    combinatorial_factor *= 2;
  }

  double energy_factor = 1.0;
  for (const SamplerParticle &p : out) {
    energy_factor *= p.momentum.Dot(hypersurface.cell(p.cell_index).dsigma);
  }
  for (const SamplerParticle &p : in) {
    energy_factor /= p.momentum.Dot(hypersurface.cell(p.cell_index).dsigma);
  }

  double cells_factor = compute_cells_factor(in, out, hypersurface);

  if (elastic) {
    assert(phase_space_factor == 1.0);
    assert(spin_factor == 1.0);
    assert(combinatorial_factor == 1.0);
  }
  const double w = phase_space_factor * spin_factor * combinatorial_factor *
                   energy_factor * cells_factor;

  if (debug_printout_ == 1) {
    std::cout << "Trying 2->2 (" << in << "->" << out << ") with w = " << w
              << "; ph. sp. factor * 1e6 = " << phase_space_factor * 1.e6
              << "; spin_factor = " << spin_factor
              << "; comb. factor = " << combinatorial_factor
              << "; energy factor = " << energy_factor
              << "; cell factor = " << cells_factor << "; srts = " << srts
              << std::endl;
  }

  // No need to do it: w = std::min(1.0, w);
  if (smash::random::canonical() < w) {
    particles[in_index1] = out[0];
    particles[in_index2] = out[1];
    if (debug_printout_ == 1) {
      std::cout << in << "->" << out << ";" << srts << std::endl;
    }
  }
}


void MicrocanonicalSampler::renormalize_momenta(
    const smash::FourVector &required_total_momentum,
    SamplerParticleList &particles) {

  // Centralize momenta
  QuantumNumbers conserved = QuantumNumbers(particles);
  std::cout << particles.size() << " particles. " << std::endl;
  std::cout << "Required 4-momentum: " << required_total_momentum << std::endl;
  std::cout << "Sampled 4-momentum: " << conserved.momentum << std::endl;
  const smash::ThreeVector mom_to_add =
      (required_total_momentum.threevec() - conserved.momentum.threevec()) /
      particles.size();
  std::cout << "Adjusting momenta by " << mom_to_add << std::endl;
  for (auto &particle : particles) {
    const double m = particle.type->mass();
    const smash::ThreeVector mom3 = particle.momentum.threevec() + mom_to_add;
    const double Energy = std::sqrt(m * m + mom3.sqr());
    particle.momentum = smash::FourVector(Energy, mom3);
  }

  // Boost every particle to the common center of mass frame
  conserved = QuantumNumbers(particles);
  const smash::ThreeVector beta_CM_generated = conserved.momentum.velocity();
  const smash::ThreeVector beta_CM_required =
      required_total_momentum.velocity();

  double E = 0.0;
  double E_expected = required_total_momentum.abs();
  for (auto &particle : particles) {
    particle.momentum = particle.momentum.lorentz_boost(beta_CM_generated);
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
    for (const auto &particle : particles) {
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
  for (auto &particle : particles) {
    const double m = particle.type->mass();
    smash::ThreeVector mom3 = particle.momentum.threevec();
    mom3 *= (1 + a);
    const double Energy = std::sqrt(m * m + mom3.sqr());
    particle.momentum = smash::FourVector(Energy, mom3);
    particle.momentum = particle.momentum.lorentz_boost(-beta_CM_required);
  }
  QuantumNumbers conserved_final = QuantumNumbers(particles);
  std::cout << "Obtained total momentum: " << conserved_final.momentum
            << std::endl;
  if (conserved_final.momentum.abs() > required_total_momentum.abs() + 1.e-3) {
    throw std::runtime_error("Sum of particle masses is larger than "
      "the required patch energy.");
  }
}

void MicrocanonicalSampler::test_3body_phase_space_sampling() {
  const double sqrts = 6.0;
  ParticleType pi("pi", 0.138, 0.0, 0x111),
               N("N", 0.938, 0.0, 0x2212),
               Delta("D", 1.232, 0.115, 0x2224);
  SamplerParticle a{smash::FourVector(), &pi, 0, false},
                  b{smash::FourVector(), &N, 0, false},
                  c{smash::FourVector(), &Delta, 0, false};
  for (size_t i = 0; i < 100000; i++) {
    MicrocanonicalSampler::sample_3body_phase_space(sqrts, a, b, c);
    std::cout << a.momentum << b.momentum << c.momentum << std::endl;
  }
}

void MicrocanonicalSampler::print_rejection_stats() {
  std::cout << "2->3 acceptance: " << accepted23_count_ << "/"
            << accepted23_count_ + rejected23_count_ << std::endl;
  std::cout << "3->2 acceptance: " << accepted32_count_ << "/"
            << accepted32_count_ + rejected32_count_ << std::endl;
}

std::ostream &
operator<<(std::ostream &out,
           const MicrocanonicalSampler::SamplerParticleList &list) {
  const size_t N = list.size();
  out << list[0].type->name();
  for (size_t i = 1; i < N; i++) {
    out << "," << list[i].type->name();
  }
  return out;
}
