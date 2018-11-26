#include <gsl/gsl_sf_bessel.h>

#include "hydro_cells.h"

#include "smash/constants.h"

#include <fstream>

HyperSurfacePatch::HyperSurfacePatch(
    const std::string &input_file,
    const std::function<bool(const smash::ParticleTypePtr)> &is_sampled) {
  for (const smash::ParticleType &ptype : smash::ParticleType::list_all()) {
    if (is_sampled(&ptype)) {
      sampled_types_.push_back(&ptype);
    }
  }

  cells_.clear();
  read_from_file(input_file);
  compute_totals();
}

void HyperSurfacePatch::read_from_file(const std::string &filename) {
  std::ifstream infile(filename);
  if (!infile.good()) {
    throw std::runtime_error("Could not open file");
  }
  std::string line;
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    double ds0, ds1, ds2, ds3, v1, v2, v3, T, muB, muS, muQ;
    // clang-format off
    if (!(iss >> ds0 >> ds1 >> ds2 >> ds3
              >> v1 >> v2 >> v3
              >> T >> muB >> muS >> muQ)) {
      break;
    }
    // clang-format on
    assert(T >= 0.0);
    const double gamma = 1.0 / std::sqrt(1.0 - v1 * v1 - v2 * v2 - v3 * v3);
    smash::FourVector u(1.0, v1, v2, v3);
    u *= gamma;
    assert(std::abs(u.abs() - 1.0) < 1.e-15);
    cells_.push_back({{ds0, ds1, ds2, ds3}, u, T, muB, muS, muQ});
  }
}

void HyperSurfacePatch::compute_totals() {
  B_tot_nonint_ = 0.0;
  S_tot_nonint_ = 0.0;
  Q_tot_nonint_ = 0.0;
  pmu_tot_ = smash::FourVector();

  const double hbarc = smash::hbarc;
  const double factor = 1.0 / (2.0 * M_PI * M_PI * hbarc * hbarc * hbarc);
  for (const smash::ParticleTypePtr t : sampled_types_) {
    const double m = t->mass();
    for (const hydro_cell &cell : cells_) {
      const double mu = cell.muB * t->baryon_number() +
                        cell.muS * t->strangeness() +
                        cell.muQ * t->charge();
      // dsigma in the frame, where umu = (1, 0, 0, 0)
      // dsigma[0] in this frame should be equal to dsigma_mu u^mu in any frame
      smash::FourVector dsigma = cell.dsigma.LorentzBoost(cell.u.velocity());
      double x = (mu - m) / cell.T;
      x = std::exp(x);
      const double m_over_T = m / cell.T;
      const double k2 = gsl_sf_bessel_Kn_scaled(2, m / cell.T);
      const double k1 = gsl_sf_bessel_Kn_scaled(1, m / cell.T);
      const double number_from_cell = cell.T * m * m * x * dsigma.x0() * k2 *
                                      t->pdgcode().spin_degeneracy() * factor;
      B_tot_nonint_ += t->baryon_number() * number_from_cell;
      S_tot_nonint_ += t->strangeness() * number_from_cell;
      Q_tot_nonint_ += t->charge() * number_from_cell;
      smash::FourVector pmu_cell(dsigma.x0() * (3.0 * k2 + m_over_T * k1),
                                 -dsigma.x1() * k2,
                                 -dsigma.x2() * k2,
                                 -dsigma.x3() * k2);
      pmu_cell *= (cell.T * cell.T * m * m * x *
                   t->pdgcode().spin_degeneracy() * factor);
      pmu_cell = pmu_cell.LorentzBoost(-cell.u.velocity());
      pmu_tot_ += pmu_cell;
      // std::cout << t.name() << " number: " << number_from_cell << std::endl;
      // std::cout << t.name() << ": " << pmu_cell << std::endl;
    }
  }
  B_tot_ = static_cast<int>(std::round(B_tot_nonint_));
  S_tot_ = static_cast<int>(std::round(S_tot_nonint_));
  Q_tot_ = static_cast<int>(std::round(Q_tot_nonint_));
}

std::ostream &operator<<(std::ostream &out, const HyperSurfacePatch &patch) {
  // clang-format off
  return out << "Patch: p^mu = " << patch.pmu()
             << ", B = " << patch.B()
             << ", S = " << patch.S()
             << ", Q = " << patch.Q()
             << ", n_cells = " << patch.Ncells();
  // clang-format on
}
