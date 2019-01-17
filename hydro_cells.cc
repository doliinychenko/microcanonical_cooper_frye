#include <gsl/gsl_sf_bessel.h>

#include "hydro_cells.h"

#include "smash/constants.h"

#include <fstream>

HyperSurfacePatch::HyperSurfacePatch(
    const std::string &input_file,
    const std::function<bool(const smash::ParticleTypePtr)> &is_sampled,
    bool quantum_statistics) :
    quantum_statistics_(quantum_statistics),
    quantum_series_max_terms_(100),
    quantum_series_rel_precision_(1e-12) {
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
      const double T = cell.T;
      const double z = m / T;
      double mu_m_over_T = (mu - m) / T;
      if (mu_m_over_T > 0 and quantum_statistics_) {
        std::cout << "Warning: quantum expressions for " << t->name() <<
                     " do not converge, m < chemical potential." << std::endl;
      }
      // Compute parts of expressions, different for classical vs.
      // quantum statistics: x1 for density, x2 for energy, and x3
      // for momentum
      double x1 = 0.0, x2 = 0.0, x3 = 0.0;
      for (unsigned int k = 1; k < quantum_series_max_terms_; k++) {
        if (k > 1 and !quantum_statistics_) {
          break;
        }
        const double k1 = gsl_sf_bessel_Kn_scaled(1, z * k);
        const double k2 = gsl_sf_bessel_Kn_scaled(2, z * k);
        const double x = std::exp(mu_m_over_T * k);
        double x1_summand = z * z / k * k2 * x;
        double x2_summand = z * z / (k * k) * (3 * k2 + k * z * k1) * x;
        double x3_summand = x1_summand / k;
        // std::cout << "k = " << k
        //           << ", x1_summand*factor*1000 = " << x1_summand*factor*1000
        //           << ", x2_summand = " << x2_summand <<
        //           << ", x3_summand = " << x3_summand << std::endl;
        if (k > 1 and
            x1_summand < x1 * quantum_series_rel_precision_ and
            x2_summand < x2 * quantum_series_rel_precision_ and
            x3_summand < x3 * quantum_series_rel_precision_) {
          break;
        }
        if (k % 2 == 0 and t->is_baryon()) {
          x1_summand = -x1_summand;
          x2_summand = -x2_summand;
          x3_summand = -x3_summand;
        }
        x1 += x1_summand;
        x2 += x2_summand;
        x3 += x3_summand;
      }

      const double number_from_cell = T*T*T * x1 * dsigma.x0() *
                                      t->pdgcode().spin_degeneracy() * factor;
      B_tot_nonint_ += t->baryon_number() * number_from_cell;
      S_tot_nonint_ += t->strangeness() * number_from_cell;
      Q_tot_nonint_ += t->charge() * number_from_cell;
      smash::FourVector pmu_cell(dsigma.x0() * x2,
                                 -dsigma.x1() * x3,
                                 -dsigma.x2() * x3,
                                 -dsigma.x3() * x3);
      pmu_cell *= T*T*T*T * t->pdgcode().spin_degeneracy() * factor;
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
             << ", n_cells = " << patch.Ncells()
             << ", quantum statistics "
             << (patch.quantum_statistics() ? "ON" : "OFF");
  // clang-format on
}
