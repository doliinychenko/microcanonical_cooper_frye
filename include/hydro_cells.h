#ifndef MICROCANONICAL_SAMPLER_HYDRO_CELLS_H
#define MICROCANONICAL_SAMPLER_HYDRO_CELLS_H

#include <functional>
#include <vector>

#include "smash/fourvector.h"
#include "smash/particletype.h"

class HyperSurfacePatch {

public:
  enum class InputFormat { DimaNaiveFormat = 0, MUSIC_ASCII_3plus1D = 1 };

  struct hydro_cell {
    smash::FourVector r;
    smash::FourVector dsigma;
    smash::FourVector u;
    double T;
    double muB;
    double muS;
    double muQ;
  };

  /**
   * Construct the hypersurface patch: read the list of elements from file,
   * compute total energy, momentum, and quantum numbers.
   * \param[in] input_file A file, where each file encodes a hypersurface
   *            element and has a format
   *            ds0 ds1 ds2 ds3 v1 v2 v3 T muB muS muQ, where
   *            (ds0 ds1 ds2 ds3) - normal 4-vector
   *            (v1 v2 v3) - collective velocity
   *            T - temperature in GeV
   *            muB muS muQ - baryon, strangeness, and charge chemical
   *                          potentials in GeV
   * \param[in] is_sampled a boolean function, true for particle types to be
   *            sampled. Total conserved quantities depend on it.
   * \param[in] quantum_statistics If true, then quantum statistics
   *            is taken into account for energy density and density
   *            calculation.
   */
  HyperSurfacePatch(
      const std::string &input_file, InputFormat read_in_format,
      const std::function<bool(const smash::ParticleTypePtr)> &is_sampled,
      bool quantum_statistics);
  /// Cannot be copied
  HyperSurfacePatch(const HyperSurfacePatch &) = delete;
  /// Cannot be copied
  HyperSurfacePatch &operator=(const HyperSurfacePatch &) = delete;

  int B() const { return B_tot_; }
  int S() const { return S_tot_; }
  int Q() const { return Q_tot_; }
  bool quantum_statistics() const { return quantum_statistics_; }
  smash::FourVector pmu() const { return pmu_tot_; }
  const std::vector<hydro_cell> &cells() const { return cells_; }
  const hydro_cell &cell(size_t index) const { return cells_[index]; }
  size_t Ncells() const { return cells_.size(); }

private:
  /// Read in cells from file of MUSIC format
  void read_from_MUSIC_file(const std::string &filename);

  /// Read in cells from file
  void read_from_file(const std::string &filename);

  /// Compute total 4-momentum, baryon number, strangeness, and charge
  void compute_totals();

  std::vector<hydro_cell> cells_;
  std::vector<smash::ParticleTypePtr> sampled_types_;
  smash::FourVector pmu_tot_;
  double B_tot_nonint_, S_tot_nonint_, Q_tot_nonint_;
  int B_tot_, S_tot_, Q_tot_;
  bool quantum_statistics_;
  /// Maximal number of terms in the series for quantum formulas
  const unsigned int quantum_series_max_terms_;
  /// Relative precision, at which quantum series summation stops
  const double quantum_series_rel_precision_;
};

/// Convenient printout of the HyperSurfacePatch class
std::ostream &operator<<(std::ostream &out, const HyperSurfacePatch &patch);

#endif // MICROCANONICAL_SAMPLER_HYDRO_CELLS_H
