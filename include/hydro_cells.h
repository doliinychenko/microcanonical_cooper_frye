#ifndef MICROCANONICAL_SAMPLER_HYDRO_CELLS_H
#define MICROCANONICAL_SAMPLER_HYDRO_CELLS_H

#include <vector>

#include "smash/fourvector.h"

class HyperSurfacePatch {

public:
 struct hydro_cell {
   smash::FourVector dsigma;
   smash::FourVector u;
   double T;
   double muB;
   double muS;
   double muQ;
 };

 HyperSurfacePatch(const std::string& input_file){
   cells_.clear();
   read_from_file(input_file);
   compute_totals();
 };
 /// Cannot be copied
 HyperSurfacePatch(const HyperSurfacePatch &) = delete;
 /// Cannot be copied
 HyperSurfacePatch &operator=(const HyperSurfacePatch &) = delete;


 int B() const { return B_tot_; }
 int S() const { return S_tot_; }
 int Q() const { return Q_tot_; }
 const std::vector<hydro_cell>& cells() const { return cells_; }
 int Ncells() const { return cells_.size(); }

 private:
  /// Read out cells from file
  void read_from_file(const std::string& filename);

  /// Compute total 4-momentum, baryon number, strangeness, and charge
  void compute_totals();

  std::vector<hydro_cell> cells_;
  smash::FourVector pmu_tot_;
  double B_tot_nonint_, S_tot_nonint_, Q_tot_nonint_;
  int B_tot_, S_tot_, Q_tot_;
};

#endif  // MICROCANONICAL_SAMPLER_HYDRO_CELLS_H
