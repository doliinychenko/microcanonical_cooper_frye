#include "hydro_cells.h"

#include "smash/hadgas_eos.h"

#include <fstream>

void HyperSurfacePatch::read_from_file(const std::string& filename) {
  std::ifstream infile(filename);
  if (!infile.good()) {
    throw std::runtime_error("Could not open file");
  }
  std::string line;
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    double ds0, ds1, ds2, ds3, u0, u1, u2, u3, T, muB, muS, muQ;
    if (!(iss >> ds0 >> ds1 >> ds2 >> ds3
              >> u0 >> u1 >> u2 >> u3
              >> T >> muB >> muS >> muQ)) {
      break;
    }
    assert(T >= 0.0);
    smash::FourVector u(u0, u1, u2, u3);
    assert(std::abs(u.abs() - 1.0) < 1.e-15);
    cells_.push_back({{ds0, ds1, ds2, ds3}, u, T, muB, muS, muQ});
  }
}

void HyperSurfacePatch::compute_totals() {
  return;
}
