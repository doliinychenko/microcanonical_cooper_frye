/*
 * Purpose:   Converts Steinheimer format to MUSIC format.
 * Compiling: g++ -Wall -Wextra -std=c++11 -O3 -o converter steinheimer_to_music_converter.cc
 * Running:   ./converter input_file output_file
 */

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cout << " Error: 2 arguments expected." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const std::string inputfile_name(argv[1]), outputfile_name(argv[2]);
  std::ifstream inputfile;
  inputfile.open(inputfile_name);
  std::ofstream outputfile;
  outputfile.open(outputfile_name);
  std::cout << "Converting " << inputfile_name << " to " << outputfile_name
            << std::endl;

  static constexpr double hbarc = 0.197327053;

  std::string line;
  size_t line_counter = 0;
  while (std::getline(inputfile, line)) {
    line_counter++;
    std::istringstream iss(line);
    double t, x, y, z, T, muB, muS, muQ, vx, vy, vz, ds0, ds1, ds2, ds3;
    if (!(iss >> x >> y >> z
              >> T >> muB >> muS
              >> vx >> vy >> vz
              >> ds0 >> ds1 >> ds2 >> ds3)) {
      break;
    }
    const double gamma = 1.0 / std::sqrt(1.0 - vx * vx - vy * vy - vz * vz);
    muQ = 0.0;
    t = 20.0;
    if (T <= 1E-4) {
      std::cout << "Suspiciously low temperature T = "
                << T << " GeV." << std::endl;
    }

    if (line_counter % 100000 == 0) {
      std::cout << "Cell " << line_counter << std::endl;
    }

    /**
     * In iSS Milne coordinates are used and u_mu with lower index is read in:
     * dsigma_dot_u = tau*(da0*gammaT + ux*da1 + uy*da2 + uz*da3/tau);
     */
    const double tau = std::sqrt(t*t - z*z);
    const double eta = std::atanh(z / t), cheta = std::cosh(eta),
                                          sheta = std::sinh(eta);
    const double ds0_Milne = (ds0 * cheta - ds3 * sheta) / tau;
    const double ds1_Milne = ds1 / tau, ds2_Milne = ds2 / tau;
    const double ds3_Milne = (ds0 * sheta - ds3 * cheta);
    const double u0_Milne = gamma * (cheta - vz * sheta);
    const double u1_Milne = gamma * vx;
    const double u2_Milne = gamma * vy;
    const double u3_Milne = gamma * (vz * cheta - sheta);
    const double dummy_e = 1.0;  // It's only for viscous corrections in iSS
    outputfile << tau << " " << x << " " << y << " " << eta << " "
               << ds0_Milne << " " << ds1_Milne << " "
               << ds2_Milne << " " << ds3_Milne << " "
               << u0_Milne << " " << u1_Milne << " "
               << u2_Milne << " " << u3_Milne << " "
               << dummy_e / hbarc << " " << T / hbarc << " "
               << muB / hbarc << " " << muS / hbarc << " "
               << muQ / hbarc << " 0 0 0 0 0 0 0 0 0 0 0" << std::endl;
    // Check that the umu*dsigmamu remains invariant after conversion
    const double umu_dsigmamu = gamma * (ds0 - ds1 * vx - ds2 * vy - ds3 * vz);
    const double umu_dsigmamu_music = tau * (ds0_Milne * u0_Milne +
                                             ds1_Milne * u1_Milne +
                                             ds2_Milne * u2_Milne +
                                             ds3_Milne * u3_Milne / tau);
    if (std::abs(umu_dsigmamu - umu_dsigmamu_music) > 1.e-12) {
      std::cout << "u^mu * dsigma_mu should be invariant: "
                << umu_dsigmamu << " == " << umu_dsigmamu_music << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
}
