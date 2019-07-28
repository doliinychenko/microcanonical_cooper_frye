#include <gsl/gsl_sf_bessel.h>

#include "hydro_cells.h"

#include "smash/constants.h"

#include <fstream>

HyperSurfacePatch::HyperSurfacePatch(
    const std::string &input_file, InputFormat read_in_format,
    const std::array<double, 3> &eta_for_2Dhydro,
    const std::function<bool(const smash::ParticleTypePtr)> &is_sampled,
    bool q_stat)
    : quantum_statistics_(q_stat) {
  for (const smash::ParticleType &ptype : smash::ParticleType::list_all()) {
    if (is_sampled(&ptype)) {
      sampled_types_.push_back(&ptype);
    }
  }

  cells_.clear();
  switch (read_in_format) {
  case InputFormat::MUSIC_ASCII_3plus1D:
    read_from_MUSIC_file(input_file);
    break;
  case InputFormat::Steinheimer:
    read_from_Steinheimer_file(input_file);
    break;
  case InputFormat::DimaNaiveFormat:
    read_from_file(input_file);
    break;
  case InputFormat::VISH_2files:
    // In this case input_file is the name of a folder
    read_from_VISH_2files(input_file, eta_for_2Dhydro);
    break;
  default:
    throw std::runtime_error("Unknown input file format");
  }
  compute_totals();
}

HyperSurfacePatch::HyperSurfacePatch(const HyperSurfacePatch &big_patch,
   std::vector<hydro_cell>::iterator patch_begin,
   std::vector<hydro_cell>::iterator patch_end) {
  quantum_statistics_ = big_patch.quantum_statistics();
  sampled_types_ = big_patch.sampled_types();
  cells_.clear();
  cells_.resize(std::distance(patch_begin, patch_end));
  std::copy(patch_begin, patch_end, cells_.begin());
  sum_up_totals_from_cells();
}

void HyperSurfacePatch::sum_up_totals_from_cells() {
  B_tot_nonint_ = 0.0;
  S_tot_nonint_ = 0.0;
  Q_tot_nonint_ = 0.0;
  pmu_tot_ = smash::FourVector();
  for (hydro_cell &a_cell : cells_) {
    pmu_tot_ += a_cell.pmu;
    B_tot_nonint_ += a_cell.B;
    S_tot_nonint_ += a_cell.S;
    Q_tot_nonint_ += a_cell.Q;
  }
  B_tot_ = static_cast<int>(std::round(B_tot_nonint_));
  S_tot_ = static_cast<int>(std::round(S_tot_nonint_));
  Q_tot_ = static_cast<int>(std::round(Q_tot_nonint_));
}


void HyperSurfacePatch::read_from_MUSIC_file(const std::string &filename) {
  std::ifstream infile(filename);
  if (!infile.good()) {
    throw std::runtime_error("Could not open file " + filename);
  }
  std::string line;
  size_t line_counter = 0;
  while (std::getline(infile, line)) {
    line_counter++;
    std::istringstream iss(line);
    double ds0, ds1, ds2, ds3, u0, u1, u2, u3, En, T, muB, muS, muQ, tau, x, y,
        eta;
    // clang-format off
    if (!(iss >> tau >> x >> y >> eta
              >> ds0 >> ds1 >> ds2 >> ds3
              >> u0 >> u1 >> u2 >> u3
              >> En >> T >> muB)) {
      break;
    }
    // clang-format on
    smash::FourVector u_Milne(u0, u1, u2, u3);
    assert(T >= 0.0);
    smash::FourVector u_Milne_test(u0, u1, u2, u3);
    if (std::abs(u_Milne_test.sqr() - 1.0) > 1.e-9) {
      std::cout << "Warning at reading from MUSIC output (line "
                << line_counter << "): "
                << "u_Milne (u_eta multiplied by tau) = " << u_Milne_test
                << ", u^2 == 1 is not fulfilled with error "
                << std::abs(u_Milne_test.sqr() - 1.0) << std::endl;
    }
    En *= smash::hbarc;
    T *= smash::hbarc;
    muB *= smash::hbarc;
    muS = 0.0;
    muQ = 0.0;
    const double umu_dsigmamu_music = tau * (ds0 * u0 + ds1 * u1 +
                                             ds2 * u2 + ds3 * u3 / tau);
    // Transforming from Milne to Cartesian
    const double ch_eta = std::cosh(eta);
    const double sh_eta = std::sinh(eta);
    const double t = tau * ch_eta, z = tau * sh_eta;
    smash::FourVector u(u0 * ch_eta + u3 * sh_eta, u1, u2,
                        u0 * sh_eta + u3 * ch_eta);
    if (std::abs(u.sqr() - 1.0) > 1.e-3) {
      std::cout << "u*u should be 1, u*u = " << u.sqr() << std::endl;
    }
    smash::FourVector ds(tau * ch_eta * ds0 - ds3 * sh_eta, - tau * ds1,
                         - tau * ds2, tau * sh_eta * ds0 - ds3 * ch_eta);
    if (ds.sqr() < 0) {
      std::cout << "dsigma^2 < 0, dsigma = " << ds
             << ", original dsigma = " << smash::FourVector(ds0, ds1, ds2, ds3)
             << ", T = " << T << ", muB = " << muB << ", cell "
             << line_counter << std::endl;
    }

    // Check that the umu*dsigmamu remains invariant after conversion
    const double umu_dsigmamu = ds.Dot(u);
    if (std::abs(umu_dsigmamu - umu_dsigmamu_music) > 1.e-4) {
      std::cout << "u^mu * dsigma_mu should be invariant: "
                << umu_dsigmamu << " == " << umu_dsigmamu_music << std::endl;
    }

    cells_.push_back({{t, x, y, z}, ds, u, {0.0, 0.0, 0.0, 0.0},
                      T, muB, muS, muQ, 0.0, 0.0, 0.0});
  }
  std::cout << cells_.size() << " cells read from the file " << filename
            << std::endl;
}

void HyperSurfacePatch::read_from_Steinheimer_file(const std::string &fname) {
  std::ifstream infile(fname);
  if (!infile.good()) {
    throw std::runtime_error("Could not open file " + fname);
  }
  std::cout << "Reading cells from file " << fname << std::endl;
  std::string line;
  size_t line_counter = 0;
  while (std::getline(infile, line)) {
    line_counter++;
    std::istringstream iss(line);
    double t, x, y, z, T, muB, muS, muQ, vx, vy, vz, ds0, ds1, ds2, ds3;
    // clang-format off
    if (!(iss >> x >> y >> z
              >> T >> muB >> muS
              >> vx >> vy >> vz
              >> ds0 >> ds1 >> ds2 >> ds3)) {
      break;
    }
    // clang-format on
    const double gamma = 1.0 / std::sqrt(1.0 - vx * vx - vy * vy - vz * vz);
    smash::FourVector u(1, vx, vy, vz), ds{ds0, ds1, ds2, ds3};
    u *= gamma;
    if (ds.sqr() < 0) {
      std::cout << "dsigma^2 < 0, dsigma = " << ds
                << ", T = " << T << ", muB = " << muB << ", cell "
                << line_counter << std::endl;
    }
    if (std::abs(u.abs() - 1.0) > 1.e-6) {
      std::cout << u << " norm should be 1!" << std::endl;
      throw std::runtime_error("Wrong u^mu.");
    }
    assert(T >= 0.0);
    muQ = 0.0;
    t = 0.0;
    cells_.push_back({{t, x, y, z}, ds, u, {0.0, 0.0, 0.0, 0.0},
                     T, muB, muS, muQ, 0.0, 0.0, 0.0});
    if (line_counter % 100000 == 0) {
      std::cout << "Cell " << line_counter << std::endl;
    }
  }
}

void HyperSurfacePatch::read_from_file(const std::string &filename) {
  cells_.clear();
  std::ifstream infile(filename);
  if (!infile.good()) {
    throw std::runtime_error("Could not open file " + filename);
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
    cells_.push_back({{0, 0, 0, 0}, {ds0, ds1, ds2, ds3}, u,
                      {0.0, 0.0, 0.0, 0.0}, T, muB, muS, muQ, 0.0, 0.0, 0.0});
  }
}

void HyperSurfacePatch::read_from_VISH_2files(const std::string &folder_name,
                                const std::array<double, 3> &eta_for_2Dhydro) {
  const std::string infile1_name = folder_name + "/surface.dat";
  const std::string infile2_name = folder_name + "/decdat2.dat";
  const std::string infile3_name = folder_name + "/surface_mu.dat";
  std::ifstream infile1(infile1_name),
                infile2(infile2_name),
                infile3(infile3_name);
  if (!infile1.good()) {
    throw std::runtime_error("Could not open file " + infile1_name);
  }
  if (!infile2.good()) {
    throw std::runtime_error("Could not open file " + infile2_name);
  }
  if (!infile3.good()) {
    throw std::runtime_error("Could not open file " + infile3_name);
  }
  std::cout << "Reading midrapidity cells from " << infile1_name << ", "
            << infile2_name << ", and " << infile3_name << std::endl;
  std::vector<hydro_cell> midrapidity_cells;
  midrapidity_cells.clear();
  std::string line1, line2, line3;
  int line_counter = 0;
  while (true) {
    bool read_line1 = static_cast<bool>(std::getline(infile1, line1));
    bool read_line2 = static_cast<bool>(std::getline(infile2, line2));
    bool read_line3 = static_cast<bool>(std::getline(infile3, line3));
    if ((read_line1 != read_line2) || (read_line1 != read_line3)) {
      throw std::runtime_error("Line " + std::to_string(line_counter) +
                   ": files " + infile1_name + ", " + infile2_name + " and "
                    + infile3_name + " have different number of lines.");
    }
    if (!read_line1 && !read_line2 && !read_line3) {
      break;
    }
    std::istringstream iss1(line1), iss2(line2), iss3(line3);
    double tau_0, tau_fo, x_fo, y_fo;
    iss1 >> tau_0 >> tau_fo >> x_fo >> y_fo;
    double tau, da0, da1, da2, vx, vy, vz, Edec, Bn, Tdec, muB, muS, muQ;
    iss2 >> tau >> da0 >> da1 >> da2 >> vx >> vy
         >> Edec >> Bn >> Tdec >> muB >> muS;
    muQ = 0.0;
    // Add additional chemical potentials from surface_mu.dat
    double mu_u, mu_d, mu_s;
    iss3 >> mu_u >> mu_d >> mu_s;
    muB += (mu_u + 2.0*mu_d);
    muQ += (mu_u - mu_d);
    muS += (mu_d - mu_s);
    vz = 0.0;
    // tau is present in both files, should be the same
    assert(std::abs(tau_fo - tau) < 1e-5);
    assert(Tdec > 0);
    const double gamma = 1.0 / std::sqrt(1.0 - vx * vx - vy * vy - vz * vz);
    smash::FourVector u(1.0, vx, vy, vz);
    u *= gamma;
    assert(std::abs(u.abs() - 1.0) < 1.e-15);
    // For eta = 0, there is no difference between Milne and Cartesian
    midrapidity_cells.push_back({{tau_fo, x_fo, y_fo, 0},
                   {da0 * tau, -da1 * tau, -da2 * tau, 0.0}, u,
                   {0.0, 0.0, 0.0, 0.0}, Tdec, muB, muS, muQ, 0.0, 0.0, 0.0});
    line_counter++;
  }
  std::cout << "Finished reading, read " << line_counter
            << " cells." << std::endl;
  const double eta_min = eta_for_2Dhydro[0], eta_max = eta_for_2Dhydro[1],
               deta = eta_for_2Dhydro[2];
  std::cout << "Assuming boost invariance, extending midrapidity cells to the "
            << "eta range (" << eta_min << ", "
            << eta_max << ") with d_eta = " << deta << std::endl;
  assert(eta_min < eta_max);
  assert(deta > 0.01);
  double eta = eta_min;
  cells_.clear();
  while (eta <= eta_max + 1.e-6) {
    std::cout << "eta = " << eta << std::endl;
    const double ch_eta = std::cosh(eta),
                 sh_eta = std::sinh(eta);
    for (const auto &cell : midrapidity_cells) {
      const smash::FourVector r = {cell.r.x0() * ch_eta,
                                   cell.r.x1(),
                                   cell.r.x2(),
                                   cell.r.x0() * sh_eta},
                              u = {cell.u.x0() * ch_eta,
                                   cell.u.x1(),
                                   cell.u.x2(),
                                   cell.u.x0() * sh_eta},
                           dsig = {cell.dsigma.x0() * ch_eta * deta,
                                   cell.dsigma.x1() * deta,
                                   cell.dsigma.x2() * deta,
                                   cell.dsigma.x0() * sh_eta * deta};
      assert(std::abs(u.Dot(dsig) - cell.u.Dot(cell.dsigma) * deta) < 1.e-9);
      cells_.push_back({r, dsig, u, {0.0, 0.0, 0.0, 0.0},
                       cell.T, cell.muB, cell.muS, cell.muQ,
                       0.0, 0.0, 0.0});
    }
    eta += deta;
  }
}


void HyperSurfacePatch::compute_totals() {
  std::cout << "Computing 4-momentum and charges in cells" << std::endl;
  const double hbarc = smash::hbarc;
  const double factor = 1.0 / (2.0 * M_PI * M_PI * hbarc * hbarc * hbarc);
  unsigned int cell_counter = 0;
  #pragma omp parallel for
  for (auto it = cells_.begin(); it < cells_.end(); it++) {
    hydro_cell& this_cell = *it;
    #pragma omp atomic
    ++cell_counter;
    if (cell_counter % 100000 == 0) {
      printf("Cell %d\n", cell_counter);
    }
    this_cell.pmu = smash::FourVector();
    this_cell.B = 0.0;
    this_cell.S = 0.0;
    this_cell.Q = 0.0;
    const double T = this_cell.T;
    for (const smash::ParticleTypePtr t : sampled_types_) {
      const double m = t->mass();
      const double mu = this_cell.muB * t->baryon_number() +
                        this_cell.muS * t->strangeness() +
                        this_cell.muQ * t->charge();
      // dsigma in the frame, where umu = (1, 0, 0, 0)
      // dsigma[0] in this frame should be equal to dsigma_mu u^mu in any frame
      smash::FourVector dsigma = this_cell.dsigma.LorentzBoost(
                                                       this_cell.u.velocity());
      const double z = m / T;
      double mu_m_over_T = (mu - m) / T;
      if (mu_m_over_T > 0 and quantum_statistics_) {
        std::cout << "Warning: quantum expressions for " << t->name()
                  << " do not converge, m < chemical potential." << std::endl;
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
        if (k > 1 and x1_summand < x1 * quantum_series_rel_precision_ and
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

      const double number_from_cell = x1 * dsigma.x0() *
                                      t->pdgcode().spin_degeneracy();
      this_cell.B += t->baryon_number() * number_from_cell;
      this_cell.S += t->strangeness() * number_from_cell;
      this_cell.Q += t->charge() * number_from_cell;
      smash::FourVector pmu_cell(dsigma.x0() * x2, -dsigma.x1() * x3,
                                 -dsigma.x2() * x3, -dsigma.x3() * x3);
      pmu_cell *= t->pdgcode().spin_degeneracy();
      pmu_cell = pmu_cell.LorentzBoost(-this_cell.u.velocity());
      this_cell.pmu += pmu_cell;
      // std::cout << t.name() << " number: " << number_from_cell << std::endl;
      // std::cout << t.name() << ": " << pmu_cell << std::endl;
    }
    const double a = T * T * T * factor;
    this_cell.pmu *= T * a;
    this_cell.B *= a;
    this_cell.S *= a;
    this_cell.Q *= a;
  }
  this->sum_up_totals_from_cells();
}

std::vector<HyperSurfacePatch> HyperSurfacePatch::split(double E_patch_max) {
  // Compute variances of temperature and muB on the hypersurface
  // They are needed for clustering metric
  double mean_T = 0.0, mean_muB = 0.0,
         mean_T_sqr = 0.0, mean_muB_sqr = 0.0;
  for (const HyperSurfacePatch::hydro_cell &this_cell : cells_) {
    mean_T += this_cell.T;
    mean_muB += this_cell.muB;
    mean_T_sqr += this_cell.T * this_cell.T;
    mean_muB_sqr += this_cell.muB * this_cell.muB;
  }
  mean_T /= Ncells();
  mean_muB /= Ncells();
  mean_T_sqr /= Ncells();
  mean_muB_sqr /= Ncells();
  double mean_dT_sqr = mean_T_sqr - mean_T * mean_T,
         mean_dmuB_sqr = mean_muB_sqr - mean_muB * mean_muB;
  if (mean_dT_sqr < 1.e-3) {
    mean_dT_sqr = 0.0;
  }
  if (mean_dmuB_sqr < 1.e-3) {
    mean_dmuB_sqr = 0.0;
  }
  double inv_dT_sqr = (mean_dT_sqr < 1.e-3) ? 0.0 : 1.0 / mean_dT_sqr,
         inv_dmuB_sqr = (mean_dmuB_sqr < 1.e-3) ? 0.0 : 1.0 / mean_dmuB_sqr;
  std::cout << "Hypersurface temperature: " << mean_T << "±"
            << std::sqrt(mean_dT_sqr) << "  GeV." << std::endl;;
  std::cout << "Hypersurface baryochemical potential: " << mean_muB << "±"
            << std::sqrt(mean_dmuB_sqr) << "  GeV." << std::endl;
  /**
   * My home-made patch splitting algorithm:
   * 1) Start from the non-clustered cell with largest energy
   * 2) Sort cells by distance d^2 =
   *       (dx^2 + dy^2 + dz^2)/d0^2 + dT^2/sigT^2 + dmuB^2/sigmuB^2),
   *       where d0 [fm] is a heuristically defined constant.
   * 3) Pick up next cells to the patch until total energy is more than E_patch
   *    or no cells left. Then start over from 1).
   *
   * Main advantage of this algorithm is that one gets patches of similar total
   * energy, while the variation of temperature and muB is not too large within
   * a patch.
   */
  std::vector<HyperSurfacePatch> patches;
  std::vector<hydro_cell>::iterator nonclustered_begin = cells_.begin();
  size_t patch_counter = 0;
  while (nonclustered_begin != cells_.end()) {
    // (1) Find cell with maximum energy from non-clustered cells
    hydro_cell* max_energy_cell_ref = &(*nonclustered_begin);
    double maximum_energy = 0.0;
    for (auto it = nonclustered_begin; it < cells_.end(); it++) {
      if (it->pmu.x0() > maximum_energy) {
        maximum_energy = it->pmu.x0();
        max_energy_cell_ref = &(*it);
      }
    }
    // Now copy the cell, because soon the cells_ vector will be sorted
    // and references will mix up.
    hydro_cell max_energy_cell = *max_energy_cell_ref;
    std::cout << "Max energy at cell with T [GeV] = " << max_energy_cell.T
              << ", muB [GeV] = " << max_energy_cell.muB
              << ", r = " << max_energy_cell.r << std::endl;

    // (2) Sort non-clustered cells from smallest to largest d^2
    constexpr double d0 = 2.0;  // fm
    constexpr double inv_d0_sqr = 1.0 / (d0 * d0);
    std::sort(nonclustered_begin, cells_.end(),
              [&](const hydro_cell& a, const hydro_cell& b) {
      double da2 = (a.r.threevec() - max_energy_cell.r.threevec()).sqr();
      da2 *= inv_d0_sqr;
      double tmp = (a.T - max_energy_cell.T);
      da2 += tmp * tmp * inv_dT_sqr;
      tmp = (a.muB - max_energy_cell.muB);
      da2 += tmp * tmp * inv_dmuB_sqr;

      double db2 = (b.r.threevec() - max_energy_cell.r.threevec()).sqr();
      db2 *= inv_d0_sqr;
      tmp = (b.T - max_energy_cell.T);
      db2 += tmp * tmp * inv_dT_sqr;
      tmp = (b.muB - max_energy_cell.muB);
      db2 += tmp * tmp * inv_dmuB_sqr;

      return da2 < db2;
    });

    // (3) Collect cells to patch until energy is enough or no cells left
    smash::FourVector pmu_patch = smash::FourVector();
    std::vector<hydro_cell>::iterator cluster_start = nonclustered_begin;
    while (pmu_patch.sqr() < E_patch_max * E_patch_max &&
           nonclustered_begin != cells_.end()) {
      pmu_patch += nonclustered_begin->pmu;
      std::advance(nonclustered_begin, 1);
    }
/*
    double dmax = 0.0;
    for (auto it = cluster_start; it < nonclustered_begin; it++) {
      const double d = (it->r.threevec() - max_energy_cell.r.threevec()).abs();
      if (d > dmax) {
        dmax = d;
      }
    }
    // Emergency printout
    if (patch_counter == 0) {
      for (auto it = cluster_start; it < nonclustered_begin; it++) {
        const auto vel = it->u.velocity();
        std::cout << it->r.x1() << " " << it->r.x2() << " " << it->r.x3() << " "
              << it->T << " " << it->muB << " " << it->muS << " "
              << vel.x1() << " " << vel.x2() << " " << vel.x3() << " "
              << it->dsigma.x0() << " " << it->dsigma.x1() << " "
              << it->dsigma.x2() << " " << it->dsigma.x3() << std::endl;
      }
    }
    std::cout << "Max geometrical distance to cluster center [fm]: " << dmax << std::endl;
*/
    patches.push_back(HyperSurfacePatch(*this, cluster_start, nonclustered_begin));
    patch_counter++;
    std::cout << "Patch " << patch_counter << ". " << patches.back() << std::endl;
  }
  return patches;
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
