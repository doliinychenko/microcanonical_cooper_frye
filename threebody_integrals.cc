#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <streambuf>
#include <string>

#include "threebody_integrals.h"

#include <gsl/gsl_sf_ellint.h>

#include "smash/integrate.h"
#include "smash/kinematics.h"
#include "smash/pow.h"

ThreeBodyIntegrals::ThreeBodyIntegrals() {
  // lmsder, lmder, lmniel
  const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
  // allocate gsl fitter
  gsl_fitter_ =
      gsl_multifit_fdfsolver_alloc(T, n_tabulation_points, n_fit_param);
}

ThreeBodyIntegrals::~ThreeBodyIntegrals() {
  gsl_multifit_fdfsolver_free(gsl_fitter_);
}

int ThreeBodyIntegrals::gslfit_f(const gsl_vector *param, void *data,
                                 gsl_vector *f) {
  struct integral_parametrization p = {
      {gsl_vector_get(param, 0), gsl_vector_get(param, 1),
       gsl_vector_get(param, 2), gsl_vector_get(param, 3),
       gsl_vector_get(param, 4)},
      0.0};
  std::vector<std::pair<double, double>> *tabulated =
      reinterpret_cast<std::vector<std::pair<double, double>> *>(data);

  const size_t n = tabulated->size();
  for (size_t i = 0; i < n; i++) {
    const auto xy = tabulated->at(i);
    gsl_vector_set(f, i, (xy.second - integral_fit_function(xy.first, p)));
  }

  return GSL_SUCCESS;
}

int ThreeBodyIntegrals::gslfit_df(const gsl_vector *param, void *data,
                                  gsl_matrix *Jac) {
  struct integral_parametrization p = {
      {gsl_vector_get(param, 0), gsl_vector_get(param, 1),
       gsl_vector_get(param, 2), gsl_vector_get(param, 3),
       gsl_vector_get(param, 4)},
      0.0};
  struct integral_parametrization p_tmp = {
      {p.p[0], p.p[1], p.p[2], p.p[3], p.p[4]}, 0.0};
  const auto a = p.p;
  std::vector<std::pair<double, double>> *tabulated =
      reinterpret_cast<std::vector<std::pair<double, double>> *>(data);

  const size_t n = tabulated->size();

  for (size_t i = 0; i < n; i++) {
    const double x = tabulated->at(i).first;
    /* Jacobian matrix J(i,j) = df_i / dparam_j */
    const double t1 = 1.0 + a[1] * x, t2 = 1.0 + a[3] * x,
                 u1 = std::pow(t1, -a[2]), u2 = std::pow(t2, -a[4]);
    gsl_matrix_set(Jac, i, 0, u1 - u2);
    gsl_matrix_set(Jac, i, 1, -a[0] * a[2] * x * u1 / t1);
    gsl_matrix_set(Jac, i, 2, -a[0] * std::log(t1) * u1);
    gsl_matrix_set(Jac, i, 3, a[0] * a[4] * x * u2 / t2);
    gsl_matrix_set(Jac, i, 4, a[0] * std::log(t2) * u2);
    for (size_t j = 0; j < 5; ++j) {
      p_tmp.p[j] += 0.001;
      const double num_der =
          (integral_fit_function(x, p_tmp) - integral_fit_function(x, p)) /
          0.001;
      const double diff = (gsl_matrix_get(Jac, i, j) - num_der) / num_der;
      if (std::abs(diff) > 0.01) {
        std::cout << "An. vs num.: " << i << " " << gsl_matrix_get(Jac, i, j)
                  << " " << num_der << std::endl;
      }
      p_tmp.p[j] = p.p[j];
    }
  }
  return GSL_SUCCESS;
}
int ThreeBodyIntegrals::gslfit_fdf(const gsl_vector *param, void *data,
                                   gsl_vector *f, gsl_matrix *Jac) {
  gslfit_f(param, data, f);
  gslfit_df(param, data, Jac);

  return GSL_SUCCESS;
}
void ThreeBodyIntegrals::print_gsl_fitter_state(size_t iter,
                                                gsl_multifit_fdfsolver *s) {
  printf("iter: %3lu par = %15.8f %15.8f %15.8f %15.8f %15.8f"
         " |f(x)| = %g\n",
         iter, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1),
         gsl_vector_get(s->x, 2), gsl_vector_get(s->x, 3),
         gsl_vector_get(s->x, 4), gsl_blas_dnrm2(s->f));
}

void ThreeBodyIntegrals::add(double m1, double m2, double m3) {
  masses_3body masses = {m1, m2, m3};
  std::sort(masses.begin(), masses.end());
  auto it = saved_integrals_.find(masses);
  if (it != saved_integrals_.end()) {
    // std::cout << "Masses " << m1 << " " << m2 << " " << m3
    //          << " are already in the map." << std::endl;
    return; // no need to recalculate
  }

  std::vector<std::pair<double, double>> tabulated;
  integral_tabulate(m1, m2, m3, tabulated);

  gsl_multifit_function_fdf fdf;

  /* define the function to be minimized */
  fdf.f = &gslfit_f;
  fdf.df = NULL;  // &gslfit_df;  // set to NULL for finite-difference Jacobian
  fdf.fdf = NULL; //&gslfit_fdf; // Sometimes computing function AND Jacobian is
                  //faster
  fdf.n = n_tabulation_points;
  fdf.p = n_fit_param;
  fdf.params = &tabulated;

  double param_init[n_fit_param] = {1.0, 0.3, 0.6, 0.4, 1.3};
  gsl_vector_view gsl_param_init =
      gsl_vector_view_array(param_init, n_fit_param);
  gsl_multifit_fdfsolver_set(gsl_fitter_, &fdf, &gsl_param_init.vector);

  int status;
  size_t iter = 0;

  // print_gsl_fitter_state(iter, gsl_fitter_);

  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(gsl_fitter_);
    // printf("status = %s\n", gsl_strerror(status));
    // print_gsl_fitter_state(iter, gsl_fitter_);
    if (status) {
      break;
    }
    status =
        gsl_multifit_test_delta(gsl_fitter_->dx, gsl_fitter_->x, 1e-4, 1e-4);
  } while (status == GSL_CONTINUE && iter < 300);

  /*
  gsl_matrix *covar = gsl_matrix_alloc(n_fit_param, n_fit_param);
  gsl_matrix *Jac = gsl_matrix_alloc(n_tabulation_points, n_fit_param);
  gsl_multifit_fdfsolver_jac(gsl_fitter_, Jac);
  gsl_matrix_fprintf(stdout, Jac, "%15.8f");
  gsl_multifit_covar(Jac, 0.0, covar);
  gsl_matrix_free(covar);
  gsl_matrix_free(Jac);
  */

  // printf("status = %s\n", gsl_strerror(status));

  struct integral_parametrization a;
  for (size_t i = 0; i < n_fit_param; ++i) {
    a.p[i] = gsl_vector_get(gsl_fitter_->x, i);
  }

  double max_rel_diff = 0.0; //, x_of_max = -1.0;
  for (const auto x : tabulated) {
    const double parametrized = integral_fit_function(x.first, a);
    const double rel_diff = (x.second - parametrized) / parametrized;
    if (std::abs(rel_diff) > max_rel_diff) {
      max_rel_diff = std::abs(rel_diff);
      // x_of_max = x.first;
    }
    // std::cout << x.first << " " << x.second << " " << parametrized
    //          << " rel. diff. = " << (x.second - parametrized) / parametrized
    //          << std::endl;
  }
  // std::cout << "Max. rel. diff. = " << max_rel_diff << std::endl;
  // std::cout << "x of Max. rel. diff. = " << x_of_max << std::endl;
  a.max_rel_diff = max_rel_diff;
  saved_integrals_[masses] = a;
}

ThreeBodyIntegrals::integral_parametrization
ThreeBodyIntegrals::retrieve(double m1, double m2, double m3) {
  masses_3body masses = {m1, m2, m3};
  std::sort(masses.begin(), masses.end());
  auto it = saved_integrals_.find(masses);
  if (it != saved_integrals_.end()) {
    return it->second;
  } else {
    std::cout << "Integral for masses " << m1 << " " << m2 << " " << m3
              << " is missing in the table." << std::endl;
    return {{0.0, 0.0, 0.0, 0.0, 0.0}, 0.0};
  }
}

double ThreeBodyIntegrals::integral_fit_function(
    double x, const struct integral_parametrization &a) {
  return a.p[0] * (std::pow(1.0 + a.p[1] * x, -a.p[2]) -
                   std::pow(1.0 + a.p[3] * x, -a.p[4]));
}

void ThreeBodyIntegrals::integral_tabulate(
    double m1, double m2, double m3,
    std::vector<std::pair<double, double>> &tabulated) {
  tabulated.clear();
  const double m_sum = m1 + m2 + m3;
  const double a = 4.0 * M_PI * std::sqrt(m1 * m2 * m3) / std::pow(m_sum, 1.5);
  const double x_max = 5.0; // 0.1 * std::max(1.0, m_sum);
  const double dx = x_max / n_tabulation_points;
  for (unsigned int i = 0; i < n_tabulation_points; i++) {
    // Take tabulation points non-uniformely from 0 to infinity, so that
    // the whole region is covered and the largest density of points is near 0,
    // where best precision is needed. At small x I don't trust the integration
    // itself, because it results in really small numbers. While the
    // parametrization is made to reproduce the right asympotics at 0.
    const double x = 0.015 + dx * i;
    // const double x = std::tan(0.02 + (0.4 * M_PI) * i / n_tabulation_points);
    const double srts = x + m_sum;
    const auto result = integrate_(m1 + m2, srts - m3, [&](double m12) {
      return smash::pCM(srts, m12, m3) * smash::pCM(m12, m1, m2);
    });
    const double g =
        result.value() / (srts * (1.0 + (a - 1) / (1 + x)) * x * x) * 16.0 -
        1.0;
    // std::cout << "Tabulating: " << x << " " << g << std::endl;
    tabulated.push_back(std::make_pair(x, g));
  }
}

double ThreeBodyIntegrals::value(double srts, double m1, double m2, double m3) {
  const double m_sum = m1 + m2 + m3;
  if (m_sum > srts) {
    return 0.0;
  }
  const double a = 4.0 * M_PI * std::sqrt(m1 * m2 * m3) / std::pow(m_sum, 1.5);
  const double x = srts - m_sum;
  const double g = integral_fit_function(x, retrieve(m1, m2, m3));
  return (g + 1.0) * (1.0 + (a - 1) / (1 + x)) * x * x /
         (256.0 * M_PI * M_PI * M_PI);
}

void ThreeBodyIntegrals::save_to_file(std::string filename) {
  std::cout << "Saving tabulated 3-body integrals: " << filename << std::endl;
  std::ofstream output_file;
  output_file.open(filename);
  for (const auto &integral : saved_integrals_) {
    output_file << std::fixed << std::setprecision(16);
    ThreeBodyIntegrals::masses_3body m = integral.first;
    ThreeBodyIntegrals::integral_parametrization p = integral.second;
    output_file << m[0] << " " << m[1] << " " << m[2] << " ";
    for (size_t i = 0; i < n_fit_param; i++) {
      output_file << p.p[i] << " ";
    }
    output_file << p.max_rel_diff << std::endl;
  }
}

void ThreeBodyIntegrals::get_from_file(std::string filename) {
  std::cout << "Reading 3-body integrals from file " << filename << std::endl;
  std::ifstream input_file;
  input_file.open(filename, std::ios::in);
  std::string line;
  while (std::getline(input_file, line)) {
    std::istringstream iss(line);
    double m1, m2, m3, p0, p1, p2, p3, p4, perr;
    if (!(iss >> m1 >> m2 >> m3 >> p0 >> p1 >> p2 >> p3 >> p4 >> perr)) {
      std::cout << "Suspicious line: " << line << std::endl;
      break;
    }
    ThreeBodyIntegrals::masses_3body m{m1, m2, m3};
    std::sort(m.begin(), m.end());
    ThreeBodyIntegrals::integral_parametrization p{{p0, p1, p2, p3, p4}, perr};
    // If it is already in the map, prefer freshly computed integrals to
    // the saved ones.
    auto it = saved_integrals_.find(m);
    if (it == saved_integrals_.end()) {
      saved_integrals_[m] = p;
    }
  }
  std::cout << "Total integrals in the map:  " << saved_integrals_.size()
            << std::endl;
}

double ThreeBodyIntegrals::analytical_value(double srts, double m1, double m2,
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

void test_3body_integrals_precision() {
  ThreeBodyIntegrals test_integrals;
  test_integrals.add(1.2, 1.3, 1.4);
  test_integrals.add(0.2, 0.3, 0.4);
  std::cout << "Numerical integral values: "
            << test_integrals.value(5.0, 1.2, 1.3, 1.4) * 1E8 << " "
            << test_integrals.value(1.0, 0.2, 0.3, 0.4) * 1E8 << std::endl;
  std::cout << "Analytical integral values: "
            << ThreeBodyIntegrals::analytical_value(5.0, 1.2, 1.3, 1.4) * 1E8
            << " "
            << ThreeBodyIntegrals::analytical_value(1.0, 0.2, 0.3, 0.4) * 1E8
            << std::endl;
}

std::ostream &
operator<<(std::ostream &out,
           const ThreeBodyIntegrals::integral_parametrization &p) {
  out << "par: (" << p.p[0] << ", " << p.p[1] << ", " << p.p[2] << ", "
      << p.p[3] << ", " << p.p[4] << "), max. rel. diff = " << p.max_rel_diff;
  return out;
}
