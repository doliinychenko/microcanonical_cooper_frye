#ifndef THREEBODY_INTEGRALS_H
#define THREEBODY_INTEGRALS_H

#include <map>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_vector.h>

#include "smash/forwarddeclarations.h"
#include "smash/integrate.h"
#include "smash/particles.h"
#include "smash/quantumnumbers.h"

/**
 * This class is used to compute relativistic 3-body phase-space integrals
 * given by eq. (10) of arXiv:1710.00665 [hep-ph].
 *
 * The integrals are tabulated numerically, asymptotics are divided out to
 * improve precision. The remaining function after dividing out is
 * fit by a specific parametrization as a function of sqrt{s} - (m1 + m2 +m3)
 * for every set of m1, m2, m3.
 *
 * There is a functionality to save and retrieve these parametrizations
 * from file. Finally, it turned out that there is an analytical expression,
 * so it is compared to tabulated integrals too.
 */
class ThreeBodyIntegrals {
public:
  static constexpr int n_tabulation_points = 100;
  static constexpr int n_fit_param = 5;
  ThreeBodyIntegrals();
  ~ThreeBodyIntegrals();

  typedef std::array<double, 3> masses_3body;
  struct integral_parametrization {
    std::array<double, n_fit_param> p;
    double max_rel_diff;
  };
  void add(double m1, double m2, double m3);
  integral_parametrization retrieve(double m1, double m2, double m3);
  double value(double x, double m1, double m2, double m3);
  void save_to_file(std::string filename);
  void get_from_file(std::string filename);
  size_t number_of_integrals() const { return saved_integrals_.size(); }
  /* I have found rather late, that the 3-body phase space integrals, that are
   * tabulated and parametrized in this class actually have analytical
   * expression, which uses elliptic integrals. The expression is given
   * in Eqs. (54-58) of [hep-ph/9406404]. I have checked that the
   * coincidence with the values computed here numerically coincide with
   * this analytical formula with relative precision better than 10^-3.
   * After the integral is parametrized, computing analytical expression
   * turns out to be around 2.5 times slower than accessing the
   * parametrization. */
  static double analytical_value(double srts, double m1, double m2, double m3);

private:
  static double integral_fit_function(double x,
                                      const integral_parametrization &a);
  void integral_tabulate(double m1, double m2, double m3,
                         std::vector<std::pair<double, double>> &tabulated);
  static void print_gsl_fitter_state(size_t iter, gsl_multifit_fdfsolver *s);
  static int gslfit_df(const gsl_vector *x, void *data, gsl_matrix *Jac);
  static int gslfit_f(const gsl_vector *x, void *data, gsl_vector *f);
  static int gslfit_fdf(const gsl_vector *param, void *data, gsl_vector *f,
                        gsl_matrix *Jac);
  /** Mapping sorted masses (m1, m2, m3) to the parametrized
   *  3-body phase space integral.
   */
  std::map<masses_3body, integral_parametrization> saved_integrals_;
  smash::Integrator integrate_;
  gsl_multifit_fdfsolver *gsl_fitter_;
};

void test_3body_integrals_precision();

std::ostream &operator<<(std::ostream &out,
                         const ThreeBodyIntegrals::integral_parametrization &p);

#endif // THREEBODY_INTEGRALS_H
