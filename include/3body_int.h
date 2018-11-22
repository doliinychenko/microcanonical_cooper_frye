#ifndef THREE_BODY_INT_H
#define THREE_BODY_INT_H

#include <map>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "smash/integrate.h"
#include "smash/particles.h"
#include "smash/forwarddeclarations.h"
#include "smash/quantumnumbers.h"

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

 private:
  static double integral_fit_function(double x, const integral_parametrization &a);
  void integral_tabulate(double m1, double m2, double m3,
    std::vector<std::pair<double, double>> &tabulated);
  static void print_gsl_fitter_state(size_t iter, gsl_multifit_fdfsolver *s);
  static int gslfit_df(const gsl_vector *x, void *data, gsl_matrix *Jac);
  static int gslfit_f(const gsl_vector *x, void *data, gsl_vector *f);
  static int gslfit_fdf(const gsl_vector *param, void *data,
                        gsl_vector *f, gsl_matrix *Jac);
  /** Mapping sorted masses (m1, m2, m3) to the parametrized
   *  3-body phase space integral.
   */
  std::map<masses_3body, integral_parametrization> saved_integrals_;
  smash::Integrator integrate_;
  gsl_multifit_fdfsolver *gsl_fitter_;
};

std::ostream &operator<<(std::ostream &out,
    const ThreeBodyIntegrals::integral_parametrization &p);

double compute_R2(double srts, double m1, double m2);
double compute_sum_R3(const smash::ParticleTypePtrList& sampled_types,
                      ThreeBodyIntegrals& three_body_int,
                      const smash::QuantumNumbers &cons);
double compute_sum_R2(const smash::ParticleTypePtrList& sampled_types,
                      const smash::QuantumNumbers &cons);
int type_count(const smash::ParticleList &particles,
               const smash::ParticleTypePtr t);
void sample_3body_phase_space(double srts,
                              smash::ParticleData &a,
                              smash::ParticleData &b,
                              smash::ParticleData &c);
bool quantum_numbers_match(const smash::ParticleTypePtrList& a,
                           const smash::QuantumNumbers& qn);
void random_two_to_three(smash::ParticleList &particles,
                         const smash::ParticleTypePtrList& sampled_types,
                         ThreeBodyIntegrals& three_body_int,
                         double V);
void random_three_to_two(smash::ParticleList &particles,
                         const smash::ParticleTypePtrList& sampled_types,
                         ThreeBodyIntegrals& three_body_int,
                         double V);

void renormalize_momenta(smash::ParticleList& plist,
                         smash::FourVector required_4mom);
void initialize_random_number_generator();
void load_smash_particles();

#endif  // THREE_BODY_INT_H

