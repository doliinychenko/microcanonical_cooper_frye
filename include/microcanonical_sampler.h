#ifndef MICROCANONICAL_SAMPLER_H
#define MICROCANONICAL_SAMPLER_H

#include <array>
#include <map>
#include <vector>

#include "smash/fourvector.h"
#include "smash/particletype.h"

#include "hydro_cells.h"
#include "threebody_integrals.h"

class MicrocanonicalSampler {
 public:
  struct SamplerParticle {
    smash::FourVector momentum;
    smash::ParticleTypePtr type;
    size_t cell_index;
  };
  typedef std::vector<SamplerParticle> SamplerParticleList;

  struct QuantumNumbers {
    smash::FourVector momentum;
    int B;
    int S;
    int Q;
    QuantumNumbers(const SamplerParticleList& plist) {
      momentum = smash::FourVector();
      B = 0;
      S = 0;
      Q = 0;
      for (const SamplerParticle& particle : plist) {
        momentum += particle.momentum;
        B += particle.type->baryon_number();
        S += particle.type->strangeness();
        Q += particle.type->charge();
      }
    }
  };

  MicrocanonicalSampler(const std::function<bool(const smash::ParticleTypePtr)>
                        &is_sampled,
                        int debug_printout,
                        bool quantum_statistics);
  /// Cannot be copied
  MicrocanonicalSampler(const MicrocanonicalSampler &) = delete;
  /// Cannot be copied
  MicrocanonicalSampler &operator=(const MicrocanonicalSampler &) = delete;

  void initialize(const HyperSurfacePatch& hypersurface);
  void one_markov_chain_step(const HyperSurfacePatch& hypersurface);

  /// Set the level of debug printout
  void set_debug_printout(int a) { debug_printout_ = a; }
  /// Set quantum statistics
  void set_quantum_statistics(bool qs) { quantum_statistics_ = qs; }
  /// Gives read-only access to particles
  const SamplerParticleList& particles() const { return particles_; }
  void print_rejection_stats();
 private:
  /**
   */
  double compute_R2(double srts, double m1, double m2);
  size_t N_available_channels3(std::array<int, 3>& BSQ, double srts);
  size_t N_available_channels2(std::array<int, 3>& BSQ, double srts);
  void sample_3body_phase_space(double srts,
                                SamplerParticle &a,
                                SamplerParticle &b,
                                SamplerParticle &c);
  double mu_minus_E_over_T(const SamplerParticle& p,
                           const HyperSurfacePatch& hypersurface);
  double compute_cells_factor(const SamplerParticleList& in,
                              const SamplerParticleList& out,
                              const HyperSurfacePatch& hypersurface);
  double compute_spin_factor(const SamplerParticleList& in,
                             const SamplerParticleList& out);

  void random_two_to_three(const HyperSurfacePatch& hypersurface);
  void random_three_to_two(const HyperSurfacePatch& hypersurface);
  void renormalize_momenta(const smash::FourVector& required_4mom);

  SamplerParticleList particles_;
  std::vector<smash::ParticleTypePtr> sampled_types_;
  /** All possible combinations of 3 particle species with given quantum
   *  numbers B, S, Q defined by the key std::array<int,3>. The key maps
   *  to a vector of sorted triplets. The vector is sorted by the sum
   *  of pole masses. This is to simplify the search of channels below
   *  given sqrt (total energy).
   */
  std::map<std::array<int,3>,
           std::vector<std::array<smash::ParticleTypePtr,3>>> channels3_;
  /// Same for 2-specie combinations
  std::map<std::array<int,3>,
           std::vector<std::array<smash::ParticleTypePtr,2>>> channels2_;
  /** Pre-computed sum of the masses for 3-body channels and t2-body channels
   *  in the same order that in channels3_ and channels2_
   */
  std::map<std::array<int,3>, std::vector<double>> thresholds3_, thresholds2_;
  ThreeBodyIntegrals three_body_int_;
  int debug_printout_;
  bool quantum_statistics_;
  int accepted23_count_ = 0, accepted32_count_ = 0,
      rejected23_count_ = 0, rejected32_count_ = 0;
};

/// Convenient printout of the SamplerParticleList
std::ostream &operator<<(std::ostream &out,
  const MicrocanonicalSampler::SamplerParticleList &list);

#endif  // MICROCANONICAL_SAMPLER_H
