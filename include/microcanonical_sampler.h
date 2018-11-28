#ifndef MICROCANONICAL_SAMPLER_H
#define MICROCANONICAL_SAMPLER_H

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
                        int debug_printout);
  /// Cannot be copied
  MicrocanonicalSampler(const MicrocanonicalSampler &) = delete;
  /// Cannot be copied
  MicrocanonicalSampler &operator=(const MicrocanonicalSampler &) = delete;

  void initialize(const HyperSurfacePatch& hypersurface);
  void one_markov_chain_step(const HyperSurfacePatch& hypersurface);

  /// Set the level of debug printout
  void set_debug_printout(int a) { debug_printout_ = a; }
  /// Gives read-only access to particles
  const SamplerParticleList& particles() const { return particles_; }
 private:
  /**
   */
  double compute_R2(double srts, double m1, double m2);
  /// Compute sum of 3-body integrals for given quantum numbers
  double compute_sum_R3(const QuantumNumbers &cons);
  /// Compute sum of 2-body integrals for given quantum numbers
  double compute_sum_R2(const QuantumNumbers &cons);
  void sample_3body_phase_space(double srts,
                                SamplerParticle &a,
                                SamplerParticle &b,
                                SamplerParticle &c);
  bool quantum_numbers_match(const smash::ParticleTypePtrList& a,
                             const QuantumNumbers& qn);
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
  ThreeBodyIntegrals three_body_int_;
  int debug_printout_;
};

/// Convenient printout of the SamplerParticleList
std::ostream &operator<<(std::ostream &out,
  const MicrocanonicalSampler::SamplerParticleList &list);

#endif  // MICROCANONICAL_SAMPLER_H
