#include <fstream>

#include "microcanonical_sampler/statistics_summary.h"

#include "smash/pow.h"

Statistics::Statistics() {
  n_events_ = 0;
  B_midy_.clear();
  S_midy_.clear();
  Q_midy_.clear();
  E_midy_.clear();
}

void Statistics::add_event(const std::vector<MicrocanonicalSampler::SamplerParticleList>& particles) {
  int B_midy = 0, S_midy = 0, Q_midy = 0;
  double E_midy = 0.0;
  for (const MicrocanonicalSampler::SamplerParticleList& plist : particles) {
    for (const MicrocanonicalSampler::SamplerParticle& particle : plist) {
      const smash::FourVector pmu = particle.momentum;
      const double y = 0.5 * std::log((pmu.x0() - pmu.x3()) / (pmu.x0() + pmu.x3()));
      if (std::abs(y) > 1.0) {
        continue;
      }
      B_midy += particle.type->baryon_number();
      S_midy += particle.type->strangeness();
      Q_midy += particle.type->charge();
      E_midy += particle.momentum.x0();
    }
  }
  B_midy_.push_back(B_midy);
  S_midy_.push_back(S_midy);
  Q_midy_.push_back(Q_midy);
  E_midy_.push_back(E_midy);
  n_events_++;
}

void Statistics::printout(const std::string output_filename) {
  // Compute central moments from 1 to 4 of B, S, Q
  double B_midy_mean = 0, S_midy_mean = 0, Q_midy_mean = 0, E_midy_mean = 0,
         B2c = 0, S2c = 0, Q2c = 0, E2c = 0,
         B3c = 0, S3c = 0, Q3c = 0, E3c = 0,
         B4c = 0, S4c = 0, Q4c = 0, E4c = 0;
  for (int i = 0; i < n_events_; i++) {
    B_midy_mean += B_midy_[i];
    S_midy_mean += S_midy_[i];
    Q_midy_mean += Q_midy_[i];
    E_midy_mean += E_midy_[i];
  }
  B_midy_mean /= n_events_;
  S_midy_mean /= n_events_;
  Q_midy_mean /= n_events_;
  E_midy_mean /= n_events_;

  for (int i = 0; i < n_events_; i++) {
    B2c += smash::pow_int(B_midy_[i] - B_midy_mean, 2);
    B3c += smash::pow_int(B_midy_[i] - B_midy_mean, 3);
    B4c += smash::pow_int(B_midy_[i] - B_midy_mean, 4);

    S2c += smash::pow_int(S_midy_[i] - S_midy_mean, 2);
    S3c += smash::pow_int(S_midy_[i] - S_midy_mean, 3);
    S4c += smash::pow_int(S_midy_[i] - S_midy_mean, 4);

    Q2c += smash::pow_int(Q_midy_[i] - Q_midy_mean, 2);
    Q3c += smash::pow_int(Q_midy_[i] - Q_midy_mean, 3);
    Q4c += smash::pow_int(Q_midy_[i] - Q_midy_mean, 4);

    E2c += smash::pow_int(E_midy_[i] - E_midy_mean, 2);
    E3c += smash::pow_int(E_midy_[i] - E_midy_mean, 3);
    E4c += smash::pow_int(E_midy_[i] - E_midy_mean, 4);
  }
  B2c /= n_events_;
  B3c /= n_events_; 
  B4c /= n_events_; 

  S2c /= n_events_; 
  S3c /= n_events_; 
  S4c /= n_events_; 

  Q2c /= n_events_; 
  Q3c /= n_events_; 
  Q4c /= n_events_; 

  E2c /= n_events_; 
  E3c /= n_events_; 
  E4c /= n_events_; 

  // Cumulants 1-3 coicide with central moments, so compute 4th cumulants
  const double B4k = B4c - 3.0 * B2c * B2c,
               S4k = S4c - 3.0 * S2c * S2c,
               Q4k = Q4c - 3.0 * Q2c * Q2c,
               E4k = E4c - 3.0 * E2c * E2c;
  std::ofstream stat_output_file(output_filename, std::ios::out);
  stat_output_file << "# E, net B, net S, net Q at |y| < 1 \n"
            << E_midy_mean << " "
            << E2c / E_midy_mean << " "
            << E3c / E2c << " "
            << E4k / E2c << " "
            << B_midy_mean << " "
            << B2c / B_midy_mean << " "
            << B3c / B2c << " "
            << B4k / B2c << " "
            << S_midy_mean << " "
            << S2c / S_midy_mean << " "
            << S3c / S2c << " "
            << S4k / S2c << " "
            << Q_midy_mean << " "
            << Q2c / Q_midy_mean << " "
            << Q3c / Q2c << " "
            << Q4k / Q2c << " "
            << std::endl; 
}
