#ifndef STATISTICS_SUMMARY_H
#define STATISTICS_SUMMARY_H

#include "microcanonical_sampler.h"

class Statistics {
 public:
  Statistics();

  void add_event(const std::vector<MicrocanonicalSampler::SamplerParticleList>& particles);

  void printout(std::string output_filename);
 private:
  int n_events_;
  std::vector<int> B_midy_, S_midy_, Q_midy_, E_midy_;
};
#endif  // STATISTICS_SUMMARY_H
