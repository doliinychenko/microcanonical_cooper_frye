#ifndef SAMPLER_PARTICLETYPE_LIST_H
#define SAMPLER_PARTICLETYPE_LIST_H

#include <vector>

#include "smash/pdgcode.h"

namespace sampler {

class ParticleType;
class ParticleTypePtr;

class ParticleType {
 public:
  ParticleType(std::string n, double m, double w, smash::PdgCode pdg) :
    name_(n), mass_(m), width_(w), pdg_(pdg) {}
  const std::string &name() const { return name_; };
  double mass() const { return mass_; }
  double width() const { return width_; }
  bool is_stable() const { return width_ < 1.e-5; }
  smash::PdgCode pdgcode() const { return pdg_; }
  bool is_hadron() const { return pdg_.is_hadron(); }
  bool is_baryon() const { return pdg_.is_baryon(); }
  bool is_meson() const { return pdg_.is_meson(); }
  int baryon_number() const { return pdg_.baryon_number(); }
  int strangeness() const { return pdg_.strangeness(); }
  int charge() const { return pdg_.charge(); }
  int charmness() const { return pdg_.charmness(); }
  int spin_degeneracy() const { return pdg_.spin_degeneracy(); }

  static const std::vector<ParticleType> &list_all();
  ParticleTypePtr operator&() const;

 private:
  std::string name_;
  double mass_;
  double width_;
  smash::PdgCode pdg_;
};

class ParticleTypePtr {
 public:
  const ParticleType &operator*() const { return lookup(); }
  const ParticleType *operator->() const { return std::addressof(lookup()); }
  ParticleTypePtr() = default;

  bool operator==(const ParticleTypePtr &rhs) const {
    return index_ == rhs.index_;
  }
  bool operator!=(const ParticleTypePtr &rhs) const {
    return index_ != rhs.index_;
  }
  bool operator<(const ParticleTypePtr &rhs) const {
    return index_ < rhs.index_;
  }
  operator bool() const { return index_ != 0xffff; }
 private:
  friend ParticleTypePtr ParticleType::operator&() const;
  explicit ParticleTypePtr(std::uint16_t i) : index_(i) {}
  const ParticleType &lookup() const {
    assert(index_ != 0xffff);
    return ParticleType::list_all()[index_];
  }
  std::uint16_t index_ = 0xffff;
};


enum class ParticleListFormat {
  SMASH,
  iSS,
};

void read_particle_list(std::string inputfile_name,
                        ParticleListFormat format);

}  // namespace sampler

#endif // SAMPLER_PARTICLETYPE_LIST_H
