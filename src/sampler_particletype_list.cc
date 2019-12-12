#include "microcanonical_sampler/sampler_particletype_list.h"
#include "smash/stringfunctions.h"

#include <fstream>

namespace {
  const std::vector<sampler::ParticleType> *all_particle_types;
}  // unnamed namespace

namespace sampler {

const std::vector<ParticleType> &ParticleType::list_all() {
  assert(all_particle_types);
  return *all_particle_types;
}

ParticleTypePtr ParticleType::operator&() const {
  const auto offset = this - std::addressof(list_all()[0]);
  return ParticleTypePtr(static_cast<uint16_t>(offset));
}

/**
 * Construct an antiparticle name-string from the given name-string for the
 * particle and its PDG code.
 *
 * \param[in] name the name-string of the particle to convert
 * \param[in] code the pdgcode of the particle to convert
 * \return the name-string of the converted antiparticle
 */
static std::string antiname(const std::string &name, smash::PdgCode code) {
  std::string basename, charge;

  if (name.find("⁺⁺") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁺⁺") + 1);
    charge = "⁻⁻";
  } else if (name.find("⁺") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁺") + 1);
    charge = "⁻";
  } else if (name.find("⁻⁻") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁻⁻") + 1);
    charge = "⁺⁺";
  } else if (name.find("⁻") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁻") + 1);
    charge = "⁺";
  } else if (name.find("⁰") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁰") + 1);
    charge = "⁰";
  } else {
    basename = name;
    charge = "";
  }

  // baryons & strange mesons: insert a bar
  if (code.baryon_number() != 0 || code.strangeness() != 0) {
    constexpr char bar[] = "\u0305";
    basename.insert(smash::utf8::sequence_length(basename.begin()), bar);
  }

  return basename + charge;
}

/**
 * Construct a charge string, given the charge as integer.
 *
 * \param[in] charge charge of a particle
 * \return the corresponding string to write out this charge
 * \throw runtime_error if the charge is not an integer between -2 and 2
 */
static std::string chargestr(int charge) {
  switch (charge) {
    case 2:
      return "⁺⁺";
    case 1:
      return "⁺";
    case 0:
      return "⁰";
    case -1:
      return "⁻";
    case -2:
      return "⁻⁻";
    default:
      throw std::runtime_error("Invalid charge " + std::to_string(charge));
  }
}


void read_particle_list(std::string inputfile_name,
                        ParticleListFormat format) {
  static std::vector<ParticleType> full_particletype_list;
  full_particletype_list.clear();
  std::ifstream infile(inputfile_name);
  if (!infile.good()) {
    throw std::runtime_error("Could not open file " + inputfile_name);
  }

  if (format == ParticleListFormat::SMASH) {
    while (true) {
      std::string line;
      if (!std::getline(infile, line)) {
        break;
      }
      const auto hash_pos = line.find('#');
      if (hash_pos != std::string::npos) {
        // Found a comment, remove it from the line and look further
        line = line.substr(0, hash_pos);
      }
      if (line.find_first_not_of(" \t") == std::string::npos) {
        // Only whitespace (or nothing) on this line. Next, please.
        continue;
      }
      line = smash::trim(line);
      std::istringstream lineinput(line);
      std::string name, parity_string;
      double mass, width;
      lineinput >> name >> mass >> width >> parity_string;
      const std::string failure_message = "While loading the ParticleType"
          " data:\nFailed to convert the input string \"" + line +
          "\" to the expected data types.";
      if (lineinput.fail()) {
        throw std::runtime_error(failure_message);
      }
      std::vector<std::string> pdgcode_strings;
      // We expect at most 4 PDG codes per multiplet.
      pdgcode_strings.reserve(4);
      // read additional PDG codes (if present)
      while (!lineinput.eof()) {
        pdgcode_strings.push_back("");
        lineinput >> pdgcode_strings.back();
        if (lineinput.fail()) {
          throw std::runtime_error(failure_message);
        }
      }
      if (pdgcode_strings.size() < 1) {
        throw std::runtime_error(failure_message);
      }
      std::vector<smash::PdgCode> pdgcode;
      pdgcode.resize(pdgcode_strings.size());
      std::transform(pdgcode_strings.begin(), pdgcode_strings.end(),
                     pdgcode.begin(),
                     [](const std::string &s) { return smash::PdgCode(s); });

      // add all states to type list
      for (size_t i = 0; i < pdgcode.size(); i++) {
        std::string full_name = name;
        if (pdgcode.size() > 1) {
          // for multiplets: add charge string to name
          full_name += chargestr(pdgcode[i].charge());
        }
        full_particletype_list.emplace_back(full_name, mass, width, pdgcode[i]);
        // std::cout << "Setting particle type: "
        //          << full_particletype_list.back().name() << std::endl;
        if (pdgcode[i].has_antiparticle()) {
          smash::PdgCode anti = pdgcode[i].get_antiparticle();
          full_name = antiname(full_name, pdgcode[i]);
          full_particletype_list.emplace_back(full_name, mass, width, anti);
          // std::cout << "Setting antiparticle type: "
          //          << full_particletype_list.back().name() << std::endl;
        }
      }
    }
  } else if (format == ParticleListFormat::iSS) {
    while (true) {
      std::string line;
      if (!std::getline(infile, line)) {
        break;
      }
      const auto hash_pos = line.find('#');
      if (hash_pos != std::string::npos) {
        // Found a comment, remove it from the line and look further
        line = line.substr(0, hash_pos);
      }
      line = smash::trim(line);
      if (line == "") {
        continue;
      }
      std::istringstream lineinput(line);
      std::string name;
      double mass, width;
      int bar, str, charm, bottom, charge, pdg;
      unsigned int spin_degeneracy, isospin_degeneracy, number_of_decays;
      lineinput >> pdg >> name >> mass >> width >> spin_degeneracy
                >> bar >> str >> charm >> bottom >> isospin_degeneracy
                >> charge >> number_of_decays;
      const std::string failure_message = "While loading the ParticleType"
          " data:\nFailed to convert the input string \"" + line +
          "\" to the expected data types.";
      if (lineinput.fail()) {
        throw std::runtime_error(failure_message);
      }
      /* Some baryons do not fit the PDG standard numbering scheme.
         This causes a mismatch between SMASH pdg treatment and iSS table.
         the following fixes this mismatch on a case-by-case basis. */
      if (pdg == 2110 || pdg == 2210 || pdg == 12110 || pdg == 12210) {
        pdg += 19920009;
      }
      smash::PdgCode pdgcode = smash::PdgCode::from_decimal(pdg);
      if (mass > 0) {
        // SMASH spin_degeneracy is bugged for photon, so testing only for m > 0
        assert(pdgcode.spin_degeneracy() == spin_degeneracy);
      }
      assert(pdgcode.baryon_number() == bar);
      assert(pdgcode.strangeness() == str);
      assert(pdgcode.charmness() == charm);
      assert(pdgcode.bottomness() == bottom);
      assert(pdgcode.charge() == charge);
      full_particletype_list.emplace_back(name, mass, width, pdgcode);
      // std::cout << full_particletype_list.back().name()
      //          << " " << mass << " " << pdgcode << std::endl;
      if (pdgcode.baryon_number() > 0) {
       full_particletype_list.emplace_back("anti-" + name,
                                      mass, width, pdgcode.get_antiparticle());
         // std::cout << full_particletype_list.back().name() << " "
         //          << mass << " " << pdgcode.get_antiparticle() << std::endl;
      }
      if (bar < 0) {
        throw std::runtime_error("Unexpected baryon number " +
                                 std::to_string(bar));
      }
      // Skip decays as they are not necessary for the sampler
      for (unsigned int j = 0; j < number_of_decays; j++) {
        std::getline(infile, line);
      }
    }
  } else {
    throw std::runtime_error("Unknown particle list file format");
  }
  all_particle_types = &full_particletype_list;
  std::cout << "Particle list initialized with "
            << all_particle_types->size()
            << " particle species (this includes charge states) from "
            << inputfile_name << std::endl;
}

}  // namespace sampler
