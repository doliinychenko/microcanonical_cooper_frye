#ifndef MICROCANONICAL_SAMPLING_MAIN_H
#define MICROCANONICAL_SAMPLING_MAIN_H

#include "smash/particletype.h"

#include <vector>

#include "microcanonical_sampler.h"
#include "hydro_cells.h"

#include "smash/particletype.h"
#include "smash/isoparticletype.h"
#include "smash/decaymodes.h"

/**
 *  Initializes random seed. It is a copy from SMASH 1.5 code. In SMASH 1.6
 *  this function is already present as a callable, so as soon as
 *  one switches to SMASH 1.6, this function can be substituted by calling
 *  SMASH function.
 */
void initialize_random_number_generator();

/**
 * Loads SMASH default particles and their decay modes. With SMASH 1.6
 * this function should be substituted by a SMASH callable.
 */
void load_smash_particles();

/**
 * Count, how many particles of type t are in the cell cell_number.
 */
int type_count(const MicrocanonicalSampler::SamplerParticleList &particles,
               const smash::ParticleTypePtr t, size_t cell_number);

/**
 * A function, which defines, which species will be sampled. For
 * a given species, if it returns true, the species will be sampled.
 */
bool is_sampled_type(const smash::ParticleTypePtr t);

/**
 * Main sampling routine. Reads the list of cells (hypersurface elements)
 * from hypersurface_input_file of format hypersurface_file_format.
 * Then prepares the sampler and samples.
 * N_warmup - number of warmup steps of Markov chain. These steps are
 *            done before printout. The number of steps should be large
 *            enough for Markov chain to reach stationary state.
 * N_decorrelate - number of Markov chain steps before next printout.
 * N_printout - number of printouts.
 * The total number of Markov chain steps is therefore
 * N_warmup + N_decorrelate * N_printout.
 */
void sample(std::string hypersurface_input_file,
            HyperSurfacePatch::InputFormat hypersurface_file_format,
            std::vector<smash::ParticleTypePtr>printout_types,
            int N_warmup, int N_decorrelate, int N_printout);
#endif  // MICROCANONICAL_SAMPLING_MAIN_H

