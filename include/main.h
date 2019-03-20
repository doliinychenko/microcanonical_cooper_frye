#ifndef MICROCANONICAL_SAMPLING_MAIN_H
#define MICROCANONICAL_SAMPLING_MAIN_H

#include <vector>

#include "hydro_cells.h"
#include "microcanonical_sampler.h"

#include "smash/decaymodes.h"
#include "smash/isoparticletype.h"
#include "smash/particletype.h"

/// Provides samples to reproduce results of arXiv:1902.09775
void reproduce_arxiv_1902_09775();

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
 * max_mass - maximal mass of hadrons to be sampled [GeV]. It influences
 * both sampling and computation of total energy and charged on
 * the hypersurface.
 * Epatch - max. energy in the patch
 */
void sample(std::string hypersurface_input_file,
            HyperSurfacePatch::InputFormat hypersurface_file_format,
            std::vector<smash::ParticleTypePtr> printout_types, int N_warmup,
            int N_decorrelate, int N_printout,
            double max_mass, double Epatch);
#endif // MICROCANONICAL_SAMPLING_MAIN_H
