#ifndef MICROCANONICAL_SAMPLING_MAIN_H
#define MICROCANONICAL_SAMPLING_MAIN_H

#include <vector>

#include "hydro_cells.h"
#include "microcanonical_sampler.h"

/// Generates proper random seed
int64_t generate_63bit_seed();

/// Provides samples to reproduce results of arXiv:1902.09775
void reproduce_arxiv_1902_09775();

/**
 * Main sampling routine. Reads the list of cells (hypersurface elements)
 * from hypersurface_input_file of format hypersurface_file_format.
 * Then prepares the sampler and samples.
 * The sampled particles are written to output_file.
 * eta_for_2Dhydro - eta_min, eta_max, d_eta for 2+1D hydro
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
void sample(const std::string hypersurface_input_file,
            HyperSurfacePatch::InputFormat hypersurface_file_format,
            const std::array<double, 3> &eta_for_2Dhydro,
            const std::string output_file,
            const std::string patches_output_filename, size_t N_warmup,
            size_t N_decorrelate, size_t N_printout,
            double max_mass, double Epatch);

/// Prints out information about all command-line options
void usage(const int rc, const std::string &progname);
#endif // MICROCANONICAL_SAMPLING_MAIN_H
