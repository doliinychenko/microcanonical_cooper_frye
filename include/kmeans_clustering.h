#ifndef KMEANS_CLUSTERING
#define KMEANS_CLUSTERING

#include <vector>

#include "hydro_cells.h"

/**
 * Adaptation of k_means clustering algorithm from
 * http://www.goldsborough.me/c++/python/cuda/
 * 2017/09/10/20-32-46-exploring_k-means_in_python,_c++_and_cuda/
 *
 * Clusters cells by collective velocity. Produces exactly <number_of_clusters>
 * clusters. Vector assignments has the same size as the vector of cells and
 * for each cell contains the index of cluster, to which it belongs.
 *
 * The function returns the vector with centers of clusters.
 *
 * The number of iterations is a detail of the algorithm, which tells, how many
 * adjustments of cluster centers should be made at maximum.
 */
std::vector<smash::ThreeVector>
k_means(const std::vector<HyperSurfacePatch::hydro_cell> &cells,
        size_t number_of_clusters, size_t number_of_iterations,
        std::vector<size_t> &assignments);

/// Used to test the k_means clustering function above
void test_clustering();
#endif // KMEANS_CLUSTERING
