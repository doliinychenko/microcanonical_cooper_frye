#ifndef KMEANS_CLUSTERING
#define KMEANS_CLUSTERING

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <random>
#include <vector>

#include "smash/random.h"
#include "smash/threevector.h"

#include "hydro_cells.h"

// Adaptation from http://www.goldsborough.me/c++/python/cuda/2017/09/10/20-32-46-exploring_k-means_in_python,_c++_and_cuda/
std::vector<smash::ThreeVector> k_means(
                  const std::vector<HyperSurfacePatch::hydro_cell>& cells,
                  size_t number_of_clusters,
                  size_t number_of_iterations,
                  std::vector<size_t>& assignments) {

  assert(assignments.size() == cells.size());
  // Pick centroids as random points from the dataset.
  std::vector<smash::ThreeVector> means(number_of_clusters);
  for (auto& cluster : means) {
    size_t i_cell = smash::random::uniform_int<size_t>(0, cells.size() - 1);
    cluster = cells[i_cell].u.velocity();
  }

  for (size_t iteration = 0; iteration < number_of_iterations; ++iteration) {
    // Find assignments.
    for (size_t cell = 0; cell < cells.size(); ++cell) {
      double best_distance = std::numeric_limits<double>::max();
      size_t best_cluster = 0;
      for (size_t cluster = 0; cluster < number_of_clusters; ++cluster) {
        const double distance = (cells[cell].u.velocity() -
                                 means[cluster]).sqr();
        if (distance < best_distance) {
          best_distance = distance;
          best_cluster = cluster;
        }
      }
      assignments[cell] = best_cluster;
    }

    // Sum up and count cells for each cluster.
    std::vector<smash::ThreeVector> new_means(number_of_clusters);
    std::vector<size_t> counts(number_of_clusters, 0);
    for (size_t cell = 0; cell < cells.size(); ++cell) {
      const auto cluster = assignments[cell];
      new_means[cluster] += cells[cell].u.velocity();
      counts[cluster] += 1;
    }

    // Divide sums by counts to get new centroids.
    for (size_t cluster = 0; cluster < number_of_clusters; ++cluster) {
      // Turn 0/0 into 0/1 to avoid zero division.
      const auto count = std::max<size_t>(1, counts[cluster]);
      means[cluster] = new_means[cluster] / count;
    }
  }

  return means;
}

#endif // KMEANS_CLUSTERING
