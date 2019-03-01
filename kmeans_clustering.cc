#include "include/kmeans_clustering.h"

#include <algorithm>

#include "smash/angles.h"
#include "smash/random.h"
#include "smash/threevector.h"

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

void test_clustering() {
  std::vector<HyperSurfacePatch::hydro_cell>cells;
  for (int i = 0; i < 100000; i++) {
    smash::Angles phitheta;
    phitheta.distribute_isotropically();
    const smash::ThreeVector v = smash::random::canonical() *
                                 phitheta.threevec();
    const double gamma = 1.0 / std::sqrt(1.0 - v.sqr());
    const smash::FourVector u(gamma, gamma * v);
    cells.push_back({smash::FourVector(),
                     smash::FourVector(),
                     u, 0.0, 0.0, 0.0, 0.0});
  }
  std::vector<size_t> assignments;
  assignments.resize(cells.size());
  k_means(cells, 10, 50, assignments);
  for (size_t i = 0; i < cells.size(); i++) {
    std::cout << cells[i].u.velocity() << " " << assignments[i] << std::endl;
  }

  std::array<int, 20> patch_counter{};
  for (size_t i = 0; i < cells.size(); i++) {
    patch_counter[assignments[i]]++;
  }
  for (size_t i = 0; i < 20; i++) {
    std::cout << patch_counter[i] << " ";
  }
  std::cout << std::endl;
}

