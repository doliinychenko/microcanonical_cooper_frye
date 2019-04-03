cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(microcanonical_sampling)

# Tell cmake where to find modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

add_executable(microcanonical test_sampler.cc microcanonical_sampler.cc hydro_cells.cc)
include_directories(include)

# Set the relevant generic compiler flags (optimisation + warnings)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -fopenmp -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Werror=switch -std=c++11 -mfpmath=sse")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

find_package(SMASH)
find_package(GSL 1.15 REQUIRED)

message("GSL includes: ${GSL_INCLUDE_DIR}")
include_directories(${GSL_INCLUDE_DIR})
target_link_libraries(microcanonical ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})

if(${SMASH_FOUND})
  include_directories(${SMASH_INCLUDE_DIR})
  target_link_libraries(microcanonical ${SMASH_LIBRARIES})
endif(${SMASH_FOUND})
