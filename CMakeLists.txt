cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(microcanonical_sampling)

find_package(Git REQUIRED)
find_package(GSL 1.15 REQUIRED)

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
ExternalProject_Add(
    smash
    PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/
    GIT_REPOSITORY https://github.com/smash-transport/smash.git
    TIMEOUT 10
    UPDATE_COMMAND ${GIT_EXECUTABLE} pull
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
)

add_executable(microcanonical test_sampler.cc microcanonical_sampler.cc hydro_cells.cc sampler_particletype_list.cc
               src/smash/src/pdgcode.cc
               src/smash/src/threevector.cc
               src/smash/src/fourvector.cc
               src/smash/src/stringfunctions.cc )
include_directories(include)
include_directories(src/smash/src/include)

# Set the relevant generic compiler flags (optimisation + warnings)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -fopenmp -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -std=c++11 -mfpmath=sse")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Werror=switch")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

message("GSL includes: ${GSL_INCLUDE_DIR}")
include_directories(${GSL_INCLUDE_DIR})
target_link_libraries(microcanonical ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})

add_dependencies(microcanonical smash)
