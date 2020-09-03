cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(microcanonical_sampling)

find_package(Git REQUIRED)
find_package(GSL 1.15 REQUIRED)
#find_package(OpenMP REQUIRED)

execute_process(COMMAND git clone  https://github.com/smash-transport/smash.git
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
set(SOURCE_FILES
    src/microcanonical_sampler.cc
    src/hydro_cells.cc
    src/sampler_particletype_list.cc
    src/statistics_summary.cc

    smash/src/pdgcode.cc
    smash/src/threevector.cc
    smash/src/fourvector.cc
    smash/src/stringfunctions.cc
)
add_executable(microcanonical src/test_sampler.cc ${SOURCE_FILES})
add_library(microcanonical_sampling ${SOURCE_FILES})
#target_link_libraries(microcanonical_sampling PRIVATE OpenMP::OpenMP_CXX)

include_directories(include)
include_directories(smash/src/include)

# Set the relevant generic compiler flags (optimisation + warnings)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3") #-fopenmp -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -std=c++11 -mfpmath=sse")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Werror=switch")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

message("GSL includes: ${GSL_INCLUDE_DIR}")
include_directories(${GSL_INCLUDE_DIR})
target_link_libraries(microcanonical ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
