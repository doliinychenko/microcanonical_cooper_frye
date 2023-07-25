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

# Copypasted from SMASH cmake
include(TestBigEndian)
TEST_BIG_ENDIAN(IS_BIG_ENDIAN)
if(IS_BIG_ENDIAN)
    message(STATUS "Big endian architecture detected.")
    option(FORCE_USE_ON_BIG_ENDIAN "Force cmake to setup and compile on big endian machine." OFF)
    if(NOT FORCE_USE_ON_BIG_ENDIAN)
        message(FATAL_ERROR " \n" # See e.g. https://stackoverflow.com/a/51035045 for formatting
                            " Using SMASH on big endian machines is not officially supported.\n"
                            " Please, open an issue at https://github.com/smash-transport/smash/issues\n"
                            " if you are interested in SMASH team supporting big endian architectures.\n"
                            " At your own risk you can pass -DFORCE_USE_ON_BIG_ENDIAN=TRUE to cmake in\n"
                            " order to setup and compile SMASH on a big endian machine.")
    endif()
    add_definitions("-DBIG_ENDIAN_ARCHITECTURE")
else()
    message(STATUS "Little endian architecture detected.")
    add_definitions("-DLITTLE_ENDIAN_ARCHITECTURE")
endif()

add_executable(microcanonical src/test_sampler.cc ${SOURCE_FILES})
add_library(microcanonical_sampling ${SOURCE_FILES})
#target_link_libraries(microcanonical_sampling PRIVATE OpenMP::OpenMP_CXX)

include_directories(include)
include_directories(smash/src/include)

# Set the relevant generic compiler flags (optimisation + warnings)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3") #-fopenmp -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -std=c++17 -mfpmath=sse")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Werror=switch")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

message("GSL includes: ${GSL_INCLUDE_DIR}")
include_directories(${GSL_INCLUDE_DIR})
target_link_libraries(microcanonical ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
