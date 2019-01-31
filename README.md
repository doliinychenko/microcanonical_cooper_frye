# Sampling particles on a hypersurface with account of energy, momentum, baryon number, strangeness and charge conservation.

Turning relativistic hydrodynamics into
particles including local event-by-event conservation laws: energy, momentum, baryon number,
strangeness, electric charge. Isospin projection conservation follows from the
last three.

The idea is to apply Metropolis algorithm on a particlization hypersurface.
In general, Metropolis algorithm always follows these steps:
1. Propose a new state (set of particles in my case) given the previous one
2. Accept or reject it with properly computed probability, which depends
   on the desired stationary distribution and on the proposal function.

In this project the proposal function is very similar to 2<->3 collisions,
which ensures all the conservation laws. One can imagine the outcome as a Markov chain,
that travels in a special subspace of {cell(i), p_i}, where i varies from 1 to N, and N
is variable), where the total energy, momentum and quantum numbers are the same.

## Inputs

A file with a list of hypersurface elements. Each element is defined by
T    - temperature
muB  - baryon chemical potential
muS  - strangeness chemical potential
muQ  - electric charge chemical potential
v    - 3-velocity of the hypersurface element
dsigma^mu - normal 4-vector to the hypersurface element, note that the index is upper
       For definition see Eq. (2) in https://arxiv.org/pdf/1206.3371.pdf

## Outputs

Sampled particles, characterized by
- type
- momentum
- number of cell

## Current state

The project is currently at the research and development stage. It
is tested in a number of simple cases and already optimized to run fast enough
for practical usage. The code still needs more documentation and clean up,
but this is postponed until the paper is written.

The sampled local distribution is limited to Boltzmann, but it is extendable
to viscous-corrected distribution and quantum statistics.


## Prequisites

1. Eigen3 library (eigen.tuxfamily.org)
2. GNU Scientific Library (GSL) >= 1.15
3. SMASH library (smash-transport.github.io), which requires the first two anyway.

## Compiling

- Set the SMASH_DIR environment variable by executing

      export SMASH_DIR=[...]/smash

- Copy the cmake files to find the libraries for the project

      cp -r $SMASH_DIR/cmake ..

- Run cmake and make

      mkdir build && cd build
      cmake ..
      make

- If Pythia is not found, try

      cmake .. -DPythia_CONFIG_EXECUTABLE=[...]/pythia8235/bin/pythia8-config
