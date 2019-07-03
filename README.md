# Sampling particles on a hypersurface with local event-by-event energy, momentum, baryon number, strangeness, and charge conservation.

## Table of Contents

- [What is this code doing?](#purpose)
- [How to install it?](#installing)
- [How to use it?](#usage)
  - [Running and command line options](#running)
  - [Inputs](#inputs)
  - [Outputs](#outputs)
- [Testing physics](#tests)
- [Troubleshooting and feature requests](#bugs)

## What is this code doing? <a name = "purpose"></a>

**Spoiler: Turning relativistic hydrodynamics into particles including local
event-by-event conservation laws: energy, momentum, baryon number, strangeness,
and electric charge.**

Also, under particular configuration this code can be applied as a fast
microcanonical sampler (see [usage](#usage)). It turns into a microcanonical sampler, if
the provided hypersurface is particularly simple --- a static box.
Microcanonical samplers were implemented before (see [a paper by Werner and
Aichelin](https://arxiv.org/pdf/nucl-th/9503021.pdf) and [a paper by Becattini
and Ferroni](https://arxiv.org/pdf/hep-ph/0407117.pdf)), but these
implementations used different methods and the codes are proprietary.


#### Background
One of the tools to simulate relativistic heavy ion collisions is relativistic
hydrodynamics.  To compare its outcome with experiments, eventually
hydrodynamical variables (pressure, energy density, etc) have to be transformed
into particles. This is called particlization. Usually particlization is done
in a way that does not conserve energy, momentum, and quantum numbers (baryon
number, strangeness, charge) in each single event (event is a single run of a
hydroynamical simulation). Instead, only conservation on the event average is
fulfilled. For more detailes on the usual particlization see [a paper by P.
Huovinen and H. Petersen](https://arxiv.org/pdf/1206.3371.pdf), mainly Sections
3 and 4.

#### Physical applications

Event-by event local conservation laws are not always important, most studies,
which address particle rapidity, transverse momentum, or radial angle distributions,
don't really need them to be local and even-by-event.

However, there are cases, where they play an important role:

- Studying correlations and fluctuations of particle numbers (in contrast to means,
   which are typically a more popular subject of research)
- Switching from fluctuating / stochastic hydrodynamics to particles
- Studying small collision systems (in my field Au+Au nuclei collision is still
  a large system, haha; small system is, for example, d+Au or C+C)
- Studying rare particles, such as Omega(sss) baryon
- Studying the background for chiral magnetic effect


#### Implementation idea

The hypersurface is split into patches. Conservation is enforced in every patch.
This provides localness. The "event-by-event" property is much harder to achieve.
One has to sample a particular distribution with restrictions. This is usually hard
and the only way out is to apply Metropolis algorithm.

So I apply a Metropolis (or a Markov chain Monte-Carlo, it's the same) on a
particlization hypersurface.  In general, Metropolis algorithm always follows
these steps:
1. Propose a new state (set of particles in my case) given the previous one
2. Accept or reject it with properly computed probability, which depends
   on the desired stationary distribution and on the proposal function.

In this project the proposal function is very similar to 2<->3 collisions,
which ensures all the conservation laws. One can imagine the outcome as a Markov chain,
that travels in a special subspace of {cell(i), p_i}, where i varies from 1 to N, and N
is variable), where the total energy, momentum and quantum numbers are the same.
If this is not very clear, don't worry, read [my preprint](https://arxiv.org/pdf/1902.09775.pdf).

#### Current state

The project is already at the stage, when I know that it works, because
it is tested in a number of simple cases and already optimized to run fast enough
for practical usage. For the derivation, explanations, and first results see
[this preprint](https://arxiv.org/pdf/1902.09775.pdf). A more detailed paper
is under preparation. However, some particular features (such as viscous corrections)
are still missing and the further testing is ongoing.



## How to install <a name = "installing"></a>

This project relies heavily on the [SMASH](https://smash-transport.github.io/) library.
It uses particles and their properties from SMASH, as well as some nice SMASH
infrastructure.


### Prequisites

Before compiling install these:

1. Eigen3 library (eigen.tuxfamily.org)
2. GNU Scientific Library (GSL) >= 1.15
3. SMASH library (smash-transport.github.io), which requires the first two anyway.

### Compiling

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

## How to use it? <a name = "usage"></a>


### Inputs <a name = "input"></a>

A file with a list of hypersurface elements. Each element is defined by
T    - temperature
muB  - baryon chemical potential
muS  - strangeness chemical potential
muQ  - electric charge chemical potential
v    - 3-velocity of the hypersurface element
dsigma^mu - normal 4-vector to the hypersurface element, note that the index is upper
       For definition see Eq. (2) in https://arxiv.org/pdf/1206.3371.pdf

### Outputs  <a name = "output"></a>

Sampled particles, characterized by
- type
- momentum
- number of cell


## Testing physics <a name = "testing"></a>

### Reproducing the results shown in Fig. 1 of arxiv.org/pdf/1902.09775.pdf

- Run

      cd build
      ./microcanonical -r > sampled_output.txt
      python ../scripts/mult_and_corr.py sampled_output.txt

  It takes around 2 minutes to generate 10^5 samples.

## Troubleshooting and feature requests   <a name = "bugs"></a>

  If you are interested in some feature, you can either create an issue right
here, at github, or write me an e-mail (doliinychenko [at] lbl.gov). Same thing
if you find a bug.
