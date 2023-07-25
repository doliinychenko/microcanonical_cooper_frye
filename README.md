[![Build Status](https://travis-ci.org/doliinychenko/microcanonical_cooper_frye.svg?branch=master)](https://travis-ci.org/doliinychenko/microcanonical_cooper_frye)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

# Sampling particles on a hypersurface with local event-by-event energy, momentum, baryon number, strangeness, and charge conservation.

## Table of Contents

- [What is this code doing?](#purpose)
- [How to install it?](#installing)
- [How to use it?](#usage)
  - [Running and command line options](#running)
  - [Inputs](#inputs)
  - [Outputs](#outputs)
  - [Using as a microcanonical sampler](#microcanonical)
- [Testing physics](#tests)
- [Troubleshooting and feature requests](#bugs)

## What is this code doing? <a name = "purpose"></a>

**Spoiler: Turning relativistic hydrodynamics into particles including local
event-by-event conservation laws: energy, momentum, baryon number, strangeness,
and electric charge.**

Also, under particular configuration this code can be applied as a fast
microcanonical sampler (see [microcanonical](#microcanonical)). It turns into a microcanonical sampler, if
the provided hypersurface is particularly simple --- a static box.
Microcanonical samplers were implemented before (see [a paper by Werner and
Aichelin](https://arxiv.org/pdf/nucl-th/9503021.pdf) and [a paper by Becattini
and Ferroni](https://arxiv.org/pdf/hep-ph/0407117.pdf)), but these
implementations used different methods, the codes are proprietary,
and I believe that my method is faster.


#### Background
One of the tools to simulate relativistic heavy ion collisions is relativistic
hydrodynamics.  To compare its outcome with experiments, eventually
hydrodynamical variables (pressure, energy density, etc) have to be transformed
into particles. This is called particlization. Usually particlization is done
in a way that does not conserve energy, momentum, and quantum numbers (baryon
number, strangeness, charge) in each single event (event is a single run of a
hydrodynamical simulation). Instead, only conservation on the event average is
fulfilled. For more detailes on the usual particlization see [a paper by P.
Huovinen and H. Petersen](https://arxiv.org/pdf/1206.3371.pdf), mainly Sections
3 and 4.

#### Physical applications

Event-by event local conservation laws are not always important, most studies,
which address particle rapidity, transverse momentum, or radial angle distributions,
don't really need conservation laws to be local and even-by-event.

However, there are cases, where these play an important role:

- Studying correlations and fluctuations of particle numbers (in contrast to means,
   which are typically a more popular subject of research)
- Switching from fluctuating / stochastic hydrodynamics to particles
- Studying small collision systems (in my field Au+Au nuclei collision is still
  a large system, even though it may be 20 just femtometers large;
  an actual small system is, for example, d+Au or C+C)
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
which ensures all the conservation laws. One can imagine the outcome as a
Markov chain, that travels in a special subspace of {cell(i), p_i}, where i
varies from 1 to N, and N is variable), where the total energy, momentum and
quantum numbers are the same.  If this is not very clear, don't worry, read [my
preprint](https://arxiv.org/pdf/1902.09775.pdf).

#### Current state

The project is already at the stage, when I know that it works, because
it is tested in a number of simple cases and already optimized to run fast enough
for practical usage. For the derivation, explanations, and first results see
[this preprint](https://arxiv.org/pdf/1902.09775.pdf). A more detailed paper
is [also available](https://arxiv.org/pdf/2001.08176). Some particular
features (such as viscous corrections)
are still missing and the further testing is ongoing.



## How to install <a name = "installing"></a>

This project uses some code from the [SMASH](https://smash-transport.github.io/),
mostly particles, their properties, and some nice SMASH infrastructure.


### Prequisites

Before compiling install these:

1. cmake
2. GNU Scientific Library (GSL) >= 1.15
3. git

### Compiling

- In the main project directory run

      mkdir build && cd build
      cmake ..
      make

## How to use it? <a name = "usage"></a>

If you are using this code, please cite [this publication](http://inspirehep.net/record/1722034/export/hlxu).


### Running and command line options  <a name = "running"></a>

Most of the features are controlled through the command-line options.
Run

      ./microcanonical --help

to see a list of these options.

### Inputs <a name = "input"></a>

**Hypersurface** is the first thing that the sampler needs. It is a list
of computational hydrodynamic cells.  Eventually, the code needs
to know the following quantities for every cell:
- space-time position x
- temperature T
- baryon chemical potential muB
- strangeness chemical potential muS
- electric charge chemical potential muQ
- 3-velocity v
- normal 4-vector to the hypersurface element,
  for definition see Eq. (2) in [Particlization in hybrid models](https://arxiv.org/pdf/1206.3371.pdf).

The sampler can read hypersurface files in different formats.
All formats are multiline files, where each line represents one hydrodynamic computational cell.

One can specify the input file and format with a command-line option, for example:

      ./microcanonical --surface  ../MUSIC_hypersurfaces_200GeV/hyper5.dat,MUSIC_format
      ./microcanonical --surface  ../UrQMD_hyper/test_hyper.dat,Steinheimer_format
      ./microcanonical --surface  ../from_VISHNU/,VISH_format

Format specifications: To be described here.

**Particle list** to be sampled is the second thing that the sampler needs.  It
is the list of possible particles to sample. The list itself is given in a file
in [SMASH format](http://theory.gsi.de/~smash/doc/1.6/inputparticles.html)
or in iSS format.
Normally, it is expected that the users just take the full SMASH particles.txt
file and remove particles until they get the desired particle list.

By default the default SMASH particle list is used. Here is an example of sampling
a custom particle list:

      ./microcanonical --particles  ../my_edited_particles.txt,SMASH
      ./microcanonical --particles  only_pi0.txt,SMASH

Parameters like patch size or number of events can be specified via command-line options:

      ./microcanonical --energy_patch 15.5
      ./microcanonical --nevents 10000

### Outputs  <a name = "output"></a>

**Sampled particles** are by default printed out into sampled_particles.dat.
This can be changed, as in this example

      ./microcanonical --output_file analyze_next_week/sampled_on_my_lovely_hypersurface134.out

The output format is

      # sample 0
      t x y z E px py pz pdgid
      t x y z E px py pz pdgid
      t x y z E px py pz pdgid
      ...
      # sample 1
      t x y z E px py pz pdgid
      t x y z E px py pz pdgid
      t x y z E px py pz pdgid
      ...

Here (t, x, y, z) [fm/c, fm, fm, fm] is a Cartesian space-time position of the particle,
(E, px, py, pz) [GeV, GeV/c, GeV/c, GeV/c] is its 4-momentum,
pdgid is a [standard Particle Data Group code defining particle identity](http://pdg.lbl.gov/2002/montecarlorpp.pdf).

**Hypersurface partitioning into patches** is not dumped by default, but can be optionally printed out:

      ./microcanonical --labeled_hyper_output  my_hypersurface_partitioning.txt

### Using as a microcanonical sampler  <a name = "microcanonical"></a>

Microcanonical sampling (see [a paper by Werner and
Aichelin](https://arxiv.org/pdf/nucl-th/9503021.pdf) and [a paper by Becattini
and Ferroni](https://arxiv.org/pdf/hep-ph/0407117.pdf)) is a special case of
this sampler. If you wish to use it in the microcanonical sampler mode, you can
take advantage of a script intended specifically for this purpose:

      ../scripts/microcanonical_sampler_box.sh  -V <Volume in fm^3> -T <Temperature in GeV>
      
It may be confusing that a temperature has to be specified for a microcanonical sampler instead
of the energy. This is done for a convenient comparison with grand-canonical samplers.
The temperature is used to compute total energy E in a grand-canonical case, and then this
energy E is used for microcanonical sampling.

## Testing physics <a name = "testing"></a>

### Reproducing the results shown in Fig. 1 of arxiv.org/pdf/1902.09775.pdf

- Run

      cd build
      ./microcanonical -r > sampled_output.txt
      python ../scripts/mult_and_corr.py sampled_output.txt

  It takes around 2 minutes to generate 10^5 samples.

## Troubleshooting and feature requests   <a name = "bugs"></a>

  If you are interested in some feature, you can either create an issue right
here, at github, or write me an e-mail (dmytrooliinychenko [at] gmail). Same thing
if you find a bug.
