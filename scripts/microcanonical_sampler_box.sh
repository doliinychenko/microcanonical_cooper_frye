#!/bin/bash

programname=$0

function usage {
    echo "Microcanonical sampling of particles in a static box"
    echo
    echo "usage: $programname -V Volume -T temperature [-p particles_file] [-s path/to/sampler]"
    echo
    echo "  -V Volume:  volume of the box in fm^3"
    echo "  -T Temperature:  temperature of particles in the box in GeV"
    echo "  -p particles.txt"
    echo "     files with hadron list and their decays in SMASH format"
    echo "     (see http://theory.gsi.de/~smash/doc/1.6/inputparticles.html)"
    echo "     By default full SMASH hadron list with masses < 2.5 GeV is used."
    echo "     Note: sampler may not work for an arbitrary list of hadrons."
    echo "     For example, sampling only pi+ and pi- gas does not work,"
    echo "     because 2<->3 collisions are not enough for them, one needs"
    echo "     pi^+pi- <-> pi^+pi-pi^+pi-. But if one adds pi^0, then sampling works."
    echo "  -s path to the sampler binary"
    echo "     default is ./microcanonical"
    echo
   exit 1
}

sampler="./microcanonical"
particles_decaymodes=""

while getopts ":V:T:p:s:h" opt; do
  case $opt in
    V)
      V="$OPTARG"
      ;;
    T)
      T="$OPTARG"
      ;;
    p)
      particles="-p $OPTARG"
      ;;
    s)
      sampler="$OPTARG"
      ;;
    h)
      usage
      exit 1
      ;;
    \?)
      echo "Invalid option: -$OPTARG"
      usage
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument."
      usage
      exit 1
      ;;
  esac
done

if [ "$#" -lt 2 ] || [ "$#" -gt 8 ]; then
  usage
  exit 1
fi

echo "Volume $V fm^3, Temperature $T GeV"
if [ ! -z $particles ]; then
  echo "Taking particle list from ${particles}"
fi

hyper_file="hydro_cells_static_V${V}_T${T}"
echo "${V} 0.0 0.0 0.0  0.0 0.0 0.0 ${T} 0.0 0.0 0.0" > ${hyper_file}
echo "Running ${sampler}"
time ${sampler} ${particles},SMASH -s ${hyper_file},Dima_format -o static_box_V${V}_T${T}.dat -n 1000000
echo "Finished running ${sampler}"
