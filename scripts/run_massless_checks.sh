#!/bin/bash

Vlist="8 16 24 32 48 64 125 250 500 900 1000 1100 1300 1500 2000 2200"
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

echo "π0            0.0001   0.0    -      111" > massless_plist

for V in ${Vlist}
do
    ${SCRIPTPATH}/microcanonical_sampler_box.sh -p massless_plist,SMASH -V ${V} -T 0.15 -n 1000000
done

touch massless_sampling_results.txt
for V in ${Vlist}
do
    python ${SCRIPTPATH}/n_mean_var.py static_box_V${V}_T0.15.dat >> massless_sampling_results.txt
done
