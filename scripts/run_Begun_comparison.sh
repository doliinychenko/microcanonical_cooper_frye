#/bin/bash

Vlist=`python -c "import numpy as np; a = np.logspace(np.log(0.69), np.log(20), num = 200, base = np.exp(1.0)); print ' '.join(map(str, a))"`

for V in ${Vlist}
do
    ../scripts/microcanonical_sampler_box.sh -V $V -T 0.16 -p "../../particle_table_iSS_SMASH/smash_particles_iSS_format.dat,iSS" &> sampler_printout
    python ../scripts/av_mult_counter.py static_box_V${V}_T0.16.dat sampler_printout
done
