import sys

input_file = sys.argv[1]
sampler_printout = sys.argv[2]

with open(sampler_printout, 'r') as f:
    while True:
        line = f.readline()
        if 'Obtained total momentum:' in line:
            energy = float(line.split('(')[1].split()[0])
            break


pi0 = 0
rho0 = 0
eta = 0
with open(input_file, 'r') as f:
    while True:
        line = f.readline()
        if (not line): break
        if ('#' in line): continue
        if (' 111 ' in line): pi0 += 1
        if (' 113 ' in line): rho0 += 1
        if (' 221 ' in line): eta += 1
print energy, pi0, rho0, eta
