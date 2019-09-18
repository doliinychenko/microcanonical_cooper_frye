import sys
import numpy as np
input_file = sys.argv[1]
V = float(input_file.split('_V')[1].split('_T')[0])
T = float(input_file.split('_T')[1].split('.dat')[0])
mean_gce = V * T*T*T / (3.141592635**2 * 0.197327053**3)

nevents = 0
mult = []
with open(input_file, 'r') as f:
    # header
    f.readline()
    f.readline()
    while(True):
        line = f.readline()
        if (not line):
            nevents -= 1
            break
        if ('#' in line):
            if (not 'end' in line):
                nevents += 1
            elif (nevents > 1):
                mult.append(mul)
            mul = 0
        else:
            mul += 1
x = np.array(mult)

mean = x.mean()
dx = x - mean
mu2 = (dx*dx).mean()
mu3 = (dx*dx*dx).mean()
mu4 = (dx*dx*dx*dx).mean()
mu5 = (dx*dx*dx*dx*dx).mean()
mu6 = (dx*dx*dx*dx*dx*dx).mean()
mu7 = (dx*dx*dx*dx*dx*dx*dx).mean()
mu8 = (dx*dx*dx*dx*dx*dx*dx*dx).mean()
sigma = np.sqrt(mu2)
m3 = mu3 / sigma**3
m4 = mu4 / sigma**4
m5 = mu5 / sigma**5
m6 = mu6 / sigma**6
m7 = mu7 / sigma**7
m8 = mu8 / sigma**8

k2 = mu2
k3 = mu3
k4 = mu4 - 3 * mu2 * mu2
c2 = mu2 - mean
c3 = 2*mean - 3*k2 + k3
c4 = -6*mean + 11*k2 - 6*k3 + k4

mean_err = sigma / np.sqrt(nevents)
sigma2_over_mean_err = np.sqrt(mu4 - sigma**4) / mean / np.sqrt(nevents)
k3_k2_err = np.sqrt(6 * sigma**2 / nevents) # np.sqrt(9 - 6*m4 + m3**2 * (35 + 9 * m4) / 4 - 3 * m3 * m5 + m6) / np.sqrt(nevents)
k4_k2_err = 2 * sigma * k3_k2_err #np.sqrt(-m4**2 + 4*m4**3 + 16 * m3**2 * (1 + m4) - 8 * m3 * m5 - 4 * m4 * m6 + m8) / np.sqrt(nevents)

assert(x.size == nevents)
print '%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f' % (mean_gce, mean, mean_err,  k2 / mean, sigma2_over_mean_err, k3/k2, k3_k2_err, k4/k2, k4_k2_err)
