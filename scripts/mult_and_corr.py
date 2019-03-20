# coding=UTF-8
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys

inputfile = sys.argv[1]
print inputfile
with open(inputfile, 'r') as f:
    for _ in xrange(19): f.readline()
    labels = f.readline().decode("utf-8").split()
    mult_size = len(f.readline().strip().split())
    Ncells = mult_size / len(labels)

print "# Ncells: ", Ncells
mult_all = np.loadtxt(inputfile, dtype = int, skiprows = 20)
Nsp = len(labels)
mult_tot = np.copy(mult_all[:, 0:Nsp])
print "# Hadron, means in cells, variance of total"
for i in xrange(1, Ncells):
    mult_tot += mult_all[:, i*Nsp:(i+1)*Nsp]
types_of_interest = [u'π⁺', u'π⁻', u'K⁺', u'K̅⁻', u'N⁺', u'Λ', u'η', u'Ξ⁻', u'Ω⁻']
for i in types_of_interest:
    ind = labels.index(i)
    print i,
    for j in xrange(Ncells):
        x = mult_all[:, j*Nsp + ind]
        print x.mean(),
    x = mult_tot[:, ind]
    print x.std()**2 / x.mean()


print "Correlations: (p,K)/(K,K), (pi,K)/(K,K),",\
      "(pi,p)/(pi,pi), (La,p)/(p,p), (p,K)/(K,K), (pi,p)/(pi,pi)"
ipipl, ipimin, iKpl, iKmin, ip, iLa = tuple([labels.index(i) for i in types_of_interest][:6])
dp = mult_tot[:, ip] - mult_tot[:, ip].mean()
dKp = mult_tot[:, iKpl] - mult_tot[:, iKpl].mean()
dpip = mult_tot[:, ipipl] - mult_tot[:, ipipl].mean()
dKmi = mult_tot[:, iKmin] - mult_tot[:, iKmin].mean()
dpimi = mult_tot[:, ipimin] - mult_tot[:, ipimin].mean()
dLa = mult_tot[:, iLa] - mult_tot[:, iLa].mean()
print (dp*dKp).mean() / mult_tot[:, iKpl].std()**2,\
      (dpip*dKp).mean() / mult_tot[:, iKpl].std()**2,\
      (dpip*dp).mean() / mult_tot[:, ipipl].std()**2,\
      (dLa*dp).mean() / mult_tot[:, ip].std()**2,\
      (dp*dKmi).mean() / mult_tot[:, iKmin].std()**2,\
      (dpimi*dp).mean() / mult_tot[:, ipimin].std()**2
