import sys
sys.path.append("navin3.7/lib/python3.7/site-packages")
import numpy as np
import hod as h
from math import log10

p = h.cosmo()
q = h.hodpars()
p.Om0 = 0.307115
p.w0 = -1
p.wa = 0
p.Omk = 0.0
p.hval = 0.6777
p.Omb = 0.048206
p.th = 2.726
p.s8 = 0.8228
p.nspec = 0.96
p.ximax = log10(8.0)
p.cfac = 1.0
q.Mmin = 13.0
q.siglogM = 0.5
q.Msat = 14.0
q.alpsat = 1.0
q.Mcut = 13.5
q.csbycdm = 1.0
q.fac = 1.0

a = h.hod(p, q)

def getdblarr(r):
    temp=h.doubleArray(r.size)
    for i in range(r.size):
        temp[i]=r[i]
    return temp

mass = np.linspace(9.0, 16.0, 100)
ncen = mass*1.0
nsat = mass*1.0

a.init_Nc_spl(getdblarr(mass), getdblarr(ncen), ncen.size)
a.init_Ns_spl(getdblarr(mass), getdblarr(nsat), nsat.size)
