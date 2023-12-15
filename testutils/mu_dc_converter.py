#!/env /usr/bin/python

import numpy as np
import os,sys
import h5py as hdf5
import re

try:
    file=sys.argv[1]
    what=sys.argv[2]
except:
    print("Usage: mu_dc_converter.py #FILE #mu/dc (optional Wterm Vterm)")

try:
    w_term=sys.argv[3]
    v_term=sys.argv[4]
except:
    w_term=47.12
    v_term=48.33
    
f = hdf5.File(file, 'r')
regex_dmft=re.compile(".*dmft")
iters = list(filter(regex_dmft.match, f))
if len(iters)==0:
    regex_stat=re.compile(".*stat")
    iters = list(filter(regex_stat.match, f))
    niter=iters[-2]
else:
    niter=iters[-3]


xmu = f[str(niter)+'/mu/value'][()]
dc=f[niter+"/ineq-001/dc/value"][()]
f_occ = f[niter+"/ineq-001/occ/value"][()]
gdensold = f[niter+"/gdensold/value"][()]

nuf=0.0
nuf_occ=0.0
nuc_type1=0.0
nuc_type2=0.0

for i in range(4):
    nuf_occ = nuf_occ + f_occ[i,0,i,0] + f_occ[i,1,i,1]
for i in [4,5,8,9]:
    nuc_type1 = nuc_type1 + gdensold[i,0,i,0] + gdensold[i,1,i,1]
for i in [6,7,10,11]:
    nuc_type2 = nuc_type2 + gdensold[i,0,i,0] + gdensold[i,1,i,1]

if what=="mu":
    shifted_mu = np.real(xmu) + w_term*(np.real(nuf_occ)-4) + v_term*(np.real(nuc_type1+nuc_type2)-8)
    print(shifted_mu)
if what=="dc":
    shifted_dc=np.zeros(12)
    for i in range(0,4):
        shifted_dc[i] = dc[i,0] + w_term*(np.real(nuc_type1+nuc_type2)-8)
    for i in range(4,12):
        shifted_dc[i] = w_term*(np.real(nuc_type1+nuc_type2)-8)
    print(*shifted_dc, sep=", ")


