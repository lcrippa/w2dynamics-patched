#!/usr/bin/env python3
# imports
import h5py as hdf5
import numpy as np
import scipy.integrate as integrate
import os
import glob
import itertools
import sys
### here we need to add the path to the w2dynamics main directory
#auxdir="$HOME/opt/w2dynamics/lib/python3.*/site-packages/w2dyn"
try:
    auxdir="/home/hpc/b158cb/b158cb10/opt/w2dynamics/lib/python3.9/site-packages/w2dyn"
except:
    auxdir="/home/lorenzo/opt/w2dynamics/lib/python3.10/site-packages/w2dyn"
sys.path.insert(0,auxdir)
### in order to import the functions available in input.py
import auxiliaries.input as io
import auxiliaries.compound_index as indexobject

### load hamiltonian; spin_orbit=True, because the Hamiltonian
### we read, is spin-dependent

N_f = 4
N_c = 8
N_orb = N_f+N_c
N_spin = 2
u1 = 57.95
w_term = 47.12 
v_term = 48.33 
hfshift = -3.5*u1
which_iter='last'

f=hdf5.File(sys.argv[1],"r")
try:
    which_iter = sys.argv[2]
except:
    which_iter = 'last'


if which_iter == 'all':
    Listiter=list(f.keys())[4:-3]
else:
    Listiter=[list(f.keys())[-4]]

for Niter in Listiter:
    xmu = f[Niter+'/mu/value'][()]
    gdensold=np.zeros(shape=(N_orb,N_spin,N_orb,N_spin),dtype=complex)
    f_occ=np.zeros(shape=(N_f,N_spin,N_f,N_spin),dtype=complex)

    f_occ = f[Niter+"/ineq-001/occ/value"]
    gdensold = f[Niter+"/gdensold/value"]

    nuf=0.0
    nuf_occ=0.0
    nuc_type1=0.0
    nuc_type2=0.0

    for i in range(4):
        nuf = nuf + gdensold[i,0,i,0] + gdensold[i,1,i,1]
        nuf_occ = nuf_occ + f_occ[i,0,i,0] + f_occ[i,1,i,1]
    for i in [4,5,8,9]:
        nuc_type1 = nuc_type1 + gdensold[i,0,i,0] + gdensold[i,1,i,1]
    for i in [6,7,10,11]:
        nuc_type2 = nuc_type2 + gdensold[i,0,i,0] + gdensold[i,1,i,1]

    real_xmu=np.real(xmu) + w_term*(np.real(nuf_occ)-4) + v_term*(np.real(nuc_type1+nuc_type2)-8)
    print(np.real(xmu),real_xmu,np.real(nuf_occ),np.real(nuc_type1),np.real(nuc_type2),np.real(nuf_occ)+np.real(nuc_type1)+np.real(nuc_type2))

