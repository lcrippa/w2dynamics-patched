# imports
import h5py as hdf5
import numpy as np

import os
import glob
import sys

### load stuff
filename="sc-2017-02-23-Thu-17-19-06.hdf5"
f=hdf5.File(filename,"r")
siw = f["dmft-010"]["ineq-001"]["siw-full"]["value"].value
iw = f[".axes"]["iw"].value
Niw=siw.shape[-1]

### make (2*Nb,2*Nb,:) out of (Nb,2,Nb,2,:)
siw=siw.reshape(2,2,-1)

### unitary matrix to rotate
sx=np.zeros(shape=(2,2))
sx[0,1]=1
sx[1,0]=1
sz=np.zeros(shape=(2,2))
sz[0,0]=1
sz[1,1]=-1

### rotation angle
b=0.7
A=b*sx+np.sqrt(1-b**2)*sz

siw_rot=np.zeros_like(siw)

### rotate
for w in range(0,Niw):
    for i in range(0,2):
        for j in range(0,2):
            for ip in range(0,2):
                for jp in range(0,2):
                    siw_rot[i,j,w]=siw_rot[i,j,w]+A[i,ip]*A[j,jp]*siw[ip,jp,w]


### now write out the rotated self energy
np.set_printoptions(suppress=True)
np.set_printoptions(precision=9)


data=np.column_stack((iw,np.real(siw_rot[0,0,:]),np.imag(siw_rot[0,0,:])))
np.savetxt("siw_rot00.dat",data)

data=np.column_stack((iw,np.real(siw_rot[0,1,:]),np.imag(siw_rot[0,1,:])))
np.savetxt("siw_rot01.dat",data)

data=np.column_stack((iw,np.real(siw_rot[1,0,:]),np.imag(siw_rot[1,0,:])))
np.savetxt("siw_rot10.dat",data)

data=np.column_stack((iw,np.real(siw_rot[1,1,:]),np.imag(siw_rot[1,1,:])))
np.savetxt("siw_rot11.dat",data)
