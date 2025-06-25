# imports
import h5py as hdf5
import numpy as np

import os
import glob
import sys

### read out stuff
filename="new_sc_highstat-2017-03-04-Sat-20-53-34.hdf5"
f=hdf5.File(filename,"r")
siw = f["dmft-004"]["ineq-001"]["siw-full"]["value"].value
iw = f[".axes"]["iw"].value
Niw=siw.shape[-1]

### rotation matrix
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

### rotate siw in spin space
for w in range(0,Niw):
   for i in range(0,2):
      for j in range(0,2):
         for ip in range(0,2):
            for jp in range(0,2):
               siw_rot[0,i,0,j,w]+=A[i,ip]*A[j,jp]*siw[0,ip,0,jp,w]
               siw_rot[1,i,1,j,w]+=A[i,ip]*A[j,jp]*siw[1,ip,1,jp,w]


### now write the new self-energy
np.set_printoptions(suppress=True)
np.set_printoptions(precision=9)

siw_rot=siw_rot.reshape(4,4,-1)

data=np.column_stack((iw,np.real(siw_rot[0,0,:]),np.imag(siw_rot[0,0,:])))
np.savetxt("siw_rot00.dat",data)
data=np.column_stack((iw,np.real(siw_rot[0,1,:]),np.imag(siw_rot[0,1,:])))
np.savetxt("siw_rot01.dat",data)
data=np.column_stack((iw,np.real(siw_rot[0,2,:]),np.imag(siw_rot[0,2,:])))
np.savetxt("siw_rot02.dat",data)
data=np.column_stack((iw,np.real(siw_rot[0,3,:]),np.imag(siw_rot[0,3,:])))
np.savetxt("siw_rot03.dat",data)

data=np.column_stack((iw,np.real(siw_rot[1,0,:]),np.imag(siw_rot[1,0,:])))
np.savetxt("siw_rot10.dat",data)
data=np.column_stack((iw,np.real(siw_rot[1,1,:]),np.imag(siw_rot[1,1,:])))
np.savetxt("siw_rot11.dat",data)
data=np.column_stack((iw,np.real(siw_rot[1,2,:]),np.imag(siw_rot[1,2,:])))
np.savetxt("siw_rot12.dat",data)
data=np.column_stack((iw,np.real(siw_rot[1,3,:]),np.imag(siw_rot[1,3,:])))
np.savetxt("siw_rot13.dat",data)

data=np.column_stack((iw,np.real(siw_rot[2,0,:]),np.imag(siw_rot[2,0,:])))
np.savetxt("siw_rot20.dat",data)
data=np.column_stack((iw,np.real(siw_rot[2,1,:]),np.imag(siw_rot[2,1,:])))
np.savetxt("siw_rot21.dat",data)
data=np.column_stack((iw,np.real(siw_rot[2,2,:]),np.imag(siw_rot[2,2,:])))
np.savetxt("siw_rot22.dat",data)
data=np.column_stack((iw,np.real(siw_rot[2,3,:]),np.imag(siw_rot[2,3,:])))
np.savetxt("siw_rot23.dat",data)

data=np.column_stack((iw,np.real(siw_rot[3,0,:]),np.imag(siw_rot[3,0,:])))
np.savetxt("siw_rot30.dat",data)
data=np.column_stack((iw,np.real(siw_rot[3,1,:]),np.imag(siw_rot[3,1,:])))
np.savetxt("siw_rot31.dat",data)
data=np.column_stack((iw,np.real(siw_rot[3,2,:]),np.imag(siw_rot[3,2,:])))
np.savetxt("siw_rot32.dat",data)
data=np.column_stack((iw,np.real(siw_rot[3,3,:]),np.imag(siw_rot[3,3,:])))
np.savetxt("siw_rot33.dat",data)
