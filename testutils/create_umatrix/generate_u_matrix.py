#!env /usr/bin/python\
import numpy as np
import os,sys

lbd=1
alpha22=1
NBANDS=4

alpha22=alpha22**2
interactionmatrix=np.zeros([NBANDS,2,NBANDS,2,NBANDS,2,NBANDS,2])
list_of_indices=np.array([[1,2,4,3],\
                          [2,1,3,4],\
                          [1,3,4,2],\
                          [3,1,2,4],\
                          [2,3,3,2],\
                          [3,2,2,3],\
                          [1,4,4,1],\
                          [4,1,1,4],\
                          [2,4,3,1],\
                          [4,2,1,3],\
                          [3,4,2,1],\
                          [4,3,1,2]])
x_coeffs=-1*np.array([-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])
y_coeffs=-1*np.array([+1,+1,-1,-1,-1,-1,-1,-1,-1,-1,+1,+1])

for iband in range(4):
    for jband in range(4):
        for kband in range(4):
            for lband in range(4):
                orbvector=[iband+1,jband+1,kband+1,lband+1]
                for itest in range(np.shape(list_of_indices)[0]):
                    
                    if np.array_equal(orbvector,list_of_indices[itest]):
                        interactionmatrix[iband,0,jband,0,kband,0,lband,0] = -1*lbd*alpha22*(x_coeffs[itest]+(y_coeffs[itest]))
                        interactionmatrix[iband,1,jband,1,kband,1,lband,1] = -1*lbd*alpha22*(x_coeffs[itest]+(y_coeffs[itest]))
                        break

kanamori_umatrix = np.genfromtxt("u_matrix_template.dat",dtype='str')
np.shape(kanamori_umatrix)
kanamori_umatrix[0:3,:]

for iline in range(np.shape(kanamori_umatrix)[0]):
    uvalue=float(kanamori_umatrix[iline,8])
    iband=int(kanamori_umatrix[iline,0])-1
    jband=int(kanamori_umatrix[iline,2])-1
    kband=int(kanamori_umatrix[iline,4])-1
    lband=int(kanamori_umatrix[iline,6])-1
    if kanamori_umatrix[iline,1]=="u":
        ispin=0
    else:
        ispin=1
    if kanamori_umatrix[iline,3]=="u":
        jspin=0
    else:
        jspin=1
    if kanamori_umatrix[iline,5]=="u":
        kspin=0
    else:
        kspin=1
    if kanamori_umatrix[iline,6]=="u":
        lspin=0
    else:
        lspin=1
    interactionmatrix[iband,ispin,jband,jspin,kband,kspin,lband,lspin] += uvalue

def sud(ispin):
    if(ispin)==0:
        spinstring="u"
    else:
        spinstring="d"
    return spinstring


fname="u_matrix.dat"
f=open(fname,'w')
f.write(' '.join([str(NBANDS),"BANDS"]))
f.write("\n")
for iband in range(4):
    for ispin in range(2):
        for jband in range(4):
            for jspin in range(2):
                for kband in range(4):
                    for kspin in range(2):
                        for lband in range(4):
                            for jspin in range(2):
                                uvalue=interactionmatrix[iband,ispin,jband,jspin,kband,kspin,lband,lspin]
                                if abs(uvalue) > 1e-15:
                                    orbvector=[str(iband+1),sud(ispin),str(jband+1),sud(jspin),str(kband+1),sud(kspin),str(lband+1),sud(kspin)]
                                    for icomp in range(8):
                                        f.write(orbvector[icomp])
                                        f.write(" ")
                                    f.write(str(uvalue))
                                    f.write("\n")
f.close()
