#!/usr/bin/env python3

import numpy as np
import os,sys
sys.path.insert(0,"../../")
import w2dyn.dmft.mf_phonons as phonons

densities=np.ones([12,2,12,2])
tmp_dc_full=np.zeros([12,2,12,2])*1j

a11=1.0
a22=1.0
beta=1.0
ll=1.0

def orbvalley_index(orbindex,valleyindex,elec_type):
    if elec_type=="f":
        Norb=2
    if elec_type=="c":
        Norb=4
    orbvalley = 4*(elec_type=="c") + Norb * valleyindex + orbindex
    return orbvalley

ph=phonons.phonon_1bo_params(a11,beta)

# f+f terms without density (coming from anticommutators in u_matrix), x and y
for phononterm in [[ph.VfX,ph.VfX],[ph.VfY,ph.VfY]]:
    firstop=phonons.onebody_op("f",phononterm[0])
    firstop.create_1bo()
    secondop=phonons.onebody_op("f",phononterm[1])
    secondop.create_1bo()
    twobodyterm=phonons.twobody_op(firstop,secondop,-1*ll*a22**2)
    twobodyterm.create_2bo()
    twobodyterm.normalorder()
    twobodyterm.decouple()
    dc_coeff_array,dc_op_array,dc_dens_array = twobodyterm.print_to_dc_decoupled()
    for iterm in range(len(dc_coeff_array)):
        if dc_dens_array[iterm][0]=="NODENS":
            for ispin in range(2):
                iindex = orbvalley_index(dc_op_array[iterm][0],dc_op_array[iterm][1],dc_op_array[iterm][2])
                jindex = orbvalley_index(dc_op_array[iterm][3],dc_op_array[iterm][4],dc_op_array[iterm][5])
                tmp_dc_full[iindex,ispin,jindex,ispin] += dc_coeff_array[iterm]
            
# f+c terms, x and y
for phononterm in [[ph.VfX,ph.VcX],[ph.VfY,ph.VcY]]:
    firstop=phonons.onebody_op("f",phononterm[0])
    firstop.create_1bo()
    secondop=phonons.onebody_op("c",phononterm[1])
    secondop.create_1bo()
    twobodyterm=phonons.twobody_op(firstop,secondop,-1*ll*a22)
    twobodyterm.create_2bo()
    twobodyterm.normalorder()
    twobodyterm.decouple()
    dc_coeff_array,dc_op_array,dc_dens_array = twobodyterm.print_to_dc_decoupled()
    for iterm in range(len(dc_coeff_array)):
        if dc_dens_array[iterm][0]=="NODENS":
            for ispin in range(2):
                iindex = orbvalley_index(dc_op_array[iterm][0],dc_op_array[iterm][1],dc_op_array[iterm][2])
                jindex = orbvalley_index(dc_op_array[iterm][3],dc_op_array[iterm][4],dc_op_array[iterm][5])
                tmp_dc_full[iindex,ispin,jindex,ispin] += dc_coeff_array[iterm]
        else:
            for ispin in range(2):
                iindex = orbvalley_index(dc_op_array[iterm][0],dc_op_array[iterm][1],dc_op_array[iterm][2])
                jindex = orbvalley_index(dc_op_array[iterm][3],dc_op_array[iterm][4],dc_op_array[iterm][5])
                dens_iindex = orbvalley_index(dc_dens_array[iterm][0],dc_dens_array[iterm][1],dc_dens_array[iterm][2])
                dens_jindex = orbvalley_index(dc_dens_array[iterm][3],dc_dens_array[iterm][4],dc_dens_array[iterm][5])
                tmp_dc_full[iindex,ispin,jindex,ispin] += dc_coeff_array[iterm] * densities[dens_iindex,ispin,dens_jindex,ispin]
                
# c+c terms, x and y
for phononterm in [[ph.VcX,ph.VcX],[ph.VcY,ph.VcY]]:
    firstop=phonons.onebody_op("c",phononterm[0])
    firstop.create_1bo()
    secondop=phonons.onebody_op("c",phononterm[1])
    secondop.create_1bo()
    twobodyterm=phonons.twobody_op(firstop,secondop,-1*ll)
    twobodyterm.create_2bo()
    twobodyterm.normalorder()
    twobodyterm.decouple()
    dc_coeff_array,dc_op_array,dc_dens_array = twobodyterm.print_to_dc_decoupled()
    for iterm in range(len(dc_coeff_array)):
        if dc_dens_array[iterm][0]=="NODENS":
            for ispin in range(2):
                iindex = orbvalley_index(dc_op_array[iterm][0],dc_op_array[iterm][1],dc_op_array[iterm][2])
                jindex = orbvalley_index(dc_op_array[iterm][3],dc_op_array[iterm][4],dc_op_array[iterm][5])
                tmp_dc_full[iindex,ispin,jindex,ispin] += dc_coeff_array[iterm]
        else:
            for ispin in range(2):
                iindex = orbvalley_index(dc_op_array[iterm][0],dc_op_array[iterm][1],dc_op_array[iterm][2])
                jindex = orbvalley_index(dc_op_array[iterm][3],dc_op_array[iterm][4],dc_op_array[iterm][5])
                dens_iindex = orbvalley_index(dc_dens_array[iterm][0],dc_dens_array[iterm][1],dc_dens_array[iterm][2])
                dens_jindex = orbvalley_index(dc_dens_array[iterm][3],dc_dens_array[iterm][4],dc_dens_array[iterm][5])
                tmp_dc_full[iindex,ispin,jindex,ispin] += dc_coeff_array[iterm] * densities[dens_iindex,ispin,dens_jindex,ispin]
   
tmp_dc_full_bak=np.copy(tmp_dc_full[:,0,:,0])   
   
if False:       
    tmp_dc_full=np.real(tmp_dc_full)       
    for i in range(12):
        print(str(tmp_dc_full[i,0,0,0])+"  "+str(tmp_dc_full[i,0,1,0])+"  "+str(tmp_dc_full[i,0,2,0])+"  "+str(tmp_dc_full[i,0,3,0])+"  "+str(tmp_dc_full[i,0,4,0])+"  "+str(tmp_dc_full[i,0,5,0])+"  "+str(tmp_dc_full[i,0,6,0])+"  "+str(tmp_dc_full[i,0,7,0])+"  "+str(tmp_dc_full[i,0,8,0])+"  "+str(tmp_dc_full[i,0,9,0])+"  "+str(tmp_dc_full[i,0,10,0])+"  "+str(tmp_dc_full[i,0,11,0]))

if False: 
    print(" ")
    tmp_dc_full=tmp_dc_full[:,0,:,0]-np.transpose(tmp_dc_full[:,0,:,0])      
    for i in range(12):
        print(str(tmp_dc_full[i,0])+"  "+str(tmp_dc_full[i,1])+"  "+str(tmp_dc_full[i,2])+"  "+str(tmp_dc_full[i,3])+"  "+str(tmp_dc_full[i,4])+"  "+str(tmp_dc_full[i,5])+"  "+str(tmp_dc_full[i,6])+"  "+str(tmp_dc_full[i,7])+"  "+str(tmp_dc_full[i,8])+"  "+str(tmp_dc_full[i,9])+"  "+str(tmp_dc_full[i,10])+"  "+str(tmp_dc_full[i,11]))

if True: 
    print(" ")
    tmp_dc_full=tmp_dc_full[:,0,:,0]-np.transpose(tmp_dc_full[:,0,:,0]) 
    for i in range(12):
        for j in range(12):
            if tmp_dc_full_bak[i,j] != 0:
                tmp_dc_full[i,j]=1
            else:
                tmp_dc_full[i,j]=0
    tmp_dc_full=np.real(tmp_dc_full)            
    for i in range(12):
        print(str(tmp_dc_full[i,0])+"  "+str(tmp_dc_full[i,1])+"  "+str(tmp_dc_full[i,2])+"  "+str(tmp_dc_full[i,3])+"  "+str(tmp_dc_full[i,4])+"  "+str(tmp_dc_full[i,5])+"  "+str(tmp_dc_full[i,6])+"  "+str(tmp_dc_full[i,7])+"  "+str(tmp_dc_full[i,8])+"  "+str(tmp_dc_full[i,9])+"  "+str(tmp_dc_full[i,10])+"  "+str(tmp_dc_full[i,11]))





if False:

    firstop=phonons.onebody_op("f",ph.VfX)
    firstop.create_1bo()
    secondop=phonons.onebody_op("f",ph.VfX)
    secondop.create_1bo()
    twobodyterm=phonons.twobody_op(firstop,secondop,-1*ll)
    twobodyterm.create_2bo()
    twobodyterm.normalorder()
    twobodyterm.prettyprint()
    twobodyterm.decouple()
    print(np.imag(ph.VcY))
    #twobodyterm.prettyprint_decoupled()
