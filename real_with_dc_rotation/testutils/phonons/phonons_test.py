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
            # orb1 val1 spin1 type1 orb2 val2 spin2 type2
            iindex = orbvalley_index(dc_op_array[iterm][0],dc_op_array[iterm][1],dc_op_array[iterm][3])
            ispin  = dc_op_array[iterm][2]
            jindex = orbvalley_index(dc_op_array[iterm][4],dc_op_array[iterm][5],dc_op_array[iterm][7])
            jspin  = dc_op_array[iterm][6]
            tmp_dc_full[iindex,ispin,jindex,jspin] += dc_coeff_array[iterm]
            
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
            # orb1 val1 spin1 type1 orb2 val2 spin2 type2
            iindex = orbvalley_index(dc_op_array[iterm][0],dc_op_array[iterm][1],dc_op_array[iterm][3])
            ispin  = dc_op_array[iterm][2]
            jindex = orbvalley_index(dc_op_array[iterm][4],dc_op_array[iterm][5],dc_op_array[iterm][7])
            jspin  = dc_op_array[iterm][6]
            tmp_dc_full[iindex,ispin,jindex,jspin] += dc_coeff_array[iterm]
        else:
            iindex = orbvalley_index(dc_op_array[iterm][0],dc_op_array[iterm][1],dc_op_array[iterm][3])
            ispin  = dc_op_array[iterm][2]
            jindex = orbvalley_index(dc_op_array[iterm][4],dc_op_array[iterm][5],dc_op_array[iterm][7])
            jspin  = dc_op_array[iterm][6]
            
            dens_iindex = orbvalley_index(dc_dens_array[iterm][0],dc_dens_array[iterm][1],dc_dens_array[iterm][3])
            dens_ispin  = dc_dens_array[iterm][2]
            dens_jindex = orbvalley_index(dc_dens_array[iterm][4],dc_dens_array[iterm][5],dc_dens_array[iterm][7])
            dens_jspin  = dc_dens_array[iterm][6]
            tmp_dc_full[iindex,ispin,jindex,jspin] += dc_coeff_array[iterm] * densities[dens_iindex,dens_ispin,dens_jindex,dens_jspin]
                
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
            iindex = orbvalley_index(dc_op_array[iterm][0],dc_op_array[iterm][1],dc_op_array[iterm][3])
            ispin  = dc_op_array[iterm][2]
            jindex = orbvalley_index(dc_op_array[iterm][4],dc_op_array[iterm][5],dc_op_array[iterm][7])
            jspin  = dc_op_array[iterm][6]
            tmp_dc_full[iindex,ispin,jindex,jspin] += dc_coeff_array[iterm]
        else:
            iindex = orbvalley_index(dc_op_array[iterm][0],dc_op_array[iterm][1],dc_op_array[iterm][3])
            ispin  = dc_op_array[iterm][2]
            jindex = orbvalley_index(dc_op_array[iterm][4],dc_op_array[iterm][5],dc_op_array[iterm][7])
            jspin  = dc_op_array[iterm][6]
            
            dens_iindex = orbvalley_index(dc_dens_array[iterm][0],dc_dens_array[iterm][1],dc_dens_array[iterm][3])
            dens_ispin  = dc_op_array[iterm][2]
            dens_jindex = orbvalley_index(dc_dens_array[iterm][4],dc_dens_array[iterm][5],dc_dens_array[iterm][7])
            dens_jspin  = dc_op_array[iterm][6]
            
            tmp_dc_full[iindex,ispin,jindex,jspin] += dc_coeff_array[iterm] * densities[dens_iindex,dens_ispin,dens_jindex,dens_jspin]
   
tmp_dc_full_bak=np.copy(tmp_dc_full[:,0,:,0])   
   
print("iband jband ispin jspin element")

for iband in range(12):
    for jband in range(12):
        for ispin in range(2):
            for jspin in range(2):
                if ispin == 0:
                    ispin_="↑"
                else:
                    ispin_="↓"
                if jspin == 0:
                    jspin_="↑"
                else:
                    jspin_="↓"
                print(iband,jband,ispin_,jspin_,tmp_dc_full[iband,ispin,jband,jspin])




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
