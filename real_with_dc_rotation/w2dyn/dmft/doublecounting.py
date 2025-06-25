"""Classes for double-counting correction"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
import w2dyn.dmft.orbspin as orbspin
import w2dyn.dmft.atoms as atoms
import sys
import w2dyn.dmft.mf_phonons as phonons
import pickle

# import w2dyn.dmft.lattice as lattice
# import w2dyn.auxiliaries.transform as tf

class DoubleCounting:
    """Base class for Double counting correction"""
    def __init__(self):
        self.dc_value = None
    def get(self, impurity_results=None):
        """Gets the value of the double counting correction.

        This method is called *after* the solution to the impurity problems in
        iteration `iter_no` (zero-based) and is expected to yield a correction
        to the self-energy.
        """
        return self.dc_value
    def set(self,dc_value):
        """Sets the value of the double counting correction internally."""
        self.dc_value = dc_value

class Zero(DoubleCounting):
    """Trivial double counting that always yields zero"""
    def __init__(self, nbands, nspins):
        self.self_cons = False
        self.dc_value = np.zeros((nbands, nspins, nbands, nspins))

class UserSupplied(DoubleCounting):
    """Double counting as supplied by user"""
    def __init__(self, dc_user_supplied):
        DoubleCounting.__init__(self)
        self.self_cons = False
        self.dc_value = dc_user_supplied

class ShellAveraged(DoubleCounting):
    r"""Base class for l-shell-averaged double counting.

    The most common double-counting schemes work with fillings and interactions
    that are averaged over the orbitals of l-subshells of a single atom or
    ligand, e.g. averaged over 5 d-bands or 3 p-bands [1].

    Therefore, this class provides some subshell-averaged quantities:

      - `nshells`:   total number of shells over every atom
      - `slices`:    list of slices for selection of the different shells
      - `dens_sum`:  total filling of a shell in LDA
      - `dens_mean`: average filling of an orbital in a shell in LDA
      - `ubar`:      the mean value of the interaction between shells

    Averaged interaction
    --------------------
    The code computes the mean U and J for every sub-shell. According to Czyzyk
    and Sawatzky [1] and Parragh [2], these are given as:

        Ubar        = 1/N^2    \sum_{mn} U_{mnmn}
        Ubar - Jbar = 1/N(N-1) \sum_{mn} (U_{mnmn} - U_{mnnm})

    (the second line is divided by N(N-1) because the term for m=n is zero by
    construction). Comparing this with the formula for the density-density part
    of the U matrix, which is spin-dependent (s,t):

        U_{ms,nt} = U_{mnmn} - U_{mnnm} \delta_{st}

    we see that for one shell we can write this as (spin-dependent) mean U:

        ubar' = [ Ubar-Jbar Ubar      ] = sum U_{ms,nt} (/) [ N(N-1) N^2    ]
                [ Ubar      Ubar-Jbar ]    mn               [ N^2    N(N-1) ]

    where the bracketed slash `(/)` denotes element-wise division, which is
    omitted for N = 1 in order to avoid 0/0).

    [1]: Czyzyk and Sawatzky, Phys. Rev. B. 49, 14211 (1994)
    [2]: Nicolaus Parragh, Ph.D. thesis
    """
    def __init__(self, densities, atom_list, u_full):
        self.self_cons = False
        if atom_list is None:
            ndorb = u_full.shape[0]
            atom_list = [atoms.Atom(ndorb, 0, ndorb, None)]

        self.atom_list = atom_list

        # split problem into sub-shells.
        self.slices = []
        self.sizes = []
        for atom in self.atom_list:
            self.slices.append(atom.dslice)
            self.sizes.append(atom.nd)
            if atom.np > 0:
                self.slices.extend(atom.ligslices)
                self.sizes.extend([atom.npplig] * atom.nlig)

        self.nshells = len(self.slices)
        self.sizes = np.array(self.sizes)
        self.norbitals, self.nspins = u_full.shape[:2]

        # compute mean density for the different manifolds
        # spin-dependent LDA occupations averaged over the orbitals in the
        # subspace.
        self.dens_mean = np.array([densities[sl].mean() for sl in self.slices])
        self.dens_sum = np.array([densities[sl].sum() for sl in self.slices])
 
        # Compute the mean U and J for every sub-shell.
        ubar = np.zeros((self.nshells, self.nshells))
        jbar = np.zeros((self.nshells, self.nshells))
        for ishell1,shell1 in enumerate(self.slices):
            for ishell2,shell2 in enumerate(self.slices):
                ubar[ishell1, ishell2] = u_full[shell1,0,shell2,1].mean(1).mean(0)
                jbar[ishell1, ishell2] = u_full[shell1,0,shell2,0].sum(1).sum(0)

        # jbar normalization avoiding 0 division.
        for ipos, isize in enumerate(self.sizes):
            for jpos, jsize in enumerate(self.sizes):
                jbar[ipos, jpos] /= isize
                if ipos == jpos:
                    if isize == 1:
                        jbar[ipos, jpos] = 0
                    else:
                        jbar[ipos, jpos] /= isize - 1
                else:
                    jbar[ipos, jpos] /= jsize

        # We actually computed Ubar - Jbar as Jbar
        jbar = ubar - jbar
        self.ubar = ubar
        self.jbar = jbar

        self.dc_shell = None
        self.dc_value = None

    def from_shell(self, dc_shell):
        dc_diag = np.zeros((self.norbitals, self.nspins))
        for isl, sl in enumerate(self.slices):
            dc_diag[sl] = dc_shell[isl, np.newaxis]

        # note the minus sign!
        self.dc_shell = dc_shell
        self.dc_value = -orbspin.promote_diagonal(dc_diag)

class FullyLocalisedLimit(ShellAveraged):
    r"""Double counting in the fully localised limit.

    The fully localised limit double counting (also known as Anisimov's double
    counting) assumes the correlation included in the non-interacting problem
    to represent the atomic limit of the lattice problem:

        DC_d = Ubar_{dd} (N_d - 1/2) + Ubar_{dp} N_p
        DC_p = Ubar_{pp} (N_p - 1/2) + Ubar_{dp} N_d

    When the indices I, J run over "subshells" corresponding either to d-atom
    or a p-atom/ligand, one can use the definition of `ubar` (see documentation
    of `ShellAveraged`) to recast the formula into:

        DC_{Is} = \sum_{Jt} Ubar_{Is,Jt} N_{Jt} - 1/2 \sum_t Ubar_{Is,It}

    Note that the subtraction of 1/2 only occurs for sub-shell local terms.
    """
    def __init__(self, *args):
        ShellAveraged.__init__(self, *args)

        dc_shell = self.ubar.dot(self.dens_sum) - self.jbar.dot(self.dens_sum/2)
        for I in range(self.nshells):
            dc_shell[I] -= self.ubar[I, I] * 0.5
            dc_shell[I] += self.jbar[I, I] * 0.5

        self.from_shell(dc_shell)

class AroundMeanField(ShellAveraged):
    r"""Double counting in the fully localised limit.

    The around mean field double counting can be seen as a correction to the
    fully localised limit for systems away from half filling:

        DC_d = Ubar_{dd} (N_d - n_d) + Ubar_{dp} N_p
        DC_p = Ubar_{pp} (N_p - n_d) + Ubar_{dp} N_d

    When the indices I, J run over "subshells" corresponding either to d-atom
    or a p-atom/ligand, one can use the definition of `ubar` (see documentation
    of `ShellAveraged`) to recast the formula into:

        DC_{Is} = \sum_{Jt} Ubar_{Is,Jt} N_{Jt} - \sum_t Ubar_{Is,It} nbar_{It}

    Note that the subtraction of `nbar` only occurs for sub-shell local terms.
    """
    def __init__(self, *args):
        ShellAveraged.__init__(self, *args)
        
        dc_shell = self.ubar.dot(self.dens_sum) - self.jbar.dot(self.dens_sum/2)
        for I in range(self.nshells):
            dc_shell[I] -= self.ubar[I, I] * self.dens_mean[I]
            dc_shell[I] += self.jbar[I, I] * self.dens_mean[I]

        self.from_shell(dc_shell)

class SelfConsistent(DoubleCounting):
    def __init__(self, ineq_list, densities, atom_list, u_full):
        self.ineq_list = ineq_list
        self.shell_av = {}
        self.shell_av['densities'] = densities
        self.shell_av['atom_list'] = atom_list
        self.shell_av['u_full'] = u_full
        if atom_list is None:
            ndorb = u_full.shape[0]
            atom_list = [atoms.Atom(ndorb, 0, ndorb, None)]

        self.atom_list = atom_list
        # split problem into sub-shells.
        # we only want the slices for the d manifold of each atom
        self.d_slices = []
        for atom in self.atom_list:
            self.d_slices.append(atom.dslice)

        self.norbitals, self.nspins = u_full.shape[:2]

        # compute mean density for the different d manifolds       
        # spin-dependent LDA occupations averaged over the orbitals in the
        # subspace.
        self.d_dens_sum = np.array([densities[sl].sum() for sl in self.d_slices])
        self.dc_value = np.zeros((self.norbitals, self.nspins, self.norbitals, self.nspins))
        self.self_cons = True

class Wterm(SelfConsistent):
    def __init__(self, u2, w, v, j, ll, a22, a11, beta, shifts, rotationmatrix, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.sav = ShellAveraged(self.shell_av['densities'],
                                 self.shell_av['atom_list'],
                                 self.shell_av['u_full'])
        self.u2 = u2
        self.w = w
        self.v = v
        self.jj = j
        self.shifts = shifts
        self.ll = ll
        self.a22 = a22
        self.a11 = a11
        self.beta = beta
        if rotationmatrix is not None:
            try:
                with open(rotationmatrix, "rb") as f:
                    rotmat = pickle.load(f)
                    self.rotationmatrix = np.kron(np.eye(2), rotmat)
            except:
                print("No valid rotation matrix provided!")
                sys.exit()
        else:
            self.rotationmatrix = np.eye(24)

    def rotate_to_bandbasis(self, mat):
        mat_t = mat.transpose(1, 0, 3, 2)
        mat_t = mat_t.reshape(24,24)

        mat_t = np.dot(np.linalg.inv(self.rotationmatrix),np.matmul(mat_t,self.rotationmatrix))

        mat_t = mat_t.reshape(2, 12, 2, 12)
        mat_t = mat_t.transpose(1, 0, 3, 2)
        return mat_t
        
    def rotate_to_orbbasis(self, mat):
        mat_t = mat.transpose(1, 0, 3, 2)
        mat_t = mat_t.reshape(24,24)

        mat_t = np.dot(self.rotationmatrix,np.matmul(mat_t,np.linalg.inv(self.rotationmatrix)))

        mat_t = mat_t.reshape(2, 12, 2, 12)
        mat_t = mat_t.transpose(1, 0, 3, 2)
        return mat_t
    
    def orbvalley_index_j(self,orbindex,valleyindex,elec_type):
        if elec_type=="f":
            orbvalley= 2 * valleyindex + orbindex
        if elec_type=="c":
            orbvalley= 4 + 4 * valleyindex + orbindex + 2
        return orbvalley
        
    def orbvalley_index(self,orbindex,valleyindex,elec_type):
        if elec_type=="f":
            Norb=2
        if elec_type=="c":
            Norb=4
        orbvalley = 4*(elec_type=="c") + Norb * valleyindex + orbindex
        return orbvalley
        

    def get(self, siws=None, smoms=None, giws=None, occs=None, densities=None):
        if densities is None: #This just sets up a diagonal density matrix if none is provided. Shouldn't be needed
            #sys.exit("Fix, we need a proper density matrix")
            densities = np.zeros((self.norbitals, self.nspins, self.norbitals, self.nspins))
            for iorb in range(12):
                for ispin in range(2):
                    densities[iorb,ispin,iorb,ispin] = self.shell_av['densities'][iorb,ispin]

        #rotate to orbital basis: this DC works in the original basis!

        densities = self.rotate_to_orbbasis(densities)

        #densities = np.real(densities)
        diagonal_densities = orbspin.extract_diagonal(densities)
        shifts = np.array(self.shifts)
        tmp_dc_full = np.zeros((self.norbitals, self.nspins, self.norbitals, self.nspins))*1j
        
        #################
        # Hartree shift #
        #################
        
        for ispin in range(2):
            for iorb in range(0,4):
                tmp_dc_full[iorb,ispin,iorb,ispin] += shifts[0]
            for iorb in range(4,12):
                tmp_dc_full[iorb,ispin,iorb,ispin] += shifts[1]
                
        #################        
        # U2 terms #
        #################
                    
        for ispin in range(2):
            for iorb_f in range(4):
                tmp_dc_full[iorb_f,ispin,iorb_f,ispin] +=  6.0 * self.u2 * np.sum(diagonal_densities[slice(0,4,None), :] - np.full_like(diagonal_densities, 0.5)[slice(0,4,None), :])
        
        #################        
        # W and V terms #
        #################
        
        if False: #this is the current convention, the one that is not exploding in the loop. This is copypasted from the original "alex" branch. We will eventually need to move away from this
            shifts = np.array(self.shifts)
            diag_dc = []
            for atom in self.atom_list:
                nu_corr = np.sum(diagonal_densities[atom.dslice, :] - np.full_like(diagonal_densities, 0.5)[atom.dslice, :])
                nu_uncorr = sum(np.sum(diagonal_densities[sl, :] - np.full_like(diagonal_densities, 0.5)[sl, :]) for sl in atom.ligslices)
                atslices = [atom.dslice, *atom.ligslices]
                atshifts = shifts[:len(atslices)]
                for i, shift in enumerate(atshifts):
                    if i == 0:
                        diag_dc.append(self.w[0] * (nu_uncorr - nu_corr) - self.v * nu_uncorr)
                    else:
                        diag_dc.append(0.0)
                        
            #double minus trick
            self.sav.from_shell(np.array(diag_dc))
            tmp_dc_full += -self.sav.dc_value
        
        else: #this is the most general expression of the mf-decoupled term, but currently it explodes in the loop. It's split between W and V terms
            #W term 1 (S284)
            for ispin in range(2):
                for ival_c in range(2):
                    for iorb_c in range(4):
                        index_c = self.orbvalley_index(iorb_c,ival_c,"c")
                        tmp_dc_full[index_c,ispin,index_c,ispin] += self.w[iorb_c] * np.sum(diagonal_densities[slice(0,4,None), :] - np.full_like(diagonal_densities, 0.5)[slice(0,4,None), :])

            #W term 2 (S284)
            for ispin in range(2):
                for iorb_f in range(4):
                    for val_c in range(2):
                        for iorb_c in range(4):
                            index_c = self.orbvalley_index(iorb_c,ival_c,"c")
                            tmp_dc_full[iorb_f,ispin,iorb_f,ispin] += self.w[iorb_c] * np.sum(diagonal_densities[index_c, :] - np.full_like(diagonal_densities, 0.5)[index_c, :])

            #W term 3 (S284) This is the Fock term, which for now we don't consider
            if True:
                for ispin in range(2):
                    for jspin in range(2):
                        for iorb_f in range(4):
                            for ival_c in range(2):
                                for iorb_c in range(4):
                                    index_c = self.orbvalley_index(iorb_c,ival_c,"c")
                                    tmp_dc_full[iorb_f,ispin,index_c,jspin] += -1.0 * self.w[iorb_c] * densities[index_c,jspin,iorb_f,ispin]
                                    tmp_dc_full[index_c,jspin,iorb_f,ispin] += -1.0 * self.w[iorb_c] * densities[index_c,jspin,iorb_f,ispin]

        
            #V term (S292)
            for ispin in range(2):
                for iorb_c in range(4,12):
                    tmp_dc_full[iorb_c,ispin,iorb_c,ispin] += self.v * np.sum(diagonal_densities[slice(4,12,None), :] - np.full_like(diagonal_densities, 0.5)[slice(4,12,None), :])   
                    
        ##########
        # J term # 
        ##########
        
        #(S288) of the PRL
        
        if False:
            for iorb in range(2):
                for ivalley in range(2):
                    f_index_1 = self.orbvalley_index_j(iorb,ivalley,"f")
                    c_index_1 = self.orbvalley_index_j(iorb,ivalley,"c")
                    f_index_2 = self.orbvalley_index_j(not(iorb),not(ivalley),"f")
                    c_index_2 = self.orbvalley_index_j(not(iorb),not(ivalley),"c")
                    for s1 in range(2):
                        for s2 in range(2):
                            #term 1 (S288)
                            tmp_dc_full[f_index_1,s1,f_index_1,s2] += -1.0 * self.jj * (densities[c_index_1,s2,c_index_1,s1] - 0.5*(s1 == s2))
                            #term 2 (S288)
                            tmp_dc_full[c_index_1,s2,c_index_1,s1] += -1.0 * self.jj * (densities[f_index_1,s1,f_index_1,s2] - 0.5*(s1 == s2))
                            #term 3 (S288)
                            tmp_dc_full[c_index_2,s2,c_index_1,s1] += self.jj * densities[f_index_1,s1,f_index_2,s2]
                            #term 4 (S288)
                            tmp_dc_full[f_index_1,s1,f_index_2,s2] += self.jj * densities[c_index_2,s2,c_index_1,s1]
                            #term 5 (S288)
                            tmp_dc_full[f_index_1,s1,c_index_1,s1] += self.jj * densities[c_index_1,s2,f_index_1,s2]
                            tmp_dc_full[c_index_1,s1,f_index_1,s1] += self.jj * densities[c_index_1,s2,f_index_1,s2]
                            #term 6 (S288)
                            tmp_dc_full[f_index_1,s1,c_index_1,s1] += -1.0 * self.jj * densities[c_index_2,s2,f_index_2,s2]
                            tmp_dc_full[c_index_1,s1,f_index_1,s1] += -1.0 * self.jj * densities[c_index_2,s2,f_index_2,s2]    
        else:
            #Here we try in the way we wrote it in the PRX supplementary or equivalently Haoyu's PRL
            Jmatrix = np.zeros_like(tmp_dc_full)
            
            #First term of S16
            for iorb in range(2):
                for jorb in range(2):
                    for ival in range(2):
                        for jval in range(2):
                            for ispin in range(2):
                                for jspin in range(2):
                                    ival_sign = (-1)**ival #valley 1 is plus valley, valley 2 is minus valley
                                    jval_sign = (-1)**jval
                                    prefactor = ( ival_sign*jval_sign + (-1)**(iorb+1+jorb+1) ) #+1 because python indexes start from 0, but it makes no difference
                                    
                                    i_f = self.orbvalley_index_j(iorb,ival,"f")
                                    j_f = self.orbvalley_index_j(jorb,jval,"f")
                                    i_c = self.orbvalley_index_j(iorb,ival,"c")
                                    j_c = self.orbvalley_index_j(jorb,jval,"c")
                                                                        
                                    Jmatrix[i_f,ispin,j_f,jspin] += 0.5 * self.jj * prefactor * (densities[j_c,jspin,i_c,ispin] - 0.5 * ((ispin == jspin) and (i_c == j_c)) )

            #Second term of S16
            for iorb in range(2):
                for jorb in range(2):
                    for ival in range(2):
                        for jval in range(2):
                            for ispin in range(2):
                                for jspin in range(2):
                                    ival_sign = (-1)**ival
                                    jval_sign = (-1)**jval
                                    prefactor = ( ival_sign*jval_sign + (-1)**(iorb+1+jorb+1) )
                                    
                                    i_f = self.orbvalley_index_j(iorb,ival,"f")
                                    j_f = self.orbvalley_index_j(jorb,jval,"f")
                                    i_c = self.orbvalley_index_j(iorb,ival,"c")
                                    j_c = self.orbvalley_index_j(jorb,jval,"c")
                                    
                                    Jmatrix[i_c,ispin,j_c,jspin] += 0.5 * self.jj * prefactor * (densities[j_f,jspin,i_f,ispin] - 0.5 * ((ispin==jspin) and (i_f == j_f)) )

            #########################################
            ## THE FOCK TERM                       ##
            #########################################
            
            without_etas = False
            S290 = True
            
            if not S290: 
                #Define V3 and V4, (S17) and (S18):
                V3 = 0.0
                for iorb in range(2):
                    for ival in range(2):
                        for ispin in range(2):
                            signval = (-1) ** ival
                            deltaval = signval * (-1) ** (iorb + 2) #not +1, because python index starts from 0
                            i_f = self.orbvalley_index_j(iorb,ival,"f")
                            i_c = self.orbvalley_index_j(iorb,ival,"c")
                            if without_etas:
                                V3 += densities[i_f,ispin,i_c,ispin] * (1 == deltaval)
                            else:
                                V3 += densities[i_f,ispin,i_c,ispin] * (1 == deltaval) * signval

                V4 = 0.0
                for iorb in range(2):
                    for ival in range(2):
                        for ispin in range(2):
                            signval = (-1) ** ival
                            deltaval = signval * (-1) ** (iorb + 2) #not +1, because python index starts from 0
                            i_f = self.orbvalley_index_j(iorb,ival,"f")
                            i_c = self.orbvalley_index_j(iorb,ival,"c")
                            if without_etas:
                                V4 += densities[i_f,ispin,i_c,ispin] * (-1 == deltaval)
                            else:
                                V4 += densities[i_f,ispin,i_c,ispin] * (-1 == deltaval) * signval

                #(S19)
                for iorb in range(2):
                    for ival in range(2):
                        for ispin in range(2):
                            signval = (-1) ** ival
                            deltaval = signval * (-1) ** (iorb + 2) #not +1, because python index starts from 0
                            i_f = self.orbvalley_index_j(iorb,ival,"f")
                            i_c = self.orbvalley_index_j(iorb,ival,"c")
                            if without_etas:
                                #firt term
                                Jmatrix[i_f,ispin,i_c,ispin] += self.jj * np.conj(V3) * (1 == deltaval) 
                                #second term
                                Jmatrix[i_f,ispin,i_c,ispin] += self.jj * np.conj(V4) * (-1 == deltaval)
                                # h.c. of first term
                                Jmatrix[i_c,ispin,i_f,ispin] += self.jj * V3 * (1 == deltaval)
                                #second term
                                Jmatrix[i_c,ispin,i_f,ispin] += self.jj * V4 * (-1 == deltaval)
                            else:
                                #firt term
                                Jmatrix[i_f,ispin,i_c,ispin] += self.jj * np.conj(V3) * (1 == deltaval) * signval 
                                #second term
                                Jmatrix[i_f,ispin,i_c,ispin] += self.jj * np.conj(V4) * (-1 == deltaval) *signval
                                # h.c. of first term
                                Jmatrix[i_c,ispin,i_f,ispin] += self.jj * V3 * (1 == deltaval) *signval
                                #second term
                                Jmatrix[i_c,ispin,i_f,ispin] += self.jj * V4 * (-1 == deltaval) *signval
                if without_etas:
                    print("Fock term from V3 and V4, NO eta")
                else:
                    print("Fock term from V3 and V4, YES eta")
            else:
                #Fock term from (S290) Song-Bernevig
                for iorb in range(2):
                    for jorb in range(2):
                        for ival in range(2):
                            for jval in range(2):
                                for ispin in range(2):
                                    for jspin in range(2):
                                        sign_ival = (-1) ** ival
                                        sign_jval = (-1) ** jval
                                        deltaval = sign_ival * sign_jval  + (-1)**(iorb + 1 + jorb + 1)
                                        i_f = self.orbvalley_index_j(iorb,ival,"f")
                                        i_c = self.orbvalley_index_j(iorb,ival,"c")
                                        j_f = self.orbvalley_index_j(jorb,jval,"f")
                                        j_c = self.orbvalley_index_j(jorb,jval,"c")
                                        Jmatrix[i_f,ispin,i_c,ispin] += 0.5 * self.jj * deltaval * densities[j_c,jspin,j_f,jspin]
                                        Jmatrix[i_c,ispin,i_f,ispin] += 0.5 * self.jj * deltaval * densities[j_c,jspin,j_f,jspin]
                print("Fock term from S290")
            
            #Sum the J matrix to the DC
            print("Summing J term")
            tmp_dc_full += Jmatrix
            
        ################################
        # Phonons mean-field decoupled #
        ################################
        
        ph=phonons.phonon_1bo_params(self.a11,self.beta)
        
        # f+f terms without density (coming from anticommutators in u_matrix), x and y
        for phononterm in [[ph.VfX,ph.VfX],[ph.VfY,ph.VfY]]:
            firstop=phonons.onebody_op("f",phononterm[0])
            firstop.create_1bo()
            secondop=phonons.onebody_op("f",phononterm[1])
            secondop.create_1bo()
            twobodyterm=phonons.twobody_op(firstop,secondop,-1*self.ll*self.a22**2)
            twobodyterm.create_2bo()
            twobodyterm.normalorder()
            twobodyterm.decouple()
            dc_coeff_array,dc_op_array,dc_dens_array = twobodyterm.print_to_dc_decoupled()
            for iterm in range(len(dc_coeff_array)):
                if dc_dens_array[iterm][0]=="NODENS":
                    for ispin in range(2):
                        iindex = self.orbvalley_index(dc_op_array[iterm][0],dc_op_array[iterm][1],dc_op_array[iterm][2])
                        jindex = self.orbvalley_index(dc_op_array[iterm][3],dc_op_array[iterm][4],dc_op_array[iterm][5])
                        tmp_dc_full[iindex,ispin,jindex,ispin] += dc_coeff_array[iterm]
                    
        # f+c terms, x and y
        for phononterm in [[ph.VfX,ph.VcX],[ph.VfY,ph.VcY]]:
            firstop=phonons.onebody_op("f",phononterm[0])
            firstop.create_1bo()
            secondop=phonons.onebody_op("c",phononterm[1])
            secondop.create_1bo()
            twobodyterm=phonons.twobody_op(firstop,secondop,-1*self.ll*self.a22)
            twobodyterm.create_2bo()
            twobodyterm.normalorder()
            twobodyterm.decouple()
            dc_coeff_array,dc_op_array,dc_dens_array = twobodyterm.print_to_dc_decoupled()
            for iterm in range(len(dc_coeff_array)):
                if dc_dens_array[iterm][0]=="NODENS":
                    for ispin in range(2):
                        iindex = self.orbvalley_index(dc_op_array[iterm][0],dc_op_array[iterm][1],dc_op_array[iterm][2])
                        jindex = self.orbvalley_index(dc_op_array[iterm][3],dc_op_array[iterm][4],dc_op_array[iterm][5])
                        tmp_dc_full[iindex,ispin,jindex,ispin] += dc_coeff_array[iterm]
                else:
                    for ispin in range(2):
                        iindex = self.orbvalley_index(dc_op_array[iterm][0],dc_op_array[iterm][1],dc_op_array[iterm][2])
                        jindex = self.orbvalley_index(dc_op_array[iterm][3],dc_op_array[iterm][4],dc_op_array[iterm][5])
                        dens_iindex = self.orbvalley_index(dc_dens_array[iterm][0],dc_dens_array[iterm][1],dc_dens_array[iterm][2])
                        dens_jindex = self.orbvalley_index(dc_dens_array[iterm][3],dc_dens_array[iterm][4],dc_dens_array[iterm][5])
                        tmp_dc_full[iindex,ispin,jindex,ispin] += dc_coeff_array[iterm] * densities[dens_iindex,ispin,dens_jindex,ispin]
                        
        # c+c terms, x and y
        for phononterm in [[ph.VcX,ph.VcX],[ph.VcY,ph.VcY]]:
            firstop=phonons.onebody_op("c",phononterm[0])
            firstop.create_1bo()
            secondop=phonons.onebody_op("c",phononterm[1])
            secondop.create_1bo()
            twobodyterm=phonons.twobody_op(firstop,secondop,-1*self.ll)
            twobodyterm.create_2bo()
            twobodyterm.normalorder()
            twobodyterm.decouple()
            dc_coeff_array,dc_op_array,dc_dens_array = twobodyterm.print_to_dc_decoupled()
            for iterm in range(len(dc_coeff_array)):
                if dc_dens_array[iterm][0]=="NODENS":
                    for ispin in range(2):
                        iindex = self.orbvalley_index(dc_op_array[iterm][0],dc_op_array[iterm][1],dc_op_array[iterm][2])
                        jindex = self.orbvalley_index(dc_op_array[iterm][3],dc_op_array[iterm][4],dc_op_array[iterm][5])
                        tmp_dc_full[iindex,ispin,jindex,ispin] += dc_coeff_array[iterm]
                else:
                    for ispin in range(2):
                        iindex = self.orbvalley_index(dc_op_array[iterm][0],dc_op_array[iterm][1],dc_op_array[iterm][2])
                        jindex = self.orbvalley_index(dc_op_array[iterm][3],dc_op_array[iterm][4],dc_op_array[iterm][5])
                        dens_iindex = self.orbvalley_index(dc_dens_array[iterm][0],dc_dens_array[iterm][1],dc_dens_array[iterm][2])
                        dens_jindex = self.orbvalley_index(dc_dens_array[iterm][3],dc_dens_array[iterm][4],dc_dens_array[iterm][5])
                        tmp_dc_full[iindex,ispin,jindex,ispin] += dc_coeff_array[iterm] * densities[dens_iindex,ispin,dens_jindex,ispin]
        
        
        ######################
        # Sum all and return #
        ######################
        
        #rotate the DC to the band basis! this is where we do the calculations

        tmp_dc_full = self.rotate_to_bandbasis(tmp_dc_full)
        
        #print("      ")
        #for row in tmp_dc_full[:,0,:,0]:
        #    print(" ".join(f"{elem.real: >7.4f}{'+' if elem.imag >= 0 else '-'}{abs(elem.imag):<7.4f}i" for elem in row))
        #print("      ")

        #if np.any(np.iscomplex(tmp_dc_full)):
        #    sys.exit("DC has finite imaginary part!")
        #else:
        #    tmp_dc_full=np.real(tmp_dc_full)
        
        #set numerically small elements to 0
        tmp_dc_full[abs(tmp_dc_full) < 1e-10] = 0
        
        self.dc_value = tmp_dc_full

        #here there was a self.sav.from_shell: this promoted the dc array to the diagonal of a matrix, with a minus sign. Then there was another minus sign
        #when assigning the return value of this function. We create directly the matrix, so we don't need to add the minus sign to compensate the one coming
        #from sav.from_shell
        #print("DIAGONAL DC = ",orbspin.extract_diagonal(self.dc_value)[:,0])
        return self.dc_value

class SigZeroBar(SelfConsistent):
    def get(self, siws=None, smoms=None, giws = None, occs = None):
        self.dc_value[...] = 0
        if siws is None:
            return self.dc_value
        for ineq, siw in zip(self.ineq_list, siws):
            niw = siw.shape[0]
            szero = siw[niw//2] # lowest Matsubara point of self-energy, 
                                # since ReS not very w-dependent. might not 
                                # be good enough for high T. Should do 
                                # (simple) extrapolation instead... 

            szerobar = szero.trace(0,1,3).trace().real
            szerobar /= ineq.nd * self.nspins
            szerobar = szerobar * np.eye(ineq.nd * self.nspins).reshape(
                                    ineq.nd, self.nspins, ineq.nd, self.nspins)
            ineq.d_setpart(self.dc_value, -szerobar)
        return self.dc_value

class SigInfBar(SelfConsistent):
    def get(self, siws=None, smoms=None, giws = None, occs = None):
        self.dc_value[...] = 0
        if smoms is None:
            return self.dc_value
            
        for ineq, smom in zip(self.ineq_list, smoms):
            sinf = smom[0]     # 0-th moment
            sinfbar = sinf.trace(0,1,3).trace()
            sinfbar /= ineq.nd * self.nspins
            sinfbar = sinfbar * np.eye(ineq.nd * self.nspins).reshape(
                                    ineq.nd, self.nspins, ineq.nd, self.nspins)
            ineq.d_setpart(self.dc_value, -sinfbar)
        return self.dc_value

class Trace(SelfConsistent):
    def get(self, siws=None, smoms=None, giws = None, occs = None):

        if giws is None:
            # initializing with FLL double counting... better than 0
            dc_fll = FullyLocalisedLimit(self.shell_av['densities'],
                                         self.shell_av['atom_list'],
                                         self.shell_av['u_full'])
            self.dc_value = dc_fll.get()
            return self.dc_value

        for ineq, giw, occ, occ0 in zip(self.ineq_list, giws, occs, self.d_dens_sum):
            # DOS as "approximated" by first Matsubara frequency
            iw0 = giw.shape[0]//2
            dos_fermi_level = -1/np.pi * orbspin.trace(giw[iw0].imag)

            # impurity density from double occupations.
            # TODO: this is tied to diagonal hybr
            dens_impurity = orbspin.promote_diagonal(
                        orbspin.extract_diagonal(occ))
            dens_impurity = np.array(dens_impurity).astype(np.double)
            dens_impurity = np.sum(dens_impurity)

            # we do not use the non interacting density from Weiss field
            # but that from the DFT Hamiltonian.
            # that is a lot more stable

            # impurity non-interacting density (from Weiss field)
            #iw = tf.matfreq(beta, 'fermi', imp_result.problem.niwf)
            #eye = np.eye(imp_result.problem.norbitals * imp_result.problem.nspins) \
            #        .reshape(imp_result.problem.norbitals, imp_result.problem.nspins,
            #                 imp_result.problem.norbitals, imp_result.problem.nspins)
            #g0iw = orbspin.invert(imp_result.problem.g0iwinv)
            #g0iw_model = orbspin.invert(1j * iw * eye + imp_result.problem.muimp)
            #dens_model = lattice.fermi(-imp_result.problem.muimp) # FIXME offdiag
            #dens0_impurity = (g0iw - g0iw_model).real.sum(0) / imp_result.problem.beta \
            #                 + dens_model
            if ineq.occ_dc_number is not None:
                dens0_impurity = ineq.occ_dc_number
            else:
                dens0_impurity = occ0
            
            # If the DOS at the Fermi level is large, we do not want to adjust
            # too much for stability reasons (fudging)
            if dos_fermi_level > 1.:
                delta_dc = (dens0_impurity - dens_impurity)/dos_fermi_level
            else:
                delta_dc = dens0_impurity - dens_impurity

            # Fudging factor
            delta_dc *= 0.2

            delta_dc = delta_dc * np.eye(ineq.nd * self.nspins).reshape(
                                    ineq.nd, self.nspins, ineq.nd, self.nspins)
            ineq.d_setpart(self.dc_value, ineq.d_downfold(self.dc_value) - delta_dc )
        return self.dc_value



# Experimental non-benchmarked DC Method!
class Fixed_dp_Distance():
  def get(self, dc, siw_dd, smom_dd, iwf, natoms, fixed_orbitals=None):
    print("    *** FIXING t2g-p ORBITAL DISTANCE ***")
    
    if fixed_orbitals == None:
      print('    ==> Fixing standard orbitals: 1, 3, 5')
      fixed_orbitals = np.array([0,2,4])
    
    # SINGLE ORBITAL SELECTION
    if np.isscalar(fixed_orbitals):
      print('    Fixing to orbital:', fixed_orbitals)
      no_iwn = np.asarray(siw_dd).shape[1]
      siw_t2g_real_extrapolated = 0
      #print '    ==> Using EXTRAPOLATION!'
      #siw_t2g_real_extrapolated = np.real(siw_dd[0][int(no_iwn/2),fixed_orbitals,1,fixed_orbitals,1]) \
            #- ( np.real(siw_dd[0][int(no_iwn/2)+1,fixed_orbitals,1,fixed_orbitals,1]) - np.real(siw_dd[0][int(no_iwn/2),fixed_orbitals,1,fixed_orbitals,1]) ) \
            #/ ( iwf[int(no_iwn/2)+1] - iwf[int(no_iwn/2)]) * iwf[int(no_iwn/2)]
      #siw_t2g_real_avg = siw_t2g_real_extrapolated

      # Alternative p+d DC Method
      print('    ==> Using first Sigma value!')
      siw_t2g_real_avg = np.real(siw_dd[0][no_iwn//2,
                                           fixed_orbitals, 1,
                                           fixed_orbitals, 1])

      ## Alternative p+d DC Method
      #print '    ==> Using Sigma oo value!'
      #siw_t2g_real_avg = np.real(siw_dd[0][int(0),fixed_orbitals,1,fixed_orbitals,1])
      
    # MULTIPLE ORBITAL SELECTION
    else:
      print('    Fixing to orbitals:', fixed_orbitals)
      no_iwn = np.asarray(siw_dd).shape[1]

      # Variant 1
      #print '    ==> Using EXTRAPOLATION!'
      #siw_t2g_real_extrapolated = np.zeros(np.asarray(fixed_orbitals).shape[0])
      #for b in range(0,np.asarray(fixed_orbitals).shape[0]):
        #siw_t2g_real_extrapolated[b] = np.real(siw_dd[0][int(no_iwn/2),fixed_orbitals[b],1,fixed_orbitals[b],1]) \
              #- ( np.real(siw_dd[0][int(no_iwn/2)+1,fixed_orbitals[b],1,fixed_orbitals[b],1]) - np.real(siw_dd[0][int(no_iwn/2),fixed_orbitals[b],1,fixed_orbitals[b],1]) ) \
              #/ ( iwf[int(no_iwn/2)+1] - iwf[int(no_iwn/2)]) * iwf[int(no_iwn/2)]
      #siw_t2g_real_avg = 0
      #for b in range(0,np.asarray(fixed_orbitals).shape[0]):
        #siw_t2g_real_avg = siw_t2g_real_avg + siw_t2g_real_extrapolated[b]
      #siw_t2g_real_avg = siw_t2g_real_avg/np.asarray(fixed_orbitals).shape[0]    

      # Alternative p+d DC Method
      siw_t2g_real_avg = 0
      for b in range(0, np.asarray(fixed_orbitals).shape[0]):
          print('    ==> Using first Sigma value!')
          siw_t2g_real_avg = (siw_t2g_real_avg
                              + np.real(siw_dd[0][no_iwn//2,
                                                  fixed_orbitals[b], 1,
                                                  fixed_orbitals[b], 1]))
          siw_t2g_real_avg = siw_t2g_real_avg/np.asarray(fixed_orbitals).shape[0]

      ## Alternative p+d DC Method
      #siw_t2g_real_avg = 0
      #print '    ==> Using Sigma oo value!'
      #for b in range(0,np.asarray(fixed_orbitals).shape[0]):
       #siw_t2g_real_avg = siw_t2g_real_avg + np.real(siw_dd[0][int(0),fixed_orbitals[b],1,fixed_orbitals[b],1])
      #siw_t2g_real_avg = siw_t2g_real_avg/np.asarray(fixed_orbitals).shape[0]

      # Z Factor Method
      # siw_t2g_real_avg = 0
      # for b in range(0, np.asarray(fixed_orbitals).shape[0]):
      #  print('    ==> Using first Sigma value times Z-factor!')
      # siw_t2g_real_avg = (siw_t2g_real_avg
      #                     + np.real(siw_dd[0][no_iwn//2,
      #                                         fixed_orbitals[b], 1,
      #                                         fixed_orbitals[b], 1]))
      # siw_t2g_real_avg = siw_t2g_real_avg/np.asarray(fixed_orbitals).shape[0]

      # for b in range(0, np.asarray(fixed_orbitals).shape[0]):
      #     alpha = (-(np.real(siw_dd[0][no_iwn//2 + 1,
      #                                  fixed_orbitals[b], 1,
      #                                  fixed_orbitals[b], 1])
      #                - np.real(siw_dd[0][no_iwn//2,
      #                                    fixed_orbitals[b], 1,
      #                                    fixed_orbitals[b], 1]))
      #              / (iwf[no_iwn//2+1] - iwf[no_iwn//2]))
      #   Z = 1/(1+alpha)
      #   print('    ==> alpha:', alpha, ', Z:', Z)
      # siw_t2g_real_avg = Z * siw_t2g_real_avg


    a=0
    for a in range(0,natoms):
      for b in range(0,np.asarray(siw_dd).shape[2]):
        idx = np.int(b + a*dc.shape[0]/natoms)
        #print '              DC Bands:', idx, a
        dc[idx,0,idx,0] = siw_t2g_real_avg
        dc[idx,1,idx,1] = siw_t2g_real_avg
    
    
    return dc
