""" @package postprocessing

Provides methods for annotation and post-processing of QMC output scripts.
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
import inspect
from functools import wraps
import math

def divide(enum, denom, enum_err=None, denom_err=None):
    """divides and propagates the uncertainty, assuming no covariance"""
    has_enum_err = enum_err is not None
    has_denom_err = denom_err is not None
    if has_enum_err or has_denom_err:
        if not has_enum_err: enum_err = np.zeros_like(denom_err)
        if not has_denom_err: denom_err = np.zeros_like(enum_err)
        ratio_err = np.abs(enum_err/enum)**2 + np.abs(denom_err/denom)**2
        return enum/denom, ratio_err
    else:
        return enum/denom

def multiply(left, right, left_err=None, right_err=None):
    """multiplies and propagates the uncertainty, assuming no covariance"""
    has_left_err = left_err is not None
    has_right_err = right_err is not None
    if has_left_err or has_right_err:
        if not has_left_err: left_err = np.zeros_like(right_err)
        if not has_right_err: right_err = np.zeros_like(left_err)
        product_err = 0
        return left*right, product_err
    else:
        return left*right

def fslice(x, minval, maxval, stride=None):
    """returns a slice of a sorted float interval"""
    if minval is not None: minval = np.searchsorted(x, minval)
    if maxval is not None: maxval = np.searchsorted(x, maxval)
    if stride is not None: stride = np.round(stride/(x[1]-x[0]))
    return slice(minval, maxval, stride)

def slice_center(subsize, size):
    """Returns a slice for the center elements of a bigger array"""
    if subsize < 0 or size < 0:
        raise ValueError("sizes must be bigger than zero")
    if subsize > size:
        raise ValueError("subsize > size")
    if subsize % 2 != size % 2:
        raise ValueError("size and subsize are different MOD 2")
    if subsize == size:
        return slice(None)
    else:
        startidx = (size - subsize)//2
        cslice = slice(startidx, startidx + subsize)
        return cslice

def unwind_dim(dim):
    """Allows treatment of inequivalent atoms and such"""
    def unwind_decorator(func):
        @wraps(func)
        def unwind_wrapper(*args):
            adim = args[0].ndim
            if adim == dim - 1:
                return func(*args)
            if adim != dim:
                raise ValueError("expected %d or %d dimensions, got shape: %s"
                                 % (dim - 1, dim, str(args[0].shape)))
            dimsize = args[0].shape[0]
            unwound_args = []
            for arg in args:
                if not np.ndim(arg):
                    arg = (arg,) * dimsize
                unwound_args.append(arg)

            res = [func(*fargs) for fargs in zip(*unwound_args)]
            # now res contains tuples or arrays
            if isinstance(res[0], tuple):
                nfields = len(res[0])
                res = tuple(np.asarray([elem[i] for elem in res])
                            for i in range(nfields))
                #print "shapes:", [elem.shape for elem in res]
                return res
            else:
                res = np.asarray(res)
                #print "shape:", np.asarray(res).shape
                return res
        return unwind_wrapper
    return unwind_decorator

def dia_slice(A, *dims):
    """ Returns a slicing object to be used for A that extracts diagonals """
    ndim = len(dims)//2
    dims = np.reshape(dims, (ndim, 2))
    slices = [slice(None)] * A.ndim
    for idim, (dim1, dim2) in enumerate(dims):
        shape = int(A.shape[dim1])
        if shape != A.shape[dim2]:
            raise NotImplementedError("not supported yet: %d/%d" % (dim1,dim2)
                                      + str(A.shape))
        extshape = [1]*ndim
        extshape[idim] = -1
        indices = np.arange(shape).reshape(extshape)
        slices[dim1] = slices[dim2] = indices
    return tuple(slices)

# ----- derived quantities -----

def impurity_density(occ, occ_err=None):
    """density of the impurity model"""
    def get_diags(val):
        return np.diagonal(np.diagonal(val, 0, -4, -2), 0, -3, -2)
    if occ_err is not None:
        return get_diags(occ), get_diags(occ_err)
    else:
        return get_diags(occ)

def impurity_electrons(occ, occ_err=None):
    """total number of electrons in the impurity model"""
    def get_traces(val):
        return np.trace(np.trace(val, 0, -3, -1), 0, -2, -1)
    if occ_err is not None:
        return get_traces(occ), np.sqrt(get_traces(np.abs(occ_err)**2))
    else:
        return get_traces(occ)

def lattice_electrons(gdensnew, gdensnew_err=None):
    """number of electrons in the lattice model, total and per spin"""
    # (..., orb, sp, orb, sp) -> (..., sp)
    gdensnew = np.diagonal(np.trace(np.real(gdensnew), 0, -4, -2), 0, -2, -1)
    return (gdensnew[..., 0], gdensnew[..., 1], np.sum(gdensnew, axis=-1))

def lattice_density(gdensnew, gdensnew_err=None):
    """number of electrons in the lattice model per band and spin"""
    # (..., orb, sp, orb, sp) -> (orb, sp)
    return np.diagonal(np.diagonal(np.real(gdensnew), 0, -4, -2), 0, -3, -2)

def sigmaiw_improved(gsigmaiw, giw, gsigmaiw_err=None, giw_err=None):
    """self energy in Matsubara expansion from improved estimators"""
    return divide(gsigmaiw, giw, gsigmaiw_err, giw_err)

@unwind_dim(5)
def moment(occ, occ_err=None):
    """spontaneuous local moment"""
    mu = occ[:,0,:,0] - occ[:,0,:,1] - occ[:,1,:,0] + occ[:,1,:,1]
    if occ_err is not None:
        mu_err = np.sqrt((occ_err**2).sum(-1).sum(1))
        return mu, mu_err
    else:
        return mu

def sigmaiw_short(iw, sigmaiw):
    """self energy in Matsubara expansion in the range 0 ... 20"""
    sl = fslice(iw, 0, 20)
    return iw[sl], sigmaiw[...,sl]

def get_tauint(data):
    """Bias-corrected logarithmic binning analysis"""
    # See arXiv:1810.05079v2, Sec. III.B
    ndat, nm = data.shape
    m = 2**np.arange(nm)
    lbap = m[None,:-1] * (4 * data[:,1:] - data[:,:-1]) / data[:,:1]
    tauint = (lbap - 1)/2
    return tauint

def get_spectrum(data, tau=None):
    """Perform logarithmic spectral analysis"""
    import scipy.optimize as sp_opt

    # See arXiv:1810.05079v2, Sec. IV.B
    def corr_wedge_conv(tau, m):
        alpha = np.exp(-1/tau)
        return alpha/(1 - alpha)**2 * (1 - alpha**m)**2/m

    # Cut away infinite variance points
    trusted = np.isfinite(data).all(0).nonzero()[0].max() + 1
    data = data[:,:trusted]
    # get grid: Eq. (33)
    ndat, nm = data.shape
    m = 2**np.arange(nm)
    if tau is None: tau = m[:-1].copy()
    # set up GLS: Eqs. (29) & (30)
    A = corr_wedge_conv(tau[None,:], m[:-1,None])
    b = m[None,:-1] * (2 * data[:,1:] - data[:,:-1])
    sigma_invsqrt = (m[:-1])**(-0.5)
    # GLS to OLS: Eq. (36)
    A = sigma_invsqrt[:,None] * A
    b = sigma_invsqrt[None,:] * b
    # Use NNLS: Eq. (36)
    x = np.array([sp_opt.nnls(A, bi)[0] for bi in b])
    # Divide by norm to pass from autocov to autocorr spectra
    x /= x.sum(1)[:,None]
    return x

def spec_tau_est(x, tau=None):
    """Compute tauint estimate from spectrum"""
    def corr_sum(tau):
        alpha = np.exp(-1/tau)
        return (1 + alpha)/(1 - alpha)

    ndat, ntau = x.shape
    if tau is None: tau = 2**np.arange(ntau)
    tauint_comp = corr_sum(tau)
    tauint = x.dot(tauint_comp)
    tauint = (tauint - 1) / 2
    return tauint

def gtau_lbap(gtau_blocks, gtau_blocks_err=None):
    """Bias-corrected logarithmic binning analysis for G(tau)"""
    # See arXiv:1810.05079v2, Sec. III.B
    nm = gtau_blocks.shape[-1]
    gtau_flat = gtau_blocks.reshape(-1, nm)
    lbap = get_tauint(gtau_flat)
    return lbap.reshape(gtau_blocks.shape[:-1] + (nm-1,))

def gtau_spec(gtau_blocks, gtau_blocks_err=None):
    """Logarithmic autocorrelation spectrum for G(tau)"""
    nm = gtau_blocks.shape[-1]
    gtau_flat = gtau_blocks.reshape(-1, nm)
    spectrum = get_spectrum(gtau_flat)
    return spectrum.reshape(gtau_blocks.shape[:-1] + spectrum.shape[-1:])

def gtau_tauint(gtau_blocks, gtau_blocks_err=None):
    """Integrated autocorrelation time estimate for G(tau)"""
    nm = gtau_blocks.shape[-1]
    gtau_flat = gtau_blocks.reshape(-1, nm)
    tauint = spec_tau_est(get_spectrum(gtau_flat))
    tauint = tauint.reshape(gtau_blocks.shape[:-1])
    return tauint

def meta_from_func(func, axes, fields=["value"], base=None):
    """Returns a meta dictionary based on a function"""
    meta = {}
    meta["axes"] = axes
    meta["fields"] = fields
    meta["func"] = func
    meta["desc"] = inspect.getdoc(func)
    if base is None:
        base = []
        error_mode = False
        for arg in inspect.getfullargspec(func).args:
            if arg.endswith("_err"):
                error_mode = True
                continue
            if error_mode: raise ValueError("errors have to come last")
            base.append(arg.replace("_","-"))
    meta["base"] = tuple(base)
    return meta

@unwind_dim(4)
def get_ggstraight_ph(giw, niw4f):
    """Helper function for getting the straight part from GG

    The "straight" disconnected part of the two-particle Green's function is
    given by:

       GG_AB(iv, iv', iw) = G_A(iv) G_B(iv') delta(iw,0)

    and is returned as six-dimensional array GG(A,B,iv,iv'), omitting the
    bosonic frequency for brevity.
    """
    nband, nspin, niw = giw.shape
    iw4f_slice = slice_center(niw4f, niw)
    giw = giw.reshape(-1, niw)[:, iw4f_slice]
    gg_straight = np.tensordot(giw, giw, ((),()))  # i,iv,j,iv'
    gg_straight = gg_straight.transpose(0, 2, 1, 3)
    return gg_straight.reshape(nband, nspin, nband, nspin, niw4f, niw4f)

@unwind_dim(4)
def get_ggcross_ph(giw, niw4b, niw4f=None):
    """Helper function for getting the cross part from GG

    The "cross" disconnected part of the two-particle Green's function is
    given by:

       GG_AB(iv, iv', iw) = G_A(iv) G_A(iv+iw) delta(iv,iv') delta(A,B)

    and is returned as four-dimensional array GG(A,iv,iw), omitting the index
    and fermionic frequency for the second particle for brevity.
    """
    def move_slice(sl, offset):  # Moves slice_center() by some offset
        return slice(sl.start + offset, sl.stop + offset)

    nband, nspin, niw = giw.shape
    if niw4f is None: niw4f = niw - niw4b - 1
    iw4f_slice = slice_center(niw4f, niw)
    iw4b_range = range(-(niw4b//2), niw4b//2+1)

    giw = giw.reshape(-1, niw)
    gg_cross =  np.dstack(giw[:, iw4f_slice, None] *
                          giw[:, move_slice(iw4f_slice, iwbos), None]
                          for iwbos in iw4b_range)  # i,iv,iw
    return gg_cross.reshape(nband, nspin, niw4f, niw4b)

@unwind_dim(4)
def get_g4iw_disconn_ph(giw, niw4b, niw4f):
    nband, nspin, niw = giw.shape
    c = np.zeros(shape=(nband, nspin, nband, nspin, niw4f, niw4f,niw4b),dtype=complex)
    iw4b0 = niw4b//2
    c[...,iw4b0] = get_ggstraight_ph(giw, niw4f)
    diasl = dia_slice(c, 0, 2, 1, 3, 4, 5) # b,s,b,s,iv,iv',iw
    c[diasl] -= get_ggcross_ph(giw, niw4b, niw4f)
    return c

@unwind_dim(4)
def get_chi_ph(giw, g4iw, giw_err=None, g4iw_err=None):
    """generalised susceptibility (particle-hole channel)"""
    iw4b0 = g4iw.shape[-1]//2

    chi = g4iw.copy()
    chi[..., iw4b0] -= get_ggstraight_ph(giw, g4iw.shape[-2])
    return chi

@unwind_dim(4)
def get_g4iw_conn_ph(giw, g4iw, giw_err=None, g4iw_err=None):
    """connected part of susceptiblity (particle-hole channel)"""
    conn = get_chi_ph(giw, g4iw, giw_err, g4iw_err)
    diasl = dia_slice(conn, 0, 2, 1, 3, 4, 5) # b,s,b,s,iv,iv',iw
    conn[diasl] += get_ggcross_ph(giw, g4iw.shape[-1], g4iw.shape[-2])
    return conn

@unwind_dim(4)
def get_ggstraight_pp(giw, g4iw_pp_shape):
    """Computes GG = G(iv)G(iv')delta(iw',-iv-iv')"""
    #print giw.shape, g4iw_pp_shape
    assert giw.shape[-3:] == g4iw_pp_shape[-5:-2], "Shape mismatch"
    dotnot = ((),())
    nneq = g4iw_pp_shape[0]
    N = g4iw_pp_shape[-3]
    K = g4iw_pp_shape[-1]
    KhpNm1 = + K//2 - N + 1
    # TODO slow
    chi0_pp = np.zeros(shape=g4iw_pp_shape, dtype=complex)
    #chi0_pp[...] = np.nan
    for m in range(N):
        for n in range(N):
            ktarg = KhpNm1 + m + n
            if 0 <= ktarg < K:
                chi0_pp[...,m,n,ktarg] = \
                           np.tensordot(giw[...,m], giw[...,n], dotnot)
    #
    return chi0_pp

@unwind_dim(4)
def get_chi_pp(giw, g4iw_pp, giw_err=None, g4iw_pp_err=None):
    """generalised susceptibility (particle-particle channel)"""
    iw4st = (giw.shape[-1] - g4iw_pp.shape[-2])//2
    iw4sl = slice(iw4st, -iw4st)
    chi0_pp = get_ggstraight_pp(giw[...,iw4sl], g4iw_pp.shape)
    return g4iw_pp - chi0_pp


def get_siw_mom(u_matrix, rho1, rho2):
    """Calculates the high frequency moments of the self-energy"""

    # reshape the arrays to contracted band and spin indices
    nbands = int(u_matrix.shape[0]/2)

    rho1_contracted=rho1.reshape((nbands*2),(nbands*2))
    rho2_contracted=rho2.reshape((nbands*2),(nbands*2),(nbands*2),(nbands*2))

    # Use antisymmetrisized U-matrix
    u_antisym=-0.5*(u_matrix-np.swapaxes(u_matrix,0,1)-np.swapaxes(u_matrix,2,3)+np.swapaxes(np.swapaxes(u_matrix,0,1),2,3))

    sigma_inf=np.tensordot(u_antisym,rho1_contracted,axes=([1,2],[0,1]))

    ## For formulas see Wang, Dang, Millis PRB 84, 073104
    #L=-0.5*(np.tensordot(u_antisym,u_antisym,axes=([2,3],[0,1])))

    #first_term=0.25*np.tensordot(u_antisym,u_antisym,axes=(1,2))
    #second_term=-np.tensordot(u_antisym,u_antisym,axes=(3,0))

    #k1=np.tensordot(first_term,rho2_contracted,axes=([1,2,3,4],[2,3,0,1]))
    #k2=np.tensordot(second_term,rho2_contracted,axes=([1,2,3,4],[0,2,1,3]))

    #K=np.add(k1,k2)

    #sigma_1=-np.add(\
            #np.add(K,\
            #np.tensordot(L,rho1_contracted,axes=([1,2],[0,1]))),\
            #-np.tensordot(sigma_inf,sigma_inf,axes=(1,0)))

    sigma_inf=sigma_inf.reshape(nbands,2,nbands,2)

    # collect computed moments
    ### TODOcompl: rho2 not yet available;
    ### setting incomplete sigma_1 to zero here
    sigma_1=np.zeros(shape=(nbands,2,nbands,2))
    smoms = np.asarray([sigma_inf, sigma_1])
    return smoms


def get_sztau_sz0_diag(ntau_n0, ntau_n0_err=None):
    "< S_z(tau) S_z(0) > diagonal in orbitals"
    nbands=ntau_n0.shape[0]
    ntau=ntau_n0.shape[-1]
    sztau_sz0=np.zeros((ntau))

    for b1 in range(0,nbands):
        # this is not the formula that Werner wrote in his spin-freezing PRL,
        # but used to produce the plot
        sztau_sz0[:] = sztau_sz0[:] + \
            ntau_n0[b1,0,b1,0,:] + ntau_n0[b1,1,b1,1,:] - \
            ntau_n0[b1,1,b1,0,:] - ntau_n0[b1,0,b1,1,:]

    if ntau_n0_err is not None:
        sztau_sz0_err = np.zeros_like(sztau_sz0)
        for b1 in range(0, nbands):
            sztau_sz0_err[:] = sztau_sz0_err[:] + \
                ntau_n0_err[b1,0,b1,0,:]**2 + ntau_n0_err[b1,1,b1,1,:]**2 + \
                ntau_n0_err[b1,1,b1,0,:]**2 + ntau_n0_err[b1,0,b1,1,:]**2
        sztau_sz0_err = np.sqrt(sztau_sz0_err)
        return sztau_sz0, sztau_sz0_err

    return sztau_sz0

def get_sztau_sz0(ntau_n0, ntau_n0_err=None):
    "< S_z(tau) S_z(0) > full"
    nbands=ntau_n0.shape[0]
    ntau=ntau_n0.shape[-1]
    sztau_sz0=np.zeros((ntau))

    for b1 in range(0,nbands):
        for b2 in range(0,nbands):

            sztau_sz0[:]=sztau_sz0[:]+(ntau_n0[b1,0,b2,0,:]+ntau_n0[b1,1,b2,1,:]-ntau_n0[b1,1,b2,0,:]-ntau_n0[b1,0,b2,1,:])

    if ntau_n0_err is not None:
        sztau_sz0_err = np.zeros_like(sztau_sz0)
        for b1 in range(0, nbands):
            for b2 in range(0, nbands):
                sztau_sz0_err[:] = sztau_sz0_err[:] + (ntau_n0_err[b1,0,b2,0,:]**2
                                                       + ntau_n0_err[b1,1,b2,1,:]**2
                                                       + ntau_n0_err[b1,1,b2,0,:]**2
                                                       + ntau_n0_err[b1,0,b2,1,:]**2)
        sztau_sz0_err = np.sqrt(sztau_sz0_err)
        return sztau_sz0, sztau_sz0_err

    return sztau_sz0

def get_Ntau_N0(ntau_n0, ntau_n0_err=None):
    "< N(tau) N(0) > full"
    if ntau_n0_err is not None:
        return np.sum(ntau_n0, axis=(0, 1, 2, 3)), \
               np.sqrt(np.sum(ntau_n0_err**2, axis=(0, 1, 2, 3)))
    return np.sum(ntau_n0, axis=(0, 1, 2, 3))

def get_sztau_sz0_orb_resolved(ntau_n0, ntau_n0_err=None):
    "< S_z(tau) S_z(0) > orbital resolved"
    nbands=ntau_n0.shape[0]
    ntau=ntau_n0.shape[-1]
    sztau_sz0=np.zeros((nbands,ntau))

    for b1 in range(0,nbands):
        sztau_sz0[b1,:]=sztau_sz0[b1,:]+(ntau_n0[b1,0,b1,0,:]+ntau_n0[b1,1,b1,1,:]-ntau_n0[b1,1,b1,0,:]-ntau_n0[b1,0,b1,1,:])

    if ntau_n0_err is not None:
        sztau_sz0_err = np.zeros_like(sztau_sz0)
        for b1 in range(0, nbands):
            sztau_sz0_err[b1, :] = sztau_sz0_err[b1, :] + (ntau_n0_err[b1,0,b1,0,:]**2
                                                            + ntau_n0_err[b1,1,b1,1,:]**2
                                                            + ntau_n0_err[b1,1,b1,0,:]**2
                                                            + ntau_n0_err[b1,0,b1,1,:]**2)
        sztau_sz0_err = np.sqrt(sztau_sz0_err)
        return sztau_sz0, sztau_sz0_err

    return sztau_sz0

def get_sztau_sz0_offdiag(ntau_n0, ntau_n0_err=None):
    "< S_z(tau) S_z(0) > orbital offdiagonal components"
    nbands=ntau_n0.shape[0]
    ntau=ntau_n0.shape[-1]
    sztau_sz0=np.zeros((ntau))

    for b1 in range(0,nbands):
        for b2 in range(0,nbands):
            if b1!=b2:
                sztau_sz0[:]=sztau_sz0[:]+(ntau_n0[b1,0,b2,0,:]+ntau_n0[b1,1,b2,1,:]-ntau_n0[b1,1,b2,0,:]-ntau_n0[b1,0,b2,1,:])

    if ntau_n0_err is not None:
        sztau_sz0_err = np.zeros_like(sztau_sz0)
        for b1 in range(0, nbands):
            for b2 in range(0, nbands):
                if b1 != b2:
                    sztau_sz0_err[:] = sztau_sz0_err[:] + (ntau_n0_err[b1,0,b2,0,:]**2
                                                            + ntau_n0_err[b1,1,b2,1,:]**2
                                                            + ntau_n0_err[b1,1,b2,0,:]**2
                                                            + ntau_n0_err[b1,0,b2,1,:]**2)
        sztau_sz0_err = np.sqrt(sztau_sz0_err)
        return sztau_sz0, sztau_sz0_err

    return sztau_sz0

# Derived quantities for the iterations
derived_quantities = {
    "imp-density": meta_from_func(
                        impurity_density,
                        ["ineq","band","spin"], ["value", "error"], ["occ"]
                        ),
    "imp-electrons": meta_from_func(
                        impurity_electrons,
                        ["ineq"], ["value", "error"]
                        ),
    "moment": meta_from_func(
                        moment,
                        ["ineq","band1","band2"], ["value","error"],
                        ["occ"]
                        ),
    "latt-density": meta_from_func(lattice_density,
                        ["lda-band", "spin"], ["value"]
                        ),
    "latt-electrons": meta_from_func(lattice_electrons,
                        [], ["spin-up", "spin-dn", "total"]
                        ),
    "chi-ph": meta_from_func(get_chi_ph,
                        ['ineq','band1','spin1','band2','spin2','iwf-g4','iwf-g4','iwb-g4'],
                        ['value'],
                        ['giw', 'g4iw']
                        ),
    "chi-pp": meta_from_func(get_chi_pp,
                        ['ineq','band1','spin1','band2','spin2','iwf-g4','iwf-g4','iwb-g4'],
                        ['value'],
                        ['giw', 'g4iw-pp'],
                        ),
    "g4iw-conn-ph": meta_from_func(get_g4iw_conn_ph,
                        ['ineq','band1','spin1','band2','spin2','iwf-g4','iwf-g4','iwb-g4'],
                        ['value'],
                        ['giw', 'g4iw']
                        ),
    "sztau-sz0": meta_from_func(get_sztau_sz0,
                        ['ineq', 'tausus'],
                        ['value', "error"],
                        ['ntau-n0']
                        ),
    "sztau-sz0-diag": meta_from_func(get_sztau_sz0_diag,
                        ['ineq', 'tausus'],
                        ['value', "error"],
                        ['ntau-n0']
                        ),
    "sztau-sz0-orb-resolved": meta_from_func(get_sztau_sz0_orb_resolved,
                        ['ineq', 'band', 'tausus'],
                        ['value', "error"],
                        ['ntau-n0']
                        ),
    "sztau-sz0-offdiag": meta_from_func(get_sztau_sz0_offdiag,
                        ['ineq', 'tausus',],
                        ['value', "error"],
                        ['ntau-n0']
                        ),
    "denstau-dens0": meta_from_func(get_Ntau_N0,
                        ['ineq', 'tausus'],
                        ['value', "error"],
                        ['ntau-n0']
                        ),
    "gtau-tauint": meta_from_func(
                        gtau_tauint,
                        ["ineq","band","spin","taubin"],
                        ["value"],
                        ["gtau-blocks"]
                        ),
    "gtau-spec": meta_from_func(
                        gtau_spec,
                        ["ineq","band","spin","taubin","m"],
                        ["value"],
                        ["gtau-blocks"]
                        ),
    "gtau-lbap": meta_from_func(
                        gtau_lbap,
                        ["ineq","band","spin","taubin","m"],
                        ["value"],
                        ["gtau-blocks"]
                        ),
    }
