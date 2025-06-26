from __future__ import division, print_function
import pathlib
import ctypes as _c
from typing import Sequence, Literal
from enum import IntEnum
import numpy as np

# ----- Accumulators

_lib = _c.CDLL(pathlib.Path(__file__).parent / "libCTQMC.so")

def _annotate_accumulators(lib):
    lib.accumulator_create.argtypes = [_c.c_bool, _c.c_int64, _c.c_int64, _c.c_int64]
    lib.accumulator_create.restype = _c.c_void_p
    lib.accumulator_n.argtypes = [_c.c_void_p]
    lib.accumulator_n.restype = _c.c_int64
    lib.accumulator_add.argtypes = [_c.c_void_p]
    lib.accumulator_add.restype = None
    lib.accumulator_buffer.argtypes = [_c.c_void_p]
    lib.accumulator_buffer.restype = _c.POINTER(_c.c_double)
    lib.accumulator_reset.argtypes = [_c.c_void_p]
    lib.accumulator_reset.restype = None
    lib.accumulator_delete.argtypes = [_c.c_void_p]
    lib.accumulator_delete.restype = None
    lib.accumulator_mean.argtypes = [_c.c_void_p, _c.POINTER(_c.c_double)]
    lib.accumulator_mean.restype = _c.c_bool

_annotate_accumulators(_lib)

class Accumulator:
    def __new__(cls, *args, **kwargs):
        instance = super().__new__(cls)
        instance._acc = None
        return instance

    def __init__(self, iscomplex, typ, ncomp, nlevels=0):
        self.iscomplex = iscomplex
        self.typ = typ
        self.ncomp = ncomp
        self.nlevels = nlevels

        # Initialize accumulator
        self._acc = _lib.accumulator_create(iscomplex, typ, ncomp, nlevels)
        if self._acc is None:
            raise ValueError("invalid arguments")

        # Get buffer
        bufptr = _lib.accumulator_buffer(self._acc)
        if iscomplex:
            dbuffer = np.ctypeslib.as_array(bufptr, (2 * ncomp,))
            self._buffer = dbuffer.view(np.complex128)
        else:
            self._buffer = np.ctypeslib.as_array(bufptr, (ncomp,))

    def reset(self):
        _lib.accumulator_reset(self._acc)

    @property
    def n(self):
        return _lib.accumulator_n(self._acc)

    def add(self, val):
        self._buffer += val
        _lib.accumulator_add(self._acc)

    @property
    def mean(self):
        if self.iscomplex:
            out = np.empty(self.ncomp, np.complex128)
        else:
            out = np.empty(self.ncomp, np.double)
        out_ptr = out.ctypes.data_as(_c.POINTER(_c.c_double))
        if not _lib.accumulator_mean(self._acc, out_ptr):
            raise RuntimeError("Accumulator does not store mean")
        return out

    def __del__(self):
        if self._acc is not None:
            _lib.accumulator_delete(self._acc)


class MeanAcc(Accumulator):
    def __init__(self, ncomp, iscomplex):
        Accumulator.__init__(self, iscomplex, 1, ncomp)


# --- Samplers

STAT_BOSE = 1
STAT_FERMI = -1

def _annotate_samplers(lib):
    lib.rhist_create.argtypes = [_c.c_int64, _c.c_double,
                                 _c.c_int64, _c.c_int64,
                                 _c.c_double]
    lib.rhist_create.restype = _c.c_void_p
    lib.rhist2d_create.argtypes = [_c.c_int64, _c.c_double,
                                   _c.c_int64, _c.c_int64,
                                   _c.c_int64, _c.c_int64,
                                   _c.c_double, _c.c_double,
                                   _c.c_double, _c.c_double]
    lib.rhist2d_create.restype = _c.c_void_p
    lib.rhist3d_create.argtypes = [_c.c_int64, _c.c_double,
                                   _c.c_int64, _c.c_int64,
                                   _c.c_int64, _c.c_int64,
                                   _c.c_int64, _c.c_int64,
                                   _c.c_double, _c.c_double,
                                   _c.c_double, _c.c_double,
                                   _c.c_double, _c.c_double,
                                   _c.c_double, _c.c_double,
                                   _c.c_double]
    lib.rhist3d_create.restype = _c.c_void_p
    lib.nfourier_create.argtypes = [_c.c_int64, _c.c_double,
                                    _c.c_int64, _c.c_int64, _c.POINTER(_c.c_int64),
                                    _c.c_double, _c.c_bool]
    lib.nfourier_create.restype = _c.c_void_p
    lib.nfourier2d_create.argtypes = [_c.c_int64, _c.c_double,
                                      _c.c_int64, _c.c_int64,
                                      _c.c_int64, _c.c_int64,
                                      _c.c_double, _c.c_double,
                                      _c.c_double, _c.c_double]
    lib.nfourier2d_create.restype = _c.c_void_p
    lib.nfft_create.argtypes = [_c.c_int64, _c.c_double,
                                _c.c_int64, _c.c_int64, _c.POINTER(_c.c_int64),
                                _c.c_double, _c.c_bool, _c.c_bool]
    lib.nfft_create.restype = _c.c_void_p
    lib.nfft2d_create.argtypes = [_c.c_int64, _c.c_double,
                                  _c.c_int64, _c.c_int64,
                                  _c.c_int64, _c.c_int64,
                                  _c.c_double, _c.c_double,
                                  _c.c_double, _c.c_double,
                                  _c.c_bool]
    lib.nfft2d_create.restype = _c.c_void_p
    lib.nfft3d_create.argtypes = [_c.c_int64, _c.c_double,
                                  _c.c_int64, _c.c_int64,
                                  _c.c_int64, _c.c_int64,
                                  _c.c_int64, _c.c_int64,
                                  _c.c_double, _c.c_double,
                                  _c.c_double, _c.c_double,
                                  _c.c_double, _c.c_double,
                                  _c.c_double, _c.c_double,
                                  _c.c_double,
                                  _c.c_bool]
    lib.nfft3d_create.restype = _c.c_void_p
    lib.legendre_create.argtypes = [_c.c_int64, _c.c_double,
                                    _c.c_int64, _c.c_int64,
                                    _c.c_double]
    lib.legendre_create.restype = _c.c_void_p
    lib.sampler_ncomp.argtypes = [_c.c_void_p]
    lib.sampler_ncomp.restype = _c.c_int64
    lib.sampler_nacc.argtypes = [_c.c_void_p]
    lib.sampler_nacc.restype = _c.c_int64
    lib.sampler_npost.argtypes = [_c.c_void_p]
    lib.sampler_npost.restype = _c.c_int64
    lib.sampler_statsign.argtypes = [_c.c_void_p]
    lib.sampler_statsign.restype = _c.c_int64
    lib.sampler_beta.argtypes = [_c.c_void_p]
    lib.sampler_beta.restype = _c.c_double
    lib.sampler_add.argtypes = [
                    _c.c_void_p, _c.c_int64, _c.POINTER(_c.c_int64),
                    _c.POINTER(_c.c_double), _c.POINTER(_c.c_double),
                    _c.POINTER(_c.c_double)
                    ]
    lib.sampler_add.restype = None
    lib.sampler_postprocess.argtypes = [_c.c_void_p, _c.POINTER(_c.c_double),
                                        _c.POINTER(_c.c_double)]
    lib.sampler_postprocess.restype = None
    lib.sampler_delete.argtypes = [_c.c_void_p]
    lib.sampler_delete.restype = None
    lib.sampler_register_register.argtypes = [_c.c_void_p, _c.c_void_p,
                                              _c.c_void_p]
    lib.sampler_register_register.restype = None
    lib.sampler2d_ncomp.argtypes = [_c.c_void_p]
    lib.sampler2d_ncomp.restype = _c.c_int64
    lib.sampler2d_nacc.argtypes = [_c.c_void_p]
    lib.sampler2d_nacc.restype = _c.c_int64
    lib.sampler2d_npost.argtypes = [_c.c_void_p]
    lib.sampler2d_npost.restype = _c.c_int64
    lib.sampler2d_statsign1.argtypes = [_c.c_void_p]
    lib.sampler2d_statsign1.restype = _c.c_int64
    lib.sampler2d_statsign2.argtypes = [_c.c_void_p]
    lib.sampler2d_statsign2.restype = _c.c_int64
    lib.sampler2d_beta.argtypes = [_c.c_void_p]
    lib.sampler2d_beta.restype = _c.c_double
    lib.sampler2d_add.argtypes = [
                    _c.c_void_p, _c.c_int64, _c.POINTER(_c.c_int64),
                    _c.POINTER(_c.c_double), _c.POINTER(_c.c_double),
                    _c.POINTER(_c.c_double), _c.POINTER(_c.c_double)
                    ]
    lib.sampler2d_add.restype = None
    lib.sampler2d_postprocess.argtypes = [_c.c_void_p, _c.POINTER(_c.c_double),
                                        _c.POINTER(_c.c_double)]
    lib.sampler2d_postprocess.restype = None
    lib.sampler2d_delete.argtypes = [_c.c_void_p]
    lib.sampler2d_delete.restype = None
    lib.sampler2d_register_register.argtypes = [_c.c_void_p, _c.c_void_p,
                                                _c.c_void_p]
    lib.sampler2d_register_register.restype = None
    lib.sampler3d_ncomp.argtypes = [_c.c_void_p]
    lib.sampler3d_ncomp.restype = _c.c_int64
    lib.sampler3d_nacc.argtypes = [_c.c_void_p]
    lib.sampler3d_nacc.restype = _c.c_int64
    lib.sampler3d_npost.argtypes = [_c.c_void_p]
    lib.sampler3d_npost.restype = _c.c_int64
    lib.sampler3d_statsign1.argtypes = [_c.c_void_p]
    lib.sampler3d_statsign1.restype = _c.c_int64
    lib.sampler3d_statsign2.argtypes = [_c.c_void_p]
    lib.sampler3d_statsign2.restype = _c.c_int64
    lib.sampler3d_statsign3.argtypes = [_c.c_void_p]
    lib.sampler3d_statsign3.restype = _c.c_int64
    lib.sampler3d_beta.argtypes = [_c.c_void_p]
    lib.sampler3d_beta.restype = _c.c_double
    lib.sampler3d_add.argtypes = [
                    _c.c_void_p, _c.c_int64, _c.POINTER(_c.c_int64),
                    _c.POINTER(_c.c_double), _c.POINTER(_c.c_double),
                    _c.POINTER(_c.c_double), _c.POINTER(_c.c_double),
                    _c.POINTER(_c.c_double)
                    ]
    lib.sampler3d_add.restype = None
    lib.sampler3d_postprocess.argtypes = [_c.c_void_p, _c.POINTER(_c.c_double),
                                        _c.POINTER(_c.c_double)]
    lib.sampler3d_postprocess.restype = None
    lib.sampler3d_delete.argtypes = [_c.c_void_p]
    lib.sampler3d_delete.restype = None
    lib.sampler3d_register_register.argtypes = [_c.c_void_p, _c.c_void_p,
                                              _c.c_void_p]
    lib.sampler3d_register_register.restype = None

_annotate_samplers(_lib)


class _SamplerBase(object):
    def __init__(self, handle):
        self._handle = handle
        self.ncomp = _lib.sampler_ncomp(handle)
        self.nacc = _lib.sampler_nacc(handle)
        self.npost = _lib.sampler_npost(handle)
        self.beta = _lib.sampler_beta(handle)
        self.statsign = _lib.sampler_statsign(handle)

    def mkacc(self):
        return np.zeros((self.ncomp, self.nacc), np.double)

    def add(self, i, x, fx, acc):
        i = np.ascontiguousarray(i, np.int64)
        x = np.ascontiguousarray(x, np.double)
        fx = np.ascontiguousarray(fx, np.double)
        acc = np.ascontiguousarray(acc, np.double)
        n = i.size
        if x.size != n or fx.size != n:
            raise ValueError("Arrays must be of same size")
        if acc.shape != (self.ncomp, self.nacc):
            raise ValueError("Accumulator array of wrong size")

        c_double_pointer = _c.POINTER(_c.c_double)
        _lib.sampler_add(self._handle, n,
                         i.ctypes.data_as(_c.POINTER(_c.c_int64)),
                         x.ctypes.data_as(c_double_pointer),
                         fx.ctypes.data_as(c_double_pointer),
                         acc.ctypes.data_as(c_double_pointer)
                         )

    def postprocess(self, acc):
        acc = np.ascontiguousarray(acc)

        # XXX make this generic
        if not np.issubdtype(acc.dtype, np.complexfloating):
            raise RuntimeError("meh")

        #print('acc.shape', acc.shape)
        #print('self.ncomp', self.ncomp)
        #print('self.nacc', self.nacc)
        if acc.size != self.ncomp * self.nacc:
            raise ValueError("Accumulator array of wrong size")

        out = np.zeros((self.ncomp, self.npost), np.complex128)
        c_double_pointer = _c.POINTER(_c.c_double)
        _lib.sampler_postprocess(self._handle,
                                 acc.ctypes.data_as(c_double_pointer),
                                 out.ctypes.data_as(c_double_pointer)
                                )

        return out

    def __del__(self):
        _lib.sampler_delete(self._handle)


class _Sampler2DBase:
    def __init__(self, handle, npostshape):
        self._handle = handle
        self.ncomp = _lib.sampler2d_ncomp(handle)
        self.nacc = _lib.sampler2d_nacc(handle)
        self.npost = _lib.sampler2d_npost(handle)
        self.beta = _lib.sampler2d_beta(handle)
        self.statsign = [_lib.sampler2d_statsign1(handle),
                         _lib.sampler2d_statsign2(handle)]
        self.npostshape = npostshape
        assert len(npostshape) == 2, "Wrong length of npostshape"
        assert npostshape[0] * npostshape[1] == self.npost, \
            "Inconsistent sizes in npostshape"

    def mkacc(self):
        return np.zeros((self.ncomp, self.nacc), np.double)

    def add(self, i, x, y, fx, acc):
        i = np.ascontiguousarray(i, np.int64)
        x = np.ascontiguousarray(x, np.double)
        y = np.ascontiguousarray(y, np.double)
        fx = np.ascontiguousarray(fx, np.double)
        acc = np.ascontiguousarray(acc, np.double)
        n = i.size
        if x.size != n or y.size != n or fx.size != n:
            raise ValueError("Arrays must be of same size")
        if acc.shape != (self.ncomp, self.nacc):
            raise ValueError("Accumulator array of wrong size")

        c_double_pointer = _c.POINTER(_c.c_double)
        _lib.sampler2d_add(self._handle, n,
                           i.ctypes.data_as(_c.POINTER(_c.c_int64)),
                           x.ctypes.data_as(c_double_pointer),
                           y.ctypes.data_as(c_double_pointer),
                           fx.ctypes.data_as(c_double_pointer),
                           acc.ctypes.data_as(c_double_pointer)
                           )

    def postprocess(self, acc):
        acc = np.ascontiguousarray(acc)

        # XXX make this generic
        if not np.issubdtype(acc.dtype, np.complexfloating):
            raise RuntimeError("meh")

        #print('acc.shape', acc.shape)
        #print('self.ncomp', self.ncomp)
        #print('self.nacc', self.nacc)
        if acc.size != self.ncomp * self.nacc:
            raise ValueError("Accumulator array of wrong size")

        out = np.zeros((self.ncomp, self.npost), np.complex128)
        c_double_pointer = _c.POINTER(_c.c_double)
        _lib.sampler2d_postprocess(self._handle,
                                   acc.ctypes.data_as(c_double_pointer),
                                   out.ctypes.data_as(c_double_pointer)
                                   )
        out = np.reshape(out, (self.ncomp, *self.npostshape[::-1]))
        out = np.transpose(out, (0, 2, 1))

        return out

    def __del__(self):
        _lib.sampler2d_delete(self._handle)


class _Sampler3DBase:
    def __init__(self, handle, npostshape):
        self._handle = handle
        self.ncomp = _lib.sampler3d_ncomp(handle)
        self.nacc = _lib.sampler3d_nacc(handle)
        self.npost = _lib.sampler3d_npost(handle)
        self.beta = _lib.sampler3d_beta(handle)
        self.statsign = [_lib.sampler3d_statsign1(handle),
                         _lib.sampler3d_statsign2(handle),
                         _lib.sampler3d_statsign3(handle)]
        self.npostshape = npostshape
        assert len(npostshape) == 3, "Wrong length of npostshape"
        assert npostshape[0] * npostshape[1] * npostshape[2] == self.npost, \
            "Inconsistent sizes in npostshape"

    def mkacc(self):
        return np.zeros((self.ncomp, self.nacc), np.double)

    def add(self, i, x, y, z, fx, acc):
        i = np.ascontiguousarray(i, np.int64)
        x = np.ascontiguousarray(x, np.double)
        y = np.ascontiguousarray(y, np.double)
        z = np.ascontiguousarray(y, np.double)
        fx = np.ascontiguousarray(fx, np.double)
        acc = np.ascontiguousarray(acc, np.double)
        n = i.size
        if x.size != n or y.size != n or z.size != n or fx.size != n:
            raise ValueError("Arrays must be of same size")
        if acc.shape != (self.ncomp, self.nacc):
            raise ValueError("Accumulator array of wrong size")

        c_double_pointer = _c.POINTER(_c.c_double)
        _lib.sampler3d_add(self._handle, n,
                           i.ctypes.data_as(_c.POINTER(_c.c_int64)),
                           x.ctypes.data_as(c_double_pointer),
                           y.ctypes.data_as(c_double_pointer),
                           z.ctypes.data_as(c_double_pointer),
                           fx.ctypes.data_as(c_double_pointer),
                           acc.ctypes.data_as(c_double_pointer)
                           )

    def postprocess(self, acc):
        acc = np.ascontiguousarray(acc)

        # XXX make this generic
        if not np.issubdtype(acc.dtype, np.complexfloating):
            raise RuntimeError("meh")

        #print('acc.shape', acc.shape)
        #print('self.ncomp', self.ncomp)
        #print('self.nacc', self.nacc)
        if acc.size != self.ncomp * self.nacc:
            raise ValueError("Accumulator array of wrong size")

        out = np.zeros((self.ncomp, self.npost), np.complex128)
        c_double_pointer = _c.POINTER(_c.c_double)
        _lib.sampler3d_postprocess(self._handle,
                                   acc.ctypes.data_as(c_double_pointer),
                                   out.ctypes.data_as(c_double_pointer)
                                   )
        out = np.reshape(out, (self.ncomp, *self.npostshape[::-1]))
        out = np.transpose(out, (0, 3, 2, 1))

        return out

    def __del__(self):
        _lib.sampler3d_delete(self._handle)


class RHistSampler(_SamplerBase):
    def __init__(self, ncomp, beta, statsign, nbins, coeff):
        if ncomp < 0:
            raise ValueError("number of components must be nonneg.")
        if beta <= 0:
            raise ValueError("beta must be positive")
        if statsign not in (-1, 1):
            raise ValueError("statistics sign is wrong")
        if nbins < 0:
            raise ValueError("number of bins must be nonneg.")

        handle = _lib.rhist_create(ncomp, beta, statsign, nbins, coeff)
        super().__init__(handle)


class RHist2DSampler(_Sampler2DBase):
    def __init__(self, ncomp, beta, statsigns, nbins, coeffs):
        if ncomp < 0:
            raise ValueError("number of components must be nonneg.")
        if beta <= 0:
            raise ValueError("beta must be positive")
        if any(sign not in (-1, 1) for sign in statsigns):
            raise ValueError("statistics sign is wrong")
        if any(nbin < 0 for nbin in nbins):
            raise ValueError("number of bins must be nonneg.")

        handle = _lib.rhist2d_create(ncomp, beta,
                                     statsigns[0], statsigns[1],
                                     nbins[0], nbins[1],
                                     coeffs[0, 0], coeffs[1, 0],
                                     coeffs[0, 1], coeffs[1, 1])
        super().__init__(handle, nbins)


class RHist3DSampler(_Sampler3DBase):
    def __init__(self, ncomp, beta, statsigns, nbins, coeffs):
        if ncomp < 0:
            raise ValueError("number of components must be nonneg.")
        if beta <= 0:
            raise ValueError("beta must be positive")
        if any(sign not in (-1, 1) for sign in statsigns):
            raise ValueError("statistics sign is wrong")
        if any(nbin < 0 for nbin in nbins):
            raise ValueError("number of bins must be nonneg.")

        handle = _lib.rhist3d_create(ncomp, beta,
                                     statsigns[0], statsigns[1], statsigns[2],
                                     nbins[0], nbins[1], nbins[2],
                                     coeffs[0, 0], coeffs[1, 0], coeffs[2, 0],
                                     coeffs[0, 1], coeffs[1, 1], coeffs[2, 1],
                                     coeffs[0, 2], coeffs[1, 2], coeffs[2, 2])
        super().__init__(handle, nbins)


class FourierSampler(_SamplerBase):
    def __init__(self, ncomp, beta, statsign, niw, coeff, use_nfft=True, sparse_frequencies=[]):
        if ncomp < 0:
            raise ValueError("number of components must be nonneg.")
        if beta <= 0:
            raise ValueError("beta must be positive")
        if statsign not in (-1, 1):
            raise ValueError("statistics sign is wrong")
        if niw < 0:
            raise ValueError("number of frequencies must be nonneg.")
        
        if sparse_frequencies:
            niw = len(sparse_frequencies)
            use_sparse = True
        else:
            use_sparse = False

        sparse_frequencies = np.ascontiguousarray(sparse_frequencies, np.int64)
        sparse_frequencies_arg = sparse_frequencies.ctypes.data_as(_c.POINTER(_c.c_int64))

        if use_nfft:
            handle = _lib.nfft_create(ncomp, beta, statsign, niw, sparse_frequencies_arg, coeff, True, use_sparse)
        else:
            handle = _lib.nfourier_create(ncomp, beta, statsign, niw, sparse_frequencies_arg, coeff, use_sparse)
        
        super().__init__(handle)
        #print("TEST", hasattr(self, '_handle'))


class Fourier2DSampler(_Sampler2DBase):
    def __init__(self, ncomp, beta, statsigns, niws, coeffs, use_nfft=True):
        if ncomp < 0:
            raise ValueError("number of components must be nonneg.")
        if beta <= 0:
            raise ValueError("beta must be positive")
        if any(sign not in (-1, 1) for sign in statsigns):
            raise ValueError("statistics sign is wrong")
        if any(niw < 0 for niw in niws):
            raise ValueError("number of frequencies must be nonneg.")

        if use_nfft:
            handle = _lib.nfft2d_create(ncomp, beta,
                                        statsigns[0], statsigns[1],
                                        niws[0], niws[1],
                                        coeffs[0, 0], coeffs[1, 0],
                                        coeffs[0, 1], coeffs[1, 1],
                                        True)
        else:
            #raise NotImplementedError("For 2D Fourier sampler, "
            #                          "only use_nfft=True supported")
            # FIXME
            handle = _lib.nfourier2d_create(ncomp, beta, statsigns[0], statsigns[1], niws[0], niws[1], coeffs[0, 0], coeffs[1, 0],
                                        coeffs[0, 1], coeffs[1, 1])

        super().__init__(handle, [2 * niw if s == STAT_FERMI else 2 * niw - 1
                                  for niw, s in zip(niws, statsigns)])


class Fourier3DSampler(_Sampler3DBase):
    def __init__(self, ncomp, beta, statsigns, niws, coeffs, use_nfft=True):
        if ncomp < 0:
            raise ValueError("number of components must be nonneg.")
        if beta <= 0:
            raise ValueError("beta must be positive")
        if any(sign not in (-1, 1) for sign in statsigns):
            raise ValueError("statistics sign is wrong")
        if any(niw < 0 for niw in niws):
            raise ValueError("number of frequencies must be nonneg.")

        if use_nfft:
            handle = _lib.nfft3d_create(ncomp, beta,
                                        statsigns[0], statsigns[1], statsigns[2],
                                        niws[0], niws[1], niws[2],
                                        coeffs[0, 0], coeffs[1, 0], coeffs[2, 0],
                                        coeffs[0, 1], coeffs[1, 1], coeffs[2, 1],
                                        coeffs[0, 2], coeffs[1, 2], coeffs[2, 2],
                                        True)
        else:
            raise NotImplementedError("For 3D Fourier sampler, "
                                      "only use_nfft=True supported")
            # FIXME
            # handle = _lib.nfourier2d_create(ncomp, beta, statsigns, niws)

        super().__init__(handle, [2 * niw if s == STAT_FERMI else 2 * niw - 1
                                  for niw, s in zip(niws, statsigns)])


class LegendreSampler(_SamplerBase):
    def __init__(self, ncomp, beta, statsign, nleg, coeff):
        if ncomp < 0:
            raise ValueError("number of components must be nonneg.")
        if beta <= 0:
            raise ValueError("beta must be positive")
        if statsign not in (-1, 1):
            raise ValueError("statistics sign is wrong")
        if nleg < 0:
            raise ValueError("number of coefficients must be nonneg.")

        handle = _lib.legendre_create(ncomp, beta, statsign, nleg, coeff)
        super().__init__(handle)


# worm operator types
class OP_W(IntEnum):
    CREA = 2
    ANNH = -2


def _annotate_ctqmc(lib):
    # setup / breakdown
    lib.init_qmc_rng.argtypes = [_c.c_int32]
    lib.init_qmc_rng.restype = None
    lib.reset_qmc_counters.argtypes = None
    lib.reset_qmc_counters.restype = None
    lib.set_fancy_progress_bar.argtypes = [_c.c_bool]
    lib.set_fancy_progress_bar.restype = None
    lib.set_qmc_simid.argtypes = [_c.c_int32]
    lib.set_qmc_simid.restype = None
    lib.set_sector_eta.argtypes = [_c.c_int32, _c.c_double]
    lib.set_sector_eta.restype = None
    lib.get_sector_eta.argtypes = [_c.c_int32]
    lib.get_sector_eta.restype = _c.c_double
    lib.set_qmc_taudiffmax.argtypes = [_c.c_double]
    lib.set_qmc_taudiffmax.restype = None
    lib.init_qmc_solver.argtypes = [_c.c_int32,
                                    _c.c_int32,
                                    _c.c_void_p,
                                    _c.c_void_p,
                                    _c.c_void_p,
                                    _c.c_void_p]
    lib.init_qmc_solver.restype = None
    lib.get_qmc_aborted.argtypes = None
    lib.get_qmc_aborted.restype = _c.c_int32
    lib.init_qmc_parameters.argtypes = [_c.c_int32,
                                        _c.c_char_p]
    lib.init_qmc_parameters.restype = None
    lib.dest_qmc_paras.argtypes = None
    lib.dest_qmc_paras.restype = None
    lib.dest_ctqmc.argtypes = None
    lib.dest_ctqmc.restype = None
    lib.get_qmc_config_size.argtypes = None
    lib.get_qmc_config_size.restype = _c.c_int32
    lib.get_qmc_config.argtypes = [_c.c_int32,
                                   _c.POINTER(_c.c_double),
                                   _c.POINTER(_c.c_int32),
                                   _c.POINTER(_c.c_int32),
                                   _c.POINTER(_c.c_int32),
                                   _c.POINTER(_c.c_int32)]
    lib.get_qmc_config.restype = _c.c_int32
    lib.set_qmc_config.argtypes = [_c.c_int32,
                                   _c.POINTER(_c.c_double),
                                   _c.POINTER(_c.c_int32),
                                   _c.POINTER(_c.c_int32),
                                   _c.POINTER(_c.c_int32),
                                   _c.POINTER(_c.c_int32),
                                   _c.c_int32]
    lib.set_qmc_config.restype = None
    lib.set_empty_qmc_config.argtypes = [_c.c_int32]
    lib.set_empty_qmc_config.restype = None
    # driver routines
    lib.ctqmc_calibrate.argtypes = [_c.c_bool,
                                    _c.c_int64, _c.c_int64,
                                    _c.c_int64, _c.c_int64,
                                    _c.c_void_p]
    lib.ctqmc_calibrate.restype = None
    lib.ctqmc_warmup.argtypes = [_c.c_int64]
    lib.ctqmc_warmup.restype = None
    lib.ctqmc_worm_warmup.argtypes = [_c.c_int64, _c.c_int32, _c.c_int32]
    lib.ctqmc_worm_warmup.restype = None
    lib.count_qmc_sector_steps.argtypes = [_c.c_int32, _c.c_int32,
                                           _c.c_int64, _c.c_int64,
                                           _c.POINTER(_c.c_int64),
                                           _c.POINTER(_c.c_int64)]
    lib.count_qmc_sector_steps.restype = None
    lib.ctqmc_measure.argtypes = [_c.c_int32, _c.c_int32, _c.c_double]
    lib.ctqmc_measure.restype = None
    # measurements
    lib.get_leftovers_into_ndarrays_FIXME.argtypes = [_c.c_int32, _c.c_char_p,
                                                      _c.c_int32, _c.c_void_p]
    lib.get_leftovers_into_ndarrays_FIXME.restype = None
    lib.get_sector_meas_count.argtypes = [_c.c_int32]
    lib.get_sector_meas_count.restype = _c.c_int64
    lib.get_sector_sample_count.argtypes = [_c.c_int32]
    lib.get_sector_sample_count.restype = _c.c_int64
    lib.get_global_meastargets_FIXME.argtypes = None
    lib.get_global_meastargets_FIXME.restype = _c.c_void_p
    lib.add_sampler_G.argtypes = [_c.c_void_p, _c.c_void_p, _c.c_void_p]
    lib.add_sampler_G.restype = None
    lib.add_sampler_G2iw.argtypes = [_c.c_void_p, _c.c_void_p, _c.c_void_p]
    lib.add_sampler_G2iw.restype = None
    lib.add_sampler_wormP2ph.argtypes = [_c.c_void_p, _c.c_void_p, _c.c_void_p]
    lib.add_sampler_wormP2ph.restype = None
    lib.add_sampler_wormP2pp.argtypes = [_c.c_void_p, _c.c_void_p, _c.c_void_p]
    lib.add_sampler_wormP2pp.restype = None
    lib.add_sampler_wormQQ.argtypes = [_c.c_void_p, _c.c_void_p, _c.c_void_p]
    lib.add_sampler_wormQQ.restype = None
    lib.add_sampler2d_wormRaman.argtypes = [_c.c_void_p, _c.c_void_p, _c.c_void_p]
    lib.add_sampler2d_wormRaman.restype = None
    lib.add_sampler1d_wormCustom.argtypes = [_c.c_void_p,
                                             _c.c_void_p,
                                             _c.c_void_p]
    lib.add_sampler1d_wormCustom.restype = None
    lib.add_sampler2d_wormCustom.argtypes = [_c.c_void_p,
                                             _c.c_void_p,
                                             _c.c_void_p]
    lib.add_sampler2d_wormCustom.restype = None
    lib.add_sampler3d_wormCustom.argtypes = [_c.c_void_p,
                                             _c.c_void_p,
                                             _c.c_void_p]
    lib.add_sampler3d_wormCustom.restype = None
    lib.set_accumulator_G4iw_full.argtypes = [_c.c_void_p, _c.c_void_p, _c.POINTER(_c.c_int64), _c.c_int64, _c.c_int64, _c.c_double]
    lib.set_accumulator_G4iw_full.restype = None
    lib.set_accumulator_sign.argtypes = [_c.c_void_p, _c.c_void_p]
    lib.set_accumulator_sign.restype = None
    lib.set_accumulator_dm.argtypes = [_c.c_void_p, _c.c_void_p]
    lib.set_accumulator_dm.restype = None
    lib.set_accumulator_expord.argtypes = [_c.c_void_p, _c.c_void_p]
    lib.set_accumulator_expord.restype = None
    lib.set_accumulator_ntaun0.argtypes = [_c.c_void_p, _c.c_void_p]
    lib.set_accumulator_ntaun0.restype = None
    # custom worm
    lib.get_wormstate.argtypes = None
    lib.get_wormstate.restype = _c.c_void_p
    lib.set_custom_worm.argtypes = [_c.c_void_p,
                                    _c.c_int,
                                    _c.c_void_p,
                                    _c.c_void_p]
    lib.set_custom_worm.restype = None


_annotate_ctqmc(_lib)


def init_qmc_rng(seed):
    _lib.init_qmc_rng(_c.c_int32(seed))


def reset_qmc_counters():
    _lib.reset_qmc_counters()


def set_fancy_qmc_progress_bar(fancy):
    _lib.set_fancy_progress_bar(_c.c_bool(fancy))


def set_qmc_simid(id):
    _lib.set_qmc_simid(_c.c_int32(id))


def set_wormeta(eta, sector=None):
    _lib.set_sector_eta(_c.c_int32(sector)
                        if sector is not None
                        else _c.c_int32(-1),
                        _c.c_double(eta))


def get_wormeta(sector):
    return _lib.get_sector_eta(_c.c_int32(sector))


def set_qmc_taudiffmax(taudiffmax):
    _lib.set_qmc_taudiffmax(_c.c_double(taudiffmax))


def init_qmc_solver(umat, ftau, muimp, screening):
    nbands = umat.shape[0] // 2
    nftau = ftau.shape[-1]

    assert umat.shape == (2 * nbands, 2 * nbands, 2 * nbands, 2 * nbands)
    assert ftau.shape == (nbands, 2, nbands, 2, nftau)
    assert muimp.shape == (nbands, 2, nbands, 2)
    assert screening.shape == (nbands, 2, nbands, 2, nftau)

    _lib.init_qmc_solver(_c.c_int32(nbands), _c.c_int32(nftau),
                         np.asfortranarray(umat, np.cdouble).ctypes
                         .data_as(_c.c_void_p),
                         np.asfortranarray(ftau, np.cdouble).ctypes
                         .data_as(_c.c_void_p),
                         np.asfortranarray(muimp, np.cdouble).ctypes
                         .data_as(_c.c_void_p),
                         np.asfortranarray(screening, np.double).ctypes
                         .data_as(_c.c_void_p))


def get_qmc_config():
    size = _lib.get_qmc_config_size()
    taus = np.zeros((size,), dtype=np.double)
    orbs = np.zeros((size,), dtype=np.int32)
    spins = np.zeros((size,), dtype=np.int32)
    cas = np.zeros((size,), dtype=np.int32)
    hashybs = np.zeros((size,), dtype=np.int32)
    outer_sst = _lib.get_qmc_config(_c.c_int32(size),
                                    taus.ctypes
                                    .data_as(_c.POINTER(_c.c_double)),
                                    orbs.ctypes
                                    .data_as(_c.POINTER(_c.c_int32)),
                                    spins.ctypes
                                    .data_as(_c.POINTER(_c.c_int32)),
                                    cas.ctypes
                                    .data_as(_c.POINTER(_c.c_int32)),
                                    hashybs.ctypes
                                    .data_as(_c.POINTER(_c.c_int32)))
    return taus, orbs, spins, cas, hashybs, outer_sst


def set_qmc_config(outer_sst, taus=None, orbs=None, spins=None, cas=None,
                   hashybs=None):
    if taus is None or taus.size == 0:
        return _lib.set_empty_qmc_config(_c.c_int32(outer_sst))

    assert taus.size == orbs.size == spins.size == cas.size == hashybs.size
    _lib.set_qmc_config(_c.c_int32(taus.size),
                        taus.ctypes
                        .data_as(_c.POINTER(_c.c_double)),
                        orbs.ctypes
                        .data_as(_c.POINTER(_c.c_int32)),
                        spins.ctypes
                        .data_as(_c.POINTER(_c.c_int32)),
                        cas.ctypes
                        .data_as(_c.POINTER(_c.c_int32)),
                        hashybs.ctypes
                        .data_as(_c.POINTER(_c.c_int32)),
                        _c.c_int32(outer_sst))


def get_qmc_aborted():
    return _lib.get_qmc_aborted()


def init_qmc_parameters(paramstring):
    pstr_bytes = paramstring.encode()
    _lib.init_qmc_parameters(_c.c_int32(len(pstr_bytes)),
                             _c.c_char_p(pstr_bytes))


def dest_qmc_paras():
    _lib.dest_qmc_paras()


def dest_ctqmc():
    _lib.dest_ctqmc()


def run_ctqmc_calibration(aslist, phase1pairnum,
                          phase2taugrid, phase2pairnum, nprogress, accpair):
    _lib.ctqmc_calibrate(_c.c_bool(aslist), _c.c_int64(phase1pairnum),
                         _c.c_int64(phase2taugrid), _c.c_int64(phase2pairnum),
                         _c.c_int64(nprogress),
                         accpair.ctypes.data_as(_c.c_void_p))

def run_ctqmc_warmup(nwarmups):
    _lib.ctqmc_warmup(_c.c_int64(nwarmups))


def run_ctqmc_worm_warmup(nwarmups, sector, component):
    _lib.ctqmc_worm_warmup(_c.c_int64(nwarmups), _c.c_int32(sector),
                           _c.c_int32(component))


def run_ctqmc_measurements(sector, component, maxtime):
    _lib.ctqmc_measure(_c.c_int32(sector), _c.c_int32(component),
                       _c.c_double(maxtime))


def count_qmc_sector_steps(sector, component, nfix, nwarmup):
    steps_worm = _c.c_int64()
    steps_z = _c.c_int64()
    _lib.count_qmc_sector_steps(_c.c_int32(sector),
                                _c.c_int32(component),
                                _c.c_int64(nfix),
                                _c.c_int64(nwarmup),
                                _c.byref(steps_worm),
                                _c.byref(steps_z))
    return steps_worm.value, steps_z.value


def get_FIXME(name, shape, dtype):
    array = np.zeros(shape, dtype=dtype, order='F')
    name_bytes = name.encode()
    _lib.get_leftovers_into_ndarrays_FIXME(_c.c_int32(len(name_bytes)),
                                           _c.c_char_p(name_bytes),
                                           _c.c_int32(array.size),
                                           array.ctypes.data_as(_c.c_void_p))
    return array


def get_sector_meas_count(sector):
    return _lib.get_sector_meas_count(_c.c_int32(sector))


def get_sector_sample_count(sector):
    return _lib.get_sector_sample_count(_c.c_int32(sector))


def get_global_meastargets_FIXME():
    return _lib.get_global_meastargets_FIXME()


def register_sampler_G(targets, samp, acc):
    _lib.add_sampler_G(targets, samp._handle, acc._acc)
    
def register_sampler_G2iw(targets, samp, acc):
    _lib.add_sampler_G2iw(targets, samp._handle, acc._acc)


def register_sampler_wormP2ph(targets, samp, acc):
    _lib.add_sampler_wormP2ph(targets, samp._handle, acc._acc)


def register_sampler_wormP2pp(targets, samp, acc):
    _lib.add_sampler_wormP2pp(targets, samp._handle, acc._acc)


def register_sampler_wormQQ(targets, samp, acc):
    _lib.add_sampler_wormQQ(targets, samp._handle, acc._acc)


def register_sampler2d_wormRaman(targets, samp, acc):
    _lib.add_sampler2d_wormRaman(targets, samp._handle, acc._acc)


def register_sampler1d_wormCustom(targets, samp, acc):
    _lib.add_sampler1d_wormCustom(targets, samp._handle, acc._acc)


def register_sampler2d_wormCustom(targets, samp, acc):
    _lib.add_sampler2d_wormCustom(targets, samp._handle, acc._acc)


def register_sampler3d_wormCustom(targets, samp, acc):
    _lib.add_sampler3d_wormCustom(targets, samp._handle, acc._acc)


def register_sampler_wormCustom(targets, samp, acc, ntaus):
    if ntaus == 2:
        register_sampler1d_wormCustom(targets, samp, acc)
    elif ntaus == 3:
        register_sampler2d_wormCustom(targets, samp, acc)
    elif ntaus == 4:
        register_sampler3d_wormCustom(targets, samp, acc)
    else:
        raise NotImplementedError("Sampler for ntaus > 4 not supported")
    

def set_accumulator_G4iw_full(targets, acc, ngrid, niw, nbands, beta):
    _lib.set_accumulator_G4iw_full(targets, acc._acc, ngrid, niw, nbands, beta)


def set_accumulator_sign(targets, acc):
    _lib.set_accumulator_sign(targets, acc._acc)


def set_accumulator_dm(targets, acc):
    _lib.set_accumulator_dm(targets, acc._acc)


def set_accumulator_expord(targets, acc):
    _lib.set_accumulator_expord(targets, acc._acc)


def set_accumulator_ntaun0(targets, acc):
    _lib.set_accumulator_ntaun0(targets, acc._acc)


def get_wormstate():
    return _lib.get_wormstate()


def set_custom_worm(wormstate,
                    optypes : Sequence[Sequence[OP_W]]):
    if not all(t in OP_W
               for o in optypes
               for t in o):
        raise ValueError("Invalid worm operator type in custom worm")

    ntaus = len(optypes)
    c_ntauops = (_c.c_int * ntaus)(*[len(o) for o in optypes])

    ntotalops = sum(len(o) for o in optypes)
    c_optypes = (_c.c_int * ntotalops)(*[int(t) for o in optypes for t in o])

    _lib.set_custom_worm(
        wormstate,
        ntaus,
        c_ntauops,
        c_optypes
    )
