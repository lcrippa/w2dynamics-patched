!> Module providing C wrappers for the samplers
module csamplers
    use MSamplerD
    use MSamplerZ
    use MSamplerRegister
    use iso_c_binding
    implicit none

    ! Boxing type that allows propagation of polymorphism information
    type csampler
        class(ZSampler), allocatable :: inner
    end type

    type csampler2d
        class(ZSampler2D), allocatable :: inner
    end type csampler2d

    type csampler3d
        class(ZSampler3D), allocatable :: inner
    end type csampler3d
contains
    function crhist_create(ncomp, beta, statsign, nbins, coeff) result(this) &
                                    bind(C, name='rhist_create')
        type(c_ptr) :: this
        integer(c_int64_t), value :: ncomp, statsign, nbins
        real(c_double), value :: beta, coeff

        type(csampler), pointer :: thisf

        allocate(thisf)
        call init_binning_sampler( &
                thisf%inner, int(ncomp), beta, int(statsign), int(nbins), coeff)
        this = c_loc(thisf)
    end function

    function crhist2d_create(ncomp, beta, statsign1, statsign2, nbins1, nbins2, &
                             coeff11, coeff21, coeff12, coeff22) &
                                    result(this) &
                                    bind(C, name='rhist2d_create')
        type(c_ptr) :: this
        integer(c_int64_t), value :: ncomp, statsign1, statsign2, nbins1, nbins2
        real(c_double), value :: beta, coeff11, coeff21, coeff12, coeff22

        type(csampler2d), pointer :: thisf

        allocate(thisf)
        call init_binning_sampler( &
                thisf%inner, int(ncomp), beta, &
                [int(statsign1), int(statsign2)], &
                [int(nbins1), int(nbins2)], &
                reshape([coeff11, coeff21, coeff12, coeff22], [2, 2]))
        this = c_loc(thisf)
    end function crhist2d_create

    function crhist3d_create(ncomp, beta, statsign1, statsign2, statsign3, &
                             nbins1, nbins2, nbins3, &
                             coeff11, coeff21, coeff31, &
                             coeff12, coeff22, coeff32, &
                             coeff13, coeff23, coeff33) &
                             result(this) &
                                    bind(C, name='rhist3d_create')
        type(c_ptr) :: this
        integer(c_int64_t), value :: ncomp, statsign1, statsign2, statsign3
        integer(c_int64_t), value :: nbins1, nbins2, nbins3
        real(c_double), value :: beta, coeff11, coeff21, coeff31
        real(c_double), value :: coeff12, coeff22, coeff32
        real(c_double), value :: coeff13, coeff23, coeff33

        type(csampler3d), pointer :: thisf

        allocate(thisf)
        call init_binning_sampler( &
                thisf%inner, int(ncomp), beta, &
                [int(statsign1), int(statsign2), int(statsign3)], &
                [int(nbins1), int(nbins2), int(nbins3)], &
                reshape([coeff11, coeff21, coeff31, &
                         coeff12, coeff22, coeff32, &
                         coeff13, coeff23, coeff33], [3, 3]))
        this = c_loc(thisf)
    end function crhist3d_create

    function cnfourier_create(ncomp, beta, statsign, niw, sparse_frequencies, coeff, use_sparse) result(this) &
                                    bind(C, name='nfourier_create')
        type(c_ptr) :: this
        integer(c_int64_t), value :: ncomp, statsign, niw
        real(c_double), value :: beta, coeff
        integer(c_int64_t), intent(in) :: sparse_frequencies(niw)
        logical(c_bool), value :: use_sparse

        type(csampler), pointer :: thisf

        allocate(thisf)
        call init_matsubara_sampler( &
                thisf%inner, int(ncomp), beta, int(statsign), int(niw), int(sparse_frequencies), coeff, logical(use_sparse))
        this = c_loc(thisf)
    end function cnfourier_create

    function cnfourier2d_create(ncomp, beta, statsign1, statsign2, niw1, niw2, coeff11, coeff21, coeff12, coeff22) result(this) &
        bind(C, name='nfourier2d_create')
        type(c_ptr) :: this
        integer(c_int64_t), value :: ncomp, statsign1, statsign2, niw1, niw2
        real(c_double), value :: beta, coeff11, coeff21, coeff12, coeff22

        type(csampler2d), pointer :: thisf

        allocate(thisf)
        call init_matsubara_sampler( &
        thisf%inner, int(ncomp), beta, [int(statsign1), int(statsign2)], &
        [int(niw1), int(niw2)], &
        reshape([coeff11, coeff21, coeff12, coeff22], [2, 2]))
        this = c_loc(thisf)
    end function cnfourier2d_create

    function cnfft_create(ncomp, beta, statsign, niw, sparse_frequencies, coeff, use_fftw, use_sparse) result(this) &
                                    bind(C, name='nfft_create')
        type(c_ptr) :: this
        integer(c_int64_t), value :: ncomp, statsign, niw
        real(c_double), value :: beta, coeff
        logical(c_bool), value :: use_fftw
        integer(c_int64_t), intent(inout) :: sparse_frequencies(niw)
        logical(c_bool), value :: use_sparse

        type(csampler), pointer :: thisf

        allocate(thisf)
        call init_fast_matsubara_sampler( &
                thisf%inner, int(ncomp), beta, int(statsign), int(niw), int(sparse_frequencies), coeff, &
                logical(use_fftw), logical(use_sparse))
        this = c_loc(thisf)
    end function cnfft_create

    function cnfft2d_create(ncomp, beta, statsign1, statsign2, niw1, niw2, &
                            coeff11, coeff21, coeff12, coeff22, &
                            use_fftw) result(this) &
                                      bind(C, name='nfft2d_create')
        type(c_ptr) :: this
        integer(c_int64_t), value :: ncomp, statsign1, statsign2, niw1, niw2
        real(c_double), value :: beta, coeff11, coeff21, coeff12, coeff22
        logical(c_bool), value :: use_fftw

        type(csampler2d), pointer :: thisf

        allocate(thisf)
        call init_fast_matsubara_sampler( &
                thisf%inner, int(ncomp), beta, &
                [int(statsign1), int(statsign2)], &
                [int(niw1), int(niw2)], &
                reshape([coeff11, coeff21, coeff12, coeff22], [2, 2]), &
                logical(use_fftw))
        this = c_loc(thisf)
    end function cnfft2d_create

    function cnfft3d_create(ncomp, beta, statsign1, statsign2, statsign3, &
                            niw1, niw2, niw3, &
                            coeff11, coeff21, coeff31, &
                            coeff12, coeff22, coeff32, &
                            coeff13, coeff23, coeff33, use_fftw) &
                            result(this) &
                            bind(C, name='nfft3d_create')
        type(c_ptr) :: this
        integer(c_int64_t), value :: ncomp, statsign1, statsign2, statsign3
        integer(c_int64_t), value :: niw1, niw2, niw3
        real(c_double), value :: beta, coeff11, coeff21, coeff31
        real(c_double), value :: coeff12, coeff22, coeff32
        real(c_double), value :: coeff13, coeff23, coeff33
        logical(c_bool), value :: use_fftw

        type(csampler3d), pointer :: thisf

        allocate(thisf)
        call init_fast_matsubara_sampler( &
                thisf%inner, int(ncomp), beta, &
                [int(statsign1), int(statsign2), int(statsign3)], &
                [int(niw1), int(niw2), int(niw3)], &
                reshape([coeff11, coeff21, coeff31, &
                         coeff12, coeff22, coeff32, &
                         coeff13, coeff23, coeff33], [3, 3]), &
                logical(use_fftw))
        this = c_loc(thisf)
    end function cnfft3d_create

    function clegendre_create(ncomp, beta, statsign, nleg, coeff) result(this) &
                                    bind(C, name='legendre_create')
        type(c_ptr) :: this
        integer(c_int64_t), value :: ncomp, statsign, nleg
        real(c_double), value :: beta, coeff

        type(csampler), pointer :: thisf

        allocate(thisf)
        call init_legendre_sampler( &
                thisf%inner, int(ncomp), beta, int(statsign), int(nleg), coeff)
        this = c_loc(thisf)
    end function

    function csampler_ncomp(this) result(ncomp) bind(C, name='sampler_ncomp')
        type(c_ptr), value :: this
        integer(c_int64_t) :: ncomp

        type(csampler), pointer :: thisf

        call c_f_pointer(this, thisf)
        ncomp = thisf%inner%ncomp
    end function

    function csampler_nacc(this) result(nacc) bind(C, name='sampler_nacc')
        type(c_ptr), value :: this
        integer(c_int64_t) :: nacc

        type(csampler), pointer :: thisf

        call c_f_pointer(this, thisf)
        nacc = thisf%inner%nacc
    end function

    function csampler_npost(this) result(npost) bind(C, name='sampler_npost')
        type(c_ptr), value :: this
        integer(c_int64_t) :: npost

        type(csampler), pointer :: thisf

        call c_f_pointer(this, thisf)
        npost = thisf%inner%npost
    end function

    function csampler_statsign(this) result(statsign) &
                                    bind(C, name='sampler_statsign')
        type(c_ptr), value :: this
        integer(c_int64_t) :: statsign

        type(csampler), pointer :: thisf

        call c_f_pointer(this, thisf)
        statsign = thisf%inner%statsign
    end function

    function csampler_beta(this) result(beta) bind(C, name='sampler_beta')
        type(c_ptr), value :: this
        real(c_double) :: beta

        type(csampler), pointer :: thisf

        call c_f_pointer(this, thisf)
        beta = thisf%inner%beta
    end function

    subroutine csampler_add(this, n, i, x, fx, acc) bind(C, name='sampler_add')
        type(c_ptr), value :: this
        integer(c_int64_t), value :: n
        integer(c_int64_t), intent(in) :: i(n)
        real(c_double), intent(in) :: x(n)
        complex(c_double_complex), intent(in) :: fx(n)
        complex(c_double_complex), target, intent(inout) :: acc(*)

        type(csampler), pointer :: thisf
        type(ZSamplePoint) :: sample(n)
        complex(c_double_complex), pointer :: acc_inner(:,:)

        call c_f_pointer(this, thisf)

        sample(:)%i = i
        sample(:)%x = x
        sample(:)%fx = fx
        acc_inner(1:thisf%inner%nacc, 1:thisf%inner%ncomp) => &
                acc(1:thisf%inner%nacc * thisf%inner%ncomp)
        call sampler_add(thisf%inner, sample, acc_inner)
    end subroutine

    subroutine csampler_delete(this) bind(C, name='sampler_delete')
        type(c_ptr), value :: this

        type(csampler), pointer :: thisf

        call c_f_pointer(this, thisf)
        call sampler_delete(thisf%inner)
        deallocate (thisf)
    end subroutine

    subroutine csampler_postprocess(this, acc, res) &
                                    bind(C, name='sampler_postprocess')
        type(c_ptr), value :: this
        complex(c_double_complex), target, intent(in) :: acc(*)
        complex(c_double_complex), target, intent(out) :: res(*)

        type(csampler), pointer :: thisf
        complex(c_double_complex), pointer :: res_inner(:,:)
        complex(c_double_complex), pointer :: acc_inner(:,:)

        call c_f_pointer(this, thisf)
        acc_inner(1:thisf%inner%nacc, 1:thisf%inner%ncomp) => &
                acc(1:thisf%inner%nacc * thisf%inner%ncomp)
        res_inner(1:thisf%inner%npost, 1:thisf%inner%ncomp) => &
                res(1:thisf%inner%npost * thisf%inner%ncomp)
        call sampler_postprocess(thisf%inner, acc_inner, res_inner)
    end subroutine

    subroutine csampler_register_register(creg, csamp, cacc) bind(C, name='sampler_register_register')
        use caccumulators
        type(c_ptr), value :: creg
        type(c_ptr), value :: csamp
        type(c_ptr), value :: cacc

        type(sampler_register), pointer :: reg
        type(csampler), pointer :: ZSampler
        type(TAccumulator), pointer :: acc

        call c_f_pointer(creg, reg)
        call c_f_pointer(csamp, ZSampler)
        call c_f_pointer(cacc, acc)

        call reg%register(ZSampler%inner, acc)
    end subroutine

    function csampler2d_ncomp(this) result(ncomp) bind(C, name='sampler2d_ncomp')
        type(c_ptr), value :: this
        integer(c_int64_t) :: ncomp

        type(csampler2d), pointer :: thisf

        call c_f_pointer(this, thisf)
        ncomp = thisf%inner%ncomp
    end function csampler2d_ncomp

    function csampler2d_nacc(this) result(nacc) bind(C, name='sampler2d_nacc')
        type(c_ptr), value :: this
        integer(c_int64_t) :: nacc

        type(csampler2d), pointer :: thisf

        call c_f_pointer(this, thisf)
        nacc = thisf%inner%nacc
    end function csampler2d_nacc

    function csampler2d_npost(this) result(npost) bind(C, name='sampler2d_npost')
        type(c_ptr), value :: this
        integer(c_int64_t) :: npost

        type(csampler2d), pointer :: thisf

        call c_f_pointer(this, thisf)
        npost = thisf%inner%npost
    end function csampler2d_npost

    function csampler2d_statsign1(this) result(statsign) &
                                    bind(C, name='sampler2d_statsign1')
        type(c_ptr), value :: this
        integer(c_int64_t) :: statsign

        type(csampler2d), pointer :: thisf

        call c_f_pointer(this, thisf)
        statsign = thisf%inner%statsign(1)
    end function csampler2d_statsign1

    function csampler2d_statsign2(this) result(statsign) &
                                    bind(C, name='sampler2d_statsign2')
        type(c_ptr), value :: this
        integer(c_int64_t) :: statsign

        type(csampler2d), pointer :: thisf

        call c_f_pointer(this, thisf)
        statsign = thisf%inner%statsign(2)
    end function csampler2d_statsign2

    function csampler2d_beta(this) result(beta) bind(C, name='sampler2d_beta')
        type(c_ptr), value :: this
        real(c_double) :: beta

        type(csampler2d), pointer :: thisf

        call c_f_pointer(this, thisf)
        beta = thisf%inner%beta
    end function csampler2d_beta

    subroutine csampler2d_add(this, n, i, x, y, fx, acc) bind(C, name='sampler2d_add')
        type(c_ptr), value :: this
        integer(c_int64_t), value :: n
        integer(c_int64_t), intent(in) :: i(n)
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in) :: y(n)
        complex(c_double_complex), intent(in) :: fx(n)
        complex(c_double_complex), target, intent(inout) :: acc(*)

        type(csampler2d), pointer :: thisf
        type(ZSamplePoint2D) :: sample(n)
        complex(c_double_complex), pointer :: acc_inner(:,:)

        call c_f_pointer(this, thisf)

        sample(:)%i = i
        sample(:)%x(1) = x
        sample(:)%x(2) = x
        sample(:)%fx = fx
        acc_inner(1:thisf%inner%nacc, 1:thisf%inner%ncomp) => &
                acc(1:thisf%inner%nacc * thisf%inner%ncomp)
        call sampler_add(thisf%inner, sample, acc_inner)
    end subroutine csampler2d_add

    subroutine csampler2d_delete(this) bind(C, name='sampler2d_delete')
        type(c_ptr), value :: this

        type(csampler2d), pointer :: thisf

        call c_f_pointer(this, thisf)
        call sampler_delete(thisf%inner)
        deallocate (thisf)
    end subroutine csampler2d_delete

    subroutine csampler2d_postprocess(this, acc, res) &
                                    bind(C, name='sampler2d_postprocess')
        type(c_ptr), value :: this
        complex(c_double_complex), target, intent(in) :: acc(*)
        complex(c_double_complex), target, intent(out) :: res(*)

        type(csampler2d), pointer :: thisf
        complex(c_double_complex), pointer :: res_inner(:,:)
        complex(c_double_complex), pointer :: acc_inner(:,:)

        call c_f_pointer(this, thisf)
        acc_inner(1:thisf%inner%nacc, 1:thisf%inner%ncomp) => &
                acc(1:thisf%inner%nacc * thisf%inner%ncomp)
        res_inner(1:thisf%inner%npost, 1:thisf%inner%ncomp) => &
                res(1:thisf%inner%npost * thisf%inner%ncomp)
        call sampler_postprocess(thisf%inner, acc_inner, res_inner)
    end subroutine csampler2d_postprocess

    subroutine csampler2d_register_register(creg, csamp, cacc) bind(C, name='sampler2d_register_register')
        use caccumulators
        type(c_ptr), value :: creg
        type(c_ptr), value :: csamp
        type(c_ptr), value :: cacc

        type(sampler2d_register), pointer :: reg
        type(csampler2d), pointer :: ZSampler2d
        type(TAccumulator), pointer :: acc

        call c_f_pointer(creg, reg)
        call c_f_pointer(csamp, ZSampler2d)
        call c_f_pointer(cacc, acc)

        call reg%register(ZSampler2d%inner, acc)
    end subroutine csampler2d_register_register

    function csampler3d_ncomp(this) result(ncomp) bind(C, name='sampler3d_ncomp')
        type(c_ptr), value :: this
        integer(c_int64_t) :: ncomp

        type(csampler3d), pointer :: thisf

        call c_f_pointer(this, thisf)
        ncomp = thisf%inner%ncomp
    end function csampler3d_ncomp

    function csampler3d_nacc(this) result(nacc) bind(C, name='sampler3d_nacc')
        type(c_ptr), value :: this
        integer(c_int64_t) :: nacc

        type(csampler3d), pointer :: thisf

        call c_f_pointer(this, thisf)
        nacc = thisf%inner%nacc
    end function csampler3d_nacc

    function csampler3d_npost(this) result(npost) bind(C, name='sampler3d_npost')
        type(c_ptr), value :: this
        integer(c_int64_t) :: npost

        type(csampler3d), pointer :: thisf

        call c_f_pointer(this, thisf)
        npost = thisf%inner%npost
    end function csampler3d_npost

    function csampler3d_statsign1(this) result(statsign) &
                                    bind(C, name='sampler3d_statsign1')
        type(c_ptr), value :: this
        integer(c_int64_t) :: statsign

        type(csampler3d), pointer :: thisf

        call c_f_pointer(this, thisf)
        statsign = thisf%inner%statsign(1)
    end function csampler3d_statsign1

    function csampler3d_statsign2(this) result(statsign) &
                                    bind(C, name='sampler3d_statsign2')
        type(c_ptr), value :: this
        integer(c_int64_t) :: statsign

        type(csampler3d), pointer :: thisf

        call c_f_pointer(this, thisf)
        statsign = thisf%inner%statsign(2)
    end function csampler3d_statsign2

    function csampler3d_statsign3(this) result(statsign) &
                                    bind(C, name='sampler3d_statsign3')
        type(c_ptr), value :: this
        integer(c_int64_t) :: statsign

        type(csampler3d), pointer :: thisf

        call c_f_pointer(this, thisf)
        statsign = thisf%inner%statsign(3)
    end function csampler3d_statsign3

    function csampler3d_beta(this) result(beta) bind(C, name='sampler3d_beta')
        type(c_ptr), value :: this
        real(c_double) :: beta

        type(csampler3d), pointer :: thisf

        call c_f_pointer(this, thisf)
        beta = thisf%inner%beta
    end function csampler3d_beta

    subroutine csampler3d_add(this, n, i, x, y, z, fx, acc) bind(C, name='sampler3d_add')
        type(c_ptr), value :: this
        integer(c_int64_t), value :: n
        integer(c_int64_t), intent(in) :: i(n)
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in) :: y(n)
        real(c_double), intent(in) :: z(n)
        complex(c_double_complex), intent(in) :: fx(n)
        complex(c_double_complex), target, intent(inout) :: acc(*)

        type(csampler3d), pointer :: thisf
        type(ZSamplePoint3D) :: sample(n)
        complex(c_double_complex), pointer :: acc_inner(:,:)

        call c_f_pointer(this, thisf)

        sample(:)%i = i
        sample(:)%x(1) = x
        sample(:)%x(2) = y
        sample(:)%x(3) = z
        sample(:)%fx = fx
        acc_inner(1:thisf%inner%nacc, 1:thisf%inner%ncomp) => &
                acc(1:thisf%inner%nacc * thisf%inner%ncomp)
        call sampler_add(thisf%inner, sample, acc_inner)
    end subroutine csampler3d_add

    subroutine csampler3d_delete(this) bind(C, name='sampler3d_delete')
        type(c_ptr), value :: this

        type(csampler3d), pointer :: thisf

        call c_f_pointer(this, thisf)
        call sampler_delete(thisf%inner)
        deallocate (thisf)
    end subroutine csampler3d_delete

    subroutine csampler3d_postprocess(this, acc, res) &
                                    bind(C, name='sampler3d_postprocess')
        type(c_ptr), value :: this
        complex(c_double_complex), target, intent(in) :: acc(*)
        complex(c_double_complex), target, intent(out) :: res(*)

        type(csampler3d), pointer :: thisf
        complex(c_double_complex), pointer :: res_inner(:,:)
        complex(c_double_complex), pointer :: acc_inner(:,:)

        call c_f_pointer(this, thisf)
        acc_inner(1:thisf%inner%nacc, 1:thisf%inner%ncomp) => &
                acc(1:thisf%inner%nacc * thisf%inner%ncomp)
        res_inner(1:thisf%inner%npost, 1:thisf%inner%ncomp) => &
                res(1:thisf%inner%npost * thisf%inner%ncomp)
        call sampler_postprocess(thisf%inner, acc_inner, res_inner)
    end subroutine csampler3d_postprocess

    subroutine csampler3d_register_register(creg, csamp, cacc) bind(C, name='sampler3d_register_register')
        use caccumulators
        type(c_ptr), value :: creg
        type(c_ptr), value :: csamp
        type(c_ptr), value :: cacc

        type(sampler3d_register), pointer :: reg
        type(csampler3d), pointer :: ZSampler3d
        type(TAccumulator), pointer :: acc

        call c_f_pointer(creg, reg)
        call c_f_pointer(csamp, ZSampler3d)
        call c_f_pointer(cacc, acc)

        call reg%register(ZSampler3d%inner, acc)
    end subroutine csampler3d_register_register
end module
