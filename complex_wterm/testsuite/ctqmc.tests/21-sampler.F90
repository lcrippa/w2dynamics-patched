program test_sampler
    use iso_c_binding
    use MSamplerD
    use MSamplerZ
    use testing
    implicit none

    class(DSampler), allocatable :: dsmpl
    type(DSampleBuffer) :: dbuf
    class(ZSampler), allocatable :: zsmpl
    type(ZSampleBuffer) :: zbuf
    real(c_double), allocatable :: dacc(:,:), dout(:,:)
    complex(c_double_complex), allocatable :: zacc(:,:), zout(:,:)
    complex(c_double_complex), parameter :: ref_mats(10,1) = reshape((/ &
            ( -0.6084989969270205d0,   1.366534782335945d0), &
            ( -0.5831798434076642d0, -0.1374581611060059d0), &
            (  -1.207106781186547d0, -0.7071067811865475d0), &
            (  0.5609429622277047d0, -0.6937957144489013d0), &
            (  0.7364980270008942d0, -0.9295187578871236d0), &
            (  0.7364980270008942d0,  0.9295187578871236d0), &
            (  0.5609429622277047d0,  0.6937957144489013d0), &
            (  -1.207106781186547d0,  0.7071067811865475d0), &
            ( -0.5831798434076642d0,  0.1374581611060059d0), &
            ( -0.6084989969270205d0,  -1.366534782335945d0)  &
            /), (/ 10, 1 /))
    complex(c_double_complex), parameter :: zref_mats(8,1) = reshape((/ &
            ( +0.2061806668402426d0, -0.06207626681383845d0), &
            (+0.07235416581219539d0, -0.08848320191611111d0), &
            ( +0.2580439722882742d0,  -0.2029171945311391d0), &
            (  +0.196342989042541d0, -0.06779905550000695d0), &
            ( +0.2088029642716221d0,  -0.2593846835674751d0), &
            ( +0.2953629641944935d0,  -0.1235542044790143d0), &
            ( +0.1343498363421691d0,  -0.2365612142727918d0), &
            ( +0.2925506585775105d0,  -0.2608219861974008d0) &
            /), (/ 8, 1/))
    complex(c_double_complex), parameter :: ref_mats_sparse(8,1) = reshape((/ &
            ( -0.6084989969270199d0,   -1.366534782335945d0), &
            ( -0.6084989969270205d0,   1.366534782335945d0), &
            (  -1.207106781186547d0, -0.7071067811865475d0), &
            (  0.7364980270008942d0, -0.9295187578871236d0), &
            (  0.7364980270008942d0,  0.9295187578871236d0), &
            (  -1.207106781186547d0,  0.7071067811865475d0), &
            ( -0.6084989969270205d0,  -1.366534782335945d0), &
            ( -0.6084989969270199d0,   1.366534782335945d0) &
            /), (/ 8, 1 /))
    complex(c_double_complex), parameter :: zref_mats_sparse(6,1) = reshape((/ &
            (-0.14828482432423218d0, +0.19673380004367572d0), &
            ( +0.2061806668402426d0, -0.06207626681383845d0), &
            (  +0.196342989042541d0, -0.06779905550000695d0), &
            ( +0.2088029642716221d0,  -0.2593846835674751d0), &
            ( +0.2925506585775105d0,  -0.2608219861974008d0), &
            ( +0.1688605817660662d0, -0.4341416183976562d0) &
            /), (/ 6, 1/))
    integer, parameter :: sparse_frequ(8) = (/ -31, -9, -5, -1, 1, 5, 9, 31/)
    integer, parameter :: z_sparse_frequ(6) = (/-31, -7, -1, 1, 7, 31/)


    ! REAL
    call prepare_buffer(dbuf, 2)
    dbuf%buffer(1) = DSamplePoint(i=1, x=0.3d0, fx=2.0d0)
    dbuf%buffer(2) = DSamplePoint(i=1, x=-0.8d0, fx=-1.0d0)

    ! Binning sampler
    call init_binning_sampler(dsmpl, 1, 2.0d0, -1, 8, 1.0d0)
    call fill_sampler_dd()
    call assert_close(dout, &
            8/2.0d0**2 * reshape((/ 0, 2, 0, 0, 1, 0, 0, 0 /), (/ 8, 1 /)))
    call sampler_delete(dsmpl)

    ! Legendre sampler
    call init_legendre_sampler(dsmpl, 1, 2.0d0, -1, 4, 1.0d0)
    call fill_sampler_dd()
    call assert_close(dout(1,1), 1.5d0)
    call sampler_delete(dsmpl)

    ! Matsubara sampler
    call init_matsubara_sampler(dsmpl, 1, 2.0d0, -1, 5, sparse_frequ, 1.0d0, .false.)
    call fill_sampler_dz()
    call assert_close(zout(:,:), ref_mats, rtol=2d-14)
    call sampler_delete(dsmpl)

    call init_fast_matsubara_sampler(dsmpl, 1, 2.0d0, -1, 5, sparse_frequ, 1.0d0, .true., .false.)
    call fill_sampler_dz()
    call assert_close(zout(:,:), ref_mats, rtol=2d-14)
    call sampler_delete(dsmpl)

    call init_fast_matsubara_sampler(dsmpl, 1, 2.0d0, -1, 5, sparse_frequ, 1.0d0, .false., .false.)
    call fill_sampler_dz()
    call assert_close(zout(:,:), ref_mats, rtol=2d-14)
    call sampler_delete(dsmpl)

    ! sparse
    call init_matsubara_sampler(dsmpl, 1, 2.0d0, -1, 8, sparse_frequ, 1.0d0, .true.)
    call fill_sampler_dz()
    call assert_close(zout(:,:), ref_mats_sparse, rtol=2d-14)
    call sampler_delete(dsmpl)

    call init_fast_matsubara_sampler(dsmpl, 1, 2.0d0, -1, 8, sparse_frequ, 1.0d0, .true., .true.)
    call fill_sampler_dz()
    call assert_close(zout(:,:), ref_mats_sparse, rtol=2d-14)
    call sampler_delete(dsmpl)

    call init_fast_matsubara_sampler(dsmpl, 1, 2.0d0, -1, 8, sparse_frequ, 1.0d0, .false., .true.)
    call fill_sampler_dz()
    call assert_close(zout(:,:), ref_mats_sparse, rtol=2d-14)
    call sampler_delete(dsmpl)

    ! COMPLEX
    call prepare_buffer(zbuf, 4)
    zbuf%buffer(1) = ZSamplePoint(i=1, x=-10.9d0, fx=(1d0, 2d0))
    zbuf%buffer(2) = ZSamplePoint(i=1, x=-0.2d0, fx=(3d0, 0.2d0))
    zbuf%buffer(3) = ZSamplePoint(i=1, x=0.8d0, fx=(0.5d0, 0d0))
    zbuf%buffer(4) = ZSamplePoint(i=1, x=4.6d0, fx=(-1d0, 0d0))

    ! Matsubara sampler
    call init_matsubara_sampler(zsmpl, 1, 11.0d0, -1, 4, sparse_frequ, 1.0d0, .false.)
    call fill_sampler_zz()
    call assert_close(zout(:,:), zref_mats, rtol=2d-14)
    call sampler_delete(zsmpl)

    call init_fast_matsubara_sampler(zsmpl, 1, 11.0d0, -1, 4, sparse_frequ, 1.0d0, .true., .false.)
    call fill_sampler_zz()
    call assert_close(zout(:,:), zref_mats, rtol=2d-14)
    call sampler_delete(zsmpl)

    call init_fast_matsubara_sampler(zsmpl, 1, 11.0d0, -1, 4, sparse_frequ, 1.0d0, .false., .false.)
    call fill_sampler_zz()
    call assert_close(zout(:,:), zref_mats, rtol=2d-14)
    call sampler_delete(zsmpl)

    ! sparse
    call init_matsubara_sampler(zsmpl, 1, 11.0d0, -1, 6, z_sparse_frequ, 1.0d0, .true.)
    call fill_sampler_zz()
    call assert_close(zout(:,:), zref_mats_sparse, rtol=2d-14)
    call sampler_delete(zsmpl)

    call init_fast_matsubara_sampler(zsmpl, 1, 11.0d0, -1, 6, z_sparse_frequ, 1.0d0, .true., .true.)
    call fill_sampler_zz()
    call assert_close(zout(:,:), zref_mats_sparse, rtol=2d-14)
    call sampler_delete(zsmpl)

    call init_fast_matsubara_sampler(zsmpl, 1, 11.0d0, -1, 6, z_sparse_frequ, 1.0d0, .false., .true.)
    call fill_sampler_zz()
    call assert_close(zout(:,:), zref_mats_sparse, rtol=2d-14)
    call sampler_delete(zsmpl)

contains

    subroutine fill_sampler_dd()
        if (allocated(dacc)) &
            deallocate(dacc, dout)
        allocate(dacc(dsmpl%nacc, dsmpl%ncomp), dout(dsmpl%npost, dsmpl%ncomp))
        dacc(:, :) = 0
        dout(:, :) = 0
        call sampler_add(dsmpl, dbuf%buffer(1:dbuf%n), dacc)
        call sampler_postprocess(dsmpl, dacc, dout)
    end subroutine

    subroutine fill_sampler_dz()
        if (allocated(zacc)) &
            deallocate(zacc, zout)
        allocate(zacc(dsmpl%nacc, dsmpl%ncomp), zout(dsmpl%npost, dsmpl%ncomp))
        zacc(:, :) = 0
        zout(:, :) = 0
        call sampler_add(dsmpl, dbuf%buffer(1:dbuf%n), zacc)
        call sampler_postprocess(dsmpl, zacc, zout)
    end subroutine

    subroutine fill_sampler_zz()
        if (allocated(zacc)) &
            deallocate(zacc, zout)
        allocate(zacc(zsmpl%nacc, zsmpl%ncomp), zout(zsmpl%npost, zsmpl%ncomp))
        zacc(:, :) = 0
        zout(:, :) = 0
        call sampler_add(zsmpl, zbuf%buffer(1:zbuf%n), zacc)
        call sampler_postprocess(zsmpl, zacc, zout)
    end subroutine


end program
