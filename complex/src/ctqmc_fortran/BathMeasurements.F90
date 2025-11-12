module MBathMeasurements
    use MBath
    use MDiagram
    use MDiagramMeasurements
    implicit none

    ! A type for the accumulating measurement of the 2P-GF depending on 4 fermionic frequencies,
    ! with the internal sampling of the 2-frequency GF
    type, abstract, extends(TAccumulatingMeasurement) :: TAccumulatingSamplingBathMeasurement
        class(ZSampler2D), allocatable :: g2_sampler
        integer(c_int64_t), allocatable :: sparse_freq(:)
    contains
        !procedure, pass :: add_g2sampler => samp_meas_add_g2sampler
        !procedure, pass :: c_add_g2sampler => samp_meas_c_add_g2sampler
        procedure(bath_to_g2_to_array_interface), deferred, nopass :: bath_to_g2_to_array
        procedure, pass :: diagram_to_array => bath_meas_diagram_to_array2
    end type TAccumulatingSamplingBathMeasurement

    ! A type that handles extracting contribution to measured
    ! quantities from the sign and bath part of Monte Carlo
    ! configurations and uses an accumulator for storage and
    ! processing.
    type, abstract, extends(TAccumulatingMeasurement) :: TAccumulatingBathMeasurement
    contains
        procedure(bath_to_array_interface), deferred, nopass :: bath_to_array
        procedure, pass :: diagram_to_array => bath_meas_diagram_to_array
    end type TAccumulatingBathMeasurement

    ! A type that handles extracting contribution to measured
    ! quantities from the sign and bath part of Monte Carlo
    ! configurations and uses samplers for processing.
    type, abstract, extends(TSamplingMeasurement) :: TSamplingBathMeasurement
    contains
        procedure(bath_to_sampbuf_interface), deferred, nopass :: bath_to_sampbuf
        procedure, pass :: diagram_to_sampbuf => bath_meas_diagram_to_sampbuf
    end type TSamplingBathMeasurement

    ! A type that handles extracting contribution to measured
    ! quantities from the sign and bath part of Monte Carlo
    ! configurations and uses a 2D sampler for processing.
    type, abstract, extends(TSamplingMeasurement2D) :: TSamplingBathMeasurement2D
    contains
        procedure(bath_to_sampbuf_interface2d), deferred, nopass :: bath_to_sampbuf2d
        procedure, pass :: diagram_to_sampbuf2d => bath_meas_diagram_to_sampbuf2d
    end type TSamplingBathMeasurement2d

    abstract interface
        ! Interface for the concrete measurement implementation that
        ! generates the correctly shaped contribution to the accumulator
        ! buffer out from the bath part b and MC sign sign of the
        ! current configuration.
        subroutine bath_to_array_interface(b, sign, out)
            import :: TBath, c_double, c_double_complex
            !type(TBath), intent(in) :: b
            type(TBath), intent(in) :: b
            !real(c_double), intent(in) :: sign
            complex(c_double_complex), intent(in) :: sign
            complex(c_double_complex), intent(inout) :: out(:)
        end subroutine bath_to_array_interface

        ! Interface for the concrete measurement implementation of
        ! the four frequency2p-GF at a set of sampling frequencies that
        ! generates the correctly shaped contribution to the accumulator
        ! buffer out from the bath part b and MC sign sign of the
        ! current configuration.
        subroutine bath_to_g2_to_array_interface(g2_sampler, sparse_freq, b, sign, out)
            use, intrinsic :: iso_c_binding
            import :: TBath, ZSampler2D
            class(ZSampler2D), allocatable :: g2_sampler
            integer(c_int64_t), intent(in) :: sparse_freq(:)
            type(TBath), intent(in) :: b
            !real(KINDR), intent(in) :: sign
            complex(c_double_complex), intent(in) :: sign
            complex(c_double_complex), intent(inout) :: out(:)
        end subroutine bath_to_g2_to_array_interface

        ! Interface for the concrete measurement implementation that
        ! generates the samples from the bath part b and MC sign sign of
        ! the current configuration and fills them into the passed
        ! buffer sampbuf.
        subroutine bath_to_sampbuf_interface(b, sign, buf)
            import TBath, c_double, c_double_complex, ZSampleBuffer
            !type(TBath), intent(in) :: b
            type(TBath), intent(in) :: b
            !real(c_double), intent(in) :: sign
            complex(c_double_complex), intent(in) :: sign
            type(ZSampleBuffer), intent(inout) :: buf
        end subroutine bath_to_sampbuf_interface

        ! Interface for the concrete measurement implementation that
        ! generates the samples from the bath part b and MC sign of
        ! the current configuration and fills them into the passed 2D
        ! buffer sampbuf.
        subroutine bath_to_sampbuf_interface2d(b, sign, buf)
            import TBath, c_double, c_double_complex, ZSampleBuffer2D
            !type(TBath), intent(in) :: b
            type(TBath), intent(in) :: b
            !real(KINDR), intent(in) :: sign
            complex(c_double_complex), intent(in) :: sign
            type(ZSampleBuffer2D), intent(inout) :: buf
        end subroutine bath_to_sampbuf_interface2d

    end interface

    ! Type for Green's function measurement from the bath part of the
    ! MC configuration using removal estimator.
    type, extends(TSamplingBathMeasurement) :: TGreenM
    contains
        procedure, nopass :: bath_to_sampbuf => measure_G
    end type TGreenM

    ! Type for two-frequency Green's function measurement
    ! from the bath part of the MC configuration.
    type, extends(TSamplingBathMeasurement2D) :: TG2iw
    contains
        procedure, nopass :: bath_to_sampbuf2d => measure_G2iw
    end type TG2iw

    ! Type for four-frequency 2P - Green's function measurement
    ! from the bath part of the MC configuration.
    type, extends(TAccumulatingSamplingBathMeasurement) :: TG4iw_full
    contains
        procedure, nopass :: bath_to_g2_to_array => measure_G4iw_full
    end type TG4iw_full

    ! Type for expansion order measurement from the bath part of the MC
    ! configuration.
    type, extends(TAccumulatingBathMeasurement) :: TExpansionOrderM
    contains
        procedure, nopass :: bath_to_array => measure_expord
    end type TExpansionOrderM
contains

   !subroutine samp_meas_add_g2sampler(meas, samp)
   !   class(TAccumulatingSamplingBathMeasurement), intent(inout) :: meas
   !   class(ZSampler2D), intent(inout), target :: samp
!
   !   meas%g2_sampler => samp
   !end subroutine samp_meas_add_g2sampler

   ! Add references to a sampler and its backing accumulator as above,
   ! but takes C void* csamp and cacc to reference wrapper types.
   !subroutine samp_meas_c_add_g2sampler(meas, csamp)
   !   class(TAccumulatingSamplingBathMeasurement), intent(inout) :: meas
   !   type(c_ptr), intent(inout) :: csamp
   !   type(ZSampler2D), pointer :: samp
!
   !   call c_f_pointer(csamp, samp)
   !   call meas%add_g2sampler(samp%inner)
   !end subroutine samp_meas_c_add_g2sampler

    ! Implementation of diagram_to_array for subtypes that need only
    ! the bath part and sign of the MC configuration.
    subroutine bath_meas_diagram_to_array(meas, diag, out)
        class(TAccumulatingBathMeasurement), intent(in) :: meas
        type(TDiagram), intent(in) :: diag
        complex(c_double_complex), intent(inout) :: out(:)

        call meas%bath_to_array(diag%bath, get_diagram_sign(diag), out)
    end subroutine bath_meas_diagram_to_array

    ! for 2P-G4iw measurement with internal sampling for 2-frequency GF
    subroutine bath_meas_diagram_to_array2(meas, diag, out)
        class(TAccumulatingSamplingBathMeasurement), intent(in) :: meas
        type(TDiagram), intent(in) :: diag
        complex(c_double_complex), intent(inout) :: out(:)

        call meas%bath_to_g2_to_array( &
                meas%g2_sampler, meas%sparse_freq, diag%bath, &
                get_diagram_sign(diag), out)
    end subroutine bath_meas_diagram_to_array2

    ! Implementation of diagram_to_sampbuf for subtypes that need only
    ! the bath part and sign of the MC configuration.
    subroutine bath_meas_diagram_to_sampbuf(meas, diag, sampbuf)
        class(TSamplingBathMeasurement), intent(in) :: meas
        type(TDiagram), intent(in) :: diag
        type(ZSampleBuffer), intent(inout) :: sampbuf

        call meas%bath_to_sampbuf(diag%bath, get_diagram_sign(diag), sampbuf)
    end subroutine bath_meas_diagram_to_sampbuf

    ! for 2d sampling measurement
    subroutine bath_meas_diagram_to_sampbuf2d(meas, diag, sampbuf)
        class(TSamplingBathMeasurement2D), intent(in) :: meas
        type(TDiagram), intent(in) :: diag
        type(ZSampleBuffer2D), intent(inout) :: sampbuf

        call meas%bath_to_sampbuf2d(diag%bath, get_diagram_sign(diag), sampbuf)
    end subroutine bath_meas_diagram_to_sampbuf2d

    ! Get removal estimator samples of the Green's function from the
    ! bath part and sign of the current MC configuration.
    subroutine measure_G(b, sign, buf)
        !type(TBath), intent(in) :: b
        type(TBath), intent(in) :: b
        !real(c_double), intent(in)           :: sign
        complex(c_double_complex), intent(in)           :: sign
        type(ZSampleBuffer), intent(inout) :: buf
        type(ZSamplePoint)                 :: p
        integer                           :: i, j
        complex(c_double_complex), allocatable       :: minv(:, :)
        type(TBathOperator), allocatable :: bath_crea(:), bath_annh(:)

        allocate(bath_crea, source=get_creators(b))
        allocate(bath_annh, source=get_annihilators(b))

        call prepare_buffer(buf, size(b)**2)
        call copyinv(b, minv)
        annhloop: do j = 1, size(b)
            crealoop: do i = 1, size(b)
                p%i = (bath_crea(i)%orb - 1) * 4 * get_nbands(b)&
                 + (bath_crea(i)%sp - 1) * 2 * get_nbands(b)&
                 + (bath_annh(j)%orb - 1) * 2&
                 + bath_annh(j)%sp
                p%x = bath_annh(j)%tau - bath_crea(i)%tau
                !p%fx = b%MInv_full%Mat(i, j) * sign
                p%fx = minv(i, j) * sign

                buf%buffer(size(b) * (i - 1) + j) = p
            end do crealoop
        end do annhloop
    end subroutine measure_G

    ! Get removal estimator samples of the two-frequency Green's function from the
    ! bath part and sign of the current MC configuration.
    subroutine measure_G2iw(b, sign, buf)
        !type(TBath), intent(in) :: b
        type(TBath), intent(in) :: b
        !real(KINDR), intent(in)           :: sign
        complex(c_double_complex), intent(in)           :: sign
        type(ZSampleBuffer2D), intent(inout) :: buf
        type(ZSamplePoint2D)                 :: p
        integer                           :: i, j
        complex(c_double_complex), allocatable       :: minv(:, :)
        type(TBathOperator), allocatable :: bath_crea(:), bath_annh(:)
        complex(c_double_complex) :: expfact1, expfact2
        real(c_double) :: beta
        real(c_double), parameter :: PI = acos(-1.0D0)

        allocate(bath_crea, source=get_creators(b))
        allocate(bath_annh, source=get_annihilators(b))

        beta = get_beta(b)

        call prepare_buffer(buf, size(b)**2)
        call copyinv(b, minv)
        !call debug_print_bathconfig_highprec(b)
        annhloop: do j = 1, size(b)
            crealoop: do i = 1, size(b)
                p%i = (bath_crea(i)%orb - 1) * 4 * get_nbands(b)&
                    + (bath_crea(i)%sp - 1) * 2 * get_nbands(b)&
                    + (bath_annh(j)%orb - 1) * 2&
                    + bath_annh(j)%sp
                p%x(1) = bath_annh(j)%tau
                p%x(2) = bath_crea(i)%tau

                !p%fx = b%MInv_full%Mat(i, j) * sign

                expfact1 = cmplx( &
                      cos(PI/beta * bath_annh(j)%tau), &
                      sin(PI/beta * bath_annh(j)%tau), &
                      c_double_complex)
                expfact2 = cmplx( &
                      cos(PI/beta * (-bath_crea(i)%tau)), &
                      sin(PI/beta * (-bath_crea(i)%tau)), &
                      c_double_complex)

                p%fx = minv(i,j) * sign !* expfact1 * expfact2


                buf%buffer(size(b) * (i - 1) + j) = p
           end do crealoop
        end do annhloop

    end subroutine measure_G2iw

    ! Measurement of 2P-GF depending on 4 frequencies, with dense NFFT
    ! of the two-frequency G2iw
    subroutine measure_G4iw_full(g2_sampler, sparse_freq, b, sign, out)
        use, intrinsic :: ieee_arithmetic
        use, intrinsic :: iso_c_binding
        class(ZSampler2D), allocatable :: g2_sampler
        integer(c_int64_t), intent(in) :: sparse_freq(:)
        type(TBath), intent(in) :: b
        complex(c_double_complex), intent(in)           :: sign
        complex(c_double_complex), intent(inout)        :: out(:)
        type(ZSampleBuffer2D)                :: buf
        type(ZSamplePoint2D)                 :: p
        integer                              :: i, j, li, o1, o2, o3, o4, s1, s2, s3, s4, w1, w2, w3, w4
        complex(c_double_complex), allocatable          :: minv(:, :)
        type(TBathOperator), allocatable     :: bath_crea(:), bath_annh(:)

        complex(c_double_complex), allocatable          :: sacc(:,:), sout(:,:)
        complex(c_double_complex), allocatable                        :: g2iw(:,:,:,:,:,:)
        complex(c_double_complex), allocatable                       :: g4iw(:,:,:,:,:,:,:,:,:)
        !class(ZSampler2D), allocatable :: g2samp
        integer(c_int64_t), allocatable :: sparse_freq_f(:,:)
        integer :: n_max, n, L, nbands, ncomp

        ! get largest fermionic frequency idx from list of sparse frequency indices
        n_max = maxval(sparse_freq)

        ! TODO: maximum grid size should not exceed niw=2000, change later for better asymptotic treatment
        if (n_max>=1999) then
            n_max=1999
        end if
        n = (n_max - 1)/2 + 1

        ! reshape sparse frequency list to a set of 4 frequency tuples
        L = size(sparse_freq)/4
        allocate(sparse_freq_f(4,L))
        sparse_freq_f = transpose(reshape(sparse_freq, [L, 4]))

        ! get ncomp for shape of g4iw-accu
        nbands = get_nbands(b)
        ncomp = 4*nbands*nbands

        ! allocate accumulator and buffer for g2-sampler
        allocate(sacc(g2_sampler%nacc, g2_sampler%ncomp))
        allocate(sout(g2_sampler%npost, g2_sampler%ncomp))

        ! fill buffer with sample points
        allocate(bath_crea, source=get_creators(b))
        allocate(bath_annh, source=get_annihilators(b))

        call prepare_buffer(buf, size(b)**2)
        call copyinv(b, minv)
        annhloop: do j = 1, size(b)
            crealoop: do i = 1, size(b)
                p%i = (bath_crea(i)%orb - 1) * 4 * get_nbands(b)&
                    + (bath_crea(i)%sp - 1) * 2 * get_nbands(b)&
                    + (bath_annh(j)%orb - 1) * 2&
                    + bath_annh(j)%sp
                p%x(1) = bath_annh(j)%tau
                p%x(2) = bath_crea(i)%tau
                p%fx = minv(i,j) * sign

                buf%buffer(size(b) * (i - 1) + j) = p
            end do crealoop
        end do annhloop

        ! sample and postprocess for G2(v,v')
        sacc(:, :) = 0
        sout(:, :) = 0

        call sampler_add(g2_sampler, buf%buffer(1:buf%n), sacc)
        call sampler_postprocess(g2_sampler, sacc, sout)

        !FIX ME: nfft writes nan into off diag GF(v,v'), replace them with zeros
        where (sout/=sout)
            sout = 0.0
        end where

        ! build G4(v,v',v'',v''') with G2 at sparse frequency set
        g2iw = reshape(sout, (/nbands, 2, nbands, 2, 2*n, 2*n/), order=(/ 6,5,4,3,2,1 /))
        allocate(g4iw(nbands, 2, nbands, 2, nbands, 2, nbands, 2, L))

        o1loop: do o1 = 1, nbands
            s1loop: do s1 = 1, 2
                o2loop: do o2 = 1, nbands
                    s2loop: do s2 = 1, 2
                        o3loop: do o3 = 1, nbands
                            s3loop: do s3 = 1, 2
                                o4loop: do o4 = 1, nbands
                                    s4loop: do s4 = 1, 2
                                        do li=1, L
                                            w1 = sparse_freq_f(1,li)
                                            w1 = (w1 - 1)/2 + 1 + n
                                            w2 = sparse_freq_f(2,li)
                                            w2 = (w2 - 1)/2 + 1 + n
                                            w3 = sparse_freq_f(3,li)
                                            w3 = (w3 - 1)/2 + 1 + n
                                            w4 = sparse_freq_f(4,li)
                                            w4 = (w4 - 1)/2 + 1 + n

                                            g4iw(o1,s1,o2,s2,o3,s3,o4,s4,li) =&
                                            g2iw(o2,s2,o1,s1,w2,w1)*g2iw(o4,s4,o3,s3,w4,w3)&
                                            - g2iw(o4,s4,o1,s1,w4,w1)*g2iw(o2,s2,o3,s3,w2,w3)
                                        end do
                                    end do s4loop
                                end do o4loop
                            end do s3loop
                        end do o3loop
                    end do s2loop
                end do o2loop
            end do s1loop
        end do o1loop

        ! flatten G4 and pass to accumulator buffer
        out = out + reshape(g4iw, (/ ncomp**2 * L /)) * sign

    end subroutine measure_G4iw_full


    ! Measurement of 2P-GF depending on 4 frequencies, with naive FT
    ! of the two-frequency G2iw on each sampling point
    subroutine measure_G4iw_full_naive(sparse_freq, b, sign, out)
        use, intrinsic :: ieee_arithmetic
        use, intrinsic :: iso_c_binding
        !class(ZSampler2D), intent(inout) :: g2_sampler
        integer(c_int64_t), intent(in) :: sparse_freq(:)
        !type(TBath), intent(in) :: b
        type(TBath), intent(in) :: b
        !real(KINDR), intent(in)           :: sign
        complex(c_double_complex), intent(in)           :: sign
        complex(c_double_complex), intent(inout)        :: out(:)
        type(ZSampleBuffer2D)                :: buf
        type(ZSamplePoint2D)                 :: p
        integer                              :: i, j, li, o1, o2, o3, o4, s1, s2, s3, s4
        complex(c_double_complex), allocatable          :: minv(:, :)
        type(TBathOperator), allocatable     :: bath_crea(:), bath_annh(:)

        complex(c_double_complex), allocatable          :: g4iw(:,:,:,:,:,:,:,:,:)
        complex(c_double_complex), allocatable          :: g21(:), g43(:), g41(:), g23(:)
        complex(c_double_complex), allocatable          :: g21_os(:,:,:,:), g43_os(:,:,:,:), g41_os(:,:,:,:), g23_os(:,:,:,:)
        class(ZSampler2D), allocatable :: g2samp
        integer(c_int64_t), allocatable :: sparse_freq_f(:,:), w1, w2, w3, w4
        real(c_double) ::  beta, coeff(2,2)
        integer :: niw(2), statsign(2), n_max, n, L, nbands, ncomp

        ! get largest fermionic frequency idx from list of sparse frequency indices
        n_max = maxval(sparse_freq)
        n = (n_max - 1)/2 + 1

        ! reshape sparse frequency list to a set of 4 frequency tuples
        L = size(sparse_freq)/4
        allocate(sparse_freq_f(4,L))
        sparse_freq_f = transpose(reshape(sparse_freq, [L, 4]))

        ! init sampler for G2(v,v')
        coeff(:,:) = 0
        coeff(1,1) = 1
        coeff(2,2) = 1
        statsign(1) = -1
        statsign(2) = -1
        niw(1) = n
        niw(2) = n
        beta = get_beta(b)
        nbands = get_nbands(b)
        ncomp = 4*nbands*nbands

        ! fill buffer with sample points
        allocate(bath_crea, source=get_creators(b))
        allocate(bath_annh, source=get_annihilators(b))

        call prepare_buffer(buf, size(b)**2)
        call copyinv(b, minv)
        annhloop: do j = 1, size(b)
            crealoop: do i = 1, size(b)
                p%i = (bath_crea(i)%orb - 1) * 4 * get_nbands(b)&
                    + (bath_crea(i)%sp - 1) * 2 * get_nbands(b)&
                    + (bath_annh(j)%orb - 1) * 2&
                    + bath_annh(j)%sp
                p%x(1) = bath_annh(j)%tau
                p%x(2) = -bath_crea(i)%tau
                !p%fx = b%MInv_full%Mat(i, j) * sign
                p%fx = minv(i,j) * sign

                buf%buffer(size(b) * (i - 1) + j) = p
            end do crealoop
        end do annhloop

        allocate(g21(ncomp))
        allocate(g43(ncomp))
        allocate(g41(ncomp))
        allocate(g23(ncomp))
        allocate(g4iw(nbands, 2, nbands, 2, nbands, 2, nbands, 2, L))

        do li=1, L
            w1 = sparse_freq_f(1,li)
            w2 = sparse_freq_f(2,li)
            w3 = sparse_freq_f(3,li)
            w4 = sparse_freq_f(4,li)

            g21 = 0
            g43 = 0
            g41 = 0
            g23 = 0

            do i=1, buf%n
                g21(buf%buffer(i)%i) = g21(buf%buffer(i)%i) +&
                                           nfourier2d(buf%buffer(i)%x(1), buf%buffer(i)%x(2), buf%buffer(i)%fx, w2, w1, beta)
                g43(buf%buffer(i)%i) = g43(buf%buffer(i)%i) +&
                                           nfourier2d(buf%buffer(i)%x(1), buf%buffer(i)%x(2), buf%buffer(i)%fx, w4, w3, beta)
                g41(buf%buffer(i)%i) = g41(buf%buffer(i)%i) +&
                                           nfourier2d(buf%buffer(i)%x(1), buf%buffer(i)%x(2), buf%buffer(i)%fx, w4, w1, beta)
                g23(buf%buffer(i)%i) = g23(buf%buffer(i)%i) +&
                                          nfourier2d(buf%buffer(i)%x(1), buf%buffer(i)%x(2), buf%buffer(i)%fx, w2, w3, beta)
            end do

            g21_os = reshape(g21, (/nbands, 2, nbands, 2/), order=(/ 4,3,2,1 /))
            g43_os = reshape(g43, (/nbands, 2, nbands, 2/), order=(/ 4,3,2,1 /))
            g41_os = reshape(g41, (/nbands, 2, nbands, 2/), order=(/ 4,3,2,1 /))
            g23_os = reshape(g23, (/nbands, 2, nbands, 2/), order=(/ 4,3,2,1 /))

            o1loop: do o1 = 1, nbands
                s1loop: do s1 = 1, 2
                    o2loop: do o2 = 1, nbands
                        s2loop: do s2 = 1, 2
                            o3loop: do o3 = 1, nbands
                                s3loop: do s3 = 1, 2
                                    o4loop: do o4 = 1, nbands
                                        s4loop: do s4 = 1, 2
                                            g4iw(o1,s1,o2,s2,o3,s3,o4,s4,li) =&
                                            g21_os(o2,s2,o1,s1)*g43_os(o4,s4,o3,s3)&
                                            - g41_os(o4,s4,o1,s1)*g23_os(o2,s2,o3,s3)
                                        end do s4loop
                                    end do o4loop
                                end do s3loop
                            end do o3loop
                        end do s2loop
                    end do o2loop
                end do s1loop
            end do o1loop

        end do

        ! flatten G4 and pass to accumulator buffer
        out = out + reshape(g4iw, (/ ncomp**2 * (2*n)**4 /)) * sign

    end subroutine measure_G4iw_full_naive


    ! naive FT implementation for G2iw measurement on sparse frequencies
    function nfourier2d(tau1, tau2, fx, w1, w2, beta) result(g2)
        real(c_double) :: tau1, tau2, beta
        complex(c_double_complex) :: fx, g2
        integer(c_int64_t) :: w1, w2
        complex(c_double_complex) :: expfact1, expfact2
        real(c_double), parameter :: PI = acos(-1.0D0)

        ! Do naive Fourier transform. This scales as O(N**2)
        expfact1 = cmplx( &
                    cos(PI/beta * tau1 * w1), &
                    sin(PI/beta * tau1 * w1), &
                    c_double_complex)
        expfact2 = cmplx( &
                    cos(PI/beta * tau2 * w2), &
                    sin(PI/beta * tau2 * w2), &
                    c_double_complex)
        g2 = expfact1 * expfact2 * fx/beta
    end function


    ! Accumulate expansion order contributions per flavor.
    ! FIXME: offdiag mis-/underaccounting with this (=old) method
    subroutine measure_expord(b, sign, out)
        !type(TBath), intent(in) :: b
        type(TBath), intent(in) :: b
        !real(c_double), intent(in)           :: sign
        complex(c_double_complex), intent(in)           :: sign
        !real(c_double), intent(inout)        :: out(:)
        complex(c_double_complex), intent(inout)      :: out(:)
        integer                           :: iB, iS, orderend, order

        orderend = size(out) / get_nbands(b) / 2
        do iB = 1, get_nbands(b)
            do iS = 1, 2
                order = min(orderend - 1, get_nosoper(b, iB, iS)/2)
                out((iB - 1) * 2 * orderend + (iS - 1) * orderend + order + 1) =&
                 out((iB - 1) * 2 * orderend + (iS - 1) * orderend + order + 1)&
                 + 1
            end do
        end do
    end subroutine measure_expord
end module MBathMeasurements
