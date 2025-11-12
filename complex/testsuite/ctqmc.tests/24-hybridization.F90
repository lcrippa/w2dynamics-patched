#define INFO line=__LINE__

program test_hybridization
    use iso_c_binding, only: c_double, c_double_complex
    use MHybridizationD
    use MRandomNumbers
    use Testing
    implicit none

    call test_prop()
    call test_hybr(0.1d0)
    call test_log()
    call test_discrete()
    call bench_discrete()

contains
    subroutine test_prop()
        real(c_double) :: gtau

        gtau = propagator_tau(1d0, 0d0, 0d0, 1000d0, .true.)
        call assert_close(gtau, -0.5d0, INFO)

        gtau = propagator_tau(1d0, 1000d0, 0d0, 1000d0, .true.)
        call assert_close(gtau, -0.5d0, INFO)

        gtau = propagator_tau(1d0, -1000d0, 0d0, 1000d0, .true.)
        call assert_close(gtau, 0.5d0, INFO)

        gtau = propagator_tau(0.7d0, 0d0, 1000d0, 2d0, .true.)
        call assert_close(gtau, -0.7d0, INFO)

        gtau = propagator_tau(80.0d0, 0d0, -1000d0, 2d0, .true.)
        call assert_close(gtau, 0.0d0, INFO)

        gtau = propagator_tau(1.0d0, 992d0, -2.25d0, 1000d0, .true.)
        call assert_close(gtau, -1.522997974471263d-8, INFO)

        gtau = propagator_tau(1.0d0, -4d0, 1.5d0, 3d0, .true.)
        if (gtau == gtau) call terminate(msg="Expect NaN", INFO)

        gtau = propagator_tau(1.0d0, 4d0, 2d0, -10d0, .true.)
        if (gtau == gtau) call terminate(msg="Expect NaN", INFO)
    end

    subroutine test_hybr(eps)
        real(c_double), intent(in) :: eps
        integer :: i, isp, jsp
        real(c_double) :: beta=7.0d0, tau, val, val_cmp
        real(c_double) :: ftau(0:1000), ftau_full(1, 2, 1, 2, 0:1000)
        class(DHybridization), allocatable :: hyb

        do i = 0,1000
            tau = i * beta / 1000
            ftau(i) = propagator_tau(1.0d0, tau, eps, beta, .true.)
        end do

        forall (isp=1:2, jsp=1:2) &
            ftau_full(1, isp, 1, jsp, :) = ftau(:)
        call init_linear_hybridization(hyb, ftau_full, beta, .true.)

        ! 1000 and 997 are coprime
        do  i = -997,997
            tau = i * beta / 997
            val = get_hybr_value(hyb, 1, 1, 1, 1, tau)
            val_cmp = propagator_tau(1.0d0, tau, eps, beta, .true.)
            call assert_close(val, val_cmp, atol=1d-7, INFO)
        end do
    end

    subroutine test_log
        class(DHybridization), allocatable :: hybr
        real(c_double) :: vmat(1, 2, 1, 2, 1), eps(1)
        real(c_double) :: tau, ftau, ftau_cmp, beta
        integer :: i, n

        eps(1) = 3.0
        vmat = 0.8
        beta = 14.0

        call init_log_hybridization(hybr, vmat, eps, beta, .true.)
        n = 10000
        do i = -n, n
            tau = 1.0d0 * i / n * beta
            ftau = get_hybr_value(hybr, 1, 1, 1, 1, tau)
            ftau_cmp = propagator_tau(vmat(1,1,1,1,1), tau, eps(1), beta, .true.)
            call assert_close(ftau, ftau_cmp, atol=1.5d-7)
            !write (*, '(ES26.16, ES26.16, ES26.16)') tau, ftau, ftau_cmp
        enddo
    end

    subroutine test_discrete
        class(DHybridization), allocatable :: log_hybr
        class(DHybridization), allocatable :: lin_hybr
        real(c_double), parameter :: eps(4) = (/ -13.2, -3.1, 0.1, 7.9 /)
        real(c_double), parameter :: v(4) = (/ 0.8, 1.2, 0.05, 3.7 /)
        real(c_double), parameter :: beta = 6.5
        logical, parameter :: fermion = .true.

        real(c_double) :: vmat(1, 2, 1, 2, 4), fmat(1, 2, 1, 2, 0:1000)
        real(c_double) :: tau, ftau_lin, ftau_log
        integer :: i, n

        vmat(:, :, :, :, :) = 0
        vmat(1, 1, 1, 1, :) = v
        vmat(1, 2, 1, 2, :) = 0.4 * v
        call init_log_hybridization(log_hybr, vmat, eps, beta, fermion)

        call discretize_propagator(vmat, eps, beta, fermion, fmat)
        call init_linear_hybridization(lin_hybr, fmat, beta, fermion)

        n = 3051
        do i = -n, n
            tau = 1.0d0 * i / n * beta
            ftau_lin = get_hybr_value(lin_hybr, 1, 1, 1, 1, tau)
            ftau_log = get_hybr_value(log_hybr, 1, 1, 1, 1, tau)
            call assert_close(ftau_lin, ftau_log, rtol=5d-4)

            ftau_lin = get_hybr_value(lin_hybr, 1, 2, 1, 2, tau)
            ftau_log = get_hybr_value(log_hybr, 1, 2, 1, 2, tau)
            call assert_close(ftau_lin, ftau_log, rtol=5d-4)
            !write (*, '(3ES26.6)') tau, ftau_lin, ftau_log
        enddo
    end

    subroutine bench_discrete
        real(c_double), parameter :: eps(4) = (/ -13.2, -3.1, 0.1, 7.9 /)
        real(c_double), parameter :: v(4) = (/ 0.8, 1.2, 0.05, 3.7 /)
        real(c_double), parameter :: beta = 6.5
        logical, parameter :: fermion = .true.
        integer, parameter :: n = 100000

        class(DHybridization), allocatable :: log_hybr
        class(DHybridization), allocatable :: lin_hybr
        class(RandomEngine), allocatable :: rng
        real(c_double) :: vmat(1, 2, 1, 2, 4), fmat(1, 2, 1, 2, 0:1000)
        real(c_double) :: t0, t1
        type(TBaseOper) :: crea(20), annh(20)
        real(c_double) :: hybmat(20, 20), hybsum(20, 20)
        integer :: i

        ! Initialize
        vmat(:, :, :, :, :) = 0
        vmat(1, 1, 1, 1, :) = v
        vmat(1, 2, 1, 2, :) = 0.4 * v
        call init_log_hybridization(log_hybr, vmat, eps, beta, fermion)
        call discretize_propagator(vmat, eps, beta, fermion, fmat)
        call init_linear_hybridization(lin_hybr, fmat, beta, fermion)

        ! RNG
        call init_xorshift_star(rng)
        do i = 1, size(crea)
            crea(i)%orb = 1
            crea(i)%sp = random_integer(rng, 1, 2)
            crea(i)%tau = random_real(rng, beta/2)
            annh(i)%orb = 1
            annh(i)%sp = random_integer(rng, 1, 2)
            annh(i)%tau = random_real(rng, beta)
        enddo

        hybsum(:, :) = 0

        call cpu_time(t0)
        do i = 1, n
            crea(:)%tau = crea(:)%tau + beta / (2*n)
            call copy_hybr_matrix(lin_hybr, crea, annh, hybmat)
            hybsum = hybsum + hybmat
        enddo
        write (0,*) sum(hybsum)
        call cpu_time(t1)
        write (0, "('lin:   ', F10.3, ' ns/invoc.')") &
                1d9/(1.0d0*n*size(hybsum)) * (t1 - t0)

        hybsum(:, :) = 0
        crea(:)%tau = crea(:)%tau - beta / 2

        call cpu_time(t0)
        do i = 1, n
            crea(:)%tau = crea(:)%tau + beta / (2*n)
            call copy_hybr_matrix(log_hybr, crea, annh, hybmat)
            hybsum = hybsum + hybmat
        enddo
        write (0,*) sum(hybsum)
        call cpu_time(t1)
        write (0, "('log:   ', F10.3, ' ns/invoc.')") &
                1d9/(1.0d0*n*size(hybsum)) * (t1 - t0)
    end subroutine

end program
