program test_exp_discrete
    use iso_c_binding, only: c_double, c_long
    use MLogDiscretization
    use MRandomNumbers
    use Testing
    implicit none

    call check_default
    call check_log2
    call check_segments
    call benchmark_log2

contains

    subroutine check_default
        type(TLogDiscretization) :: discr
        type(TLogGrid) :: ediscr

        call init_log_discretization(discr, 1.1920928955078125d-7, 3)
        call assert_close(discr%eps,  LOG_DISCR_FLOAT_3%eps)
        call assert_close(discr%alpha,  LOG_DISCR_FLOAT_3%alpha)
        call assert_close(discr%lambda,  LOG_DISCR_FLOAT_3%lambda)
        call assert_close(discr%root_of_2,  LOG_DISCR_FLOAT_3%root_of_2)

        call assert_close(discr%alpha ** discr%root_of_2, 2.0d0)

        call init_log_grid(ediscr, discr, 15d0, 320.0d0)
        call assert_equal(ediscr%t_max >= 15d0, .true.)
        call assert_equal(ediscr%w_max >= 320d0, .true.)
        call assert_equal(ediscr%t_max <= 30d0, .true.)
        call assert_equal(ediscr%w_max <= 640d0, .true.)
        write (0,*) lbound(ediscr%tn), ubound(ediscr%tn)
    end

    subroutine check_log2
        class(RandomEngine), allocatable :: rng
        real(c_double), parameter :: t_min = 1d-30, t_max = 180d0
        real(c_double) :: t
        integer :: i

        call init_xorshift_star(rng)
        do i = 1, 1000000
            t = random_real(rng, t_min, t_max)
            call assert_close(log2_sloppy(t), log(t) / log(2.0d0), atol=0.005d0)
        enddo
    end

    subroutine check_segments
        type(TLogGrid) :: ediscr
        class(RandomEngine), allocatable :: rng
        real(c_double), parameter :: t_max = 15d0
        real(c_double) :: t, dt, spread
        integer :: n, i

        call init_mersenne_twister(rng)
        call init_log_grid(ediscr, LOG_DISCR_FLOAT_3, t_max, 320.0d0)

        ! XXX Store this in discretization
        spread = 1.01 * sqrt(ediscr%spec%alpha)
        do i = 1, 100000
            t = random_real(rng, t_max)
            call segment_time(ediscr, t, n, dt)
            call assert_close(ediscr%tn(n) + dt, t)
            call assert_equal(n >= ediscr%n_min, .true.)
            call assert_equal(n <= ediscr%n_max, .true.)
            call assert_equal(t / ediscr%tn(n) <= spread, .true.)
            if (n /= ediscr%n_min) then
                call assert_equal(ediscr%tn(n) / t <= spread, .true.)
            endif
        enddo
    end

    subroutine benchmark_log2
        use MLibC
        real(c_double), parameter :: t_min = 1d-40, t_max = 1d70
        integer, parameter :: n = 50000000
        class(RandomEngine), allocatable :: rng
        integer :: i
        real(c_double) :: s, t0, t1, tref
        real(c_double), allocatable :: r(:)
        real(c_double), parameter :: INV_LOG2 = 1.0 / log(2.0d0)

        call init_xorshift_star(rng)
        allocate(r(n))
        do i = 1, n
            r(i) = random_real(rng, t_min, t_max)
        enddo
        call cpu_time(t0)
        s = 0
        do i = 1, n
            s = s + r(i)
        enddo
        write (0,*) s
        call cpu_time(t1)
        tref = t1 - t0

        call cpu_time(t0)
        s = 0
        do i = 1, n
            s = s + libc_log2(r(i))
        enddo
        write (0,*) s
        call cpu_time(t1)
        write (0, "('log2:   ', F10.3, ' ns/invoc.')") 1d9/n * (t1 - t0 - tref)

        call cpu_time(t0)
        s = 0
        do i = 1, n
            s = s + log2_sloppy(r(i))
        enddo
        write (0,*) s
        call cpu_time(t1)
        write (0, "('log2 sl:', F10.3, ' ns/invoc.')") 1d9/n * (t1 - t0 - tref)
    end

end program
