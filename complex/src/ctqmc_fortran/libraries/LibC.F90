module MLibC
    use iso_c_binding
    private

    interface
        !> libc abort() function
        subroutine libc_abort() bind(C, name='abort')
        end subroutine

        !> libc signal() function - install signal handler
        function libc_signal(sig, func) bind(C, name='signal')
            import :: c_funptr, c_int
            integer(c_int), intent(in), value :: sig
            type(c_funptr), intent(in), value :: func
            type(c_funptr) :: libc_signal
        end
    end interface
    public :: libc_abort, libc_signal

    interface libc_expm1
        !> libc expm1f() function
        pure function libc_expm1f(x) bind(C, name='expm1f')
            import :: c_float
            real(c_float), value :: x
            real(c_float) :: libc_expm1f
        end

        !> libc expm1() function
        pure function libc_expm1d(x) bind(C, name='expm1')
            import :: c_double
            real(c_double), value :: x
            real(c_double) :: libc_expm1d
        end
    end interface
    public :: libc_expm1, libc_expm1f, libc_expm1d

    interface libc_log2
        !> libc expm1f() function
        pure function libc_log2f(x) bind(C, name='log2f')
            import :: c_float
            real(c_float), value :: x
            real(c_float) :: libc_log2f
        end

        !> libc expm1f() function
        pure function libc_log2d(x) bind(C, name='log2')
            import :: c_double
            real(c_double), value :: x
            real(c_double) :: libc_log2d
        end
    end interface
    public :: libc_log2, libc_log2f, libc_log2d

    integer, parameter, public :: LIBC_SIGINT  =  2
    integer, parameter, public :: LIBC_SIGUSR1 = 10
    integer, parameter, public :: LIBC_SIGUSR2 = 12
    integer, parameter, public :: LIBC_SIGTERM = 15

    integer(c_intptr_t), parameter, public :: LIBC_SIG_ERR = -1
    integer(c_intptr_t), parameter, public :: LIBC_SIG_DFL = 0
    integer(c_intptr_t), parameter, public :: LIBC_SIG_IGN = 1

end module
