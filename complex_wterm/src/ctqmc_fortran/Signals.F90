module MSignal
    use iso_c_binding, only: c_int, c_int64_t, c_funptr, c_funloc
    use MLibC
    implicit none
    private

    integer(c_int64_t), volatile :: signals_fired = 0

    integer, parameter, public :: SIGNAL_INTERRUPT = LIBC_SIGINT
    integer, parameter, public :: SIGNAL_USER_1 = LIBC_SIGUSR1
    integer, parameter, public :: SIGNAL_USER_2 = LIBC_SIGUSR2
    integer, parameter, public :: SIGNAL_TERMINATE = LIBC_SIGTERM

    public :: register_signal_handler
    public :: clear_signal_handler
    public :: record_signal
    public :: any_signal_fired
    public :: clear_signals_fired

contains

    !> Registers a signal handler for a certain signal.
    !!
    !! Note that the usual restrictions for signal handlers apply (you should
    !! not do any I/O or other asynchronous calls).
    subroutine register_signal_handler(sig, handler)
        integer(c_int), intent(in) :: sig
        interface
            subroutine handler(sig) bind(C)
                import :: c_int
                integer(c_int), value :: sig
            end
            function libc_signal(sig, func) bind(C, name='signal')
                import :: c_funptr, c_int
                integer(c_int), intent(in), value :: sig
                type(c_funptr), intent(in), value :: func
                type(c_funptr) :: libc_signal
            end
        end interface
        type(c_funptr) :: prev_handler

        prev_handler = libc_signal(sig, c_funloc(handler))
        ! XXX check for SIG_ERR
    end

    subroutine clear_signal_handler(sig)
        integer(c_int), intent(in) :: sig
    end

    !> Predefined signal handler, which records caught signals
    !!
    !! This handler raises a flag corresponding to the signal being fired.
    !! If a signal is fired twice, the program is aborted.  One can check
    !! flags using `any_signal_fired()`.
    subroutine record_signal(sig) bind(C)
        integer(c_int), value :: sig
        integer(c_int64_t) :: bit

        write (0, 99) sig

        bit = ishft(1_c_int64_t, sig)
        if (iand(signals_fired, bit) /= 0) then
            call libc_abort
        endif
        signals_fired = ior(signals_fired, bit)

    99  format ('Signal No. ', I0, ' has fired')
    end

    !> Return true if record_signal has detected any signals fired
    logical function any_signal_fired()
        any_signal_fired = signals_fired /= 0
    end

    !> Clear all signals detected by record_signal
    subroutine clear_signals_fired()
        signals_fired = 0
    end
end module
