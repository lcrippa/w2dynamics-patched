!> Module providing functions for comparing values or arrays, useful for tests.
module testing
    use MPrinting
    use MSignal
    use iso_c_binding, only: c_double, c_double_complex
    implicit none
    private

    !> Numpy-like function interface:
    !!
    !!      isclose(a, b, atol=0, rtol=2*epsilon)
    !!
    !! which returns true iff `a` and `b` are equal within absolute tolerance
    !! `atol` (defaulting to 0) or relative tolerance `rtol` (default
    !! is twice the machine epsilon).
    interface isclose
        module procedure disclose_0
        module procedure zisclose_0
    end interface isclose
    public isclose

    !> Numpy-like function interface:
    !!
    !!      assert_close(a, b, atol=0, rtol=2*epsilon)
    !!
    !! which aborts the program with exit code 1 and an error message if `a`
    !! and `b` are not equal within absolute tolerance `atol` (defaulting to
    !! zero) and relative tolerance `rtol` (default  is twice the machine
    !! epsilon).
    interface assert_close
        module procedure dassert_close_0, dassert_close_1, dassert_close_2
        module procedure zassert_close_0, zassert_close_1, zassert_close_2
    end interface assert_close
    public assert_close

    !> Numpy-like function interface:
    !!
    !!      assert_equal(a, b)
    !!
    !! which aborts the program with exit code 1 and an error message if `a`
    !! and `b` are not equal.
    interface assert_equal
        module procedure iassert_equal_0, iassert_equal_1
        module procedure lassert_equal_0
    end interface
    public assert_equal

    public terminate
    public get_threshold

    ! XXX: These should probably reside somewhere else
    public read_ftau
    public dump_ftau
contains
    !> terminate the program with exit code 1 and printing a certain message
    subroutine terminate(msg, file, line)
        character(len=*), intent(in), optional :: msg, file
        integer, intent(in), optional :: line

        call write_error_message(msg, file, line)
        error stop 'Test failed'
    end subroutine

    !> Write error message
    subroutine write_error_message(msg, file, line)
        character(len=*), intent(in), optional :: msg, file
        integer, intent(in), optional :: line

        if (present(file) .or. present(line)) then
            if (present(file)) then
                write (0, "('ERROR in ',A)", advance='no') file
                if (present(line)) &
                    write (0, "(':',I0)", advance='no') line
            else
                write (0,"('ERROR on line ',I0)",advance='no') line
            end if
        else
            write (0, "('ERROR')", advance='no')
        endif
        if (present(msg)) then
            write (0, "(': ',/,A)") msg
        else
            write (0, *)
        endif
    end

    subroutine write_deviation(diff, magn, atol, rtol)
        real(c_double), intent(in) :: diff, magn
        real(c_double), intent(in) :: rtol, atol

        write (0, 10) diff, default(atol, 0.0d0)
        write (0, 11) diff/magn, default(rtol, 2 * epsilon(diff))

    10  format('absolute dev. =', ES11.3, ', tolerance =', ES11.3)
    11  format('relative dev. =', ES11.3, ', tolerance =', ES11.3)
    end

    subroutine check_shape(ashape, bshape, file, line)
        integer, intent(in) :: ashape(:), bshape(:)
        character(len=*), intent(in), optional :: file
        integer, intent(in), optional :: line

        if (size(ashape) == size(bshape)) then
            if (all(ashape == bshape)) then
                return
            endif
        endif

        call print_array(ashape, name='shape(a)')
        call print_array(bshape, name='shape(b)')
        call terminate('Arrays differ in shape', file, line)
    end

    subroutine iassert_equal_0(a, b, file, line)
        integer, intent(in) :: a, b
        character(len=*), intent(in), optional :: file
        integer, intent(in), optional :: line

        if (a == b) return

        call print_array(a, name='a')
        call print_array(b, name='b')
        call terminate('Difference between vectors', file, line)
    end

    subroutine iassert_equal_1(a, b, file, line)
        integer, intent(in) :: a(:), b(:)
        character(len=*), intent(in), optional :: file
        integer, intent(in), optional :: line

        call check_shape(shape(a), shape(b))
        if (all(a == b)) return

        call print_array(a, name='a')
        call print_array(b, name='b')
        call terminate('Difference between vectors', file, line)
    end

    subroutine lassert_equal_0(a, b, file, line)
        logical, intent(in) :: a, b
        character(len=*), intent(in), optional :: file
        integer, intent(in), optional :: line

        if (a .eqv. b) return

        write (0,*) 'a = ', a
        write (0,*) 'b = ', b
        call terminate('Difference between scalars', file, line)
    end

    pure function get_threshold(magn, atol, rtol) result(thr)
        real(c_double), intent(in) :: magn
        real(c_double), intent(in), optional :: atol, rtol
        real(c_double) :: thr, thr_r

        real(c_double), parameter :: default_rtol = 2 * epsilon(0.0d0)
        real(c_double), parameter :: default_atol = 0.0d0

        if (present(atol)) then
            thr = atol
        else
            thr = default_atol
        endif
        if (present(rtol)) then
            thr_r = rtol * magn
        else
            thr_r = default_rtol * magn
        endif
        if (thr_r > thr) then
            thr = thr_r
        endif
    end

    !> isclose() for scalars a, b
    pure logical function disclose_0(a, b, atol, rtol) result(result)
        real(c_double), intent(in) :: a, b
        real(c_double), intent(in), optional :: atol, rtol
        real(c_double) :: diff, magn, thr

        magn = max(abs(a), abs(b))
        diff = abs(a - b)
        thr = get_threshold(magn, atol, rtol)
        result = diff <= thr
    end function

    !> assert_close() for scalars a, b
    subroutine dassert_close_0(a, b, atol, rtol, file, line)
        real(c_double), intent(in) :: a, b
        real(c_double), intent(in), optional :: atol, rtol
        character(len=*), intent(in), optional :: file
        integer, intent(in), optional :: line

        real(c_double) :: diff, magn

        magn = max(abs(a), abs(b))
        diff = abs(a - b)
        if (diff <= get_threshold(magn, atol, rtol)) &
            return

        call print_array(a, name='a')
        call print_array(b, name='b')
        call write_deviation(diff, magn, atol, rtol)
        call terminate('significant difference between scalars', file, line)
    end subroutine

    !> assert_close() for vectors a, b
    subroutine dassert_close_1(a, b, atol, rtol, file, line)
        real(c_double), intent(in) :: a(:), b(:)
        real(c_double), intent(in), optional :: atol, rtol
        character(len=*), intent(in), optional :: file
        integer, intent(in), optional :: line

        real(c_double) :: diff, magn

        call check_shape(shape(a), shape(b))
        magn = max(maxval(abs(a)), maxval(abs(b)))
        diff = maxval(abs(a - b))
        if (diff <= get_threshold(magn, atol, rtol)) &
            return

        call print_array(a, name='a')
        call print_array(b, name='b')
        call write_deviation(diff, magn, atol, rtol)
        call terminate('significant difference between vectors', file, line)
    end subroutine

    !> assert_close() for matrices a, b
    subroutine dassert_close_2(a, b, atol, rtol, file, line)
        real(c_double), intent(in) :: a(:,:), b(:,:)
        real(c_double), intent(in), optional :: atol, rtol
        character(len=*), intent(in), optional :: file
        integer, intent(in), optional :: line

        real(c_double) :: diff, magn

        call check_shape(shape(a), shape(b))
        magn = max(maxval(abs(a)), maxval(abs(b)))
        diff = maxval(abs(a - b))
        if (diff <= get_threshold(magn, atol, rtol)) &
            return

        call print_array(a, name='a')
        call print_array(b, name='b')
        call write_deviation(diff, magn, atol, rtol)
        call terminate('significant difference between vectors', file, line)
    end subroutine

    !> isclose() for scalars a, b
    pure logical function zisclose_0(a, b, atol, rtol) result(result)
        complex(c_double_complex), intent(in) :: a, b
        real(c_double), intent(in), optional :: atol, rtol
        real(c_double) :: diff, magn, thr

        magn = max(abs(a), abs(b))
        diff = abs(a - b)
        thr = get_threshold(magn, atol, rtol)
        result = diff <= thr
    end function

    !> assert_close() for scalars a, b
    subroutine zassert_close_0(a, b, atol, rtol, file, line)
        complex(c_double_complex), intent(in) :: a, b
        real(c_double), intent(in), optional :: atol, rtol
        character(len=*), intent(in), optional :: file
        integer, intent(in), optional :: line

        real(c_double) :: diff, magn

        magn = max(abs(a), abs(b))
        diff = abs(a - b)
        if (diff <= get_threshold(magn, atol, rtol)) &
            return

        call print_array(a, name='a')
        call print_array(b, name='b')
        call write_deviation(diff, magn, atol, rtol)
        call terminate('significant difference between scalars', file, line)
    end subroutine

    !> assert_close() for vectors a, b
    subroutine zassert_close_1(a, b, atol, rtol, file, line)
        complex(c_double_complex), intent(in) :: a(:), b(:)
        real(c_double), intent(in), optional :: atol, rtol
        character(len=*), intent(in), optional :: file
        integer, intent(in), optional :: line

        real(c_double) :: diff, magn

        call check_shape(shape(a), shape(b))
        magn = max(maxval(abs(a)), maxval(abs(b)))
        diff = maxval(abs(a - b))
        if (diff <= get_threshold(magn, atol, rtol)) &
            return

        call print_array(a, name='a')
        call print_array(b, name='b')
        call write_deviation(diff, magn, atol, rtol)
        call terminate('significant difference between vectors', file, line)
    end subroutine

    !> assert_close() for matrices a, b
    subroutine zassert_close_2(a, b, atol, rtol, file, line)
        complex(c_double_complex), intent(in) :: a(:,:), b(:,:)
        real(c_double), intent(in), optional :: atol, rtol
        character(len=*), intent(in), optional :: file
        integer, intent(in), optional :: line

        real(c_double) :: diff, magn

        call check_shape(shape(a), shape(b))
        magn = max(maxval(abs(a)), maxval(abs(b)))
        diff = maxval(abs(a - b))
        if (diff <= get_threshold(magn, atol, rtol)) &
            return

        call print_array(a, name='a')
        call print_array(b, name='b')
        call write_deviation(diff, magn, atol, rtol)
        call terminate('significant difference between vectors', file, line)
    end subroutine

    subroutine read_ftau(filename, norbitals, beta, ftau_full, Nftau)
        ! input
        character(len=*), intent(in) :: filename
        integer, intent(in) :: norbitals
        real(c_double), intent(in) :: beta
        ! output
        complex(c_double_complex), allocatable, intent(out) :: ftau_full(:,:,:,:,:)
        !integer, intent(out) :: Nftau_check
        integer, intent(in) :: Nftau
        ! local
        integer :: i, j, n
        real(c_double) :: a, b, c
        integer :: iflav1, iflav2
        integer :: ib1, ib2, is1, is2

        ! calculate number of tau points from beta and the
        ! distance between first and second tau point

        write(*,*) "filename:", filename
        open (unit=66, file=filename, form="formatted", status="old")
        read(66, *)
        read(66, "(I3, I4, F13.8, F13.8)") i, j, a, b
        read(66, "(I3, I4, F13.8, F13.8)") i, j, a, b
        close(66)

        !tauf1 = a
        !Nftau_check = anint(beta / tauf1) + 1
        !write(*,*) "Nftau_check", Nftau_check

        ! read stuff and write it into array

        open (66, file=filename, form="formatted", status="old")
        read(66, *)

        allocate(ftau_full(norbitals, 2, norbitals, 2, Nftau))

        do ib1 = 1,norbitals
        do is1 = 1,2
        do ib2 = 1,norbitals
        do is2 = 1,2

            iflav1 = 2*(ib1 - 1) + is1
            iflav2 = 2*(ib2 - 1) + is2
            write(*,*) "iflav1, iflav2", iflav1, iflav2
            do n = 1,Nftau
                read(66, "(I4, I4, F11.5, F15.7, F14.7)") i, j, a, b, c
                ftau_full(ib1, is1, ib2, is2, n) = dcmplx(b, c)
                !write(*,*) "n, a, b, c", n, a, b, c
                !write(*,"(A10, I4, I4, F11.5, F15.7, F14.7)") "i,j,a,b,c ", i,j,a,b,c
            enddo
            if(.not.(iflav1.eq.norbitals*2.and.iflav2.eq.norbitals*2))then
                read(66, *)
                read(66, *)
            endif
        enddo
        enddo
        enddo
        enddo

    end subroutine read_ftau

    subroutine dump_ftau(filename, norbitals, beta, ftau_full, Nftau)
        ! input
        character(len=*), intent(in) :: filename
        integer, intent(in) :: norbitals
        real(c_double), intent(in) :: beta
        ! output
        complex(c_double_complex), allocatable, intent(in) :: ftau_full(:,:,:,:,:)
        integer, intent(in) :: Nftau
        ! local
        integer :: i, j, n, ib
        real(c_double) :: a, b
        real(c_double) :: tauf1
        integer :: iflav1, iflav2
        integer :: ib1, ib2, is1, is2
        real(c_double) :: tau

        ! calculate number of tau points from beta and the
        ! distance between first and second tau point

        write(*,*) "filename:", filename
        open (newunit=ib, file=filename, form="formatted", status="old")
        write(ib, *)
        write(ib, "(I3, I4, F13.8, F13.8)") i, j, a, b
        write(ib, "(I3, I4, F13.8, F13.8)") i, j, a, b
        close(ib)

        !tauf1 = a
        !Nftau = anint(beta / tauf1) + 1
        !write(*,*) "Nftau_check", Nftau_check
        tauf1 = beta / float(Nftau - 1)

        ! write stuff and write it into array

        open(66, file = filename)
        write(66, *)

        !allocate(ftau_full(norbitals, 2, norbitals, 2, Nftau_check))

        do ib1 = 1,norbitals
        do is1 = 1,2
        do ib2 = 1,norbitals
        do is2 = 1,2

            iflav1 = 2*(ib1 - 1) + is1
            iflav2 = 2*(ib2 - 1) + is2
            do n = 1,Nftau
            tau = (n-1)*tauf1
            write(66, "(I4, I4, F11.5, F15.7, F14.7)") iflav1, iflav2, tau, ftau_full(ib1, is1, ib2, is2, n), 0.0d0
                !ftau_full(ib1, is1, ib2, is2, n) = dcmplx(b, c)
                !write(*,"(A10, I4, I4, F11.5, F15.7, F14.7)") "i,j,a,b,c ", i,j,a,b,c
            enddo
            if(.not.(iflav1.eq.norbitals*2.and.iflav2.eq.norbitals*2))then
                write(66, *)
                write(66, *)
            endif
        enddo
        enddo
        enddo
        enddo

    end subroutine dump_ftau

end module testing
