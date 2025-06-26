module MPrinting
    use iso_c_binding, only: c_double, c_double_complex
    implicit none
    private

    !>     call print_array(obj, fmt=..., name='X', file=6)
    !!
    !! Interface for printing arrays.
    interface print_array
        module procedure print_inum, print_ivec, print_imat
        module procedure print_dnum, print_dvec, print_dmat
        module procedure print_znum, print_zvec, print_zmat
    end interface
    public print_array

    !> Return format edit descriptor for printing a matrix in python forat.
    !!
    !! Expects the dimensions of a matrix and a format edit descriptor for each
    !! individual item and returns a format edit descriptor for printing the
    !! name plus a matrix.  Given, e.g., a real matrix `a`, one can write
    !!
    !!      write (*, format_matrix(size(a,1), size(a,2), 'ES11.3')) 'A', A
    !!
    !! and this will print the matrix.
    public format_matrix

    !>     default(value, default_value)
    !!
    !! Get `value` if present, and `default_value` if not present.  This can
    !! be used to safely get defaults for optional arguments.
    interface default
        module procedure idefault, ddefault, sdefault
    end interface
    public default

contains
    subroutine print_inum(num, fmt, name, file)
        integer, intent(in) :: num
        character(len=*), intent(in), optional :: name, fmt
        integer, intent(in), optional :: file
        character(len=:), allocatable :: mfmt

        mfmt = format_matrix(1, 1, default(fmt, 'I10'))
        write (default(file, 6), mfmt) default(name, 'X'), num
    end subroutine

    subroutine print_ivec(vec, fmt, name, file)
        integer, intent(in) :: vec(:)
        character(len=*), intent(in), optional :: name, fmt
        integer, intent(in), optional :: file
        character(len=:), allocatable :: mfmt

        mfmt = format_matrix(1, size(vec,1), default(fmt, 'I10'))
        write (default(file, 6), mfmt) default(name, 'X'), vec
    end subroutine

    subroutine print_imat(mat, fmt, name, file)
        integer, intent(in) :: mat(:,:)
        character(len=*), intent(in), optional :: name, fmt
        integer, intent(in), optional :: file
        character(len=:), allocatable :: mfmt

        mfmt = format_matrix(size(mat,1), size(mat,2), default(fmt, 'I10'))
        write (default(file, 6), mfmt) default(name, 'X'), mat
    end subroutine

    subroutine print_dnum(num, fmt, name, file)
        real(c_double), intent(in) :: num
        character(len=*), intent(in), optional :: name, fmt
        integer, intent(in), optional :: file
        character(len=:), allocatable :: mfmt

        mfmt = format_matrix(1, 1, default(fmt, 'ES11.3'))
        write (default(file, 6), mfmt) default(name, 'X'), num
    end subroutine

    subroutine print_dvec(vec, fmt, name, file)
        real(c_double), intent(in) :: vec(:)
        character(len=*), intent(in), optional :: name, fmt
        integer, intent(in), optional :: file
        character(len=:), allocatable :: mfmt

        mfmt = format_matrix(1, size(vec,1), default(fmt, 'ES11.3'))
        write (default(file, 6), mfmt) default(name, 'X'), vec
    end subroutine

    subroutine print_dmat(mat, fmt, name, file)
        real(c_double), intent(in) :: mat(:,:)
        character(len=*), intent(in), optional :: name, fmt
        integer, intent(in), optional :: file
        character(len=:), allocatable :: mfmt

        mfmt = format_matrix(size(mat,1), size(mat,2), default(fmt, 'ES11.3'))
        write (default(file, 6), mfmt) default(name, 'X'), transpose(mat)
    end subroutine

    subroutine print_znum(num, fmt, name, file)
        complex(c_double_complex), intent(in) :: num
        character(len=*), intent(in), optional :: name, fmt
        integer, intent(in), optional :: file
        character(len=:), allocatable :: mfmt

        mfmt = format_matrix(1, 1, default(fmt, 'ES11.3'))
        write (default(file, 6), mfmt) default(name, 'X'), num
    end subroutine

    subroutine print_zvec(vec, fmt, name, file)
        complex(c_double_complex), intent(in) :: vec(:)
        character(len=*), intent(in), optional :: name, fmt
        integer, intent(in), optional :: file
        character(len=:), allocatable :: mfmt
        character(len=*), parameter :: def_fmt = "'(',ES11.3,',',ES11.3,')'"

        mfmt = format_matrix(1, size(vec,1), default(fmt, def_fmt))
        write (default(file, 6), mfmt) default(name, 'X'), vec
    end subroutine

    subroutine print_zmat(mat, fmt, name, file)
        complex(c_double_complex), intent(in) :: mat(:,:)
        character(len=*), intent(in), optional :: name, fmt
        integer, intent(in), optional :: file
        character(len=:), allocatable :: mfmt
        character(len=*), parameter :: def_fmt = "'(',ES11.3,',',ES11.3,')'"

        mfmt = format_matrix(size(mat,1), size(mat,2), default(fmt, def_fmt))
        write (default(file, 6), mfmt) default(name, 'X'), transpose(mat)
    end subroutine

    function format_matrix(rows, cols, item) result(fmt)
        integer, intent(in) :: rows, cols
        character(len=*), intent(in) :: item
        character(len=:), allocatable :: fmt, fmt_row

        allocate(character(len=30 + 2 * len(item)) :: fmt_row)
        if (cols <= 0) then
            write(fmt_row, 10)
        else if (cols == 1) then
            write(fmt_row, 11) item
        else
            write(fmt_row, 12) cols-1, item, item
        endif

        allocate(character(len=30 + 2 * len(fmt_row)) :: fmt)
        if (rows <= 0) then
            write(fmt, 20)
        else if (rows == 1) then
            write(fmt, 21) fmt_row
        else
            write(fmt, 22) rows - 1, fmt_row, fmt_row
        endif

    10  format("'[]'")
    11  format("'[',", A, ",']'")
    12  format("'[',", I0, "(", A, ",','), ", A, ",']'")
    20  format("(A,' = []')")
    21  format("(A,' = [',", A, ",']')")
    22  format("(A,' = [',/,", I0, "(", A, ",',',/),", A, ",']')")
    end

    pure function idefault(val, default_val) result(result)
        integer, intent(in), optional :: val
        integer, intent(in) :: default_val
        integer :: result

        if (present(val)) then
            result = val
        else
            result = default_val
        endif
    end

    pure function ddefault(val, default_val) result(result)
        real(c_double), intent(in), optional :: val
        real(c_double), intent(in) :: default_val
        real(c_double) :: result

        if (present(val)) then
            result = val
        else
            result = default_val
        endif
    end

    pure function sdefault(val, default_val) result(result)
        character(len=*), intent(in), optional :: val
        character(len=*), intent(in) :: default_val
        character(len=:), allocatable :: result

        if (present(val)) then
            result = val
        else
            result = default_val
        endif
    end

end module