module MLocalBase
    use MCommon
    implicit none
    private

    !> A type storing all information describing one impurity operator.
    type, public, extends(TBaseOper) :: TLocalOper
        !> Operator type, one of OpCrea, OpAnnh, OpCreaW, OpAnnhW
        integer(kind(OpDummy)) :: type
    end type TLocalOper

    interface operator ( == )
        module procedure isequal_local_oper
    end interface
    public operator ( == )

    interface swap
       module procedure p_swap
    end interface
    public swap

    interface dump
        module procedure p_dump
    end interface
    public dump

    public nflavcaeq
    public nflavcamatch
    public isworm

contains

    pure elemental function isequal_local_oper(op1, op2) result(test)
        type(TLocalOper), intent(in)  :: op1, op2
        logical                       :: test

        test = op1%tau == op2%tau &
                .and. op1%orb == op2%orb &
                .and. op1%sp == op2%sp &
                .and. op1%type == op2%type
    end function

    !> Test whether flavors and types (creator/annihilator) of two
    !> operators are equal.
    !>
    !> @param op1 First operator
    !> @param op2 Second operator
    pure function nflavcaeq(op1, op2) result(test)
        type(TLocalOper), intent(in) :: op1, op2
        logical :: test

        test = op1%orb == op2%orb &
                .and. op1%sp == op2%sp &
                .and. op1%type == op2%type
    end function

    !> Test whether flavors of two operators are equal, but types are
    !> opposite (density-like "matching" pair, but tau irrelevant).
    !>
    !> @param op1 First operator
    !> @param op2 Second operator
    pure function nflavcamatch(op1, op2) result(test)
        type(TLocalOper), intent(in) :: op1, op2
        logical :: test

        test = op1%orb == op2%orb &
                .and. op1%sp == op2%sp &
                .and. op1%type == invtype(op2%type)
    end function

    pure function isworm(op) result(test)
        type(TLocalOper), intent(in) :: op
        logical :: test

        test = iswormtype(op%type)
    end function

    pure subroutine p_swap(a, b)
        type(TLocalOper), intent(inout) :: a, b
        type(TLocalOper) :: tmp

        tmp = a
        a = b
        b = tmp
    end subroutine

    subroutine p_dump(op, unit)
        type(TLocalOper), intent(in) :: op(:)
        integer, optional :: unit
        integer :: i, runit
        character(len=4), parameter :: spin_str(2) = (/ 'up  ', 'down' /)
        character(len=3), parameter :: type_str(-2:2) = &
                (/ '[-]', ' - ', ' ? ', ' + ', '[+]' /)

        if (present(unit)) then
            runit = unit
        else
            runit = 0
        endif
        write (runit, 98) size(op)
        do i = 1, size(op)
            write (runit, 99) i, op(i)%tau, type_str(op(i)%type), &
                              op(i)%orb, spin_str(op(i)%sp)
        enddo
    98  format ('List of ', I0, ' operators:')
    99  format ('op(', I3, ') = ', F20.10, X, A, X, I3, X, A)
    end subroutine

end module MLocalBase
