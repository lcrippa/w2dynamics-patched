module MCommon
    use iso_c_binding, only: c_double
    use MSorting
    implicit none
    private

    ! operator types (debugging dummy, creator, annihilator,
    ! composed equal-time worms)
    enum, bind(c)
        enumerator :: OpDummy = 0, OpCrea = 1, OpAnnh = -1, OpCreaW = 2, OpAnnhW = -2
    end enum

    public :: OpDummy, OpCrea, OpAnnh, OpCreaW, OpAnnhW

    !> A type representing one operator in terms of imaginary time and flavor.
    type, public :: TBaseOper
        real(c_double) :: tau  !< Imaginary time
        integer     :: orb  !< Orbital quantum number
        integer     :: sp   !< Spin quantum number
    end type TBaseOper

    public :: nflaveq, invca, invtype
    public :: iswormtype, oldtype, binarysearch_taupos

contains

    !> Test whether flavors of two operators are equal.
    !>
    !> @param op1 First operator
    !> @param op2 Second operator
    pure logical function nflaveq(op1, op2)
        type(TBaseOper), intent(in) :: op1, op2

        nflaveq = op1%orb == op2%orb .and. op1%sp == op2%sp
    end function nflaveq

    !> Get opposite 'type' number (1 corr. to creator for 2 corr. to
    !> annihilator and vice versa).
    !> OBSOLETE: use invtype for new local types
    !>
    !> @param ca Type number (1 or 2)
    pure integer function invca(ca)
        integer, intent(in) :: ca
        invca = 3 - ca
    end function invca

    pure integer(kind(OpDummy)) function invtype(type)
        integer(kind(OpDummy)), intent(in) :: type
        invtype = -type
    end function invtype

    pure logical function iswormtype(type)
        integer(kind(OpDummy)), intent(in) :: type

        if (type == OpCrea .or. type == OpAnnh) then
            iswormtype = .false.
        else
            iswormtype = .true.
        end if
    end function iswormtype

    !> Get old type corresponding to new type
    pure integer function oldtype(type)
        integer(kind(OpDummy)), intent(in) :: type
        select case (type)
        case (OpCrea, OpCreaW)
            oldtype = 1
        case (OpAnnh, OpAnnhW)
            oldtype = 2
        case default
            oldtype = -1
        end select
    end function oldtype

    !> Return the smallest possible (one-based) index between first and
    !> last + 1 (both inclusive) where time tau would fit into the
    !> time-ordered array ops such that it remains time-ordered. Note
    !> that first is returned if it already indexes beyond the end of
    !> the array or is greater than last.
    !
    !> UNCHECKED PRECONDITIONS:
    !> - 1 <= first <= size(ops) + 1
    !> - 0 <= last <= size(ops)
    pure integer function binarysearch_taupos(ops, tau) result(inspos)
        class(TBaseOper), intent(in) :: ops(:)
        real(c_double), intent(in) :: tau
        real(c_double) :: ops_tau(size(ops))

        ops_tau(:) = ops%tau
        inspos = get_sorted_insert(ops_tau, nearest(tau, -1.0d0))
        inspos = inspos
    end function binarysearch_taupos

end module MCommon
