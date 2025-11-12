!> This module contains a type that stores some symmetry moves
module MSymmetryMoves
    use MParameters
    use MPermutation
    implicit none
    private

    type, public :: TSymmetryMoves
        !> The number of symmtry moves specified in the paramter file.
        integer                          :: NSymMove
        !> The symmetry moves specified in the paramter file.
        integer,allocatable              :: SymMoves(:,:)
    end type TSymmetryMoves

    public :: init_symmetry_moves, init_symmetry_moves_from_parameters

contains

    !> The following subroutine initialises the type trace. It has to do a lot
    !! more than the equivalent functions in other modules since we need to
    !! compute some more information about the system before we can start the
    !! simulation, e.g. find all the states we want to consider in the outer
    !! trace.
    subroutine init_symmetry_moves_from_parameters(this)
        type(TSymmetryMoves)                        :: this
        !local
        integer, allocatable                :: user_symmoves(:, :)

        allocate(user_symmoves, source=symmetry_moves_from_parameters())
        call init_symmetry_moves(this, user_symmoves)
    end subroutine

    !> The following subroutine initialises the type trace. It has to do a lot
    !! more than the equivalent functions in other modules since we need to
    !! compute some more information about the system before we can start the
    !! simulation, e.g. find all the states we want to consider in the outer
    !! trace.
    subroutine init_symmetry_moves(this, user_symmoves)
        type(TSymmetryMoves), intent(out) :: this
        integer, intent(in)       :: user_symmoves(:, :)

        ! check validity of user-specified moves (every flavor mapped to unique
        ! and valid flavor)
        allocate(this%SymMoves, source=process_symmetry_moves(user_symmoves))
        this%NSymMove = size(this%SymMoves, 1)
    end subroutine

    function symmetry_moves_from_parameters() result(symmoves)
        integer, allocatable :: symmoves(:, :)
        character(len=9) :: IDSymMove
        character(len=2) :: ISymMove
        integer :: i, nbands, nsymmoves

        nbands = get_Integer_Parameter("Nd")
        if (nbands <= 0) &
            error stop 'invalid number of bands'

        nsymmoves = get_Integer_Parameter("NSymMove")
        if (nsymmoves < 0 .or. nsymmoves > 99) &
            error stop 'invalid number of symmetry moves'

        allocate(symmoves(nsymmoves, 2 * nbands))
        do i = 1, nsymmoves
            write (ISymMove, '(I2.2)') i
            IDSymMove = "SymMove"//ISymMove
            call get_Integer_List(IDSymMove, symmoves(i, :))

            if (.not. is_permutation(symmoves(i, :))) then
                write (0, 99) symmoves(i, :)
    99          format ('WARNING: not valid symmetry move: ', 40(I0, ', '))
            endif
        end do
    end function

    function process_symmetry_moves(user_symmoves) result(symmoves)
        integer, intent(in) :: user_symmoves(:, :)
        integer, allocatable :: symmoves(:, :)
        logical :: valid(size(user_symmoves, 1))
        integer :: i, ivalid, nvalid

        do i = 1, size(user_symmoves, 1)
            valid(i) = is_permutation(user_symmoves(i, :))
        enddo

        nvalid = count(valid)
        allocate(symmoves(2 * nvalid, size(user_symmoves, 2)))

        ivalid = 0
        do i = 1, size(user_symmoves, 1)
            if (valid(i)) then
                ivalid = ivalid + 1
                symmoves(ivalid, :) = user_symmoves(i, :)
                symmoves(nvalid + ivalid, :) = get_inv_perm(user_symmoves(i, :))
            endif
        enddo
    end

end module MSymmetryMoves
