module MDiagram
    use iso_c_binding, only: c_double, c_double_complex
    use MBath
    use MCommon
    use MSorting
    use MImpurity
    use MPermutation
    use MWormState
    implicit none
    private

    ! Re-exports
    public :: TImpurity, TBath, TWormState
    public :: TLocalOper, TBathOperator
    public :: ZExtRange

    type, public :: TDiagram
        !> State of the local trace (cumulant)
        type(TImpurity)                 :: local

        !> State of the bath trace (determinant)
        type(TBath)                       :: bath

        !> State of the worm sector and position of worm operators
        type(TWormState)                  :: worm

        !> Weight of the full configuration
        type(ZExtRange)                   :: weight

        !> The "natural" order of the trace is the order that the operators
        !! appear in the strong-coupling expansion, i.e., before the
        !! time-ordering symbol permutes them.  The natural order is:
        !!
        !!  - first, worm operators (if present) in the order they appear in the
        !!    corelation function, e.g., for the Green's function <c c+>, first
        !!    c then c+
        !!
        !!  - afterwards, alternating c+ and c oparators.  The ordering of c
        !!    and c+ operators amongst themselves must match the ordering in
        !!    the bath (usually ordered first by flavor, then time)
        !!
        !! This array contains the mapping from time-ordered to natural
        integer, allocatable :: natural_idx(:)

        !> The time-ordered trace is simply all operators ordered by time.
        !! In case of a tie, which can happen for worm operators, the natural
        !! order is preserved.
        !!
        !! This array contains the mapping from natural to time-ordered
        integer, allocatable :: timeord_idx(:)

        !> Control whether to attempt inserting offdiagonal hybridization lines
        logical :: b_offdiag

        !> Stuff to cache
        real(c_double) :: beta
        integer :: nbands
    end type TDiagram

    type, public :: TAddMove
        ! FIXME: private
        type(TDiagram), pointer :: target => null()
        real(c_double)           :: prop_prob=0.0
        type(TImpurityGrow) :: local
        type(TBathAddMove)  :: bath

        type(TLocalOper), allocatable  :: newops(:)
        type(TBathOperator), allocatable :: crea(:), annh(:)
    end type TAddMove

    type, public :: TRemMove
        ! FIXME: private
        type(TDiagram), pointer :: target => null()
        real(c_double)           :: prop_prob=0.0
        type(TImpurityShrink) :: local
        type(TBathRemMove)  :: bath
        integer, allocatable :: local_remidx(:)
    end type TRemMove

    type, public :: TShiftMove
        !private
        type(TDiagram), pointer :: target => null()
        real(c_double)           :: prop_prob=0.0
        type(TImpurityShift) :: local
        type(TBathShiftMove)  :: bath
    end type TShiftMove

    type, public :: TPermuteMove
        !private
        type(TDiagram), pointer :: target => null()
        real(c_double)           :: prop_prob=0.0
        type(TImpurityReplace) :: local
        type(TBathMove)  :: bath
    end type TPermuteMove

    ! XXX remove this type
    type, public :: TWormAddMoveSpec
        integer(kind(SecDummy)) :: target_sector
        integer :: target_component
        type(TLocalOper), allocatable :: newops(:)
    end type TWormAddMoveSpec

    type, public :: TWormAddMove
        !private
        type(TDiagram), pointer :: target => null()
        real(c_double) :: prop_prob=0.0
        integer(kind(SecDummy)) :: target_sector
        integer :: target_component
        type(TImpurityGrow) :: local
    end type TWormAddMove

    type, public :: TWormRemMove
        !private
        type(TDiagram), pointer :: target => null()
        real(c_double)           :: prop_prob=0.0
        type(TImpurityShrink) :: local
    end type TWormRemMove

    ! XXX remove this type
    type, public :: TWormReplaceMoveSpec
        ! index in local operator array of the operator that gets
        ! disconnected from hybridization in the move
        integer :: local_repidx=-100
        ! index in 'type order' worm operator index array of the worm
        ! operator that gets connected to hybridization in the move
        integer :: worm_reptypepos=-100
        ! index in local operator array of the worm operator that gets
        ! connected to hybridization in the move
        integer :: local_wormidx=-100
    end type TWormReplaceMoveSpec

    type, public :: TWormReplaceMove
        !private
        type(TDiagram), pointer :: target => null()
        real(c_double)           :: prop_prob=0.0
        type(TWormReplaceMoveSpec) :: spec
        type(TImpurityUpdate) :: local
        integer, allocatable :: local_remidx(:)
        ! type / natural indices of operators in the equal-time object
        ! to be moved (if any)
        integer, allocatable :: eqtau_typeidx(:)
        type(TBathReplaceMove)  :: bath
        ! FIXME: unused? (duplicate of spec entries)
        integer :: local_repidx=-100, worm_reptypepos=-100, local_wormidx=-100
        ! number of operators at equal-time in the worm operator
        ! product that gets moved / replaced (1: single op)
        integer :: neqtauops=-100
        ! first index belonging to the equal time product operators in
        ! newops
        integer :: offset_eqtau_newops=-100
        ! index belonging to the operator to be disconnected from
        ! hybridization (initially at local_repidx) in newops
        integer :: offset_swapop_newops=-100
    end type TWormReplaceMove

    interface weight_ratio
        module procedure :: weight_ratio_addmove
        module procedure :: weight_ratio_remmove
        module procedure :: weight_ratio_shiftmove
        module procedure :: weight_ratio_permutemove
        module procedure :: weight_ratio_wormaddmove
        module procedure :: weight_ratio_wormremmove
        module procedure :: weight_ratio_wormreplacemove
    end interface weight_ratio
    public weight_ratio

    interface accept_update
        module procedure :: accept_addmove
        module procedure :: accept_remmove
        module procedure :: accept_shiftmove
        module procedure :: accept_permutemove
        module procedure :: accept_wormaddmove
        module procedure :: accept_wormremmove
        module procedure :: accept_wormreplacemove
    end interface accept_update
    public accept_update

    public :: init_diagram
    public :: verify_diagram
    public :: clear_diagram, set_diagram_from_local_ops
    public :: propose_addmove, init_addmove
    public :: propose_remmove, init_remmove
    public :: propose_shiftmove, init_shiftmove
    public :: propose_permutemove, init_permutemove
    public :: propose_wormaddmove, init_wormaddmove
    public :: propose_wormremmove, init_wormremmove
    public :: propose_wormreplacemove, init_wormreplacemove
    public :: get_diagram_sign

contains

    subroutine init_diagram(self, local, bath, worm, offdiag_hyb)
        type(TDiagram), intent(out) :: self
        type(TImpurity), intent(in) :: local
        type(TBath), intent(in) :: bath
        type(TWormState), intent(in) :: worm
        logical, intent(in) :: offdiag_hyb

        self%local = local
        self%bath = bath
        self%worm = worm

        self%b_offdiag = offdiag_hyb

        ! Initialize the lookups
        allocate(self%natural_idx(0))
        allocate(self%timeord_idx(0))
        call update_diagram_lookups(self)
        call update_diagram_weight(self)

        self%beta = get_beta(self%local)
        self%nbands = get_nbands(self%bath)
    end

    subroutine verify_mapping(self)
        type(TDiagram), intent(in) :: self

        integer :: i, ibath, n, nworm, local_outer_sst
        type(TLocalOper) :: curr
        type(TLocalOper), allocatable :: local_ops(:)
        type(TBathOperator), allocatable :: crea(:), annh(:)

        allocate(crea, source=get_creators(self%bath))
        allocate(annh, source=get_annihilators(self%bath))
        allocate(local_ops, source=get_operators(self%local))
        local_outer_sst = get_outerbl(self%local)
        n = size(local_ops)
        nworm = get_nworm(self%worm)

        do i = 1, n
            curr = local_ops(i)
            if (self%natural_idx(i) < 1 .or. self%natural_idx(i) > n) &
                error stop 'Invalid natural index'

            ibath = (self%natural_idx(i) - (nworm - 1)) / 2
            if (curr%type == OpCrea) then
                if (self%natural_idx(i) <= nworm) &
                    error stop 'Invalid bath position'
                if (mod(self%natural_idx(i) - nworm, 2) /= 1) &
                    error stop 'Invalid creator position'
                if (curr%orb /= crea(ibath)%orb .or. curr%sp /= crea(ibath)%sp &
                        .or. curr%tau /= crea(ibath)%tau) &
                    error stop 'Mismatch between local and bath operator'
            else if (curr%type == OpAnnh) then
                if (self%natural_idx(i) <= nworm) &
                    error stop 'Invalid bath position'
                if (mod(self%natural_idx(i) - nworm, 2) /= 0) &
                    error stop 'Invalid annihilator position'
                if (curr%orb /= annh(ibath)%orb .or. curr%sp /= annh(ibath)%sp &
                        .or. curr%tau /= annh(ibath)%tau) &
                    error stop 'Mismatch between local and bath operator'
            else
                if (self%natural_idx(i) > nworm) &
                    error stop 'Invalid worm position'
            endif
        enddo
    end

    subroutine verify_diagram(self, rtol)
        type(TDiagram), intent(in) :: self
        real(c_double), intent(in) :: rtol

        call verify_worm(self%worm, get_operators(self%local))
        call verify_mapping(self)
        call verify_impurity(self%local, rtol)
        call verify_bath(self%bath, rtol)

        ! Check things
        if (self%beta /= get_beta(self%local)) &
            error stop 'Wrong beta'
        if (self%nbands /= get_nbands(self%local)) &
            error stop 'Wrong nbands'
        call assert_close( &
                abs(limrange( &
                    get_weight(self%local) * get_det(self%bath) / self%weight)), &
                1.0d0)
    end subroutine

    subroutine bath_order_resort(self, perm)
        type(TDiagram) :: self
        integer, intent(inout) :: perm(:)

        real(c_double), allocatable :: sortval(:)
        type(TLocalOper), allocatable :: local_ops(:)
        type(TLocalOper) :: loc_op
        type(TBathOperator) :: bath_op
        integer :: nworm, nbath, i

        nbath = size(self%bath)
        nworm = get_nworm(self%worm)
        allocate(sortval(nbath))
        allocate(local_ops, source=get_operators(self%local))
        if (size(local_ops) /= size(perm)) &
            error stop 'Wrong size of perm array'

        !write (0,"('perm=',  100(I3,','))") perm
        do i = 1, nbath
            loc_op = local_ops(perm(nworm + 2 * i - 1))
            bath_op = TBathOperator(tau=loc_op%tau, orb=loc_op%orb, sp=loc_op%sp)
            sortval(i) = get_sortval(self%bath, bath_op)
        enddo
        call sort_modify_perm(sortval, perm(nworm+1:nworm+2*nbath-1:2))

        do i = 1, nbath
            loc_op = local_ops(perm(nworm + 2 * i))
            bath_op = TBathOperator(tau=loc_op%tau, orb=loc_op%orb, sp=loc_op%sp)
            sortval(i) = get_sortval(self%bath, bath_op)
        enddo
        call sort_modify_perm(sortval, perm(nworm+2:nworm+2*nbath:2))

        !write (0,"('perm=',  100(I3,','))") perm
    end subroutine

    subroutine update_diagram_lookups(self)
        type(TDiagram), intent(inout) :: self
        type(TLocalOper), allocatable :: ops(:)

        allocate(ops, source=get_operators(self%local))

        if (size(self%natural_idx) /= size(ops)) then
            deallocate(self%natural_idx)
            deallocate(self%timeord_idx)
            allocate(self%natural_idx(size(ops)))
            allocate(self%timeord_idx(size(ops)))
        endif

        call natural_order_worm_perm(ops, self%worm, &
                                     self%natural_idx)
        call inv_perm(self%natural_idx, self%timeord_idx)

        call bath_order_resort(self, self%timeord_idx)
        call inv_perm(self%timeord_idx, self%natural_idx)

        !call verify_diagram(self, 1d-3)
   end

    subroutine update_diagram_weight(diagram)
        type(TDiagram), intent(inout) :: diagram
        type(ZExtRange) :: bath_det, local_w
        complex(c_double_complex), parameter :: one=1.0
        logical :: diagram_parity

        bath_det = get_det(diagram%bath)
        local_w = get_weight(diagram%local)
        diagram%weight = local_w * bath_det

        ! Note that w2dynamics (for historical reasons) uses the convention
        ! that the weights sum up to -Z rather than Z.
        diagram_parity = perm_parity(diagram%natural_idx)
        if (.not. diagram_parity) &
            diagram%weight = -diagram%weight
    end subroutine

    function get_diagram_sign(diagram) result(s)
        type(TDiagram), intent(in) :: diagram
        complex(c_double_complex) :: s
        complex(c_double_complex), parameter :: ONE = 1.0d0

        s = sign(ONE, diagram%weight)
    end function

    !> clears the diagram of all operators
    subroutine clear_diagram(diagram, outer_sst)
        type(TDiagram), intent(inout)       :: diagram
        integer, intent(in), optional       :: outer_sst
        integer :: outer_sst_r

        if (present(outer_sst)) then
            outer_sst_r = outer_sst
        else
            outer_sst_r = get_outerbl(diagram%local)
        endif

        ! Clear local state, but keep outersst
        call clear_impurity(diagram%local, outer_sst_r)
        call clear_bath(diagram%bath)
        call wormstate_leave_wormsector(diagram%worm)

        call update_diagram_lookups(diagram)
        call update_diagram_weight(diagram)
    end subroutine

    !> Create Z-space configuration from data arrays.
    subroutine set_diagram_from_local_ops(diagram, ops, outer_sst)
        type(TDiagram), intent(inout) :: diagram
        type(TLocalOper), intent(in)  :: ops(:)
        integer, intent(in) :: outer_sst
        !local
        integer                                             :: ncrea, nannh
        type(TBathMove)                                     :: bathmove
        type(TBathOperator), allocatable                    :: crea(:), annh(:)
        real(c_double)                                         :: throwaway
        type(TLocalOper), allocatable                       :: local_ops(:)

        call set_impurity(diagram%local, ops, outer_sst)

        ncrea = count(ops(:)%type == OpCrea)
        nannh = count(ops(:)%type == OpAnnh)
        allocate(crea(ncrea), annh(nannh))

        ! XXX figure out if we can feed ops(:) in directly
        allocate(local_ops, source=get_operators(diagram%local))
        call convert_to_bath_arrays(local_ops, crea, annh)
        call init_bath_move(BathMove, diagram%bath)
        call propose_bath_move(bathmove, crea, annh)
        call compute_weight(bathmove)
        throwaway = weight_ratio(bathmove)
        call accept_update(bathmove)

        ! FIXME: change if and when serialized configurations contain such
        ! info and we want to use it, operator creation above would need to
        ! be adapted as well
        call wormstate_leave_wormsector(diagram%worm)
        call update_diagram_lookups(diagram)
        call update_diagram_weight(diagram)
    end subroutine

    !> Convert an array trops of local operators into arrays of creators
    !> and annihilators as used by the bath code.
    subroutine convert_to_bath_arrays(trops, crops, anops)
        type(TLocalOper), intent(in)     :: trops(:)
        type(TBathOperator), intent(out) :: crops(:), anops(:)
        type(TBathOperator) :: op
        integer :: i, j, k

        j = 0
        k = 0
        do i = 1, size(trops)
            op = TBathOperator(tau=trops(i)%tau, orb=trops(i)%orb, sp=trops(i)%sp)
            select case (trops(i)%type)
            case (OpCrea)
                j = j + 1
                if (j > size(crops)) &
                    error stop 'crea ops too small'
                crops(j) = op
            case (OpAnnh)
                k = k + 1
                if (k > size(anops)) &
                    error stop 'annh ops too small'
                anops(k) = op
            case (OpCreaW, OpAnnhW)
                continue
            case default
                error stop 'Cannot translate to bath operator'
            end select
        end do

        if (j /= size(crops) .or. k /= size(anops)) &
            error stop 'Wrong size of output arrays'
    end subroutine

    ! -------------------------------------------------------------------------
    ! ADD MOVE

    subroutine init_addmove(move, base)
        type(TAddMove), intent(out) :: move
        type(TDiagram), intent(inout), target :: base
        integer, parameter :: maxorder = 2

        if (associated(move%target)) &
              error stop 'Already initialized'

        move%target => base
        call init_grow_impurity(move%local, base%local)
        call init_bath_add_move(move%bath, base%bath)

        allocate(move%newops(2*maxorder))
        allocate(move%crea(maxorder), move%annh(maxorder))
    end

    subroutine propose_addmove(newmove, newops)
        type(TAddMove), intent(inout), target :: newmove
        type(TLocalOper), intent(in) :: newops(:)
        type(TDiagram), pointer :: base
        type(TBathOperator), pointer :: crea(:), annh(:)
        integer :: n

        base => newmove%target
        n = size(newops)/2
        if (n > size(newmove%newops)) &
            error stop 'Order of move is too large'

        ! local part
        call propose_grow_impurity(newmove%local, newops(:))
        newmove%newops(:2*n) = newops(:)

        ! bath part
        crea => newmove%crea(:n)
        annh => newmove%annh(:n)
        call convert_to_bath_arrays(newops, crea, annh)
        call propose_bath_add_move(newmove%bath, crea, annh)
     end

    real(c_double) function weight_ratio_addmove(proposal) result(ratio)
        class(TAddMove), intent(inout) :: proposal
        type(TDiagram), pointer :: base
        real(c_double) :: localratio

        base => proposal%target
        localratio = weight_ratio(proposal%local)
        if (localratio == 0) then
            ratio = 0
            return
        end if

        ratio = weight_ratio(proposal%bath) * localratio
    end function weight_ratio_addmove

    subroutine accept_addmove(proposal)
        class(TAddMove), intent(inout) :: proposal
        type(TDiagram), pointer :: tgtstate
        integer, allocatable :: new_idx(:)

        tgtstate => proposal%target
        call accept_update(proposal%local)
        call accept_update(proposal%bath)
        call get_newidx_addorder(proposal%local, new_idx)
        call update_wormstate_add_opers(tgtstate%worm, new_idx)
        call update_diagram_lookups(tgtstate)
        call update_diagram_weight(tgtstate)
        !call verify_worm(tgtstate%worm, get_operators(tgtstate%local))
    end subroutine accept_addmove

    ! -------------------------------------------------------------------------
    ! REMOVE MOVE

    subroutine init_remmove(move, base)
        type(TRemMove), intent(out) :: move
        type(TDiagram), intent(inout), target :: base

        if (associated(move%target)) &
                error stop 'Already initialized'

        move%target => base
        call init_shrink_impurity(move%local, base%local)
        call init_bath_rem_move(move%bath, base%bath)

        allocate(move%local_remidx(1))
    end

    subroutine propose_remmove(newmove, local_remidx)
        type(TRemMove), intent(inout) :: newmove
        integer, intent(in)           :: local_remidx(:)
        type(TDiagram), pointer :: base

        base => newmove%target
        deallocate(newmove%local_remidx)

        ! bath part: create move after getting bath spec from diagram (= local) spec
        call bath_remmove_from_local_remidx(newmove%bath, base, local_remidx)

        newmove%target => base
        newmove%local_remidx = local_remidx(:)

        ! for more than 1 pair, sort; single pairs are already ordered
        ! correctly (assuming that generate_remove is used)
        if (size(newmove%local_remidx) > 2) call sort(newmove%local_remidx)

        ! local part: copy spec into local spec and create move
        call propose_shrink_impurity(newmove%local, newmove%local_remidx)
    end

    subroutine bath_remmove_from_local_remidx(bathprop, from_diagram, remidx)
        type(TBathRemMove), intent(inout) :: bathprop
        type(TDiagram), intent(in)      :: from_diagram

        ! must be valid index into local ops, must not point to worms
        integer, intent(in)             :: remidx(:)

        integer :: bath_remidx(size(remidx))
        integer :: i, nworm, cur_crea, cur_annh, itype, bathidx

        type(TLocalOper), allocatable :: local_ops(:)

        nworm = get_nworm(from_diagram%worm)
        allocate(local_ops, source=get_operators(from_diagram%local))
        cur_crea = 1
        cur_annh = 2
        do i = 1, size(remidx)
            if (remidx(i) < 1 .or. remidx(i) > size(local_ops)) &
                error stop 'Inexistent index'
            if (remidx(i) > size(from_diagram%natural_idx)) &
                error stop 'Index array not properly sized'

            itype = local_ops(remidx(i))%type
            bathidx = (from_diagram%natural_idx(remidx(i)) - (nworm - 1)) / 2
            if (itype == OpCrea) then
                bath_remidx(cur_crea) = bathidx
                cur_crea = cur_crea + 2
            elseif (itype == OpAnnh) then
                bath_remidx(cur_annh) = bathidx
                cur_annh = cur_annh + 2
            else
                error stop 'Invalid operator type'
            endif
        enddo

        call propose_bath_rem_move(bathprop, bath_remidx)
    end subroutine bath_remmove_from_local_remidx

    real(c_double) function weight_ratio_remmove(proposal)
        class(TRemMove), intent(inout) :: proposal
        type(TDiagram), pointer :: base
        real(c_double) :: localratio

        base => proposal%target

        localratio = weight_ratio(proposal%local)
        if (localratio == 0) then
            weight_ratio_remmove = 0
            return
        end if

        weight_ratio_remmove = weight_ratio(proposal%bath) * localratio
    end function weight_ratio_remmove

    subroutine accept_remmove(proposal)
        class(TRemMove), intent(inout) :: proposal
        type(TDiagram), pointer :: tgtstate

        tgtstate => proposal%target
        call accept_update(proposal%local)
        call accept_update(proposal%bath)
        call wormstate_remove_hybops(tgtstate%worm, proposal%local_remidx)
        call update_diagram_lookups(tgtstate)
        call update_diagram_weight(tgtstate)
    end subroutine accept_remmove

    ! -------------------------------------------------------------------------
    ! SHIFT MOVE

    subroutine init_shiftmove(move, base)
        type(TShiftMove), intent(out) :: move
        type(TDiagram), intent(inout), target :: base

        if (associated(move%target)) &
                error stop 'Already initialized'

        move%target => base
        call init_shift_impurity(move%local, base%local)
        call init_bath_shift_move(move%bath, base%bath)
    end

    subroutine propose_shiftmove(newmove, deltatau)
        type(TShiftMove), intent(inout) :: newmove
        real(c_double) :: deltatau
        type(TDiagram), pointer :: base

        base => newmove%target

        ! local part: copy spec into local spec and create move
        call propose_shift_impurity(newmove%local, deltatau)

        ! bath part: create from spec
        call propose_bath_shift_move(newmove%bath, deltatau)
    end

    real(c_double) function weight_ratio_shiftmove(proposal) result(ratio)
        class(TShiftMove), intent(inout) :: proposal
        type(TDiagram), pointer :: base

        base => proposal%target
        ratio = weight_ratio(proposal%local)
    end function weight_ratio_shiftmove

    subroutine accept_shiftmove(proposal)
        class(TShiftMove), intent(inout) :: proposal
        type(TDiagram), pointer :: tgtstate

        tgtstate => proposal%target
        call accept_update(proposal%local)
        call accept_update(proposal%bath)
        call wormstate_taushift(tgtstate%worm, &
                     size(tgtstate%local), get_num_wrapped_ops(proposal%local))
        call update_diagram_lookups(tgtstate)
        call update_diagram_weight(tgtstate)
    end subroutine accept_shiftmove

    ! -------------------------------------------------------------------------
    ! PERMUTE MOVE

    subroutine init_permutemove(move, base)
        type(TPermuteMove), intent(out) :: move
        type(TDiagram), intent(inout), target :: base

        if (associated(move%target)) &
                error stop 'Already initialized'

        move%target => base
        call init_replace_impurity(move%local, base%local)
        call init_bath_move(move%bath, base%bath)
    end

    subroutine propose_permutemove(newmove, caperm, flavperm)
        type(TPermuteMove), intent(inout) :: newmove
        logical, intent(in) :: caperm
        integer, intent(in) :: flavperm(:)
        type(TDiagram), pointer :: base

        base => newmove%target

        ! local part: copy spec into local spec and create move
        call propose_permute_impurity_flavours(newmove%local, caperm, flavperm)

        ! bath part: create from changed local operators
        call bath_move_from_local_prop(newmove%bath, base, newmove%local)
    end

    subroutine bath_move_from_local_prop(bathprop, from_diagram, localprop)
        type(TBathMove), intent(inout) :: bathprop
        type(TDiagram), intent(inout), target :: from_diagram
        type(TImpurityReplace), intent(in) :: localprop

        type(TLocalOper), allocatable :: new_ops(:)
        type(TBathOperator), allocatable :: crea(:), annh(:)
        integer :: ncrea, nannh

        allocate(new_ops, source=get_operators(localprop))
        ncrea = count(new_ops(:)%type == OpCrea)
        nannh = count(new_ops(:)%type == OpAnnh)
        allocate(crea(ncrea), annh(nannh))

        call convert_to_bath_arrays(new_ops, crea, annh)
        call propose_bath_move(bathprop, crea, annh)
    end subroutine bath_move_from_local_prop

    real(c_double) function weight_ratio_permutemove(proposal) result(ratio)
        class(TPermuteMove), intent(inout) :: proposal
        type(TDiagram), pointer :: base
        real(c_double) :: localratio

        base => proposal%target
        localratio = weight_ratio(proposal%local)
        if (localratio == 0) then
            ratio = 0
            return
        end if

        ratio = weight_ratio(proposal%bath) * localratio
    end function weight_ratio_permutemove

    subroutine accept_permutemove(proposal)
        class(TPermuteMove), intent(inout) :: proposal
        type(TDiagram), pointer :: tgtstate

        tgtstate => proposal%target
        call accept_update(proposal%local)
        call accept_update(proposal%bath)
        call update_diagram_lookups(tgtstate)
        call update_diagram_weight(tgtstate)
    end subroutine accept_permutemove

    ! =========================================================================
    ! WORM ADD MOVE

    subroutine init_wormaddmove(move, base)
        type(TWormAddMove), intent(out) :: move
        type(TDiagram), intent(inout), target :: base

        if (associated(move%target)) &
                error stop 'Already initialized'

        move%target => base
        call init_grow_impurity(move%local, base%local)
    end

    subroutine propose_wormaddmove(newmove, propspec)
        type(TWormAddMove), intent(inout) :: newmove
        type(TWormAddMoveSpec), intent(in) :: propspec
        type(TDiagram), pointer :: base

        base => newmove%target

        newmove%target_sector = propspec%target_sector
        newmove%target_component = propspec%target_component

        ! local part: copy spec into local spec and create move
        call propose_grow_impurity(newmove%local, propspec%newops(:))
    end

    real(c_double) function weight_ratio_wormaddmove(proposal) result(ratio)
        class(TWormAddMove), intent(inout) :: proposal
        type(TDiagram), pointer :: base

        base => proposal%target
        ratio = weight_ratio(proposal%local)
    end function weight_ratio_wormaddmove

    subroutine accept_wormaddmove(proposal)
        class(TWormAddMove), intent(inout) :: proposal
        type(TDiagram), pointer :: tgtstate
        integer, allocatable :: wormpos(:)

        tgtstate => proposal%target
        call accept_update(proposal%local)

        ! Generate array wormpos of indices of worm operators in local
        ! operator list in 'type order', sector dependent order that
        ! allows using the indices contained by the array to easily
        ! access specific worm operator (or operator group) for purposes
        ! of measurement etc.
        select case (proposal%target_sector)
        case (SecP2, SecP2pp, SecQQ, SecRaman, SecCustom)
            ! operators were in desired order in propspec%newops, which
            ! was set up by set_optypes, so the indices in the new list
            ! in order of the new operator array can be used directly
            call get_newidx_addorder(proposal%local, wormpos)
        case default
            error stop "accept_wormaddmove: Sector unimplemented"
        end select

        call wormstate_enter_wormsector(tgtstate%worm, proposal%target_sector, &
                     proposal%target_component, wormpos)
        call update_diagram_lookups(tgtstate)
        call update_diagram_weight(tgtstate)
    end subroutine accept_wormaddmove

    ! -------------------------------------------------------------------------
    ! WORM REMOVE MOVE

    subroutine init_wormremmove(move, base)
        type(TWormRemMove), intent(out) :: move
        type(TDiagram), intent(inout), target :: base

        if (associated(move%target)) &
                error stop 'Already initialized'

        move%target => base
        call init_shrink_impurity(move%local, base%local)
    end

    subroutine propose_wormremmove(newmove)
        type(TWormRemMove), intent(inout) :: newmove
        type(TDiagram), pointer :: base

        base => newmove%target

        ! local part: copy spec into local spec and create move
        call propose_shrink_impurity(newmove%local, &
                     get_worm_oper_positions_time_ordered(base%worm))
    end

    real(c_double) function weight_ratio_wormremmove(proposal) result(ratio)
        class(TWormRemMove), intent(inout) :: proposal
        type(TDiagram), pointer :: base

        base => proposal%target
        ratio = weight_ratio(proposal%local)
    end function weight_ratio_wormremmove

    subroutine accept_wormremmove(proposal)
        class(TWormRemMove), intent(inout) :: proposal
        type(TDiagram), pointer :: tgtstate

        tgtstate => proposal%target
        call accept_update(proposal%local)
        call wormstate_leave_wormsector(tgtstate%worm)
        call update_diagram_lookups(tgtstate)
        call update_diagram_weight(tgtstate)
    end subroutine accept_wormremmove

    ! -------------------------------------------------------------------------
    ! WORM REPLACE MOVE

    subroutine init_wormreplacemove(move, base)
        type(TWormReplaceMove), intent(out) :: move
        type(TDiagram), intent(inout), target :: base

        if (associated(move%target)) &
                error stop 'Already initialized'

        move%target => base
        call init_update_impurity(move%local, base%local)
    end

    subroutine propose_wormreplacemove(newmove, propspec)
        type(TWormReplaceMove), intent(inout) :: newmove
        type(TWormReplaceMoveSpec), intent(in) :: propspec
        type(TDiagram), pointer :: base
        integer :: reptypepos
        integer :: i, iopend, nops_current, ntaus
        integer, allocatable          :: remidx(:)
        type(TLocalOper), allocatable :: newops(:)

        base => newmove%target
        newmove%spec = propspec
        reptypepos = newmove%spec%worm_reptypepos

        ! check whether the selected operator is part of a product of
        ! operators at equal time and if so set the number of operators
        ! in that product neqtauops > 1 and the array eqtau_typeidx to
        ! the operators type indices
        newmove%neqtauops = 1
        select case (get_current_worm_sector(base%worm))
        case (SecGSigma, SecH4, SecQUDdag)
            ! FIXME: check this after implementing these sectors
            if (reptypepos < 4) then
                newmove%neqtauops = 3
                newmove%eqtau_typeidx = [1, 2, 3]
            end if
        case (SecP2, SecP2pp, SecUccaa, SecUcaca, SecRaman)
            ! FIXME: check this after implementing the missing of these sectors
            ! (esp. SecUccaa, assumed changed order vs. old version)
            newmove%neqtauops = 2
            newmove%eqtau_typeidx = [((reptypepos - 1) / 2) * 2 + 1,&
                                         ((reptypepos - 1) / 2) * 2 + 2]
        case (SecQQ, SecQ4)
            ! FIXME: check this after implementing the missing of these sectors
            newmove%neqtauops = 3
            newmove%eqtau_typeidx = [((reptypepos - 1) / 3) * 3 + 1,&
                                         ((reptypepos - 1) / 3) * 3 + 2,&
                                         ((reptypepos - 1) / 3) * 3 + 3]
        case (SecNQQdag, SecQQdd)
            ! FIXME: check this after implementing these sectors
            ! (esp. SecQQdd, assumes changed order vs. old version)
            if (reptypepos < 7) then
                newmove%neqtauops = 3
                newmove%eqtau_typeidx = [((reptypepos - 1) / 3) * 3 + 1,&
                                            ((reptypepos - 1) / 3) * 3 + 2,&
                                            ((reptypepos - 1) / 3) * 3 + 3]
            else
                newmove%neqtauops = 2
                newmove%eqtau_typeidx = [((reptypepos - 1) / 2) * 2 + 1,&
                                            ((reptypepos - 1) / 2) * 2 + 2]
            end if
        case (SecCustom)

            iopend = 0
            ntaus = get_custom_worm_ntaus(base%worm, SecCustom)
            if (ntaus <= 0) &
                error stop 'Invalid ntaus'

            do i = 1, ntaus
                nops_current = get_custom_worm_ntauops(base%worm, &
                                                          SecCustom, &
                                                          i)
                iopend = iopend + nops_current
                if (reptypepos <= iopend) exit
            end do

            newmove%neqtauops = nops_current
            newmove%eqtau_typeidx = [ (i, i = iopend - nops_current + 1, iopend) ]

        case default
            error stop "propose_wormreplacemove: no neqtauops rule for sector"
        end select

        ! bath part: create move after getting bath spec from diagram (= local) spec
        call bath_replacemove_from_local_repidx(newmove%bath, base,&
                                                     propspec%local_repidx,&
                                                     propspec%local_wormidx)

        if (newmove%neqtauops > 1) then
            ! replacing operator in equal-time product: local change needed
            call local_replacespec_for_eqtauops(remidx, newops, base, &
                                                    propspec%local_repidx, &
                                                    propspec%local_wormidx, &
                                                    newmove%neqtauops, &
                                                    newmove%eqtau_typeidx, &
                                                    newmove%offset_eqtau_newops, &
                                                    newmove%offset_swapop_newops)
            newmove%local_remidx = remidx
            call propose_update_impurity(newmove%local, remidx, newops)
        end if
    end

    subroutine local_replacespec_for_eqtauops(remidx, newops, from_diagram, repidx, &
                                                    wormidx, neqtauops, eqtau_typeidx, &
                                                    offset_eqtau_newops, &
                                                    offset_swapop_newops)
        integer, intent(out), allocatable          :: remidx(:)
        type(TLocalOper), intent(out), allocatable :: newops(:)
        type(TDiagram), intent(in)      :: from_diagram
        integer, intent(in)             :: repidx, wormidx, neqtauops
        integer, intent(in)             :: eqtau_typeidx(:)
        integer, intent(out)            :: offset_eqtau_newops
        integer, intent(out)            :: offset_swapop_newops
        integer :: i, j, into, ieqtime, tau_new_first, tau_new_last
        real(c_double) :: tau_new
        type(TLocalOper), allocatable :: ops(:)

        allocate(ops, source=get_operators(from_diagram%local))
        if (repidx < 1 .or. repidx > size(ops)) &
             error stop 'Invalid replace index'

        ! neqtauops worm operators are present in the old configuration
        ! at the imag time tau_old of the (worm) operator at index
        ! wormidx in the local operator array, of which the neqtauops -
        ! 1 worm operators not at index wormidx need to be moved to the
        ! imag time tau_new of the (non-worm) operator at index repidx.
        ! to safely construct the worm operator product at tau_new
        ! without accidentally performing other changes even in the
        ! unlikely case that there are additional operators at tau_new,
        ! we need to remove and reinsert all operators at tau_new in
        ! addition to the neqtauops - 1 operators at tau_old.

        ! we determine the indices of the first and last operator at
        ! tau_new (most likely both repidx)

        call operators_at_same_time(ops, repidx, tau_new_first, tau_new_last)
        tau_new = ops(repidx)%tau

        ! we write indices of the operators to be reinserted (all at
        ! tau_new and neqtauops - 1 at tau_old) into the array of
        ! indices to be removed and the operators themselves into the
        ! array of operators to be inserted
        allocate(remidx(neqtauops + tau_new_last - tau_new_first))
        allocate(newops(size(remidx)))
        into = 1
        do i = tau_new_first, tau_new_last
            if (i == repidx) then
                offset_eqtau_newops = into
                ! special treatment necessary: determine indices of
                ! neqtauops - 1 operators at tau_old to be moved for
                ! remidx and write all neqtauops operators at equal time
                ! tau_new in the right order into newops
                do j = 1, size(eqtau_typeidx)
                    ieqtime = get_worm_oper_position(from_diagram%worm, &
                                                       eqtau_typeidx(j))
                    if (ieqtime == wormidx) then
                        ! the operator at wormidx does not need to be
                        ! removed and reinserted, instead we do this for the
                        ! operator at repidx that will take its place in the
                        ! equal time product
                        offset_swapop_newops = into
                        remidx(into) = repidx
                        newops(into) = ops(repidx)
                    else
                        remidx(into) = ieqtime
                        newops(into) = ops(ieqtime)
                        newops(into)%tau = tau_new
                    end if
                    into = into + 1
                end do
            else
                ! no special treatment necessary
                remidx(into) = i
                newops(into) = ops(i)
                into = into + 1
            end if
        end do

        call sort(remidx)
    end subroutine local_replacespec_for_eqtauops

    !> Given operators, find ops(ifirst:ilast)%tau == ops(i)%tau
    subroutine operators_at_same_time(ops, i, ifirst, ilast)
        type(TLocalOper), intent(in) :: ops(:)
        integer, intent(in) :: i
        integer, intent(out) :: ifirst, ilast

        if (i < 1 .or. i > size(ops)) &
            error stop 'Invalid index'

        do ifirst = i, 2, -1
            if (ops(ifirst-1)%tau /= ops(i)%tau) &
                exit
        enddo
        do ilast = i, size(ops)-1
            if (ops(ilast+1)%tau /= ops(i)%tau) &
                exit
        enddo
    end subroutine

    subroutine bath_replacemove_from_local_repidx(bathprop, from_diagram,&
                                                        repidx, wormidx)
        type(TBathReplaceMove), intent(out) :: bathprop
        type(TDiagram), intent(in)      :: from_diagram
        integer, intent(in)             :: repidx, wormidx

        integer(kind(OpDummy))          :: bath_reptype
        integer                         :: bath_repidx
        type(TLocalOper)                :: wormop

        integer :: nworm, wtype
        type(TLocalOper), allocatable :: local_ops(:)

        ! XXX this seems to largely duplicate bath_remmove_from_local_remidx
        nworm = get_nworm(from_diagram%worm)
        allocate(local_ops, source=get_operators(from_diagram%local))

        if (wormidx < 1 .or. wormidx > size(local_ops)) &
             error stop 'Invalid worm index'
        if (repidx < 1 .or. repidx > size(local_ops)) &
             error stop 'Invalid rep index'

        bath_reptype = local_ops(repidx)%type
        if (bath_reptype == OpCrea) then
            wtype = OpCreaW
        elseif (bath_reptype == OpAnnh) then
            wtype = OpAnnhW
        else
            error stop 'Invalid replacement operator type'
        endif
        if (local_ops(wormidx)%type /= wtype) &
            error stop 'Worm operator has wrong type'
        bath_repidx = (from_diagram%natural_idx(repidx) - (nworm - 1)) / 2

        wormop = local_ops(wormidx)
        call init_bathreplacemove(bathprop, from_diagram%bath)
        call propose_bathreplacemove(bathprop, bath_repidx, bath_reptype,&
                TBathOperator(tau=wormop%tau, orb=wormop%orb, sp=wormop%sp))
    end subroutine bath_replacemove_from_local_repidx

    real(c_double) function weight_ratio_wormreplacemove(proposal) result(ratio)
        class(TWormReplaceMove), intent(inout) :: proposal
        type(TDiagram), pointer :: base
        ! real(c_double) :: localratio

        base => proposal%target

        if (proposal%neqtauops > 1) then
            ratio = weight_ratio(proposal%local)
            if (ratio == 0) return
        else
            ratio = 1.0d0
        end if

        ratio = ratio * weight_ratio(proposal%bath)
    end function weight_ratio_wormreplacemove

    subroutine accept_wormreplacemove(proposal)
        class(TWormReplaceMove), intent(inout) :: proposal
        type(TDiagram), pointer :: tgtstate
        integer :: i, delta_local_wormidx
        integer, allocatable :: newidx_addorder(:)

        tgtstate => proposal%target

        if (proposal%neqtauops > 1) then
            call get_newidx_addorder(proposal%local, newidx_addorder)
            ! new indices of the moved equal time product can be taken
            ! from get_newidx_addorder, indices of the other worm
            ! operators change exactly if they are between local_repidx
            ! and local_wormidx by (neqtauops - 1) (sign depends on order
            ! of local_repidx and local_wormidx)
            call update_wormstate_replace_eqtauop( &
                    tgtstate%worm, &
                    proposal%eqtau_typeidx, &
                    newidx_addorder(proposal%offset_eqtau_newops &
                        : proposal%offset_eqtau_newops + proposal%neqtauops - 1), &
                    proposal%spec%local_wormidx, &
                    proposal%spec%local_repidx&
            )

            ! accept the local move, making tgtstate%local contain the
            ! changed list
            call accept_update(proposal%local)

            ! local_repidx and local_wormidx are no longer the
            ! appropriate indices for the changed list, so we need to fix
            ! them before swapping the types

            ! fix local_wormidx by determining the number of operators
            ! inserted or removed before it
            delta_local_wormidx = 0
            do i = 1, size(proposal%local_remidx)
                if (proposal%local_remidx(i) == proposal%spec%local_wormidx) then
                    ! reachable only by bugs etc.
                    error stop "accept_wormreplacemove: removed local_wormidx"
                else if (proposal%local_remidx(i) < proposal%spec%local_wormidx&
                        .and. proposal%spec%local_wormidx < proposal%spec%local_repidx) then
                    ! if local_wormidx is less than local_repidx, then any
                    ! operators removed before it will reduce its index by
                    ! one (and their reinsertion after it is irrelevant)
                    delta_local_wormidx = delta_local_wormidx - 1
                else if (proposal%local_remidx(i) > proposal%spec%local_wormidx&
                        .and. proposal%spec%local_wormidx > proposal%spec%local_repidx) then
                    ! if local_wormidx is greater than local_repidx, then
                    ! any operators removed after it will be reinserted
                    ! before it, increasing its index by one (removed
                    ! operators before it will also be reinserted before it
                    ! causing no change)
                    delta_local_wormidx = delta_local_wormidx + 1
                end if
            end do
            proposal%spec%local_wormidx = &
                    proposal%spec%local_wormidx + delta_local_wormidx

            ! fix local_repidx by taking the new index from newidx_addorder
            proposal%spec%local_repidx = newidx_addorder(proposal%offset_swapop_newops)
        else
            call update_wormstate_replace(tgtstate%worm, &
                    proposal%spec%worm_reptypepos, proposal%spec%local_repidx)
        end if

        ! worm and local operator trade places
        call swap_operator_types( &
                     tgtstate%local, proposal%spec%local_repidx, &
                     proposal%spec%local_wormidx)

        ! accept bath replace move
        call accept_update(proposal%bath)
        call update_diagram_lookups(tgtstate)
        call update_diagram_weight(tgtstate)
        !call verify_worm(tgtstate%worm, get_operators(tgtstate%local))
    end subroutine accept_wormreplacemove

end module MDiagram
