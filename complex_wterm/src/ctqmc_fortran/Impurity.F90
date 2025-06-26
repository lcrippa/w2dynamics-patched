module MImpurity
    use iso_c_binding
    use MExtRange
    use MImpurityD
    use MImpurityZ
    use MOperator
    use MPsi
    implicit none
    private

    public :: TLocalOper, ZExtRange

    type, public :: TImpurity
        type(DImpurity), allocatable :: dself
        type(ZImpurity), allocatable :: zself
    end type

    type, abstract, public :: TImpurityUpdateBase
        class(DImpurityUpdateBase), allocatable :: dself
        class(ZImpurityUpdateBase), allocatable :: zself
    end type

    type, extends(TImpurityUpdateBase), public :: TImpurityReplace
    end type

    type, extends(TImpurityUpdateBase), public :: TImpurityGrow
    end type

    type, extends(TImpurityUpdateBase), public :: TImpurityShrink
    end type

    type, extends(TImpurityUpdateBase), public :: TImpurityShift
    end type

    type, extends(TImpurityUpdateBase), public :: TImpurityUpdate
    end type

    !> Initialize impurity with Hamiltonian and inverse temperature
    interface init_impurity
        module procedure p_init_impurity
    end interface
    public init_impurity

    !> Set impurity to given set of operators and outer state
    interface set_impurity
        module procedure p_set_impurity
    end interface
    public set_impurity

    !> Clear impurity of operators and set to outer state
    interface clear_impurity
        module procedure :: p_clear_impurity
    end interface
    public clear_impurity

    !> Return inverse temperature
    interface get_beta
        module procedure p_get_beta
    end interface
    public get_beta

    !> Return relative weight of each superstate in an empty impurity.
    interface thermal_weights
        module procedure p_thermal_weights
    end interface
    public thermal_weights

    !> Return number of operators in impurity cumulant
    interface size
        module procedure p_size
    end interface
    public size

    !> Return total number of orbitals in the impurity
    interface get_nbands
        module procedure p_get_nbands
    end interface
    public get_nbands

    !> Return set of operators of impurity or move
    interface get_operators
        module procedure p_get_operators, p_update_get_operators
    end interface
    public get_operators

    !> Return outer superstate currently chosen for the impurity
    interface get_outerbl
        module procedure p_get_outerbl
    end interface
    public get_outerbl

    !> Return superstate in the impurity just before the given time
    interface sst_at_time
        module procedure p_sst_at_time
    end interface
    public sst_at_time

    !> Check invariants in the impurity
    interface verify_impurity
        module procedure p_verify_impurity
    end interface
    public verify_impurity

    !> Return impurity weight (trace of the local cumulant)
    interface get_weight
        module procedure p_get_weight
    end interface
    public get_weight

    !> Compute many-body density matrix
    interface compute_density_matrix
        module procedure p_compute_density_matrix
    end interface
    public compute_density_matrix

    !> Write out information about the impurity
    interface dump
        module procedure p_dump
    end interface
    public dump

    !> Set of outer states for operators which do not violate quantum numbers
    interface allowed_outer_states
        module procedure p_allowed_outer_states, p_update_allowed_outer_states
    end interface
    public allowed_outer_states

    !> Check if the given outer state is allowed
    interface is_allowed_outer_state
       module procedure p_is_allowed_outer_state, p_update_is_allowed_outer_state
    end interface
    public is_allowed_outer_state

    !> Initialize update to replace all operators in the impurity
    interface init_replace_impurity
        module procedure p_init_replace_impurity
    end interface
    public init_replace_impurity

    !> Initialize update to add operators to the impurity
    interface init_grow_impurity
        module procedure p_init_grow_impurity
    end interface
    public init_grow_impurity

    !> Initialize update to remove operators from the impurity
    interface init_shrink_impurity
        module procedure p_init_shrink_impurity
    end interface
    public init_shrink_impurity

    !> Initialize update to cyclically shift the time of impurity operators
    interface init_shift_impurity
        module procedure p_init_shift_impurity
    end interface
    public init_shift_impurity

    !> Initialize update to replace impurity operators with others
    interface init_update_impurity
        module procedure p_init_update_impurity
    end interface
    public init_update_impurity

    !> Propose a completely new set of operators for the impurity
    interface propose_replace_impurity
        module procedure p_propose_replace_impurity
    end interface
    public propose_replace_impurity

    !> Propose operators for insertion into the impurity
    interface propose_grow_impurity
        module procedure p_propose_grow_impurity
    end interface
    public propose_grow_impurity

    !> Propose operators positions for removal from the impurity
    interface propose_shrink_impurity
        module procedure p_propose_shrink_impurity
    end interface
    public propose_shrink_impurity

    !> Propose cyclically shifting the time of all operators
    interface propose_shift_impurity
        module procedure p_propose_shift_impurity
    end interface
    public propose_shift_impurity

    !> Return number of operators shifted across zero.
    interface get_num_wrapped_ops
        module procedure p_get_num_wrapped_ops
    end interface
    public get_num_wrapped_ops

    !> Propose a permutation of spin-orbitals for the impurity operators
    interface propose_permute_impurity_flavours
        module procedure p_propose_permute_impurity_flavours
    end interface
    public propose_permute_impurity_flavours

    !> Propose a set of indices of operators to be replaced by new ones
    interface propose_update_impurity
        module procedure p_propose_update_impurity
    end interface
    public propose_update_impurity

    !> Choose outer state of the move
    interface choose_outer_state
        module procedure p_choose_outer_state
    end interface
    public choose_outer_state

    !> Choose matching outer state
    interface choose_matching_outer_state
        module procedure p_choose_matching_outer_state
    end interface
    public choose_matching_outer_state

    !> Return false if one is barred from calling accept_update() on the move
    interface is_update_allowed
        module procedure p_is_update_allowed
    end interface
    public is_update_allowed

    !> Compute weight of the updated configuration
    interface compute_weight
        module procedure p_compute_weight
    end interface
    public compute_weight

    !> Return ratio of new to old weight if the move is to be accepted
    interface weight_ratio
        module procedure p_weight_ratio
    end interface
    public weight_ratio

    !> Accept the update
    interface accept_update
        module procedure p_accept_update
    end interface
    public accept_update

    !> Swap type of two operators, allowing them to trade the worm state
    interface swap_operator_types
        module procedure p_swap_operator_types
    end interface
    public swap_operator_types

    interface get_newidx_addorder
        module procedure p_get_newidx_addorder
    end interface
    public get_newidx_addorder

contains

    logical function is_complex_impurity(self)
        type(TImpurity), intent(in) :: self

        if (allocated(self%zself) .eqv. allocated(self%dself)) &
            error stop 'Invalid state of impurity'
        is_complex_impurity = allocated(self%zself)
    end

    subroutine p_init_impurity(self, beta, hloc, psi, outer_sst)
        type(TImpurity), intent(out) :: self
        real(c_double), intent(in) :: beta
        type(TOperEigen), intent(in) :: hloc
        type(TPsis), intent(in) :: psi
        integer :: outer_sst

        if (is_complex_opeigen(hloc) .neqv. is_complex_psis(psi)) &
            error stop 'Inconsistent hloc/psi'

        if (is_complex_opeigen(hloc)) then
            allocate(self%zself)
            call init_impurity(self%zself, beta, hloc%zeig, psi%zpsis, outer_sst)
        else
            allocate(self%dself)
            call init_impurity(self%dself, beta, hloc%deig, psi%dpsis, outer_sst)
        endif
    end

    subroutine p_set_impurity(self, ops, outer_sst)
        type(TImpurity), intent(inout) :: self
        type(TLocalOper), intent(in) :: ops(:)
        integer, intent(in) :: outer_sst

        if (is_complex_impurity(self)) then
            call set_impurity(self%zself, ops, outer_sst)
        else
            call set_impurity(self%dself, ops, outer_sst)
        endif
    end

    subroutine p_clear_impurity(self, outer_sst)
        type(TImpurity), intent(inout) :: self
        integer, intent(in) :: outer_sst

        if (is_complex_impurity(self)) then
            call clear_impurity(self%zself, outer_sst)
        else
            call clear_impurity(self%dself, outer_sst)
        endif
    end

    function p_get_beta(self) result(r)
        type(TImpurity), intent(in) :: self
        real(c_double) :: r

        if (is_complex_impurity(self)) then
            r = get_beta(self%zself)
        else
            r = get_beta(self%dself)
        endif
    end

    function p_size(self) result(r)
        type(TImpurity), intent(in) :: self
        integer :: r

        if (is_complex_impurity(self)) then
            r = size(self%zself)
        else
            r = size(self%dself)
        endif
    end

    function p_get_nbands(self) result(r)
        type(TImpurity), intent(in) :: self
        integer :: r

        if (is_complex_impurity(self)) then
            r = get_nbands(self%zself)
        else
            r = get_nbands(self%dself)
        endif
    end

    function p_get_weight(self) result(r)
        type(TImpurity), intent(in) :: self
        type(ZExtRange) :: r

        if (is_complex_impurity(self)) then
            r = get_weight(self%zself)
        else
            r = get_weight(self%dself)
        endif
    end

    function p_get_operators(self) result(r)
        type(TImpurity), intent(in) :: self
        type(TLocalOper), allocatable :: r(:)

        if (is_complex_impurity(self)) then
            r = get_operators(self%zself)
        else
            r = get_operators(self%dself)
        endif
    end

    function p_get_outerbl(self) result(r)
        type(TImpurity), intent(in) :: self
        integer :: r

        if (is_complex_impurity(self)) then
            r = get_outerbl(self%zself)
        else
            r = get_outerbl(self%dself)
        endif
    end

    subroutine p_verify_impurity(self, rtol)
        type(TImpurity), intent(in) :: self
        real(c_double), intent(in) :: rtol

        if (is_complex_impurity(self)) then
            call verify_impurity(self%zself, rtol)
        else
            call verify_impurity(self%dself, rtol)
        endif
    end

    function p_thermal_weights(self) result(weights)
        type(TImpurity), intent(in) :: self
        real(c_double), allocatable :: weights(:)

        if (is_complex_impurity(self)) then
            weights = thermal_weights(self%zself)
        else
            weights = thermal_weights(self%dself)
        endif
    end

    function p_allowed_outer_states(self) result(r)
        type(TImpurity), intent(in) :: self
        integer, pointer :: r(:)

        if (is_complex_impurity(self)) then
            r => allowed_outer_states(self%zself)
        else
            r => allowed_outer_states(self%dself)
        endif
    end

    function p_is_allowed_outer_state(self, sst) result(r)
        type(TImpurity), intent(in) :: self
        integer, intent(in) :: sst
        logical :: r

        if (is_complex_impurity(self)) then
            r = is_allowed_outer_state(self%zself, sst)
        else
            r = is_allowed_outer_state(self%dself, sst)
        endif
    end

    subroutine p_compute_density_matrix(self, sgn, out_op)
        type(TImpurity), intent(in) :: self
        complex(c_double_complex), intent(in) :: sgn
        type(TOperator), intent(inout) :: out_op

        if (is_complex_impurity(self)) then
            if (allocated(out_op%dop)) &
                error stop 'Operator may not be real'
            if (.not. allocated(out_op%zop)) &
                allocate(out_op%zop)
            call compute_density_matrix(self%zself, sgn, out_op%zop)
        else
            if (allocated(out_op%zop)) &
                error stop 'Operator may not be complex'
            if (.not. allocated(out_op%dop)) &
                allocate(out_op%dop)
            if (aimag(sgn) /= 0) &
                error stop 'sgn must be real'
            call compute_density_matrix(self%dself, real(sgn), out_op%dop)
        endif
    end subroutine

    function p_sst_at_time(self, tau) result(r)
        type(TImpurity), intent(in) :: self
        real(c_double), intent(in) :: tau
        integer :: r

        if (is_complex_impurity(self)) then
            r = sst_at_time(self%zself, tau)
        else
            r = sst_at_time(self%dself, tau)
        endif
    end

    subroutine p_swap_operator_types(self, idx1, idx2)
        type(TImpurity), intent(inout) :: self
        integer, intent(in) :: idx1, idx2

        if (is_complex_impurity(self)) then
            call swap_operator_types(self%zself, idx1, idx2)
        else
            call swap_operator_types(self%dself, idx1, idx2)
        endif
    end subroutine

    subroutine p_dump(self, unit)
        type(TImpurity), intent(in) :: self
        integer, intent(in), optional :: unit

        if (is_complex_impurity(self)) then
            call dump(self%zself, unit)
        else
            call dump(self%dself, unit)
        endif
    end subroutine

    ! -------------------------------------------------------------------------
    ! BASE MOVE

    logical function is_complex_impurity_update(self)
        class(TImpurityUpdateBase), intent(in) :: self

        if (allocated(self%zself) .eqv. allocated(self%dself)) &
            error stop 'Invalid state of impurity update'
        is_complex_impurity_update = allocated(self%zself)
    end

    function p_is_update_allowed(self) result(r)
        class(TImpurityUpdateBase), intent(in) :: self
        logical :: r

        if (is_complex_impurity_update(self)) then
            r = is_update_allowed(self%zself)
        else
            r = is_update_allowed(self%dself)
        endif
    end

    subroutine p_compute_weight(self)
        class(TImpurityUpdateBase), intent(in), target :: self

        ! the argument of a select type construct may not be a component, so
        ! we have to wrap this thing into a pointer
        class(ZImpurityUpdateBase), pointer :: zself
        class(DImpurityUpdateBase), pointer :: dself

        if (is_complex_impurity_update(self)) then
            zself => self%zself
            select type (zself)
            type is (ZImpurityReplace)
                call compute_weight(zself)
            type is (ZImpurityGrow)
                call compute_weight(zself)
            type is (ZImpurityShrink)
                call compute_weight(zself)
            type is (ZImpurityShift)
                call compute_weight(zself)
            type is (ZImpurityUpdate)
                call compute_weight(zself)
            end select
        else
            dself => self%dself
            select type (dself)
            type is (DImpurityReplace)
                call compute_weight(dself)
            type is (DImpurityGrow)
                call compute_weight(dself)
            type is (DImpurityShrink)
                call compute_weight(dself)
            type is (DImpurityShift)
                call compute_weight(dself)
            type is (DImpurityUpdate)
                call compute_weight(dself)
            end select
        endif
    end subroutine

    function p_weight_ratio(self) result(r)
        class(TImpurityUpdateBase), intent(in) :: self
        real(c_double) :: r

        if (is_complex_impurity_update(self)) then
            r = weight_ratio(self%zself)
        else
            r = weight_ratio(self%dself)
        endif
    end

    subroutine p_accept_update(self)
        class(TImpurityUpdateBase), intent(inout) :: self

        if (is_complex_impurity_update(self)) then
            call accept_update(self%zself)
        else
            call accept_update(self%dself)
        endif
    end

    subroutine p_get_newidx_addorder(self, out)
        class(TImpurityUpdateBase), intent(in), target :: self
        integer, allocatable, intent(out) :: out(:)

        ! the argument of a select type construct may not be a component, so
        ! we have to wrap this thing into a pointer
        class(ZImpurityUpdateBase), pointer :: zself
        class(DImpurityUpdateBase), pointer :: dself

        if (is_complex_impurity_update(self)) then
            zself => self%zself
            select type (zself)
            type is (ZImpurityGrow)
                call get_newidx_addorder(zself, out)
            type is (ZImpurityUpdate)
                call get_newidx_addorder(zself, out)
            end select
        else
            dself => self%dself
            select type (dself)
            type is (DImpurityGrow)
                call get_newidx_addorder(dself, out)
            type is (DImpurityUpdate)
                call get_newidx_addorder(dself, out)
            end select
        endif
    end subroutine

    function p_update_get_operators(self) result(r)
        class(TImpurityUpdateBase), intent(in), target :: self
        type(TLocalOper), allocatable :: r(:)

        ! the argument of a select type construct may not be a component, so
        ! we have to wrap this thing into a pointer
        class(ZImpurityUpdateBase), pointer :: zself
        class(DImpurityUpdateBase), pointer :: dself

        if (is_complex_impurity_update(self)) then
            zself => self%zself
            select type (zself)
            type is (ZImpurityReplace)
                r = get_operators(zself)
            type is (ZImpurityGrow)
                r = get_operators(zself)
            type is (ZImpurityShrink)
                r = get_operators(zself)
            type is (ZImpurityShift)
                r = get_operators(zself)
            type is (ZImpurityUpdate)
                r = get_operators(zself)
            class default
                error stop 'Invalid type'
            end select
        else
            dself => self%dself
            select type (dself)
            type is (DImpurityReplace)
                r = get_operators(dself)
            type is (DImpurityGrow)
                r = get_operators(dself)
            type is (DImpurityShrink)
                r = get_operators(dself)
            type is (DImpurityShift)
                r = get_operators(dself)
            type is (DImpurityUpdate)
                r = get_operators(dself)
            class default
                error stop 'Invalid type'
            end select
        endif
    end function

    subroutine p_choose_matching_outer_state(self, outer)
        class(TImpurityUpdateBase), intent(inout), target :: self
        logical, intent(in) :: outer

        ! the argument of a select type construct may not be a component, so
        ! we have to wrap this thing into a pointer
        class(ZImpurityUpdateBase), pointer :: zself
        class(DImpurityUpdateBase), pointer :: dself

        if (is_complex_impurity_update(self)) then
            zself => self%zself
            select type (zself)
            type is (ZImpurityGrow)
                call choose_matching_outer_state(zself, outer)
            type is (ZImpurityShrink)
                call choose_matching_outer_state(zself, outer)
            type is (ZImpurityUpdate)
                call choose_matching_outer_state(zself, outer)
            end select
        else
            dself => self%dself
            select type (dself)
            type is (DImpurityGrow)
                call choose_matching_outer_state(dself, outer)
            type is (DImpurityShrink)
                call choose_matching_outer_state(dself, outer)
            type is (DImpurityUpdate)
                call choose_matching_outer_state(dself, outer)
            end select
        endif
    end subroutine

    ! -------------------------------------------------------------------------
    ! REPLACE MOVE

    subroutine p_init_replace_impurity(self, t)
        type(TImpurityReplace), intent(out) :: self
        type(TImpurity), intent(inout), target :: t

        type(ZImpurityReplace), allocatable :: zself
        type(DImpurityReplace), allocatable :: dself

        if (is_complex_impurity(t)) then
            allocate(zself)
            call init_replace_impurity(zself, t%zself)
            call move_alloc(zself, self%zself)
        else
            allocate(dself)
            call init_replace_impurity(dself, t%dself)
            call move_alloc(dself, self%dself)
        endif
    end subroutine

    subroutine p_propose_replace_impurity(self, op, outer_sst)
        type(TImpurityReplace), intent(inout), target :: self
        type(TLocalOper), intent(in) :: op(:)
        integer, intent(in), optional :: outer_sst

        ! the argument of a select type construct may not be a component, so
        ! we have to wrap this thing into a pointer
        class(ZImpurityUpdateBase), pointer :: zself
        class(DImpurityUpdateBase), pointer :: dself

        if (is_complex_impurity_update(self)) then
            zself => self%zself
            select type (zself)
            type is (ZImpurityReplace)
                call propose_replace_impurity(zself, op, outer_sst)
            end select
        else
            dself => self%dself
            select type (dself)
            type is (DImpurityReplace)
                call propose_replace_impurity(dself, op, outer_sst)
            end select
        endif
    end subroutine

    function p_update_allowed_outer_states(self) result(r)
        class(TImpurityReplace), intent(in), target :: self
        integer, pointer :: r(:)

        class(ZImpurityUpdateBase), pointer :: zself
        class(DImpurityUpdateBase), pointer :: dself

        if (is_complex_impurity_update(self)) then
            zself => self%zself
            select type (zself)
            type is (ZImpurityReplace)
                r => allowed_outer_states(zself)
            end select
        else
            dself => self%dself
            select type (dself)
            type is (DImpurityReplace)
                r => allowed_outer_states(dself)
            end select
        endif
    end

    function p_update_is_allowed_outer_state(self, sst) result(r)
        class(TImpurityReplace), intent(in), target :: self
        integer, intent(in) :: sst
        logical :: r

        class(ZImpurityUpdateBase), pointer :: zself
        class(DImpurityUpdateBase), pointer :: dself

        r = .false.
        if (is_complex_impurity_update(self)) then
            zself => self%zself
            select type (zself)
            type is (ZImpurityReplace)
                r = is_allowed_outer_state(zself, sst)
            end select
        else
            dself => self%dself
            select type (dself)
            type is (DImpurityReplace)
                r = is_allowed_outer_state(dself, sst)
            end select
        endif
    end

    subroutine p_choose_outer_state(self, outer_sst)
        class(TImpurityReplace), intent(inout), target :: self
        integer, intent(in) :: outer_sst

        ! the argument of a select type construct may not be a component, so
        ! we have to wrap this thing into a pointer
        class(ZImpurityUpdateBase), pointer :: zself
        class(DImpurityUpdateBase), pointer :: dself

        if (is_complex_impurity_update(self)) then
            zself => self%zself
            select type (zself)
            type is (ZImpurityReplace)
                call choose_outer_state(zself, outer_sst)
            end select
        else
            dself => self%dself
            select type (dself)
            type is (DImpurityReplace)
                call choose_outer_state(dself, outer_sst)
            end select
        endif
    end

    ! -------------------------------------------------------------------------
    ! GROW MOVE

    subroutine p_init_grow_impurity(self, t)
        type(TImpurityGrow), intent(out) :: self
        type(TImpurity), intent(inout), target :: t

        type(ZImpurityGrow), allocatable :: zself
        type(DImpurityGrow), allocatable :: dself

        if (is_complex_impurity(t)) then
            allocate(zself)
            call init_grow_impurity(zself, t%zself)
            call move_alloc(zself, self%zself)
        else
            allocate(dself)
            call init_grow_impurity(dself, t%dself)
            call move_alloc(dself, self%dself)
        endif
    end subroutine

    subroutine p_propose_grow_impurity(self, op, outer_sst)
        type(TImpurityGrow), intent(inout), target :: self
        type(TLocalOper), intent(in) :: op(:)
        integer, intent(in), optional :: outer_sst

        ! the argument of a select type construct may not be a component, so
        ! we have to wrap this thing into a pointer
        class(ZImpurityUpdateBase), pointer :: zself
        class(DImpurityUpdateBase), pointer :: dself

        if (is_complex_impurity_update(self)) then
            zself => self%zself
            select type (zself)
            type is (ZImpurityGrow)
                call propose_grow_impurity(zself, op, outer_sst)
            end select
        else
            dself => self%dself
            select type (dself)
            type is (DImpurityGrow)
                call propose_grow_impurity(dself, op, outer_sst)
            end select
        endif
    end subroutine

    ! -------------------------------------------------------------------------
    ! SHRINK MOVE

    subroutine p_init_shrink_impurity(self, t)
        type(TImpurityShrink), intent(out) :: self
        type(TImpurity), intent(inout), target :: t

        type(ZImpurityShrink), allocatable :: zself
        type(DImpurityShrink), allocatable :: dself

        if (is_complex_impurity(t)) then
            allocate(zself)
            call init_shrink_impurity(zself, t%zself)
            call move_alloc(zself, self%zself)
        else
            allocate(dself)
            call init_shrink_impurity(dself, t%dself)
            call move_alloc(dself, self%dself)
        endif
    end subroutine

    subroutine p_propose_shrink_impurity(self, idx, outer_sst)
        type(TImpurityShrink), intent(inout), target :: self
        integer, intent(in) :: idx(:)
        integer, intent(in), optional :: outer_sst

        ! the argument of a select type construct may not be a component, so
        ! we have to wrap this thing into a pointer
        class(ZImpurityUpdateBase), pointer :: zself
        class(DImpurityUpdateBase), pointer :: dself

        if (is_complex_impurity_update(self)) then
            zself => self%zself
            select type (zself)
            type is (ZImpurityShrink)
                call propose_shrink_impurity(zself, idx, outer_sst)
            end select
        else
            dself => self%dself
            select type (dself)
            type is (DImpurityShrink)
                call propose_shrink_impurity(dself, idx, outer_sst)
            end select
        endif
    end subroutine

    ! -------------------------------------------------------------------------
    ! SHIFT MOVE

    subroutine p_init_shift_impurity(self, t)
        type(TImpurityShift), intent(out) :: self
        type(TImpurity), intent(inout), target :: t

        type(ZImpurityShift), allocatable :: zself
        type(DImpurityShift), allocatable :: dself

        if (is_complex_impurity(t)) then
            allocate(zself)
            call init_shift_impurity(zself, t%zself)
            call move_alloc(zself, self%zself)
        else
            allocate(dself)
            call init_shift_impurity(dself, t%dself)
            call move_alloc(dself, self%dself)
        endif
    end subroutine

    subroutine p_propose_shift_impurity(self, delta_tau)
        type(TImpurityShift), intent(inout), target :: self
        real(c_double), intent(in) :: delta_tau

        ! the argument of a select type construct may not be a component, so
        ! we have to wrap this thing into a pointer
        class(ZImpurityUpdateBase), pointer :: zself
        class(DImpurityUpdateBase), pointer :: dself

        if (is_complex_impurity_update(self)) then
            zself => self%zself
            select type (zself)
            type is (ZImpurityShift)
                call propose_shift_impurity(zself, delta_tau)
            end select
        else
            dself => self%dself
            select type (dself)
            type is (DImpurityShift)
                call propose_shift_impurity(dself, delta_tau)
            end select
        endif
    end subroutine

    function p_get_num_wrapped_ops(self) result(r)
        type(TImpurityShift), intent(in), target :: self
        integer :: r

        ! the argument of a select type construct may not be a component, so
        ! we have to wrap this thing into a pointer
        class(ZImpurityUpdateBase), pointer :: zself
        class(DImpurityUpdateBase), pointer :: dself

        if (is_complex_impurity_update(self)) then
            zself => self%zself
            select type (zself)
            type is (ZImpurityShift)
                r = get_num_wrapped_ops(zself)
            class default
                error stop 'Invalid type'
            end select
        else
            dself => self%dself
            select type (dself)
            type is (DImpurityShift)
                r = get_num_wrapped_ops(dself)
            class default
                error stop 'Invalid type'
            end select
        endif
    end function

    ! -------------------------------------------------------------------------
    ! PERMUTE MOVE

    subroutine p_propose_permute_impurity_flavours( &
                        self, ca_perm, flavour_perm, outer_sst)
        type(TImpurityReplace), intent(inout), target :: self
        logical, intent(in) :: ca_perm
        integer, intent(in) :: flavour_perm(:)
        integer, intent(in), optional :: outer_sst

        ! the argument of a select type construct may not be a component, so
        ! we have to wrap this thing into a pointer
        class(ZImpurityUpdateBase), pointer :: zself
        class(DImpurityUpdateBase), pointer :: dself

        if (is_complex_impurity_update(self)) then
            zself => self%zself
            select type (zself)
            type is (ZImpurityReplace)
                call propose_permute_impurity_flavours( &
                            zself, ca_perm, flavour_perm, outer_sst)
            end select
        else
            dself => self%dself
            select type (dself)
            type is (DImpurityReplace)
                call propose_permute_impurity_flavours( &
                            dself, ca_perm, flavour_perm, outer_sst)
            end select
        endif
    end subroutine

    ! -------------------------------------------------------------------------
    ! UPDATE MOVE

    subroutine p_init_update_impurity(self, t)
        type(TImpurityUpdate), intent(out) :: self
        type(TImpurity), intent(inout), target :: t

        type(ZImpurityUpdate), allocatable :: zself
        type(DImpurityUpdate), allocatable :: dself

        if (is_complex_impurity(t)) then
            allocate(zself)
            call init_update_impurity(zself, t%zself)
            call move_alloc(zself, self%zself)
        else
            allocate(dself)
            call init_update_impurity(dself, t%dself)
            call move_alloc(dself, self%dself)
        endif
    end subroutine

    subroutine p_propose_update_impurity(self, idx, op, outer_sst)
        type(TImpurityUpdate), intent(inout), target :: self
        integer, intent(in) :: idx(:)
        type(TLocalOper), intent(in) :: op(:)
        integer, intent(in), optional :: outer_sst

        ! the argument of a select type construct may not be a component, so
        ! we have to wrap this thing into a pointer
        class(ZImpurityUpdateBase), pointer :: zself
        class(DImpurityUpdateBase), pointer :: dself

        if (is_complex_impurity_update(self)) then
            zself => self%zself
            select type (zself)
            type is (ZImpurityUpdate)
                call propose_update_impurity(zself, idx, op, outer_sst)
            end select
        else
            dself => self%dself
            select type (dself)
            type is (DImpurityUpdate)
                call propose_update_impurity(dself, idx, op, outer_sst)
            end select
        endif
    end subroutine

end module
