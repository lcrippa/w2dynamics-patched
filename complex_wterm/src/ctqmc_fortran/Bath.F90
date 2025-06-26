module MBath
    use iso_c_binding, only: c_double, c_double_complex
    use MCommon, only: OpDummy
    use MBathBase
    use DMBath
    use ZMBath
    use MBlockMatrixD
    use MBlockMatrixZ
    use MHybridizationD
    use MHybridizationZ
    use MExtRange

    type, public :: TBath
        private
        type(DTBath), allocatable :: dbath
        type(ZTBath), allocatable :: zbath
    end type

    type, public :: TBathMove
        private
        type(DTBathMove), allocatable :: dmove
        type(ZTBathMove), allocatable :: zmove
    end type

    type, public :: TBathAddMove
        private
        type(DTBathAddMove), allocatable :: dmove
        type(ZTBathAddMove), allocatable :: zmove
    end type

    type, public :: TBathRemMove
        private
        type(DTBathRemMove), allocatable :: dmove
        type(ZTBathRemMove), allocatable :: zmove
    end type

    type, public :: TBathReplaceMove
        private
        type(DTBathReplaceMove), allocatable :: dmove
        type(ZTBathReplaceMove), allocatable :: zmove
    end type

    type, public :: TBathShiftMove
        private
        type(DTBathShiftMove), allocatable :: dmove
        type(ZTBathShiftMove), allocatable :: zmove
    end type

    interface copyinv
        module procedure dcopyinv, zcopyinv
    end interface
    public copyinv

    interface get_annihilators
        module procedure p_annihilators
    end interface
    public get_annihilators

    interface get_creators
        module procedure p_creators
    end interface
    public get_creators

    interface get_beta
        module procedure p_get_beta
    end interface
    public get_beta

    interface size
        module procedure p_size
    end interface
    public size

    interface get_det
        module procedure p_getdet
    end interface
    public get_det

    interface get_nbands
        module procedure p_getnbands
    end interface
    public get_nbands

    interface get_nosoper
        module procedure p_get_nosoper
    end interface
    public get_nosoper

    interface verify_bath
        module procedure p_verify_bath
    end Interface
    public verify_bath

    interface init_bath
        module procedure p_init_bath_real, p_init_bath_complex
    end interface
    public init_bath

    interface clear_bath
        module procedure p_clear_bath
    end interface
    public clear_bath

    interface debug_print_bathconfig_highprec
        module procedure p_debug_print_bathconfig_highprec
    end interface
    public debug_print_bathconfig_highprec

    interface propose_bath_move
        module procedure p_propose_bath_move
    end interface
    public propose_bath_move

    interface propose_bath_add_move
        module procedure p_propose_bath_add_move
    end interface
    public propose_bath_add_move

    interface propose_bath_rem_move
        module procedure p_propose_bath_rem_move
    end interface
    public propose_bath_rem_move

    interface propose_bathreplacemove
        module procedure p_propose_bathreplacemove
    end interface
    public propose_bathreplacemove

    interface propose_bath_shift_move
        module procedure p_propose_bath_shift_move
    end interface
    public propose_bath_shift_move

    interface init_bath_move
        module procedure p_init_bath_move
    end interface
    public init_bath_move

    interface init_bath_add_move
        module procedure p_init_bath_add_move
    end interface
    public init_bath_add_move

    interface init_bath_rem_move
        module procedure p_init_bath_rem_move
    end interface
    public init_bath_rem_move

    interface init_bathreplacemove
        module procedure p_init_bathreplacemove
    end interface
    public init_bathreplacemove

    interface init_bath_shift_move
        module procedure p_init_bath_shift_move
    end interface
    public init_bath_shift_move

    interface is_update_allowed
        module procedure move_is_update_allowed
        module procedure add_is_update_allowed
        module procedure rem_is_update_allowed
        module procedure replace_is_update_allowed
        module procedure shift_is_update_allowed
    end interface
    public is_update_allowed

    interface compute_weight
        module procedure move_compute_weight
        module procedure add_compute_weight
        module procedure rem_compute_weight
        module procedure repl_compute_weight
        module procedure shift_compute_weight
    end interface
    public compute_weight

    interface weight_ratio
        module procedure p_bath_absratio
        module procedure p_bath_absratio_add
        module procedure p_bath_absratio_rem
        module procedure p_bath_absratio_replace
        module procedure p_bath_absratio_shift
    end interface
    public weight_ratio

    interface accept_update
        module procedure p_bath_accept
        module procedure p_bath_accept_add
        module procedure p_bath_accept_rem
        module procedure p_bath_accept_replace
        module procedure p_bath_accept_shift
    end interface
    public accept_update

    interface get_sortval
        module procedure p_getsortval
    end interface
    public get_sortval

contains

    ! -------------------------------------------------------------------------
    ! BATH FUNCTIONS

    subroutine dcopyinv(this, inv)
        type(TBath), intent(in) :: this
        real(c_double), allocatable, intent(out) :: inv(:, :)

        if (allocated(this%zbath)) error stop 'complex matrix'
        inv = getinv(this%dbath)
    end

    subroutine zcopyinv(this, inv)
        type(TBath), intent(in) :: this
        complex(c_double_complex), allocatable, intent(out) :: inv(:, :)

        if (allocated(this%zbath)) then
            inv = getinv(this%zbath)
        else
            inv = getinv(this%dbath)
        end if
    end

    pure function p_annihilators(self) result(ops)
        type(TBath), intent(in) :: self
        type(TBathOperator), allocatable :: ops(:)

        if (allocated(self%zbath)) then
            ops = get_annihilators(self%zbath)
        else
            ops = get_annihilators(self%dbath)
        endif
    end

    pure elemental function p_getsortval(self, op) result(val)
        class(TBath), intent(in) :: self
        type(TBathOperator), intent(in) :: op
        real(c_double) :: val

        if (allocated(self%zbath)) then
            val = get_sortval(self%zbath, op)
        else
            val = get_sortval(self%dbath, op)
        endif
    end

    pure function p_creators(self) result(ops)
        type(TBath), intent(in) :: self
        type(TBathOperator), allocatable :: ops(:)

        if (allocated(self%zbath)) then
            ops = get_creators(self%zbath)
        else
            ops = get_creators(self%dbath)
        endif
    end

    pure function p_get_beta(self) result(beta)
        type(TBath), intent(in) :: self
        real(c_double) :: beta

        if (allocated(self%zbath)) then
            beta = get_beta(self%zbath)
        else
            beta = get_beta(self%dbath)
        endif
    end

    pure integer function p_size(self) result(r)
        type(TBath), intent(in) :: self

        if (allocated(self%zbath)) then
            r = size(self%zbath)
        else
            r = size(self%dbath)
        endif
    end

    pure function p_getdet(self) result(det)
        type(TBath), intent(in) :: self
        type(ZExtRange) :: det

        if (allocated(self%zbath)) then
            det = get_det(self%zbath)
        else
            det = get_det(self%dbath)
        endif
    end

    pure integer function p_getnbands(self) result(nbands)
        type(TBath), intent(in) :: self

        if (allocated(self%zbath)) then
            nbands = get_nbands(self%zbath)
        else
            nbands = get_nbands(self%dbath)
        endif
    end

    pure integer function p_get_nosoper(self, orb, sp) result(nosoper)
        type(TBath), intent(in) :: self
        integer, intent(in)     :: orb, sp

        if (allocated(self%zbath)) then
            nosoper = get_nosoper(self%zbath, orb, sp)
        else
            nosoper = get_nosoper(self%dbath, orb, sp)
        endif
    end function

    subroutine p_verify_bath(self, rtol)
        type(TBath), intent(in) :: self
        real(c_double), intent(in) :: rtol

        if (allocated(self%zbath)) then
            call verify_bath(self%zbath, rtol)
        else
            call verify_bath(self%dbath, rtol)
        endif
    end subroutine

    subroutine p_debug_print_bathconfig_highprec(self)
        type(TBath), intent(in) :: self

        if (allocated(self%zbath)) then
            call debug_print_bathconfig_highprec(self%zbath)
        else
            call debug_print_bathconfig_highprec(self%dbath)
        endif
    end subroutine

    subroutine p_init_bath_real(self, hybr)
        type(TBath), intent(out) :: self
        class(DHybridization), intent(in) :: hybr

        allocate(self%dbath)
        call init_bath(self%dbath, hybr)
    end subroutine

    subroutine p_init_bath_complex(self, hybr)
        type(TBath), intent(out) :: self
        class(ZHybridization), intent(in) :: hybr

        allocate(self%zbath)
        call init_bath(self%zbath, hybr)
    end subroutine

    subroutine p_clear_bath(self)
        type(TBath), intent(inout) :: self

        if (allocated(self%zbath)) then
            call clear_bath(self%zbath)
        else
            call clear_bath(self%dbath)
        endif
    end subroutine

    ! -------------------------------------------------------------------------
    ! BATH MOVE

    subroutine p_propose_bath_move(move, crea, annh)
        type(TBathMove), intent(inout) :: move
        type(TBathOperator), intent(in) :: crea(:), annh(:)

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            call propose_bath_move(move%zmove, crea, annh)
        else
            call propose_bath_move(move%dmove, crea, annh)
        endif
    end subroutine

    logical function move_is_update_allowed(move) result(t)
        type(TBathMove), intent(in) :: move

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            t = is_update_allowed(move%zmove)
        else
            t = is_update_allowed(move%dmove)
        endif
    end function

    real(c_double) function p_bath_absratio(move) result(ratio)
        type(TBathMove), intent(inout) :: move

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            ratio = weight_ratio(move%zmove)
        else
            ratio = weight_ratio(move%dmove)
        endif
    end function

    subroutine move_compute_weight(move)
        type(TBathMove), intent(inout) :: move

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            call compute_weight(move%zmove)
        else
            call compute_weight(move%dmove)
        endif
    end subroutine

    !> Copies stuff into MBath (=parent), which saves the recently accepted configuration
    subroutine p_bath_accept(move)
        type(TBathMove), intent(inout) :: move

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            call accept_update(move%zmove)
        else
            call accept_update(move%dmove)
        endif
        !call verify_bath(move%parent, 5d-2)
    end subroutine

    subroutine p_init_bath_move(this, parent)
        type(TBathMove), intent(out) :: this
        type(TBath), intent(in), target   :: parent

        if(allocated(parent%zbath))then
            call init_bath_move(this%zmove, parent%zbath)
        else
            call init_bath_move(this%dmove, parent%dbath)
        endif
    end subroutine

    ! -------------------------------------------------------------------------
    ! BATH ADD MOVE

    subroutine p_propose_bath_add_move(move, crea, annh)
        type(TBathAddMove) :: move
        type(TBathOperator), intent(in) :: crea(:), annh(:)

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
           call propose_bath_add_move(move%zmove, crea, annh)
        else
           call propose_bath_add_move(move%dmove, crea, annh)
       endif
     end subroutine

    logical function add_is_update_allowed(move) result(t)
        type(TBathAddMove), intent(in) :: move

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            t = is_update_allowed(move%zmove)
        else
            t = is_update_allowed(move%dmove)
        endif
    end function

    subroutine add_compute_weight(move)
        type(TBathAddMove), intent(inout) :: move

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            call compute_weight(move%zmove)
        else
            call compute_weight(move%dmove)
        endif
    end subroutine

    real(c_double) function p_bath_absratio_add(move) result(ratio)
        type(TBathAddMove), intent(inout) :: move

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            ratio = weight_ratio(move%zmove)
        else
            ratio = weight_ratio(move%dmove)
        endif
    end function

    !> Copies stuff into MBath (=parent), which saves the recently accepted
    !! configuration
    subroutine p_bath_accept_add(move)
        type(TBathAddMove), intent(inout) :: move

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            call accept_update(move%zmove)
        else
            call accept_update(move%dmove)
        endif
        !call verify_bath(move%parent, 5d-2)
    end subroutine

    subroutine p_init_bath_add_move(this, parent)
        type(TBathAddMove), intent(out) :: this
        type(TBath), intent(in), target   :: parent

        if(allocated(parent%zbath))then
            call init_bath_add_move(this%zmove, parent%zbath)
        else
            call init_bath_add_move(this%dmove, parent%dbath)
        endif
    end subroutine

    ! -------------------------------------------------------------------------
    ! BATH REMOVE MOVE

    subroutine p_propose_bath_rem_move(move, remidx)
        type(TBathRemMove), intent(inout)     :: move
        integer, intent(in) :: remidx(:)

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            call propose_bath_rem_move(move%zmove, remidx)
        else
            call propose_bath_rem_move(move%dmove, remidx)
        endif
    end subroutine

    logical function rem_is_update_allowed(move) result(t)
        type(TBathRemMove), intent(in) :: move

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            t = is_update_allowed(move%zmove)
        else
            t = is_update_allowed(move%dmove)
        endif
    end function

    subroutine rem_compute_weight(move)
        type(TBathRemMove), intent(inout) :: move

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            call compute_weight(move%zmove)
        else
            call compute_weight(move%dmove)
        endif
    end subroutine

    real(c_double) function p_bath_absratio_rem(move) result(ratio)
        type(TBathRemMove), intent(inout) :: move

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            ratio = weight_ratio(move%zmove)
        else
            ratio = weight_ratio(move%dmove)
        endif
    end function

    !> Copies stuff into MBath (=parent), which saves the recently accepted
    !! configuration
    subroutine p_bath_accept_rem(move)
        type(TBathRemMove), intent(inout) :: move

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            call accept_update(move%zmove)
        else
            call accept_update(move%dmove)
                endif
        !call verify_bath(move%parent, 5d-2)
    end subroutine

    subroutine p_init_bath_rem_move(this, parent)
        type(TBathRemMove), intent(out) :: this
        type(TBath), intent(in), target   :: parent

        if(allocated(parent%zbath))then
            call init_bath_rem_move(this%zmove, parent%zbath)
        else
            call init_bath_rem_move(this%dmove, parent%dbath)
        endif
    end subroutine

    ! -------------------------------------------------------------------------
    ! BATH REPLACE MOVE

    subroutine p_init_bathreplacemove(this, parent)
        type(TBathReplaceMove), intent(out) :: this
        type(TBath), intent(in), target   :: parent

        if(allocated(parent%zbath))then
            call init_bathreplacemove(this%zmove, parent%zbath)
        else
            call init_bathreplacemove(this%dmove, parent%dbath)
        endif
    end subroutine

    subroutine p_propose_bathreplacemove(move, origidx, origtype, newop)
        type(TBathReplaceMove), intent(inout)   :: move
        integer, intent(in)                   :: origidx
        integer(kind(OpDummy)), intent(in)     :: origtype
        type(TBathOperator), intent(in)       :: newop

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            call propose_bathreplacemove(move%zmove, origidx, origtype, newop)
        else
            call propose_bathreplacemove(move%dmove, origidx, origtype, newop)
        endif
    end subroutine

    logical function replace_is_update_allowed(move) result(t)
        type(TBathReplaceMove), intent(in) :: move

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            t = is_update_allowed(move%zmove)
        else
            t = is_update_allowed(move%dmove)
        endif
    end function

    subroutine repl_compute_weight(move)
        type(TBathReplaceMove), intent(inout) :: move

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            call compute_weight(move%zmove)
        else
            call compute_weight(move%dmove)
        endif
    end subroutine

    real(c_double) function p_bath_absratio_replace(move) result(ratio)
        type(TBathReplaceMove), intent(inout)  :: move

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            ratio = weight_ratio(move%zmove)
        else
            ratio = weight_ratio(move%dmove)
        endif
    end function

    subroutine p_bath_accept_replace(move)
        type(TBathReplaceMove), intent(inout)    :: move

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            call accept_update(move%zmove)
        else
            call accept_update(move%dmove)
        endif
    end subroutine

    ! -------------------------------------------------------------------------
    ! BATH SHIFT MOVE

    subroutine p_init_bath_shift_move(this, parent)
        type(TBathShiftMove), intent(out) :: this
        type(TBath), intent(in), target   :: parent

        if(allocated(parent%zbath))then
            call init_bath_shift_move(this%zmove, parent%zbath)
        else
            call init_bath_shift_move(this%dmove, parent%dbath)
        endif
    end subroutine

    subroutine p_propose_bath_shift_move(move, delta_tau)
        type(TBathShiftMove), intent(inout)   :: move
        real(c_double), intent(in)               :: delta_tau

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            call propose_bath_shift_move(move%zmove, delta_tau)
        else
            call propose_bath_shift_move(move%dmove, delta_tau)
        endif
    end subroutine

    logical function shift_is_update_allowed(move) result(t)
        type(TBathShiftMove), intent(in) :: move

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            t = is_update_allowed(move%zmove)
        else
            t = is_update_allowed(move%dmove)
        endif
    end function

    subroutine shift_compute_weight(move)
        type(TBathShiftMove), intent(inout) :: move

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            call compute_weight(move%zmove)
        else
            call compute_weight(move%dmove)
        endif
    end subroutine

    subroutine p_bath_accept_shift(move)
        type(TBathShiftMove), intent(inout)   :: move

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            call accept_update(move%zmove)
        else
            call accept_update(move%dmove)
        endif
    end subroutine

    function p_bath_absratio_shift(move) result(ratio)
        type(TBathShiftMove), intent(in)   :: move
        real(c_double) :: ratio

        if (allocated(move%zmove) .eqv. allocated(move%dmove)) &
            error stop 'Move is in an invalid state'

        if (allocated(move%zmove)) then
            ratio = weight_ratio(move%zmove)
        else
            ratio = weight_ratio(move%dmove)
        endif
    end function

end module
