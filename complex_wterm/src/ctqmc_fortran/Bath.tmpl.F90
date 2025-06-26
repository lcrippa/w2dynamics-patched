#if 0
! This is a templated Fortran file.
!
! By invoking `make`, several modules are generated in sub-directories
! double/ etc, for each valid value type (real, complex, etc.).
! In each generated module, the value VALUE_TYPE is replaced by the
! corresponding value type, and each uppercase derived type is replaced by
! a typed version.
#endif

module MBATH
    use iso_c_binding, only: c_double, c_double_complex
    use MBathBase
    use MBufferD
    use MBufferI
    use MBufferZ
    use MCommon
    use MEXT_RANGE
    use MHYBRIDIZATION
    use MPermutation
    use MUtilities
    use MSorting
    use MUPDATABLE_MATRIX
    use MBLOCK_MATRIX
    use testing
    implicit none
    private

    ! Re-exports
    public :: HYBRIDIZATION

    type, public :: BATH_TYPE
        private
        integer :: NOper = 0    !< Number of operators in the trace
        REAL_TYPE :: beta       !< Inverse temperature
        integer :: nbands       !< Number of impurity orbitals

        !> Specification of creators and annihilators, respectively, each
        !! ordered by tau. Note that, e.g., `crea` specifies the location of
        !! *impurity* creators, which correspond to bath annihilators.
        type(TBathOperator), allocatable :: crea(:)
        type(TBathOperator), allocatable :: annh(:)

        !> The hybridization function
        class(HYBRIDIZATION), allocatable :: hybr

        !> Indices for sorting the operators in blocks
        integer, allocatable :: sortidx(:, :), whichblock(:, :)

        !> The inverted matrix of the hybridization lines.
        type(BLOCK_MATRIX) :: mat

        !> Scaling factor for tau to be used in computing sortval
        REAL_TYPE :: sortval_tau_scaling
        REAL_TYPE, allocatable :: crea_sortval(:), annh_sortval(:)
    end type

    !> Base move for bath.
    !!
    !! Use propose_* routines to propose actual updates. You can erase the
    !! concrete update type by setting a pointer of type(BATH_MOVE_BASE)
    !! and then use weight_ratio() and accept_update().
    type, abstract, public :: BATH_MOVE_BASE
        private

        !> pointer to DBath, that saves properties of currently accepted
        !! bath-configuration
        type(BATH_TYPE), pointer :: parent => null()

        !> stores whether an update is allowed to be accepted
        logical :: is_update_allowed = .false.
    contains
        private
        procedure(base_accept_update), deferred :: accept_update
    end type

    abstract interface
        subroutine base_accept_update(move)
            import :: BATH_MOVE_BASE
            class(BATH_MOVE_BASE), intent(inout) :: move
        end
    end interface

    !> This module has functions to calculate the determinant and inverse
    !! hybridization matrix from scratch, given the arrays crea(:) and annh(:),
    !! which define the properties of a bath configuration.
    !! Upon acceptance, it can update MBath.
    type, extends(BATH_MOVE_BASE), public :: BATH_MOVE
        private
        integer :: newsize

        type(TBathOperator), allocatable :: crea(:) !> properties of new creators
        type(TBathOperator), allocatable :: annh(:) !> properties of new annihilators
        integer, allocatable :: crea_blocks(:), annh_blocks(:)
        type(BLOCK_REPLACE)     :: inner
    contains
        private
        procedure, non_overridable :: accept_update => p_bath_accept
    end type

    !> This module has functions to calculate the determinant and inverse
    !! hybridization matrix with the fast Sherman Morrison update formulas,
    !! given the properties of the added operators and the properties of the
    !! current bath configuration in MBath.
    !! Upon acceptance, it updates MBath.
    type, extends(BATH_MOVE_BASE), public :: BATH_ADD_MOVE
        private
        integer :: currsize, nadd, newsize

        type(TBathOperator), allocatable :: crea(:) !> properties of new creators
        type(TBathOperator), allocatable :: annh(:) !> properties of new annihilators
        integer, allocatable :: crea_pos(:), annh_pos(:)
        integer, allocatable :: crea_blocks(:), annh_blocks(:)
        REAL_TYPE, allocatable :: crea_sortval(:), annh_sortval(:)
        integer, allocatable :: crea_perm(:), annh_perm(:)

        VALUE_TYPE, allocatable :: q(:,:), r(:,:), s(:,:)
        type(BLOCK_GROW) :: inner
    contains
        private
        procedure, non_overridable :: accept_update => p_bath_accept_add
    end type

    !> This module has functions to calculate the determinant and inverse
    !! hybridization matrix with the fast Sherman Morrison update formulas,
    !! given the properties of the added operators and the properties of the
    !! current bath configuration in MBath.
    !! Upon acceptance, it updates MBath.
    type, extends(BATH_MOVE_BASE), public :: BATH_REM_MOVE
        private
        integer :: currsize, nrem, newsize

        type(TBathOperator), allocatable :: crea(:) !> properties of new creators
        type(TBathOperator), allocatable :: annh(:) !> properties of new annihilators
        integer, allocatable :: crea_pos(:), annh_pos(:)
        integer, allocatable :: crea_perm(:), annh_perm(:)

        type(BLOCK_SHRINK)      :: inner
    contains
        private
        procedure, non_overridable :: accept_update => p_bath_accept_rem
    end type

    !> This module has functions to calculate the determinant and inverse
    !! hybridization matrix with the fast Sherman Morrison update formulas,
    !! given the properties of the added operators and the properties of the
    !! current bath configuration in MBath.
    !! Upon acceptance, it updates MBath.
    type, extends(BATH_MOVE_BASE), public :: BATH_REPLACE_MOVE
        private
        integer :: currsize, nmove

        type(TBathOperator), allocatable :: ops(:), new_ops(:)
        integer, allocatable :: curr_pos(:), new_pos(:)
        integer, allocatable :: move_perm(:)
        integer(kind(OpDummy)) :: op_typ
        VALUE_TYPE, allocatable :: repl_rows(:,:), repl_cols(:,:)

        type(BLOCK_REPLACE_ROWS) :: inner_rows
        type(BLOCK_REPLACE_COLS) :: inner_cols
    contains
        private
        procedure, non_overridable :: accept_update => p_bath_accept_replace
    end type

    !> This module has functions to produce the
    !! hybridization matrix upon global tau-shifts.
    !! Upon acceptance, it updates MBath.
    type, extends(BATH_MOVE_BASE), public :: BATH_SHIFT_MOVE
        private
        integer :: currsize, crea_wrapped, annh_wrapped
        REAL_TYPE :: delta_tau

        type(TBathOperator), allocatable :: crea(:) !> properties of new creators
        type(TBathOperator), allocatable :: annh(:) !> properties of new annihilators
        integer, allocatable :: crea_perm(:), annh_perm(:)

        type(BLOCK_PERMUTE)     :: inner
    contains
        private
        procedure, non_overridable :: accept_update => p_bath_accept_shift
    end type

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

    interface print_bath
        module procedure p_print_bath
    end interface
    public print_bath

    interface getinv
        module procedure p_getinv
    end interface
    public getinv

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

    interface is_update_allowed
        module procedure move_is_update_allowed
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
        module procedure p_accept_update
    end interface
    public accept_update

    interface init_bath
        module procedure p_init_Bath
    end interface
    public init_bath

    interface clear_bath
        module procedure p_clear_bath
    end interface
    public clear_bath

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

    interface get_sortval
        module procedure p_get_sortval
    end interface
    public get_sortval

    interface setsize
        module procedure p_setsize
        module procedure p_setsize_move
        module procedure p_setsize_add
        module procedure p_setsize_rem
        module procedure p_setsize_repl
        module procedure p_setsize_shift
    end interface

    REAL_TYPE, parameter :: ONE_REAL = 1

 contains

    ! -------------------------------------------------------------------------
    ! BATH FUNCTIONS

    pure function p_annihilators(self) result(ops)
        type(BATH_TYPE), intent(in) :: self
        type(TBathOperator), allocatable :: ops(:)
        integer :: n

        n = self%noper / 2
        allocate(ops(n), source=self%annh(:n))
    end

    pure function p_creators(self) result(ops)
        type(BATH_TYPE), intent(in) :: self
        type(TBathOperator), allocatable :: ops(:)
        integer :: n

        n = self%noper / 2
        allocate(ops(n), source=self%crea(:n))
    end

    pure function p_get_beta(self) result(beta)
        type(BATH_TYPE), intent(in) :: self
        REAL_TYPE :: beta

        beta = self%beta
    end

    pure integer function p_size(self) result(r)
        type(BATH_TYPE), intent(in) :: self
        r = self%noper / 2
    end

    pure function p_getdet(self) result(det)
        type(BATH_TYPE), intent(in) :: self
        type(EXT_RANGE) :: det

        det = get_det(self%mat)
    end

    pure integer function p_getnbands(self) result(nbands)
        type(BATH_TYPE), intent(in) :: self
        nbands = self%nbands
    end

    subroutine p_setsize(self, new_size)
        type(BATH_TYPE), intent(inout) :: self
        integer, intent(in) :: new_size

        call reserve_opers(self%crea, new_size, .true.)
        call reserve_opers(self%annh, new_size, .true.)
        call reserve(self%crea_sortval, new_size, .true.)
        call reserve(self%annh_sortval, new_size, .true.)
        self%noper = 2 * new_size
    end subroutine

    pure elemental function p_get_sortval(self, op) result(val)
        class(BATH_TYPE), intent(in) :: self
        type(TBathOperator), intent(in) :: op
        REAL_TYPE :: val

        val = self%sortidx(op%orb, op%sp) + self%sortval_tau_scaling * op%tau
    end

    subroutine p_verify_bath(self, rtol)
        type(BATH_TYPE), intent(in) :: self
        REAL_TYPE, intent(in) :: rtol

        integer :: i, k
        REAL_TYPE :: prev_val
        type(TBathOperator) :: curr
        complex(c_double_complex) :: detratio

        VALUE_TYPE, allocatable :: mat(:, :)
        type(UPDATABLE_MATRIX) :: umat

        if (.not. (self%beta > 0)) &
            error stop 'Invalid inverse temperature'
        if (self%nbands < 1) &
            error stop 'Invalid number of bands'

        call assert_equal(allocated(self%crea), .true.)
        call assert_equal(allocated(self%annh), .true.)

        k = self%Noper / 2
        if (self%Noper /= 2 * k) &
            error stop 'Number of operators must be even'
        if (size(self%crea) < k) &
            error stop 'Bath creators too short'
        if (size(self%annh) < k) &
            error stop 'Bath annihilators too short'
        if (size(self%crea_sortval) < k) &
            error stop 'Bath creators too short'
        if (size(self%annh_sortval) < k) &
            error stop 'Bath annihilators too short'

        if (any(shape(self%sortidx) /= (/self%nbands, 2/))) &
            error stop 'Wrong shape of sortidx'
        if (any(shape(self%whichblock) /= (/self%nbands, 2/))) &
            error stop 'Wrong shape of sortidx'

        prev_val = 0
        do i = 1, k
            curr = self%crea(i)
            if (curr%tau < 0 .or. curr%tau > self%beta) &
                error stop 'Invalid time for operator'
            if (curr%sp < 1 .or. curr%sp > 2) &
                error stop 'Invalid spin for operator'
            if (curr%orb < 1 .or. curr%orb > self%nbands) &
                error stop 'Invalid orbital for operator'
            if (.not. (get_sortval(self, curr) >= prev_val)) &
                error stop 'Operators not sorted properly'
            prev_val = get_sortval(self, curr)
        enddo

        prev_val = 0
        do i = 1, k
            curr = self%annh(i)
            if (curr%tau < 0 .or. curr%tau > self%beta) &
                error stop 'Invalid time for operator'
            if (curr%sp < 1 .or. curr%sp > 2) &
                error stop 'Invalid spin for operator'
            if (curr%orb < 1 .or. curr%orb > self%nbands) &
                error stop 'Invalid orbital for operator'
            if (.not. (get_sortval(self, curr) >= prev_val)) &
                error stop 'Operators not sorted properly'
            prev_val = get_sortval(self, curr)
        enddo

        if (size(self%mat) /= k) &
            error stop 'Matrix is of wrong size'

        allocate (mat(k, k))
        call verify(self%mat, rtol)
        call copy_hybr_matrix(self%hybr, self%annh(1:k), &
                              self%crea(1:k), mat)
        call setmatrix(umat, mat)

        detratio = limrange(get_det(umat) / get_det(self))
        call assert_close(detratio, cmplx(1, 0, c_double_complex), rtol=rtol)
        call assert_close(getinv(umat), getinv(self%mat), rtol=rtol)

        call assert_close( &
                self%crea_sortval(:k), get_sortval(self, self%crea(:k)), &
                atol=0d0)
        call assert_close( &
                self%annh_sortval(:k), get_sortval(self, self%annh(:k)), &
                atol=0d0)
    end

    !> print the bath properties on the screen
    subroutine p_print_bath(this)
        type(BATH_TYPE)      :: this
        integer          :: i

        write(*,*) " "
        write(*,*) "bath creators:"
        write(*,*) "N, orb, spin, tau"
        do i=1,this%noper/2
        write(*,*) i, this%crea(i)%orb, this%crea(i)%sp, this%crea(i)%tau
        enddo

        write(*,*) " "
        write(*,*) "bath annihilators:"
        write(*,*) "N, orb, spin, tau"
        do i=1,this%noper/2
        write(*,*) i, this%annh(i)%orb, this%annh(i)%sp, this%annh(i)%tau
        enddo
    end subroutine

    !> Copy inverse matrix
    function p_getinv(this) result(inv)
        type(BATH_TYPE), intent(in) :: this
        VALUE_TYPE, allocatable :: inv(:, :)

        inv = getinv(this%mat)
    end

    !> prints bath configs in same format as local configs
    subroutine p_debug_print_bathconfig_highprec(this)
        type(BATH_TYPE), intent(in) :: this
        integer :: i

        write(*, "('bcfg with ', I5, ' ops')") this%Noper
        write(*, "('creators:')")
        write(*, "('  iop        tau orb sp')")
        do i = 1, this%Noper/2
            write(*, "(I5, '   ', F19.16, I3, I3)") i, this%crea(i)%tau, &
                this%crea(i)%orb, this%crea(i)%sp
        end do
        write(*, "('annihilators:')")
        write(*, "('  iop        tau orb sp')")
        do i = 1, this%Noper/2
            write(*, "(I5, '   ', F19.16, I3, I3)") i, this%annh(i)%tau, &
                this%annh(i)%orb, this%annh(i)%sp
        end do
    end subroutine       ! FIXME: replace or remove this for final version of new code

    pure integer function p_get_nosoper(bath, orb, sp) result(nosoper)
        type(BATH_TYPE), intent(in) :: bath
        integer, intent(in)     :: orb, sp
        integer                 :: i

        nosoper = 0
        do i = 1, bath%NOper/2
        if (bath%crea(i)%orb == orb .and. bath%crea(i)%sp == sp)&
            nosoper = nosoper + 1

        if (bath%annh(i)%orb == orb .and. bath%annh(i)%sp == sp)&
            nosoper = nosoper + 1
        end do
    end function

    !> Allocate stuff and set to initial values
    subroutine p_init_Bath(this, hybr)
        type(BATH_TYPE), intent(out) :: this
        class(HYBRIDIZATION), intent(in) :: hybr

        allocate(this%hybr, source=hybr)
        if (.not. is_fermionic(this%hybr)) &
            error stop 'Hybridization function must be fermionic'

        this%beta = get_beta(this%hybr)
        this%nbands = get_nbands(this%hybr)
        call init_bath_lookups(this)

        ! We have to make sure that any tau in [0, beta) that is scaled with
        ! this value stays strictly between [0, 1) even if we add an integer
        ! corresponding to the spinorbital.
        this%sortval_tau_scaling = &
            (1 - 2 * this%nbands * epsilon(ONE_REAL)) / this%beta

        call initialize(this%mat, maxval(this%whichblock))
        call setsize(this, 0)
    end subroutine

    subroutine init_bath_lookups(self)
        type(BATH_TYPE), intent(inout) :: self

        integer :: norb, nsp
        integer, allocatable :: whichblock(:), sortidx(:), idx_structure(:,:,:,:)
        logical, allocatable :: adj_flat(:,:)

        allocate(idx_structure, source=get_index_structure(self%hybr))
        norb = size(idx_structure, 1)
        nsp = size(idx_structure, 2)
        adj_flat = reshape(idx_structure /= 0, (/ norb*nsp, norb*nsp /))
        whichblock = connected_components(adj_flat)

        ! if p is a group permutation, then j = p(i) is the flavour that will
        ! end up in the i'th position of a sorted array. But we want just the
        ! opposite: the given a flavour j, which position i shall we place
        ! it? This is just the inverse permutation:
        sortidx = get_inv_perm(get_group_perm(whichblock))

        self%whichblock = reshape(whichblock, (/ norb, nsp /))
        self%sortidx = reshape(sortidx, (/ norb, nsp /))
    end subroutine

    subroutine p_clear_bath(self)
        type(BATH_TYPE), intent(inout) :: self
        integer :: block_idx(0)
        VALUE_TYPE :: mat(0, 0)

        self%NOper = 0
        call setmatrix(self%mat, block_idx, mat)
    end subroutine

    subroutine update_caches(self)
        type(BATH_TYPE), intent(inout) :: self
        integer :: n

        n = self%noper / 2
        self%crea_sortval(:n) = get_sortval(self, self%crea(:n))
        self%annh_sortval(:n) = get_sortval(self, self%annh(:n))
    end subroutine

    ! -------------------------------------------------------------------------
    ! BATH MOVE BASE

    pure logical function move_is_update_allowed(self) result(t)
        class(BATH_MOVE_BASE), intent(in) :: self

        t = self%is_update_allowed
    end

    subroutine p_accept_update(move)
        class(BATH_MOVE_BASE), intent(inout) :: move

        call move%accept_update()
    end subroutine

    ! ========================================================================
    ! BATH MOVE
    ! ========================================================================

    subroutine p_init_bath_move(this, parent)
        type(BATH_MOVE), allocatable :: this
        type(BATH_TYPE), intent(in), target   :: parent

        allocate(this)
        this%parent => parent
    end subroutine

    subroutine p_setsize_move(self, newsize)
        type(BATH_MOVE), intent(inout) :: self
        integer, intent(in) :: newsize

        call reserve_opers(self%crea, newsize)
        call reserve_opers(self%annh, newsize)
        call reserve(self%crea_blocks, newsize)
        call reserve(self%annh_blocks, newsize)
        self%newsize = newsize
    end subroutine

    !> Initialize move with properties of proposed configuration.
    subroutine p_propose_bath_move(move, crea, annh)
        type(BATH_MOVE), intent(inout) :: move
        type(TBathOperator), intent(in) :: crea(:), annh(:)

        integer :: crea_ins_perm(size(crea)), annh_ins_perm(size(annh))
        integer :: i, n

        n = size(crea)
        if (size(annh) /= n) &
            error stop 'Size of creators and annihilators must match'

        call setsize(move, n)
        crea_ins_perm = get_sort_perm(get_sortval(move%parent, crea))
        annh_ins_perm = get_sort_perm(get_sortval(move%parent, annh))

        move%crea(:n) = crea(crea_ins_perm)
        move%annh(:n) = annh(annh_ins_perm)

        do i = 1, n
            move%crea_blocks(i) = &
                    move%parent%whichblock(move%crea(i)%orb, move%crea(i)%sp)
            move%annh_blocks(i) = &
                    move%parent%whichblock(move%annh(i)%orb, move%annh(i)%sp)
        end do

        move%is_update_allowed = &
                all(move%crea_blocks(:n) == move%annh_blocks(:n))
    end

    subroutine move_compute_weight(move)
        type(BATH_MOVE), intent(inout) :: move
        VALUE_TYPE :: mat(move%newsize, move%newsize)

        call copy_hybr_matrix( &
                move%parent%hybr, move%annh(:move%newsize), &
                move%crea(:move%newsize), mat)

        call propose_replace_matrix( &
                move%inner, move%parent%mat, &
                move%annh_blocks(:move%newsize), &
                move%crea_blocks(:move%newsize), mat)

        call assert_equal(is_update_allowed(move%inner), is_update_allowed(move))
    end subroutine

    !> Calculated ratio of determinants and saves new M-matrix and new determinant into structure
    function p_bath_absratio(move) result(ratio)
        type(BATH_MOVE) :: move
        REAL_TYPE :: ratio
        type(EXT_RANGE) :: detrat

        if (.not. move%is_update_allowed) then
            ratio = 0
            return
        endif

        detrat = get_det(move%inner) / get_det(move%parent)
        ratio = abs(limrange(detrat))
    end function

    !> Copies stuff into MBath (=parent), which saves the recently accepted configuration
    subroutine p_bath_accept(move)
        class(BATH_MOVE), intent(inout) :: move

        call setsize(move%parent, move%newsize)

        move%parent%crea(:move%newsize) = move%crea(:move%newsize)
        move%parent%annh(:move%newsize) = move%annh(:move%newsize)

        call accept_update(move%inner)
        move%parent%noper = move%newsize * 2

        call update_caches(move%parent)

        !call verify_bath(move%parent, 5d-2)
    end subroutine

    ! -------------------------------------------------------------------------
    ! BATH ADD MOVE

    subroutine p_init_bath_add_move(this, parent)
        type(BATH_ADD_MOVE), allocatable :: this
        type(BATH_TYPE), intent(in), target   :: parent

        allocate(this)
        this%parent => parent

        allocate(this%q(1,1), this%r(1,1), this%s(1,1))
    end subroutine

    subroutine p_setsize_add(self, newsize, newk)
        type(BATH_ADD_MOVE), intent(inout) :: self
        integer, intent(in) :: newsize, newk

        call reserve_opers(self%crea, newsize)
        call reserve_opers(self%annh, newsize)
        call reserve(self%crea_pos, newk)
        call reserve(self%annh_pos, newk)
        call reserve(self%crea_sortval, newk)
        call reserve(self%annh_sortval, newk)
        call reserve(self%crea_blocks, newk)
        call reserve(self%annh_blocks, newk)
        call reserve(self%crea_perm, newsize)
        call reserve(self%annh_perm, newsize)
        call reserve(self%q, newsize, newk)
        call reserve(self%r, newk, newsize)
        call reserve(self%s, newk, newk)

        self%currsize = newsize - newk
        self%newsize = newsize
        self%nadd = newk
    end subroutine

    !> Moves new operator allocations for proposal into the move object and sets count
    subroutine p_propose_bath_add_move(move, crea, annh)
        type(BATH_ADD_MOVE)                    :: move
        type(TBathOperator), intent(in)        :: crea(:), annh(:)

        integer :: crea_ins_perm(size(crea)), annh_ins_perm(size(annh))
        integer :: i
        type(TBathOperator) :: crea_sorted(size(crea)), annh_sorted(size(annh))

        if (size(crea) /= size(annh)) &
            error stop 'Creators and annihilators must be of same size'

        call setsize(move, move%parent%noper/2 + size(crea), size(crea))

        move%crea(:move%currsize) = move%parent%crea(:move%currsize)
        move%annh(:move%currsize) = move%parent%annh(:move%currsize)

        move%crea_sortval(:move%nadd) = get_sortval(move%parent, crea)
        move%annh_sortval(:move%nadd) = get_sortval(move%parent, annh)
        call sort_perm(move%crea_sortval(:move%nadd), crea_ins_perm)
        call sort_perm(move%annh_sortval(:move%nadd), annh_ins_perm)
        crea_sorted = crea(crea_ins_perm)
        annh_sorted = annh(annh_ins_perm)

        do i = 1, move%nadd
            move%crea_blocks(i) = &
                move%parent%whichblock(crea_sorted(i)%orb, crea_sorted(i)%sp)
            move%annh_blocks(i) = &
                move%parent%whichblock(annh_sorted(i)%orb, annh_sorted(i)%sp)
        end do

        move%crea(move%currsize+1:move%newsize) = crea_sorted
        move%annh(move%currsize+1:move%newsize) = annh_sorted

        move%is_update_allowed = &
            all(move%crea_blocks(:move%nadd) == move%annh_blocks(:move%nadd))
    end

    subroutine add_compute_weight(move)
        type(BATH_ADD_MOVE), intent(inout) :: move

        call sorted_insert( &
                move%parent%crea_sortval(:move%currsize), &
                move%crea_sortval(:move%nadd), &
                move%crea_pos(:move%nadd))
        call grow_perm( &
                move%currsize, move%crea_pos(:move%nadd), &
                move%crea_perm(:move%newsize))

        call sorted_insert( &
                move%parent%annh_sortval(:move%currsize), &
                move%annh_sortval(:move%nadd), &
                move%annh_pos(:move%nadd))
        call grow_perm( &
                move%currsize, move%annh_pos(:move%nadd), &
                move%annh_perm(:move%newsize))

        call copy_hybr_matrix( &
                move%parent%Hybr, &
                move%annh(:move%currsize), &
                move%crea(move%currsize+1:move%newsize), &
                move%q(:move%currsize, :move%nadd))
        call copy_hybr_matrix( &
                move%parent%Hybr, &
                move%annh(move%currsize+1:move%newsize), &
                move%crea(:move%currsize), &
                move%r(:move%nadd, :move%currsize))
        call copy_hybr_matrix( &
                move%parent%Hybr, &
                move%annh(move%currsize+1:move%newsize), &
                move%crea(move%currsize+1:move%newsize), &
                move%s(:move%nadd, :move%nadd))

        call propose_grow_matrix( &
                move%inner, move%parent%mat, &
                move%annh_blocks(:move%nadd), move%crea_blocks(:move%nadd), &
                move%annh_pos(:move%nadd), move%crea_pos(:move%nadd), &
                move%r(:move%nadd, :move%currsize), &
                move%q(:move%currsize, :move%nadd), &
                move%s(:move%nadd, :move%nadd))

        call assert_equal(is_update_allowed(move%inner), is_update_allowed(move))
    end subroutine

    !> Calculates determinant ratio between new and old bath configuration
    function p_bath_absratio_add(move) result(ratio)
        type(BATH_ADD_MOVE), intent(inout)   :: move
        type(EXT_RANGE)                       :: detrat
        REAL_TYPE :: ratio

        if (.not. move%is_update_allowed) then
            ratio = 0
            return
        endif

        detrat = get_det(move%inner) / get_det(move%parent)
        ratio = abs(limrange(detrat))
    end function p_bath_absratio_add

    !> Calculates new M-matrix with Sherman Morrison formula and copies
    !! properties of newly accepted configuration into MBath structure
    subroutine p_bath_accept_add(move)
        class(BATH_ADD_MOVE), intent(inout)   :: move

        if (.not. move%is_update_allowed) &
            error stop 'Update disallowed to be accepted'
        move%is_update_allowed = .false.

        call setsize(move%parent, move%newsize)

        move%parent%crea(:move%newsize) = move%crea(move%crea_perm(:move%newsize))
        move%parent%annh(:move%newsize) = move%annh(move%annh_perm(:move%newsize))

        call accept_update(move%inner)
        move%parent%noper = move%newsize * 2

        call update_caches(move%parent)
        !call verify_bath(move%parent, 5d-2)
    end subroutine p_bath_accept_add

    ! -------------------------------------------------------------------------
    ! BATH REMOVE MOVE

    subroutine p_init_bath_rem_move(this, parent)
        type(BATH_REM_MOVE), allocatable :: this
        type(BATH_TYPE), intent(in), target   :: parent

        allocate(this)
        this%parent => parent

        ! XXX clean this up
        allocate(this%crea_pos(1))
        allocate(this%annh_pos(1))
        allocate(this%crea(2))
        allocate(this%annh(2))
        allocate(this%crea_perm(2))
        allocate(this%annh_perm(2))
    end subroutine

    subroutine p_setsize_rem(self, newsize, newk)
        type(BATH_REM_MOVE), intent(inout) :: self
        integer, intent(in) :: newsize, newk

        call reserve_opers(self%crea, newsize)
        call reserve_opers(self%annh, newsize)
        call reserve(self%crea_pos, newk)
        call reserve(self%annh_pos, newk)
        call reserve(self%crea_perm, newsize)
        call reserve(self%annh_perm, newsize)
    end subroutine

    !> Create proposal based on indices of operators to be removed
    !> (remidx(1): creator, remidx(2): annihilator)
    subroutine p_propose_bath_rem_move(move, remidx)
        type(BATH_REM_MOVE), intent(inout)     :: move
        integer, intent(in) :: remidx(:)
        integer :: crea_blocks(size(remidx)/2), annh_blocks(size(remidx)/2)
        integer :: i
        type(TBathOperator) :: crea_curr, annh_curr

        move%currsize = move%parent%noper/2
        move%nrem = size(remidx) / 2
        move%newsize = move%currsize - move%nrem

        call setsize(move, move%currsize, move%nrem)

        move%crea_pos(:move%nrem) = remidx(1::2)
        move%annh_pos(:move%nrem) = remidx(2::2)

        ! XXX we're duplicating some effort of the block matrix layer here:
        !     ordering indices and figuring out block indices. This is in
        !     service of figuring out as early as possible whether the move
        !     is valid.
        if (move%nrem > 1) then
            call sort(move%crea_pos(:move%nrem))
            call sort(move%annh_pos(:move%nrem))
        endif

        do i = 1, move%nrem
            crea_curr = move%parent%crea(move%crea_pos(i))
            crea_blocks(i) = &
                move%parent%whichblock(crea_curr%orb, crea_curr%sp)

            annh_curr = move%parent%crea(move%annh_pos(i))
            annh_blocks(i) = &
                move%parent%whichblock(annh_curr%orb, annh_curr%sp)
        end do

        move%is_update_allowed = all(crea_blocks == annh_blocks)
    end

    subroutine rem_compute_weight(move)
        type(BATH_REM_MOVE), intent(inout) :: move

        call shrink_perm( &
                move%currsize, move%crea_pos(:move%nrem), &
                move%crea_perm(:move%currsize))
        call shrink_perm( &
                move%currsize, move%annh_pos(:move%nrem), &
                move%annh_perm(:move%currsize))

        move%crea(:move%currsize) = &
                move%parent%crea(move%crea_perm(:move%currsize))
        move%annh(:move%currsize) = &
                move%parent%annh(move%annh_perm(:move%currsize))

        call propose_shrink_matrix( &
                    move%inner, move%parent%mat, &
                    move%annh_pos(:move%nrem), move%crea_pos(:move%nrem))

        call assert_equal(is_update_allowed(move%inner), is_update_allowed(move))
    end subroutine

    !> Calculates determinant ratio between new and old bath configuration
    function p_bath_absratio_rem(move) result(ratio)
        type(BATH_REM_MOVE), intent(inout)   :: move
        type(EXT_RANGE)                       :: detrat
        REAL_TYPE :: ratio

        if (.not. move%is_update_allowed) then
            ratio = 0
            return
        endif

        detrat = get_det(move%inner) / get_det(move%parent)
        ratio = abs(limrange(detrat))
    end function p_bath_absratio_rem

    !> Generates new M-matrix with Sherman Morrison and copies stuff into
    !! MBath structure
    subroutine p_bath_accept_rem(move)
        class(BATH_REM_MOVE), intent(inout)    :: move

        if (.not. move%is_update_allowed) &
            error stop 'Update disallowed to be accepted'
        move%is_update_allowed = .false.

        call setsize(move%parent, move%newsize)

        move%parent%crea(:move%newsize) = move%crea(:move%newsize)
        move%parent%annh(:move%newsize) = move%annh(:move%newsize)

        call accept_update(move%inner)
        move%parent%noper = move%newsize * 2

        call update_caches(move%parent)
        !call verify_bath(move%parent, 5d-2)
    end subroutine p_bath_accept_rem

    ! -------------------------------------------------------------------------
    ! BATH REPLACE ROW/COLUMN MOVE

    subroutine p_init_bathreplacemove(self, parent)
        type(BATH_REPLACE_MOVE), allocatable :: self
        type(BATH_TYPE), intent(in), target   :: parent

        allocate(self)
        self%parent => parent

        allocate(self%repl_rows(1,1), self%repl_cols(1,1))
    end subroutine

    subroutine p_setsize_repl(self, newsize, newk)
        type(BATH_REPLACE_MOVE), intent(inout) :: self
        integer, intent(in) :: newsize, newk

        call reserve_opers(self%ops, newsize)
        call reserve_opers(self%new_ops, newk)
        call reserve(self%move_perm, newsize)
        call reserve(self%curr_pos, newk)
        call reserve(self%new_pos, newk)
        call reserve(self%repl_rows, newk, newsize)
        call reserve(self%repl_cols, newsize, newk)
        self%currsize = newsize
        self%nmove = newk
    end subroutine

    ! Propose replacement of bath operators without change of total number.
    !
    ! move:       proposal
    ! origidx:    index of operator to be replaced in annihilator or creator list
    ! origtype:   type of replaced operator (1/2: C/A)
    ! newop:      new operator to be substitued for the replaced one
    subroutine p_propose_bathreplacemove(move, origidx, origtype, newop)
        type(BATH_REPLACE_MOVE), intent(inout)   :: move
        integer, intent(in)                   :: origidx
        integer(kind(OpDummy)), intent(in)     :: origtype
        type(TBathOperator), intent(in)       :: newop

        type(TBathOperator) :: origop
        integer :: origblock, newblock

        call setsize(move, move%parent%noper/2, 1)

        move%op_typ = origtype
        move%curr_pos(1) = origidx
        move%new_ops(1) = newop
        select case (move%op_typ)
        case (OpCrea)
            move%ops(:move%currsize) = move%parent%crea(:move%currsize)
            origop = move%parent%crea(origidx)
        case (OpAnnh)
            move%ops(:move%currsize) = move%parent%annh(:move%currsize)
            origop = move%parent%annh(origidx)
        case default
            error stop 'Invalid move type'
        end select

        origblock = move%parent%whichblock(origop%orb, origop%sp)
        newblock = move%parent%whichblock(newop%orb, newop%sp)
        move%is_update_allowed = origblock == newblock
    end

    subroutine repl_compute_weight(move)
        type(BATH_REPLACE_MOVE), intent(inout) :: move
        type(TBathOperator) :: newop
        integer :: origidx

        newop = move%new_ops(1)
        origidx = move%curr_pos(1)
        call sorted_insert(get_sortval(move%parent, move%ops(:move%currsize)), &
                           get_sortval(move%parent, newop), &
                           move%new_pos(1))
        if (move%new_pos(1) > move%curr_pos(1)) &
            move%new_pos(1) = move%new_pos(1) - 1

        move%ops(origidx) = newop
        call move_perm(move%currsize, move%curr_pos(:1), move%new_pos(:1), &
                        move%move_perm)

        select case(move%op_typ)
        case (OpCrea)
            call copy_hybr_matrix( &
                    move%parent%hybr, &
                    move%parent%annh(:move%currsize), &
                    (/ newop /), &
                    move%repl_cols(:move%currsize, :move%nmove))
            call propose_replace_cols( &
                    move%inner_cols, move%parent%mat, &
                    move%curr_pos(:move%nmove), move%new_pos(:move%nmove), &
                    move%repl_cols(:move%currsize, :move%nmove))

            call assert_equal( &
                    is_update_allowed(move%inner_cols), is_update_allowed(move))
        case (OpAnnh)
            call copy_hybr_matrix( &
                    move%parent%hybr, &
                    (/ newop /), &
                    move%parent%crea(:move%currsize), &
                    move%repl_rows(:move%nmove, :move%currsize))
            call propose_replace_rows( &
                    move%inner_rows, move%parent%mat, &
                    move%curr_pos(:move%nmove), move%new_pos(:move%nmove), &
                    move%repl_rows(:move%nmove, :move%currsize))

            call assert_equal( &
                    is_update_allowed(move%inner_rows), is_update_allowed(move))
        end select
    end subroutine

    function p_bath_absratio_replace(move) result(ratio)
        type(BATH_REPLACE_MOVE), intent(inout)  :: move
        type(EXT_RANGE)                         :: detrat
        REAL_TYPE :: ratio

        if (.not. move%is_update_allowed) then
            ratio = 0
            return
        endif

        select case (move%op_typ)
        case (OpCrea)
            detrat = get_det(move%inner_cols)
        case (OpAnnh)
            detrat = get_det(move%inner_rows)
        end select

        detrat = detrat / get_det(move%parent)
        ratio = abs(limrange(detrat))
    end function p_bath_absratio_replace

    !> Generates new M-matrix with Sherman Morrison and copies stuff into
    !! MBath structure
    subroutine p_bath_accept_replace(move)
        class(BATH_REPLACE_MOVE), intent(inout)    :: move

        if (.not. move%is_update_allowed) &
            error stop 'Update disallowed to be accepted'
        move%is_update_allowed = .false.

        select case(move%op_typ)
        case (OpCrea)
            move%parent%crea(:move%currsize) = &
                move%ops(move%move_perm(:move%currsize))
            call accept_update(move%inner_cols)

        case (OpAnnh)
            move%parent%annh(:move%currsize) = &
                move%ops(move%move_perm(:move%currsize))
            call accept_update(move%inner_rows)
        end select

        call update_caches(move%parent)
        !call verify_bath(move%parent, 5d-2)
    end subroutine p_bath_accept_replace

    ! -------------------------------------------------------------------------
    ! BATH SHIFT MOVE

    subroutine p_init_bath_shift_move(this, parent)
        type(BATH_SHIFT_MOVE), allocatable :: this
        type(BATH_TYPE), intent(in), target   :: parent

        allocate(this)
        this%parent => parent
    end subroutine

    subroutine p_setsize_shift(self, newsize)
        type(BATH_SHIFT_MOVE), intent(inout) :: self
        integer, intent(in) :: newsize

        call reserve_opers(self%crea, newsize)
        call reserve_opers(self%annh, newsize)
        call reserve(self%crea_perm, newsize)
        call reserve(self%annh_perm, newsize)
        self%currsize = newsize
    end subroutine

    !> Copies properties of tau-shift into structure, and counts how many
    !! operators end up at higher tau than before after the cyclic shift.
    subroutine p_propose_bath_shift_move(move, delta_tau)
        type(BATH_SHIFT_MOVE), intent(inout)   :: move
        REAL_TYPE, intent(in)               :: delta_tau

        call setsize(move, move%parent%noper/2)

        if (delta_tau < 0 .or. delta_tau > move%parent%beta) &
            error stop 'Invalid diff'
        move%delta_tau = delta_tau
        move%is_update_allowed = .true.
    end

    subroutine shift_compute_weight(move)
        type(BATH_SHIFT_MOVE), intent(inout) :: move
        integer :: i
        VALUE_TYPE :: iscale(move%parent%noper/2)
        VALUE_TYPE :: jscale(move%parent%noper/2)

        do i = 1, move%currsize
            move%crea(i) = move%parent%crea(i)
            move%crea(i)%tau = move%crea(i)%tau - move%delta_tau
            if (move%crea(i)%tau < 0) then
                move%crea(i)%tau = move%crea(i)%tau + move%parent%beta
                jscale(i) = -1
            else
                jscale(i) = +1
            endif
        end do
        move%crea_perm(:move%currsize) &
            = get_sort_perm(get_sortval(move%parent, move%crea(:move%currsize)))

        do i = 1, move%currsize
            move%annh(i) = move%parent%annh(i)
            move%annh(i)%tau = move%annh(i)%tau - move%delta_tau
            if (move%annh(i)%tau < 0) then
                move%annh(i)%tau = move%annh(i)%tau + move%parent%beta
                iscale(i) = -1
            else
                iscale(i) = +1
            endif
        end do
        move%annh_perm(:move%currsize) &
            = get_sort_perm(get_sortval(move%parent, move%annh(:move%currsize)))

        move%crea(:move%currsize) = move%crea(move%crea_perm(:move%currsize))
        move%annh(:move%currsize) = move%annh(move%annh_perm(:move%currsize))

        call propose_permute_matrix( &
                move%inner, move%parent%mat, &
                move%annh_perm(:move%currsize), move%crea_perm(:move%currsize), &
                iscale, jscale)

        move%is_update_allowed = is_update_allowed(move%inner)
    end subroutine

    function p_bath_absratio_shift(move) result(ratio)
        type(BATH_SHIFT_MOVE), intent(in)   :: move
        REAL_TYPE :: ratio

        ratio = 1.0
    end function

    !> Calculates the new inverse hybridization matrix by shifting the rows and
    !! columns. Calculates the new arrays that store the properties of the bath
    !! configuration. Calculates the new sign.
    !! Copies everything into DBath.
    subroutine p_bath_accept_shift(move)
        class(BATH_SHIFT_MOVE), intent(inout)   :: move

        if (.not. move%is_update_allowed) &
            error stop 'Update disallowed to be accepted'
        move%is_update_allowed = .false.

        move%parent%crea(:move%currsize) = move%crea(:move%currsize)
        move%parent%annh(:move%currsize) = move%annh(:move%currsize)

        call accept_update(move%inner)

        call update_caches(move%parent)
        !call verify_bath(move%parent%BATH_TYPE, 5d-2)
    end subroutine p_bath_accept_shift

    ! ========================================================================
    ! HELPER ROUTINES
    ! ========================================================================

    subroutine reserve_opers(buf, rowcap, preserve)
        type(TBathOperator), allocatable, intent(inout) :: buf(:)
        integer, intent(in) :: rowcap
        logical, intent(in), optional :: preserve
        type(TBathOperator), allocatable :: tmp(:)

        if (allocated(buf)) then
            if (rowcap <= size(buf, 1)) &
                return
            if (present(preserve)) then
                if (preserve) then
                    call move_alloc(buf, tmp)
                else
                    deallocate(buf)
                endif
            else
                deallocate(buf)
            endif
        endif
        allocate(buf(rowcap))
        if (allocated(tmp)) then
            buf(:size(tmp)) = tmp(:)
        endif
    end subroutine

end module
