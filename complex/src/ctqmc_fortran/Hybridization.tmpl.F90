#if 0
! This is a templated Fortran file.
!
! By invoking `make`, several modules are generated in sub-directories
! double/ etc, for each valid value type (real, complex, etc.).
! In each generated module, the value VALUE_TYPE is replaced by the
! corresponding value type, and each uppercase derived type is replaced by
! a typed version.
#endif

module MHYBRIDIZATION
    use iso_c_binding, only: c_double, c_double_complex, c_int64_t
    use MCommon, only: TBaseOper
    use MComplexNumbers, only: nan
    use MLogDiscretization
    implicit none
    private

    public :: TBaseOper

    !> Hybridization function.
    !!
    !! Stores the hybridization function. Given hybridization stengths `V` and
    !! bath energies `p`, and inverse temperature `beta, the hybridization
    !! function is given as follows for positive imaginary times `tau`:
    !! $$
    !!   \Delta_{i\sigma i'\sigma'j}(\tau) = \sum_{p} V_{p,i\sigma}
    !!        \frac{\exp(-\tau\epsilon_p)}{1+ \exp(-\beta\epsilon_p)}
    !!        V^*_{p,i'\sigma'}
    !! $$
    !!
    !! This class stores the function in some discretization. User codes will
    !! mostly interact with `copy_hybr_matrix(...)`.
    !!
    type, abstract, public :: HYBRIDIZATION
        private
        integer :: norb = -9999
        logical :: fermion
        REAL_TYPE :: beta
        integer, allocatable :: lookup(:,:,:,:)
    contains
        private
        procedure(hybr_matrix), deferred :: hybr_matrix
    end type

    interface
        subroutine hybr_matrix(self, iop, jop, out)
            import :: HYBRIDIZATION, TBaseOper, c_double, c_double_complex
            class(HYBRIDIZATION), intent(in) :: self
            class(TBaseOper), intent(in) :: iop(:), jop(:)
            VALUE_TYPE, intent(out) :: out(:, :)
        end subroutine
    end interface

    type, extends(HYBRIDIZATION) :: hybr_linear
        private
        REAL_TYPE :: invdtau
        VALUE_TYPE, allocatable :: mat(:,:)
    contains
        private
        procedure, non_overridable :: hybr_matrix => hybr_linear_matrix
    end type

    type, extends(HYBRIDIZATION) :: hybr_log
        private
        type(TLogGrid) :: grid
        VALUE_TYPE, allocatable :: mat(:,:,:,:)
    contains
        private
        procedure, non_overridable :: hybr_matrix => hybr_log_matrix
    end type

    interface init_linear_hybridization
        module procedure p_linear_init_bandspin_alloc
    end interface
    public :: init_linear_hybridization

    interface init_log_hybridization
        module procedure p_init_log_hybridization
    end interface
    public :: init_log_hybridization

    !> Return inverse temperature
    interface get_beta
        module procedure p_get_beta
    end interface
    public get_beta

    !> Return number of impurity orbitals
    interface get_nbands
        module procedure p_get_nbands
    end interface
    public get_nbands

    !> Return true if the statistics of the function is fermionic
    interface is_fermionic
        module procedure p_is_fermionic
    end interface
    public is_fermionic

    !> Hybridization lines between a set of impurity creators and annihilators
    interface copy_hybr_matrix
        module procedure p_copy_hybr_matrix
    end interface
    public :: copy_hybr_matrix

    !> Return value of hybridization function
    interface get_hybr_value
        module procedure p_get_hybr_value
    end interface
    public :: get_hybr_value

    !> Return the symmetry structure of the hybridization function
    interface get_index_structure
        module procedure p_index_structure
    end interface
    public :: get_index_structure

    !> Helper function for getting a propagator in tau
    interface propagator_tau
        module procedure p_propagator_tau
    end interface
    public :: propagator_tau

    !> Helper function for turning a log into linear discretization
    interface discretize_propagator
        module procedure p_discretize_propagator
    end interface
    public discretize_propagator

contains

    ! =========================================================================
    ! GENERIC FUNCTIONS

    subroutine init_base_hybridization(self, norb, beta, fermion)
        class(HYBRIDIZATION), intent(out) :: self
        integer, intent(in) :: norb
        REAL_TYPE, intent(in) :: beta
        logical, intent(in) :: fermion

        if (.not. (beta >= 0 .and. beta <= HUGE(beta))) &
            error stop 'Invalid inverse temperature'
        if (norb < 0) &
            error stop 'Number of orbitals must be non-negative'

        self%norb = norb
        self%beta = beta
        self%fermion = fermion
    end subroutine

    subroutine p_copy_hybr_matrix(self, iop, jop, out)
        class(HYBRIDIZATION), intent(in) :: self
        class(TBaseOper), intent(in) :: iop(:), jop(:)
        VALUE_TYPE, intent(out) :: out(:, :)

        ! check size consistency
        if (size(out, 1) /= size(iop)) &
            error stop 'Invalid number of rows'
        if (size(out, 2) /= size(jop)) &
            error stop 'Invalid number of columns'

        ! Perform checks
        call check_ops(self, iop)
        call check_ops(self, jop)
        call self%hybr_matrix(iop, jop, out)
    end subroutine

    subroutine check_ops(self, op)
        class(HYBRIDIZATION), intent(in) :: self
        class(TBaseOper), intent(in) :: op(:)

        integer :: i

        do i = 1, size(op)
            if (op(i)%orb < 1 .or. op(i)%orb > self%norb) &
                stop 'Invalid orbital index'
            if (op(i)%sp < 1 .or. op(i)%sp > 2) &
                stop 'Invalid spin index'
            if (.not. (op(i)%tau >= 0 .and. op(i)%tau <= self%beta)) &
                stop 'Invalid tau'
        end do
    end

    function p_get_hybr_value(self, i, si, j, sj, tau) result(ftau)
        class(HYBRIDIZATION), intent(in) :: self
        integer, intent(in) :: i, si, j, sj
        REAL_TYPE, intent(in) :: tau
        VALUE_TYPE :: ftau, out(1,1)
        type(TBaseOper) :: iop(1), jop(1)

        if (tau >= 0) then
            iop(1) = TBaseOper(orb=i, sp=si, tau=tau)
            jop(1) = TBaseOper(orb=j, sp=sj, tau=0)
        else
            iop(1) = TBaseOper(orb=i, sp=si, tau=0)
            jop(1) = TBaseOper(orb=j, sp=sj, tau=-tau)
        endif
        call copy_hybr_matrix(self, iop, jop, out)
        ftau = out(1, 1)
    end

    function p_index_structure(self) result(lookup)
        class(HYBRIDIZATION), intent(in) :: self
        integer, allocatable :: lookup(:,:,:,:)

        allocate(lookup, source=self%lookup)
    end

    pure function p_get_nbands(self) result(nbands)
        class(HYBRIDIZATION), intent(in) :: self
        integer :: nbands

        nbands = self%norb
    end

    pure function p_get_beta(self) result(beta)
        class(HYBRIDIZATION), intent(in) :: self
        REAL_TYPE :: beta

        beta = self%beta
    end

    pure function p_is_fermionic(self) result(test)
        class(HYBRIDIZATION), intent(in) :: self
        logical :: test

        test = self%fermion
    end

    ! =========================================================================
    ! LINEAR HYBRIDIZATION

    subroutine p_linear_init_bandspin(self, fmat, beta, fermion)
        type(hybr_linear), intent(out) :: self
        VALUE_TYPE, intent(in) :: fmat(:, :, :, :, :)
        REAL_TYPE, intent(in) :: beta
        logical, intent(in) :: fermion

        integer :: norb, ntau, ninv, iinv
        integer, allocatable :: invlookup(:,:)

        norb = size(fmat, 1)
        ntau = size(fmat, 5) - 1
        if (size(fmat, 3) /= norb) &
            error stop 'Number of orbitals are inconsistent'
        if (size(fmat, 2) /= 2 .or. size(fmat, 4) /= 2) &
            error stop 'Expected spin dimensions'
        if (.not. (beta >= 0 .and. beta <= HUGE(beta))) &
            error stop 'Invalid inverse temperature'
        if (ntau < 1) &
            error stop 'Must have at least 2 points in time'

        call init_base_hybridization(self, norb, beta, fermion)
        self%invdtau = ntau / beta

        ! Only store those distinct and non-zero components.  This means we
        ! have to store significantly less stuff, easing cache pressure.  The
        ! lookup also helps finding block indices.
        call find_lookup(fmat, 8 * epsilon(1.0d0), self%lookup, invlookup)
        ninv = size(invlookup, 2)

        allocate(self%mat(0:ntau, ninv))
        do iinv = 1, ninv
            self%mat(:, iinv) = fmat(invlookup(1,iinv), invlookup(2,iinv), &
                                     invlookup(3,iinv), invlookup(4,iinv), :)
        enddo
    end subroutine

    subroutine p_linear_init_bandspin_alloc(self, fmat, beta, fermion)
        class(HYBRIDIZATION), allocatable, intent(out) :: self
        VALUE_TYPE, intent(in) :: fmat(:, :, :, :, :)
        REAL_TYPE, intent(in) :: beta
        logical, intent(in) :: fermion
        type(hybr_linear), allocatable :: actual_self

        allocate(actual_self)
        call p_linear_init_bandspin(actual_self, fmat, beta, fermion)
        call move_alloc(actual_self, self)
    end subroutine

    subroutine hybr_linear_matrix(self, iop, jop, out)
        class(hybr_linear), intent(in) :: self
        class(TBaseOper), intent(in) :: iop(:), jop(:)
        VALUE_TYPE, intent(out) :: out(:, :)

        integer :: i, j, iinv
        REAL_TYPE :: dtau

        do j = 1, size(jop)
            do i = 1, size(iop)
                iinv = self%lookup(iop(i)%orb, iop(i)%sp, jop(j)%orb, jop(j)%sp)
                if (iinv == 0) then
                    out(i, j) = 0
                else
                    dtau = iop(i)%tau - jop(j)%tau
                    out(i, j) = hybr_linear_raw(self, iinv, dtau)
                endif
            end do
        end do
    end

    function hybr_linear_raw(self, iinv, tau) result(ftau)
        type(hybr_linear), intent(in) :: self
        integer, intent(in) :: iinv
        REAL_TYPE, intent(in) :: tau

        VALUE_TYPE :: ftau
        REAL_TYPE :: tausc, space
        integer :: itau

        if (tau >= 0) then
            tausc = tau
        else
            tausc = tau + self%beta
        endif
        tausc = tausc * self%invdtau
        itau = int(tausc)
        space = tausc - itau
        ftau = self%mat(itau, iinv)
        if (itau /= ubound(self%mat, 1)) then
            ftau = ftau + space * (self%mat(itau+1, iinv) - ftau)
        end if
        if (tau < 0 .and. self%fermion) &
            ftau = -ftau

        !write (*,"('F(', I3, ',', I3, ',', F12.6, ') = ', 2E13.6)") &
        !    i, j, tau, ftau
    end

    ! =========================================================================
    ! SUM OF POLES

    subroutine p_init_log_hybridization(self, vmat, eps, beta, fermion)
        class(HYBRIDIZATION), allocatable, intent(out) :: self
        VALUE_TYPE, intent(in) :: vmat(:, :, :, :, :)
        REAL_TYPE, intent(in) :: eps(:), beta
        logical, intent(in) :: fermion
        type(hybr_log), allocatable :: actual_self

        allocate(actual_self)
        call p_init_log_hybridization_i(actual_self, vmat, eps, beta, fermion)
        call move_alloc(actual_self, self)
    end subroutine


    subroutine p_init_log_hybridization_i(self, vmat, eps, beta, fermion)
        type(hybr_log), intent(out) :: self
        VALUE_TYPE, intent(in) :: vmat(:, :, :, :, :)
        REAL_TYPE, intent(in) :: eps(:), beta
        logical, intent(in) :: fermion

        integer :: norb, npole, ninv, iinv, ipole
        integer, allocatable :: invlookup(:,:)
        real(c_double) :: w_max
        VALUE_TYPE :: v_curr

        norb = size(vmat, 1)
        npole = size(vmat, 5)
        if (size(vmat, 3) /= norb) &
            error stop 'Number of orbitals are inconsistent'
        if (size(vmat, 2) /= 2 .or. size(vmat, 4) /= 2) &
            error stop 'Expected spin dimensions'
        if (size(eps) /= npole) &
            error stop 'Inconsistent number of poles'

        call init_base_hybridization(self, norb, beta, fermion)

        w_max = maxval(abs(eps))
        call init_log_grid(self%grid, LOG_DISCR_FLOAT_3, beta / 2, w_max)

        ! Only store those distinct and non-zero components.  This means we
        ! have to store significantly less stuff, easing cache pressure.  The
        ! lookup also helps finding block indices.
        call find_lookup(vmat, 8 * epsilon(1.0d0), self%lookup, invlookup)
        ninv = size(invlookup, 2)

        allocate(self%mat(0:3, self%grid%n_min:self%grid%n_max, 0:1, ninv))
        self%mat(:, :, :, :) = 0
        do iinv = 1, ninv
            do ipole = 1, npole
                v_curr = vmat(invlookup(1,iinv), invlookup(2,iinv), &
                              invlookup(3,iinv), invlookup(4,iinv), ipole)
                call add_propagator( &
                        self%mat(:, :, :, iinv), v_curr, eps(ipole), &
                        beta, self%grid)
            enddo
        enddo
    end subroutine

    subroutine add_propagator(acc, c, w, beta, grid)
        type(TLogGrid), intent(in) :: grid
        VALUE_TYPE, intent(inout) :: acc(0:grid%spec%k, grid%n_min:grid%n_max, 0:1)
        VALUE_TYPE, intent(in) :: c
        REAL_TYPE, intent(in) :: w, beta

        REAL_TYPE :: p0(0:grid%spec%k)
        VALUE_TYPE :: cpos, cneg
        integer :: n

        call exp_taylor(-w, p0)
        cpos = -c / (1 + exp(-beta * w))
        do n = grid%n_min, grid%n_max
            acc(:, n, 0) = acc(:, n, 0) + cpos * exp(-w * grid%tn(n)) * p0
        enddo

        call exp_taylor(w, p0)
        cneg = c / (1 + exp(beta * w))
        do n = grid%n_min, grid%n_max
            acc(:, n, 1) = acc(:, n, 1) + cneg * exp(w * grid%tn(n)) * p0
        enddo
    end subroutine

    pure subroutine exp_taylor(a, p)
        REAL_TYPE, intent(in) :: a
        REAL_TYPE, intent(out)  :: p(0:)
        integer :: i

        do i = 0, ubound(p, 1)
            if (i == 0) then
                p(i) = 1
            else if (i == 1) then
                p(i) = a
            else
                p(i) = a / i * p(i - 1)
            endif
        enddo
    end subroutine

    subroutine hybr_log_matrix(self, iop, jop, out)
        class(hybr_log), intent(in) :: self
        class(TBaseOper), intent(in) :: iop(:), jop(:)
        VALUE_TYPE, intent(out) :: out(:, :)

        integer :: i, j, iinv
        REAL_TYPE :: dtau

        do j = 1, size(jop)
            do i = 1, size(iop)
                iinv = self%lookup(iop(i)%orb, iop(i)%sp, jop(j)%orb, jop(j)%sp)
                if (iinv == 0) then
                    out(i, j) = 0
                else
                    dtau = iop(i)%tau - jop(j)%tau
                    out(i, j) = hybr_log_raw(self, iinv, dtau)
                endif
            end do
        end do
    end

    function hybr_log_raw(self, iinv, tau) result(ftau)
        type(hybr_log), intent(in) :: self
        integer, intent(in) :: iinv
        REAL_TYPE, value :: tau
        VALUE_TYPE :: ftau

        REAL_TYPE :: tau_abs, dtau
        logical :: flip
        integer :: n

        tau_abs = abs(tau)
        flip = 2 * tau_abs > self%beta
        if (flip) then
            tau = sign(tau_abs - self%beta, -tau)
            tau_abs = abs(tau)
        endif

        call segment_time(self%grid, tau_abs, n, dtau)
        ftau = eval_poly3(self%mat(:, n, signbit(tau), iinv), dtau)

        if (flip .and. self%fermion) &
            ftau = -ftau
    end

    pure function eval_poly3(p, x) result(r)
        VALUE_TYPE, intent(in) :: p(0:3)
        REAL_TYPE, intent(in) :: x
        REAL_TYPE :: xsq
        VALUE_TYPE :: y(0:1), r

        y(0) = p(1) * x + p(0)
        y(1) = p(3) * x + p(2)
        xsq = x * x
        r = y(1) * xsq + y(0)
    end function

    pure function signbit(x) result(b)
        real(c_double), value :: x
        integer(c_int64_t) :: b

        ! Relies on IEEE double
        b = transfer(x, b)
        b = ishft(b, -63)
    end

    ! =========================================================================
    ! HELPER FUNCTIONS

    pure elemental function p_propagator_tau( &
                        coeff, tau, eps, beta, fermion) result(gtau)
        VALUE_TYPE, value :: coeff
        REAL_TYPE, value :: tau, eps, beta
        logical, value :: fermion

        VALUE_TYPE :: gtau
        REAL_TYPE :: tau_abs
        logical :: flip
        REAL_TYPE, parameter :: ONE = 1
        REAL_TYPE, parameter :: PLUS_ONE_THRESHOLD = -log(4/epsilon(ONE))

        ! Check for invalid values
        tau_abs = abs(tau)
        if (.not. (tau_abs <= beta)) then
            gtau = nan(1.0d0)
            return
        endif

        ! Map back to [-beta/2, beta/2]
        flip = 2 * tau_abs > beta
        if (flip) then
            tau = sign(tau_abs - beta, -tau)
            tau_abs = abs(tau)
        endif

        ! For negative values, we can flip both energies and times
        if (signbit(tau) /= 0) then
            tau = -tau
            eps = -eps
            flip = .not. flip
        endif

        if (beta * eps > PLUS_ONE_THRESHOLD) then
            ! Use the definition
            gtau = exp(-tau * eps) / (1 + exp(-beta * eps))
        else
            ! Use the fact that exp arg is so negative that +1 does not matter
            gtau = exp((beta - tau) * eps)
        endif
        gtau = coeff * gtau

        if (fermion .and. .not. flip) &
            gtau = -gtau
    end

    subroutine p_discretize_propagator(vmat, eps, beta, fermion, fmat)
        VALUE_TYPE, intent(in) :: vmat(:, :, :, :, :)
        REAL_TYPE, intent(in) :: eps(:), beta
        logical, intent(in) :: fermion
        VALUE_TYPE, intent(out) :: fmat(:, :, :, :, 0:)

        integer :: i, s, j, t, itau, ntau
        REAL_TYPE :: tau
        VALUE_TYPE :: fcurr

        ntau = size(fmat, 5) - 1
        do itau = 0, ntau
            tau = 1.0d0 * itau / ntau * beta
            do t = 1, size(fmat, 4)
            do j = 1, size(fmat, 3)
            do s = 1, size(fmat, 2)
            do i = 1, size(fmat, 1)
                fcurr = sum(propagator_tau(vmat(i, s, j, t, :), tau, &
                                           eps(:), beta, fermion))
                fmat(i, s, j, t, itau) = fcurr
            enddo
            enddo
            enddo
            enddo
        enddo
    end

    subroutine find_lookup(fmat, rtol, lookup, invlookup)
        VALUE_TYPE, intent(in) :: fmat(:,:,:,:,:)
        REAL_TYPE, intent(in) :: rtol
        integer, intent(out), allocatable :: lookup(:,:,:,:)
        integer, intent(out), allocatable :: invlookup(:,:)

        integer :: inv(4, size(fmat,1)*size(fmat,2)*size(fmat,3)*size(fmat,4))
        integer :: i, j, s, t, ninv, iinv, ip, jp, sp, tp
        REAL_TYPE :: tol

        if (rtol < 0) &
            error stop 'Invalid relative tolerance'

        allocate(lookup(size(fmat,1), size(fmat,2), size(fmat,3), size(fmat,4)))
        lookup(:,:,:,:) = 0

        tol = rtol * maxval(abs(fmat))
        if (.not. (tol <= HUGE(tol))) &
            error stop 'Invalid values in the hybridization function'

        ninv = 0
        do t = 1, size(fmat, 4)
        do j = 1, size(fmat, 3)
            do s = 1, size(fmat, 2)
            comp: do i = 1, size(fmat, 1)
                ! First check if the component is actually signifacnt
                if (all(abs(fmat(i,s,j,t,:)) <= tol)) &
                    cycle comp

                ! Check against other components whether its a duplicate
                do iinv = 1, ninv
                    ip = inv(1, iinv)
                    sp = inv(2, iinv)
                    jp = inv(3, iinv)
                    tp = inv(4, iinv)
                    if (all(abs(fmat(i,s,j,t,:) - fmat(ip,sp,jp,tp,:)) <= tol)) then
                        lookup(i, s, j, t) = iinv
                        cycle comp
                    endif
                enddo

                ninv = ninv + 1
                inv(1, ninv) = i
                inv(2, ninv) = s
                inv(3, ninv) = j
                inv(4, ninv) = t
                lookup(i, s, j, t) = ninv
            enddo comp
            enddo
        enddo
        enddo

        allocate(invlookup, source=inv(:, :ninv))
    end subroutine

end module
