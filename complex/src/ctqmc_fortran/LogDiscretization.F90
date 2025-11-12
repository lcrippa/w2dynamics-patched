module MLogDiscretization
    use iso_c_binding, only: c_double, c_int64_t
    implicit none
    private

    !> Logarithmic discretization with guaranteed accuracy.
    !!
    !! Partitioning of a time axis [0, t_max] into logarithmic intervals with
    !! the segment centres:
    !!
    !!    t(n) = alpha**n,
    !!
    !! Provided that the function is a simple falling exponential exp(-w*t) and
    !! is approximated in each interval in a polynomial of degree k, then the
    !! absolute error of that is bound by some eps.
    type, public :: TLogDiscretization
        real(c_double) :: eps         !< tolerance
        integer :: k = -1             !< order of polynomial in each segment
        real(c_double) :: alpha       !< factor between segment centres
        real(c_double) :: root_of_2   !< 2**(1/root_of_2) == alpha
        real(c_double) :: lambda      !< upper bound of tau0 * wmax
    end type

    !> Logarithmic grid with guaranteed accuracy.
    !!
    !! Partitioning of a time axis [0, t_max] into logarithmic intervals with
    !! the segment centres:
    !!
    !!    t(n) = alpha**n,   n = nmin, ..., n_max,
    !!
    !! Provided that the function is a sum of simple falling exponentials,
    !! exp(-w*t), abs(w) <= w_max, and it is is approximated in each interval
    !! in a polynomial of degree k, then the absolute error of that is bound
    !! by some eps. For the lowest interval n = n_min, the approximation
    !! additionally holds for any 0 <= t <= t(n_min).
    type, public :: TLogGrid
        type(TLogDiscretization) :: spec
        real(c_double) :: w_max, t_max
        real(c_double), allocatable :: tn(:)
        integer :: n_min, n_max
    end type

    !> Discretization for 3-order polynomials maintaining single accuracy
    type(TLogDiscretization), parameter, public :: LOG_DISCR_FLOAT_3 = &
        TLogDiscretization( &
                eps=1.1920928955078125d-7, k=3, &
                alpha=1.0551172497206347d0, root_of_2=12.919341296760553d0, &
                lambda=0.04112731290476539d0)

    public :: init_log_discretization
    public :: init_log_grid
    public :: segment_time
    public :: log2_sloppy

contains

    !> New log discretization for given tolerance `eps` and poly order `k`
    subroutine init_log_discretization(self, eps, k)
        type(TLogDiscretization), intent(out) :: self
        real(c_double), intent(in) :: eps
        integer, intent(in) :: k

        if (k < 0) &
            error stop 'Invalid value of polynomial order'
        if (.not. (eps > 0)) &
            error stop 'tolerance must be positive'

        self%eps = eps
        self%k = k
        self%alpha = find_alpha(eps, k)
        self%root_of_2 = log(2.0d0) / log(self%alpha)
        self%lambda = find_lambda(eps, k)
    end subroutine

    !> New logarithmic grid for given time and decay cutoffs
    subroutine init_log_grid(self, spec, t_max, w_max)
        type(TLogGrid), intent(out) :: self
        type(TLogDiscretization), intent(in) :: spec
        real(c_double), intent(in) :: t_max, w_max
        integer :: n, n_min, n_max

        if (spec%k < 0) &
            error stop 'Discretization not initialized'
        if (.not. (t_max >= 0)) &
            error stop 't_max must be non-negative'
        if (.not. (w_max >= 0)) &
            error stop 'w_max must be non-negative'

        self%spec = spec
        call find_grid_bounds( &
                spec%alpha, spec%lambda, t_max, w_max, n_min, n_max)

        ! We need a plus-one here, since there are some inaccuracies in log2,
        ! which may catapult us into a further bin
        n_max = n_max + 1
        allocate(self%tn(n_min:n_max))
        do n = n_min, n_max
            self%tn(n) = spec%alpha**n
        enddo

        self%t_max = self%tn(n_max-1) * sqrt(spec%alpha)
        self%w_max = spec%lambda / self%tn(n_min)
        self%n_min = n_min
        self%n_max = n_max
    end subroutine

    !> Split time `t` into segment index `n` and time `dt` within segment
    subroutine segment_time(self, t, n, dt)
        type(TLogGrid), intent(in) :: self
        real(c_double), intent(in) :: t
        integer, intent(out) :: n
        real(c_double), intent(out) :: dt

        real(c_double) :: nr

        nr = self%spec%root_of_2 * log2_sloppy(t)
        n = max(floor(nr + 0.5), lbound(self%tn, 1))
        dt = t - self%tn(n)
    end subroutine

    pure function find_lambda(eps, k) result(lambda)
        real(c_double), intent(in) :: eps
        integer, intent(in) :: k
        real(c_double) :: lambda

        ! One has that:
        !   eps = lambda**(k+1) / gamma(k+2)
        !   log(eps) = (k + 1) * log(lambda) - log_gamma(k + 2)
        !
        lambda = exp((log(eps) + log_gamma(k+2.0d0)) / (k + 1))
    end function

    pure function find_alpha(eps, k) result(alpha)
        real(c_double), intent(in) :: eps
        integer, intent(in) :: k
        real(c_double) :: alpha, c

        ! One has that
        ! eps = poisson_gamma(k + 1)
        !       * ((sqrt(alpha) - 1) / (2 - sqrt(alpha)))**(k + 1)
        ! from which it follows:
        c = (eps / poisson_gamma(k + 1))**(1 / (k + 1.0d0))
        alpha = ((1 + 2 * c) / (1 + c))**2
    end function

    !> Implement the following expression:
    !!
    !!    γ(k) := k^k * exp(-k) / factorial(k) ≈ 1/sqrt(2π * k)
    !!
    pure function poisson_gamma(k) result(gamma)
        integer, intent(in) :: k
        real(c_double) :: gamma

        if (k == 0) then
            gamma = 1
        else
            gamma = exp(k * log(0.0d0 + k) - k - log_gamma(k + 1.0d0))
        endif
    end function

    pure subroutine find_grid_bounds(alpha, lambda, t_max, w_max, n_min, n_max)
        real(c_double), intent(in) :: alpha, lambda, t_max, w_max
        integer, intent(out) :: n_min, n_max
        real(c_double) :: log_alpha, t_min

        if (t_max == 0 .or. w_max == 0) then
            n_min = 0
            n_max = 0
        else
            t_min = lambda / w_max
            log_alpha = log(alpha)
            n_min = floor(log(t_min) / log_alpha)
            n_max = nint(log(t_max) / log_alpha)
        endif
    end subroutine

    !> Logarithm accurate within 0.005 and without bounds checks
    pure elemental function log2_sloppy(x) result(r)
        real(c_double), value :: x
        real(c_double) :: s, r
        integer(c_int64_t) :: e

        ! Assumes IEEE double floating point
        integer(c_int64_t), parameter :: e_bias = 1023

        ! Minimax polynomial to log(2) in [1,2), accurate to 0.005
        real(c_double), parameter :: b0 = -1.674873d0 - e_bias
        real(c_double), parameter :: b1 = +2.024658d0
        real(c_double), parameter :: b2 = -0.3448453d0

        e = exponent_part(x)
        s = significand_part(x)
        r = b2 * s + b1
        r = r * s + (b0 + e)
    end function

    pure function exponent_part(x) result(e)
        real(c_double), value :: x
        integer(c_int64_t) :: x_bits, e

        ! Assumes IEEE double floating point
        x_bits = transfer(x, x_bits)
        x_bits = ishft(x_bits, 1)
        e = ishft(x_bits, -53)
    end function

    pure function significand_part(x) result(s)
        real(c_double), value :: x
        real(c_double) :: s
        integer(c_int64_t) :: x_bits, s_bits

        ! Assumes IEEE double floating point
        integer(c_int64_t), parameter :: m_mask = 4503599627370495_c_int64_t
        integer(c_int64_t), parameter :: s_mask = 4607182418800017408_c_int64_t

        x_bits = transfer(x, x_bits)
        s_bits = ior(iand(x_bits, m_mask), s_mask)
        s = transfer(s_bits, s)
    end function

end module
