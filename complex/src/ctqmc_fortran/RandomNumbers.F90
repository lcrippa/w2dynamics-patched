module MRandomNumbers
    use iso_c_binding, only: c_int32_t, c_int64_t, c_double
    implicit none
    private

    !> Pseudorandom number generator (PRNG)
    !!
    !! Abstract buffer-backed random number generator, which yields 32-bit
    !! integers which should be approximately equidistributed over the whole
    !! range.
    type, abstract, public :: RandomEngine
        private
        integer(c_int32_t), allocatable :: buffer(:)
        integer :: cur = 0
    contains
        private
        procedure(buffer_size), deferred :: buffer_size
        procedure(seed), deferred :: seed
        procedure(refill), deferred :: refill
    end type

    interface
        !> Return buffer size of the RNG raw buffer
        integer function buffer_size(self)
            import :: RandomEngine
            class(RandomEngine), intent(in) :: self
        end function

        !> Seed the random number generator
        subroutine seed(self, newseed)
            import :: RandomEngine, c_int32_t
            class(RandomEngine), intent(inout) :: self
            integer(c_int32_t), intent(in) :: newseed
        end subroutine

        !> Refill the buffer
        subroutine refill(self)
            import :: RandomEngine
            class(RandomEngine), intent(inout) :: self
        end subroutine
    end interface

    type, extends(RandomEngine) :: mersenne_twister
        integer(c_int32_t) :: x(0:623)
    contains
        procedure, non_overridable :: buffer_size => mt_buffer_size
        procedure, non_overridable :: seed => mt_seed
        procedure, non_overridable :: refill => mt_refill
    end type

    type, extends(RandomEngine) :: xorshift_star_64
        integer(c_int64_t) :: s
    contains
        procedure, non_overridable :: buffer_size => xor64_buffer_size
        procedure, non_overridable :: seed => xor64_seed
        procedure, non_overridable :: refill => xor64_refill
    end type

    public :: seed_rng
    public :: rng_get_raw

    interface random_integer
       module procedure random_integer_base, random_integer_interval
    end interface
    public random_integer

    interface random_real
       module procedure random_real_01, random_real_scaled, random_real_interval
    end interface
    public random_real

    interface shuffle
        module procedure shuffle_int
    end interface
    public shuffle

    public :: random_logical
    public :: random_permutation, random_swaps, random_sample
    public :: init_mersenne_twister
    public :: init_xorshift_star

contains

    !> Seed the random number engine
    !!
    !! Seed the random number generator with a number, in effect selecting the
    !! pseudorandom sequence to use.
    !!
    !! NOTE
    !!   Depending on the RNG and in particular its state size, different
    !!   seeds may yield overlapping sequences of random numbers.
    subroutine seed_rng(self, newseed)
        class(RandomEngine), intent(inout) :: self
        integer(c_int32_t), intent(in), optional :: newseed
        integer(c_int32_t), parameter :: default_seed = 5489

        if (present(newseed)) then
            call self%seed(newseed)
        else
            call self%seed(default_seed)
        endif

        ! Invalidate the cache
        self%cur = ubound(self%buffer, 1) + 1
    end subroutine

    !> Get new pseudorandom 32-bit integer
    function rng_get_raw(self) result(y)
        class(RandomEngine), intent(inout) :: self
        integer(c_int32_t) :: y

        if (self%cur >= size(self%buffer)) &
            call self%refill()

        y = self%buffer(self%cur)
        self%cur = self%cur + 1
    end function

    ! -------------------------------------------------------------------------
    ! DISTRIBUTIONS

    function random_integer_base(self, s) result(r)
        class(RandomEngine), intent(inout) :: self
        integer(c_int32_t), intent(in) :: s
        integer(c_int32_t) :: r, l, x, t
        integer(c_int64_t) :: m

        if (s <= 0) &
            error stop 'Invalid range of numbers'

        ! see: D. Lemire, arXiv:1805.10941v4, Algorithm 5
        x = rng_get_raw(self)
        m = widen_unsigned(x) * s
        l = int(m, c_int32_t)

        if (l >= 0 .and. l < s) then
            t = mod(huge(t) - s + 1, s)

            do while (l >= 0 .and. l < t)
                x = rng_get_raw(self)
                m = widen_unsigned(x) * s
                l = int(m, c_int32_t)
            enddo
        endif

        r = int(ishft(m, -32), c_int32_t)
    end function

    function random_integer_interval(self, low, high) result(r)
        class(RandomEngine), intent(inout) :: self
        integer(c_int32_t), intent(in) :: low, high
        integer(c_int32_t) :: r

        r = low + random_integer_base(self, high - low + 1)
    end function

    function random_logical(self) result(r)
        class(RandomEngine), intent(inout) :: self
        logical :: r

        r = random_integer(self, 1) == 1
    end

    function random_real_base(self) result(r)
        class(RandomEngine), intent(inout) :: self
        real(c_double) :: r

        ! Multiplication with power of 2 should be faster than ldexp
        real(c_double), parameter :: radix = scale(1.0d0, 32)
        real(c_double), parameter :: invsc = scale(1.0d0, -64)

        r = real(widen_unsigned(rng_get_raw(self)), c_double)
        r = r + real(widen_unsigned(rng_get_raw(self)), c_double) * radix
        r = r * invsc
    end function

    function random_real_01(self) result(r)
        class(RandomEngine), intent(inout) :: self
        real(c_double) :: r

        ! Occasionally, this may round to one. If that's the case, simply retry
        do
            r = random_real_base(self)
            if (r /= 1) &
                return
        enddo
    end function

    function random_real_scaled(self, rscale) result(r)
        class(RandomEngine), intent(inout) :: self
        real(c_double), intent(in) :: rscale
        real(c_double) :: r

        if (.not. (rscale > 0)) &
            error stop 'scale parameter must be non-negative'

        do
            r = random_real_base(self) * rscale
            if (r < rscale) &
                return
        enddo
    end function

    function random_real_interval(self, rmin, rmax) result(r)
        class(RandomEngine), intent(inout) :: self
        real(c_double), intent(in) :: rmin, rmax
        real(c_double) :: r, rscale

        ! Occasionally, this may round to rmax
        rscale = rmax - rmin
        if (.not. (rscale > 0)) &
            error stop 'interval [rmin, rmax) must be non-empty'

        do
            r = rmin + random_real_base(self) * rscale
            if (r < rmax) &
                return
        enddo
    end function

    subroutine random_swaps(self, n, s)
        class(RandomEngine), intent(inout) :: self
        integer, intent(in) :: n
        integer, intent(out) :: s(:)
        integer :: i

        if (n /= size(s)) &
            error stop 'Invalid size of perm array'

        ! Use Fisher-Yates shuffle
        do i = 1, size(s)-1
            s(i) = random_integer(self, i, size(s))
        end do
        s(size(s)) = size(s)
    end subroutine

    subroutine shuffle_int(self, arr)
        class(RandomEngine), intent(inout) :: self
        integer, intent(inout) :: arr(:)
        integer :: i, j, tmp

        ! Use Fisher-Yates shuffle
        do i = 1, size(arr)-1
            j = random_integer(self, i, size(arr))

            tmp = arr(i)
            arr(i) = arr(j)
            arr(j) = tmp
        enddo

        ! XXX We only need size(arr)-1 random numbers. This is here to
        !     satisfy regression test
        j = random_integer(self, 1)
    end subroutine

    subroutine random_permutation(self, n, perm)
        class(RandomEngine), intent(inout) :: self
        integer, intent(in) :: n
        integer, intent(out) :: perm(:)
        integer :: i

        if (n /= size(perm)) &
            error stop 'Invalid size of perm array'

        do i = 1, n
            perm(i) = i
        enddo
        call shuffle(self, perm)
    end subroutine

    !> Sample k values from the set {1, ..., n} without replacement
    subroutine random_sample(self, n, k, sample)
        class(RandomEngine), intent(inout) :: self
        integer, intent(in) :: n, k
        integer, intent(out) :: sample(:)
        integer :: i, j, t

        if (k /= size(sample)) &
            error stop 'Invalid size of the sample array'
        if (k > n) &
            error stop 'Sample size is greater than set to sample from'

        ! Use Floyd's algorithm, Comm. ACM 30(9), 754 (1987), Algorithm F2
        ! This algorithm runs in O(k^2) time, but k is small in our case
        do i = 1, k
            j = n - k + i
            t = random_integer(self, 1, j)
            if (is_element_of(t, sample(:i-1))) then
                sample(i) = j
            else
                sample(i) = t
            endif
        enddo
    contains
        pure logical function is_element_of(x, set)
            integer, intent(in) :: x, set(:)
            integer :: i

            is_element_of = .false.
            do i = 1, size(set)
                if (x == set(i)) then
                    is_element_of = .true.
                    return
                endif
            enddo
        end
    end subroutine

    ! -------------------------------------------------------------------------
    ! HELPER FUNCTIONS

    pure function widen_unsigned(x) result(w)
        integer(c_int32_t), intent(in) :: x
        integer(c_int64_t) :: w
        integer(c_int64_t), parameter :: lower_32bits = 4294967295_c_int64_t

        ! Interpret the 32-bit integer as an unsigned value and store that
        ! in a 64-bit (signed) integer.
        w = iand(int(x, c_int64_t), lower_32bits)
    end function

    subroutine init_buffer(self)
        class(RandomEngine), intent(inout) :: self
        integer :: n

        n = self%buffer_size()
        allocate(self%buffer(0:n-1))
        self%cur = n
    end subroutine

    ! -------------------------------------------------------------------------
    ! MERSENNE TWISTER ROUTINES

    subroutine init_mersenne_twister(self, seed)
        class(RandomEngine), allocatable, intent(out) :: self
        integer(c_int32_t), intent(in), optional :: seed
        type(mersenne_twister), allocatable :: actual_self

        allocate(actual_self)
        call init_buffer(actual_self)
        call seed_rng(actual_self, seed)
        call move_alloc(actual_self, self)
    end subroutine

    function mt_buffer_size(self) result(n)
        class(mersenne_twister), intent(in) :: self
        integer :: n

        n = 624
    end function

    subroutine mt_seed(self, newseed)
        class(mersenne_twister), intent(inout) :: self
        integer(c_int32_t), intent(in) :: newseed

        call mt_seed_internal(newseed, self%x)
    end subroutine

    subroutine mt_refill(self)
        class(mersenne_twister), intent(inout) :: self
        integer :: i

        call mt_refill_internal(self%x)
        do i = 0, 623
            self%buffer(i) = mt_temper(self%x(i))
        enddo
        self%cur = 0
    end subroutine

    pure subroutine mt_seed_internal(s, x)
        integer, parameter :: n = 624
        integer(c_int32_t), intent(in) :: s
        integer(c_int32_t), intent(out) :: x(0:n-1)

        integer(c_int32_t), parameter :: f = 1812433253
        integer(c_int32_t) :: prev
        integer :: i

        x(0) = s
        do i = 1, n - 1
            prev = ieor(x(i-1), ishft(x(i-1), -30))
            x(i) = f * prev + i
        enddo
    end subroutine

    pure subroutine mt_refill_internal(x)
        integer, parameter :: n = 624
        integer(c_int32_t), intent(inout) :: x(0:n-1)

        integer, parameter :: m = 397
        integer(c_int32_t), parameter :: lower_mask = 2147483647  ! 0x7FFFFFFF
        integer(c_int32_t), parameter :: upper_mask = not(lower_mask)
        integer(c_int32_t), parameter :: a = -1727483681   ! 0x9908b0df

        integer(c_int32_t) :: xcurr, xa
        integer :: i

        do i = 0, n - 1
            xcurr = ior(iand(x(i), upper_mask), &
                        iand(x(mod(i+1, n)), lower_mask))
            xa = ishft(xcurr, -1)
            if (iand(xcurr, 1) /= 0) then
                xa = ieor(xa, a)
            endif
            x(i) = ieor(x(mod(i + m, n)), xa)
        enddo
    end subroutine

    pure function mt_temper(x) result(y)
        integer(c_int32_t), intent(in) :: x
        integer(c_int32_t) :: y

        integer, parameter :: u = 11, s = 7, t = 15, l = 18
        integer(c_int32_t), parameter :: &
                    d = -1, &                ! 0xffffffff
                    b = -1658038656, &       ! 0x9d2c5680
                    c = -272236544           ! 0xefc60000

        y = ieor(x, iand(ishft(x, -u), d))
        y = ieor(y, iand(ishft(y, s), b))
        y = ieor(y, iand(ishft(y, t), c))
        y = ieor(y, ishft(y, -l))
    end function

    ! -------------------------------------------------------------------------
    ! XOR SHIFT STAR 64

    subroutine init_xorshift_star(self, seed)
        class(RandomEngine), allocatable, intent(out) :: self
        integer(c_int32_t), intent(in), optional :: seed
        type(xorshift_star_64), allocatable :: actual_self

        allocate(actual_self)
        call init_buffer(actual_self)
        call seed_rng(actual_self, seed)
        call move_alloc(actual_self, self)
    end subroutine

    function xor64_buffer_size(self) result(n)
        class(xorshift_star_64), intent(in) :: self
        integer :: n

        n = 512
    end function

    subroutine xor64_seed(self, newseed)
        class(xorshift_star_64), intent(inout) :: self
        integer(c_int32_t), intent(in) :: newseed
        integer(c_int64_t) :: s

        ! xor-shift generators have a problem with "zero land", i.e., the
        ! phenomenon that they have to run for a while for the number of 0
        ! and 1 bits in the output to equilibrate if they are stronly
        ! imbalanced (this never happens if the state is seeded to zero).
        ! To combat this, we set the high 32 bits to the bitwise inverse of
        ! the low 32 bits.
        s = newseed
        self%s = ior(ishft(not(s), 32), s)
    end subroutine

    subroutine xor64_refill(self)
        class(xorshift_star_64), intent(inout) :: self

        ! The constant M_32 given in Table IV of [1].
        integer(c_int64_t), parameter :: m = 2685821657736338717_c_int64_t
        integer :: i

        do i = 0, ubound(self%buffer, 1)
            self%s = xor64_advance(self%s)
            self%buffer(i) = int(ishft(self%s * m, -32), c_int32_t)
        enddo
        self%cur = 0
    end subroutine

    pure function xor64_advance(state) result(x)
        integer(c_int64_t), intent(in) :: state
        integer(c_int64_t) :: x
        integer, parameter :: a = 12, b = 25, c = 27

        ! xorshift generator A1(12, 25, 27) of [1]
        ! [1] https://arxiv.org/abs/1402.6246v5
        x = state
        x = ieor(x, ishft(x, -a))
        x = ieor(x, ishft(x, b))
        x = ieor(x, ishft(x, -c))
    end function

end module
