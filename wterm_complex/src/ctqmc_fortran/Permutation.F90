module MPermutation
    implicit none
    private

    !>      perm_parity(p [, work])
    !!
    !! Returns .true. if the permutation in p is odd, and .false. if it is
    !! even.  Optionally, you can pass a logical work array of the same size
    !! as p.
    interface perm_parity
        module procedure perm_parity_alloc, perm_parity_work
    end interface
    public perm_parity

    !>      perm_to_swaps(p, s, [, work])
    interface perm_to_swaps
       module procedure perm_to_swaps_alloc, perm_to_swaps_work
    end interface
    public perm_to_swaps

    !>    call swap(x, y)
    !!
    !! Swaps two elements x and y of the same type.
    interface swap
        module procedure p_swap_int
    end interface
    public swap

    public is_permutation
    public get_identity_perm, identity_perm
    public get_cyclic_perm, cyclic_perm
    public get_group_perm, group_perm
    public get_grow_perm, grow_perm
    public get_grow_perm_front, grow_perm_front
    public get_shrink_perm, shrink_perm
    public get_move_perm, move_perm
    public get_inv_perm, inv_perm
    public get_perm_remove, perm_remove
    public swaps_to_perm
    public swap_parity
contains
    !> Trivial permutation of n elements
    function get_identity_perm(n) result(perm)
        integer, intent(in) :: n
        integer, allocatable :: perm(:)

        allocate(perm(n))
        call identity_perm(n, perm)
    end

    !> Trivial permutation of n elements
    subroutine identity_perm(n, perm)
        integer, intent(in) :: n
        integer, intent(out) :: perm(:)

        integer :: i

        if (size(perm) /= n) &
            error stop 'Invalid size'
        do i = 1, n
            perm(i) = i
        enddo
    end

    !> Cyclic permutation by k places.
    !!
    !! This permutation moves elements by k places towards the front, i.e.,
    !! for 0 <= k < n, the element x(1) gets the value of x(k+1) and so forth,
    !! while wrapping around.
    function get_cyclic_perm(n, k) result(perm)
        integer, intent(in) :: n, k
        integer, allocatable :: perm(:)

        allocate(perm(n))
        call cyclic_perm(n, k, perm)
    end

    !! Cyclic permutation by k places (in-place version of cyclic_perm)
    subroutine cyclic_perm(n, k, perm)
        integer, intent(in) :: n, k
        integer, intent(out) :: perm(:)

        integer :: i, p

        if (size(perm) /= n) &
            error stop 'Invalid size'

        p = modulo(k, n)
        do i = 1, n
            p = p + 1
            perm(i) = p
            if (p == n) &
                p = 0
        enddo
    end

    !> Permutation which groups items by group index.
    !!
    !! Expects an array `group`, where `group(i)` denotes the group
    !! index of the `i`-th item (a positive number). Return a permutation,
    !! which groups item by group index while otherwise preserving the order.
    function get_group_perm(group) result(perm)
        integer, intent(in) :: group(:)
        integer, allocatable :: perm(:), gstop(:)

        allocate(perm(size(group)), gstop(maxval(group)))
        call group_perm(group, perm, gstop)
    end

    !> Permutation which groups items by group index.
    !!
    !! Expects an array `group` of size `n`, where `group(i)` denotes the group
    !! index of the `i`-th item (a positive number). Fills `perm` with a
    !! permutation such that `group(perm)` is ordered by group number.
    !! `gstop(k)` then contains the last index of `group(perm)` which contains
    !! the group `k`.
    subroutine group_perm(group, perm, gstop)
        integer, intent(in) :: group(:)
        integer, intent(out) :: perm(:), gstop(:)
        integer :: i, curr

        if (size(group) /= size(perm)) &
            error stop 'Size mismatch'

        ! First sweep: determine group counts
        call histogram(group, gstop)

        ! Second sweep: Determine offsets
        call sizes_to_offsets(gstop)

        ! Third sweep
        do i = 1, size(perm)
            curr = group(i)
            gstop(curr) = gstop(curr) + 1
            perm(gstop(curr)) = i
        enddo
    end subroutine

    subroutine histogram(a, cnt)
        integer, intent(in) :: a(:)
        integer, intent(out) :: cnt(:)
        integer :: i, curr

        cnt(:) = 0
        do i = 1, size(a)
            curr = a(i)
            if (curr < 1 .or. curr > size(cnt)) &
                error stop 'Invalid group number'
            cnt(curr) = cnt(curr) + 1
        enddo
    end subroutine

    subroutine sizes_to_offsets(a)
        integer, intent(inout) :: a(:)
        integer :: curr, i, offset

        offset = 0
        do i = 1, size(a)
            curr = a(i)
            a(i) = offset
            offset = offset + curr
        enddo
    end subroutine

    !> Permutation for moving trailing elements to other positions.
    !!
    !! Into a set of n elements (x(1), ..., x(n)), we would like to insert k
    !! new (y(1), ..., y(k)) elements such that the come to lie at position
    !! (new(1), ..., new(k)).  A simple way to do that is to first insert
    !! these elements at the end and then apply a permutation such that
    !!
    !!     xupd(i) = x(perm(i))
    !!
    !! This function returns the permutation to achieve just that.
    function get_grow_perm(n, new) result(perm)
        integer, intent(in) :: n, new(:)
        integer, allocatable :: perm(:)

        allocate(perm(n + size(new)))
        call grow_perm(n, new, perm)
    end

    !> Get permutation for array growth (inplace version of get_grow_perm)
    subroutine grow_perm(n, new, perm)
        integer, intent(in) :: n, new(:)
        integer, intent(out) :: perm(:)

        integer :: ifrom, ito
        if (size(perm) /= n + size(new)) &
            error stop 'Perm array is of wrong size'

        call create_mask(new, perm)

        ifrom = 0
        do ito = 1, n + size(new)
            if (perm(ito) == 0) then
                ifrom = ifrom + 1
                perm(ito) = ifrom
            endif
        enddo
        do ito = 1, size(new)
            ifrom = ifrom + 1
            perm(new(ito)) = ifrom
        enddo
    end

    !> Permutation for moving trailing elements to other positions.
    !!
    !! Into a set of n elements (x(1), ..., x(n)), we would like to insert k
    !! new (y(1), ..., y(k)) elements such that the come to lie at position
    !! (new(1), ..., new(k)).  A simple way to do that is to first insert
    !! these elements at the FRONT and then apply a permutation such that
    !!
    !!     xupd(i) = x(perm(i))
    !!
    !! This function returns the permutation to achieve just that.
    function get_grow_perm_front(n, new) result(perm)
        integer, intent(in) :: n, new(:)
        integer, allocatable :: perm(:)

        allocate(perm(n + size(new)))
        call grow_perm_front(n, new, perm)
    end

    !> Get permutation for array growth (inplace version of get_grow_perm)
    subroutine grow_perm_front(n, new, perm)
        integer, intent(in) :: n, new(:)
        integer, intent(out) :: perm(:)

        integer :: ifrom, ito
        if (size(perm) /= n + size(new)) &
            error stop 'Perm array is of wrong size'

        call create_mask(new, perm)

        ifrom = size(new)
        do ito = 1, n + size(new)
            if (perm(ito) == 0) then
                ifrom = ifrom + 1
                perm(ito) = ifrom
            endif
        enddo
        do ifrom = 1, size(new)
            perm(new(ifrom)) = ifrom
        enddo
    end

    !> Permutation for moving elements to trailing positions.
    !!
    !! Into a set of n elements (x(1), ..., x(n)), we would like to remove k
    !! elements at position (rem(1), ..., rem(k)).  A simple way to do that is
    !! to first permute these elements to the back:
    !!
    !!     xupd(i) = x(perm(i))
    !!
    !! and then shrink the vector.  This function returns the permutation to
    !! achieve just that.  This is the inverse of get_grow_perm
    function get_shrink_perm(n, rem) result(perm)
        integer, intent(in) :: n, rem(:)
        integer, allocatable :: perm(:)

        allocate(perm(n))
        call shrink_perm(n, rem, perm)
    end

    !> Get permutation for array shrinkage (inplace version of get_shrink_perm)
    subroutine shrink_perm(n, rem, perm)
        integer, intent(in) :: n, rem(:)
        integer, intent(out) :: perm(:)

        integer :: ifrom, ito
        if (size(perm) /= n) &
            error stop 'Perm array is of wrong size'

        call create_mask(rem, perm)

        ito = 0
        do ifrom = 1, n
            if (perm(ifrom) == 0) then
                ito = ito + 1
                perm(ito) = ifrom
            end if
        end do
        do ifrom = 1, size(rem)
            ito = ito + 1
            perm(ito) = rem(ifrom)
        end do
    end

    subroutine create_mask(idx, perm)
        integer, intent(in) :: idx(:)
        integer, intent(out) :: perm(:)

        integer :: n, k, ik

        n = size(perm)
        k = size(idx)
        perm(:) = 0
        do ik = 1, k
            if (idx(ik) < 1 .or. idx(ik) > n) &
                error stop 'Invalid index'
            if (perm(idx(ik)) /= 0) &
                error stop 'Duplicate index'
            perm(idx(ik)) = 1
        enddo
    end

    !> Get permutation for moving elements
    !!
    !! Into a set of n elements (x(1), ..., x(n)), we would like to move
    !! elements at indices (from(1), ..., from(k)) such that they comet to lie
    !! at indices (to(1), ..., to(k)), respectively.  This function returns
    !! a permutation which achieves just that.
    function get_move_perm(n, from, to) result(perm)
        integer, intent(in) :: n, from(:), to(:)
        integer, allocatable :: perm(:)

        allocate(perm(n))
        call move_perm(n, from, to, perm)
    end

    !> Get permutation for moving (inplace version of get_move_perm)
    subroutine move_perm(n, from, to, perm)
        integer, intent(in) :: n, from(:), to(:)
        integer, intent(out) :: perm(:)

        ! FIXME: do not allocate stuff ...
        integer :: work(n,2)

        integer :: i, k

        k = size(from)
        call shrink_perm(n, from, work(:,1))
        call grow_perm(n - k, to, work(:,2))
        do i = 1, n
            perm(i) = work(work(i, 2), 1)
        enddo
    end

    !> Given a permutation p remove the values idx
    subroutine perm_remove(p, idx, q)
        integer, intent(in) :: p(:), idx(:)
        integer, intent(out) :: q(:)
        integer :: i, j

        if (size(p) /= size(q) + size(idx)) &
            error stop 'Invalid sizes'

        j = 0
        do i = 1, size(p)
            if (any(idx == i)) &
                cycle

            j = j + 1
            q(j) = p(i) - count(p(idx) < p(i))
        enddo
    end subroutine

    function get_perm_remove(p, idx) result(q)
        integer, intent(in) :: p(:), idx(:)
        integer, allocatable :: q(:)

        allocate(q(size(p) - size(idx)))
        call perm_remove(p, idx, q)
    end

    !> Return inverse permutation.
    !!
    !! Given a permutation p, return the permutation invp such that for any
    !! element i we have p(invp(i)) == pinv(p(i)) == i.
    function get_inv_perm(p) result(invp)
        integer, intent(in) :: p(:)
        integer, allocatable :: invp(:)

        allocate(invp(size(p)))
        call inv_perm(p, invp)
    end

    !> Get inverse permutation (non-allocating version of get_inv_perm)
    subroutine inv_perm(p, invp)
        integer, intent(in) :: p(:)
        integer, intent(out) :: invp(size(p))

        integer :: i, n

        n = size(p)
        invp(:) = 0
        do i = 1, n
            if (p(i) < 1 .or. p(i) > n) &
                error stop 'Invalid permutation'
            if (invp(p(i)) /= 0) &
                error stop 'Duplicate index in permutation'
            invp(p(i)) = i
        end do
    end

    subroutine perm_to_swaps_alloc(p, s)
        integer, intent(in) :: p(:)
        integer, intent(out) :: s(:)
        integer, allocatable :: work(:)

        allocate(work(size(p)))
        call perm_to_swaps_work(p, s, work)
    end

    subroutine perm_to_swaps_work(p, s, work)
        integer, intent(in) :: p(:)
        integer, intent(out) :: s(:)
        integer, intent(inout) :: work(size(p))
        integer :: i, j

        if (size(p) /= size(s)) &
            error stop 'Size mismatch'

        call identity_perm(size(p), work)

        ! first element is clear
        ! XXX this algorithm is quadratic, but I'm pretty sure there should be
        !     a linear one.
        do i = 1, size(s)-1
            do j = i, size(s)-1
                if (p(i) == work(j)) &
                    exit
            enddo
            s(i) = j
            if (i /= j) &
                call swap(work(i), work(j))
        enddo
        s(size(s)) = size(s)
    end

    subroutine swaps_to_perm(s, p)
        integer, intent(in) :: s(:)
        integer, intent(out) :: p(:)
        integer :: i, j

        call identity_perm(size(s), p)
        do i = 1, size(s)
            j = s(i)
            if (i /= j) &
                call swap(p(i), p(j))
        enddo
    end

    ! Swap integers
    pure subroutine p_swap_int(x, y)
        integer, intent(inout) :: x, y
        integer :: tmp

        tmp = x
        x = y
        y = tmp
    end

    !> Return parity of a LAPACK permutation array
    !!
    !! LAPACK represents the permutation as a set of 2-cycle swaps, i.e.,
    !! you're supposed to go through the array from i = 1, ..., n and swap
    !! element i with swap(i).  The parity computation is then significantly
    !! simpler.
    function swap_parity(s) result(signbit)
        integer, intent(in) :: s(:)
        logical :: signbit

        integer :: i, n

        n = size(s)
        signbit = .false.
        do i = 1, n
            if (s(i) < 1 .or. s(i) > n) &
                error stop 'Invalid swap array'
            if (s(i) /= i) &
                signbit = .not. signbit
        end do
    end

    function perm_parity_alloc(p) result(signbit)
        integer, intent(in) :: p(:)
        logical :: signbit, visited(size(p))

        signbit = perm_parity_work(p, visited)
    end

    function perm_parity_work(p, visited) result(signbit)
        integer, intent(in) :: p(:)
        logical, intent(out) :: visited(size(p))
        logical :: signbit

        integer :: cur, start

        ! Determine cycles
        signbit = .false.
        visited(:) = .false.

        do start = 1, size(p)-1
            if (visited(start)) &
                cycle
            visited(start) = .true.
            cur = p(start)
            do while (cur /= start)
                if (cur < 1 .or. cur > size(p)) &
                    error stop 'Invalid index'
                if (visited(cur)) &
                    error stop 'Not a permutation'
                visited(cur) = .true.
                signbit = .not. signbit
                cur = p(cur)
            end do
        end do
    end

    !> Returns true if and only if argument is a permutation
    logical function is_permutation(p) result(test)
        integer, intent(in) :: p(:)
        logical :: visited(size(p))

        integer :: cur, start

        test = .true.
        visited(:) = .false.
        do start = 1, size(p)
            if (visited(start)) &
                cycle
            visited(start) = .true.
            cur = p(start)
            do while (cur /= start)
                if (cur < 1 .or. cur > size(p)) then
                    test = .false.
                    return
                endif
                if (visited(cur)) then
                    test = .false.
                    return
                endif
                visited(cur) = .true.
                cur = p(cur)
            end do
        end do
    end

end module
