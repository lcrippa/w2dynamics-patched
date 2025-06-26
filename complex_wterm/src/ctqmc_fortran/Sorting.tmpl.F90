#if 0
! This is a templated Fortran file.
!
! By invoking `make`, several modules are generated in sub-directories
! double/ etc, for each valid value type (real, complex, etc.).
! In each generated module, the value VALUE_TYPE is replaced by the
! corresponding value type, and each uppercase derived type is replaced by
! a typed version.
#endif

module MSORTING
    use MPermutation
    use iso_c_binding, only: c_double
    implicit none
    private

    !>      issorted(arr)
    !!
    !! Return if array of values are sorted from smallest to largest. NaNs,
    !! if present, are ignored.
    interface issorted
        module procedure p_issorted
    end interface
    public issorted

    !>      call sort(a[, work])
    !!
    interface sort
        module procedure p_sort, p_sort_work
    end interface
    public sort

    !>      get_sort_perm(a)
    !!
    !! Return the permutation p, such that a(p) sorts the array from smallest
    !! to largest value.  Does not modify a.  Guaranteed to be a stable sort.
    interface get_sort_perm
        module procedure p_get_sort_perm
    end interface
    public get_sort_perm

    !>      call sort_perm(a, p)
    !!
    !! Sorts the array a from smallest to largest value in-place, i.e., on exit,
    !! a is replaced by its sorted version.  Writes an array p, such that a(p)
    !! on the original array a would also have sorted it.  Guaranteed to be a
    !! stable sort.
    interface sort_perm
        module procedure p_sort_perm_alloc
    end interface
    public sort_perm

    !>      call sort_modify_perm(a, p[, awork, pwork])
    !!
    !! Sorts the array a from smallest to largest value in-place, i.e., a
    !! is overrides with its sorted version.  Given an initial permutation p,
    !! update the permutation together with a.  Guaranteed to be a
    !! stable sort.
    interface sort_modify_perm
        module procedure p_sort_modify_perm_alloc, p_sort_modify_perm_work
    end interface
    public sort_modify_perm

    !>      get_sorted_insert(a, new)
    !!
    !! Given a sorted array `a`, suppose we insert the new element `new`.  Give
    !! the position `pos` which the new element would in a sorted array of the
    !! elements of `a` and `new`.  In the case where there is a tie, insert the
    !! new elements *after* any equivalent elements already in `a`.
    !!
    !! If `new` is an array, `pos` returned is an array which is filled with
    !! the position of each element after all elements have been inserted in
    !! order of appearance in new.
    interface get_sorted_insert
        module procedure p_get_sorted_insert0, p_get_sorted_insert1
    end interface
    public get_sorted_insert

    !>      call sorted_insert(a, new, pos)
    !!
    !! Given a sorted array `a`, suppose we insert the new element `new`.  Give
    !! the position `pos` which the new element would in a sorted array of the
    !! elements of `a` and `new`.  In the case where there is a tie, insert the
    !! new elements *after* any equivalent elements already in `a`.
    !!
    !! If `new` is an array, `pos` must be an array which is filled with the
    !! position of each element after all elements have been inserted in
    !! order of appearance in new.
    interface sorted_insert
        module procedure p_sorted_insert0, p_sorted_insert1
    end interface
    public sorted_insert

contains

    pure logical function p_issorted(x) result(answer)
        VALUE_TYPE, intent(in) :: x(:)
        integer :: i

        do i = 1, size(x) - 1
            if (x(i+1) < x(i)) then
                answer = .false.
                return
            endif
        enddo
        answer = .true.
    end

    subroutine p_sort(a)
        VALUE_TYPE, intent(inout) :: a(:)
        VALUE_TYPE :: awork(size(a))

        call p_sort_work(a, awork)
    end subroutine

    subroutine p_sort_work(a, awork)
        VALUE_TYPE, intent(inout) :: a(:), awork(:)

        if (size(awork) < size(a)) &
            stop 'work arrays are too small'

        awork(:) = a(:)
        call merge_split(awork, a)
    end subroutine

    recursive subroutine merge_split(a1, a2)
        VALUE_TYPE, intent(inout) :: a1(:), a2(:)

        integer :: mid
        if (size(a1) <= 1) &
            return

        mid = (1 + size(a1)) / 2
        call merge_split(a2(:mid), a1(:mid))
        call merge_split(a2(mid+1:), a1(mid+1:))
        call merge_sorted(a1(:mid), a1(mid+1:), a2)
    end subroutine

    subroutine merge_sorted(a, b, m)
        VALUE_TYPE, intent(in) :: a(:), b(:)
        VALUE_TYPE, intent(out) :: m(:)

        integer :: ia, ib, im

        ia = 1
        ib = 1
        do im = 1, size(m)
            if (.not. (b(ib) < a(ia))) then
                m(im) = a(ia)
                ia = ia + 1
                if (ia > size(a)) then
                    m(im+1:) = b(ib:)
                    return
                endif
            else
                m(im) = b(ib)
                ib = ib + 1
                if (ib > size(b)) then
                    m(im+1:) = a(ia:)
                    return
                endif
            endif
        enddo
    end subroutine

    function p_get_sort_perm(a) result(perm)
        VALUE_TYPE, intent(in) :: a(:)
        VALUE_TYPE :: acopy(size(a))
        integer, allocatable :: perm(:)

        allocate(perm(size(a)))

        acopy(:) = a(:)
        call sort_perm(acopy, perm)
    end

    subroutine p_sort_perm_alloc(a, p)
        VALUE_TYPE, intent(inout) :: a(:)
        integer, intent(inout) :: p(:)

        call identity_perm(size(a), p)
        call sort_modify_perm(a, p)
    end subroutine

    subroutine p_sort_modify_perm_alloc(a, p)
        VALUE_TYPE, intent(inout) :: a(:)
        integer, intent(inout) :: p(:)

        VALUE_TYPE :: awork(size(a))
        integer :: pwork(size(a))
        call sort_modify_perm(a, p, awork, pwork)
    end subroutine

    subroutine p_sort_modify_perm_work(a, p, awork, pwork)
        VALUE_TYPE, intent(inout) :: a(:), awork(:)
        integer, intent(inout) :: p(:), pwork(:)

        if (size(a) /= size(p)) &
            stop 'a and p must be of the same size'
        if (size(awork) < size(a) .or. size(pwork) < size(p)) &
            stop 'work arrays are too small'

        pwork(:) = p(:)
        awork(:) = a(:)
        call merge_split_perm(awork, pwork, a, p)
    end subroutine

    recursive subroutine merge_split_perm(a1, p1, a2, p2)
        VALUE_TYPE, intent(inout) :: a1(:), a2(:)
        integer, intent(inout) :: p1(:), p2(:)

        integer :: mid
        if (size(a1) <= 1) &
            return

        mid = (1 + size(a1)) / 2
        call merge_split_perm(a2(:mid), p2(:mid), a1(:mid), p1(:mid))
        call merge_split_perm(a2(mid+1:), p2(mid+1:), a1(mid+1:), p1(mid+1:))
        call merge_sorted_perm(a1(:mid), p1(:mid), a1(mid+1:), p1(mid+1:), a2, p2)
    end subroutine

    subroutine merge_sorted_perm(a, ap, b, bp, m, mp)
        VALUE_TYPE, intent(in) :: a(:), b(:)
        integer, intent(in) :: ap(:), bp(:)
        VALUE_TYPE, intent(out) :: m(:)
        integer, intent(out) :: mp(:)

        integer :: ia, ib, im

        ia = 1
        ib = 1
        do im = 1, size(m)
            if (.not. (b(ib) < a(ia))) then
                m(im) = a(ia)
                mp(im) = ap(ia)
                ia = ia + 1
                if (ia > size(a)) then
                    m(im+1:) = b(ib:)
                    mp(im+1:) = bp(ib:)
                    return
                endif
            else
                m(im) = b(ib)
                mp(im) = bp(ib)
                ib = ib + 1
                if (ib > size(b)) then
                    m(im+1:) = a(ia:)
                    mp(im+1:) = ap(ia:)
                    return
                endif
            endif
        enddo
    end

    pure function p_get_sorted_insert0(a, new) result(pos)
        VALUE_TYPE, intent(in) :: a(:), new
        integer :: pos

        call sorted_insert(a, new, pos)
    end

    function p_get_sorted_insert1(a, new) result(pos)
        VALUE_TYPE, intent(in) :: a(:), new(:)
        integer, allocatable :: pos(:)

        allocate(pos(size(new)))
        call sorted_insert(a, new, pos)
    end

    pure subroutine p_sorted_insert0(a, new, pos)
        VALUE_TYPE, intent(in) :: a(:), new
        integer, intent(out) :: pos
        integer :: lo, mid, hi

        lo = 1
        hi = size(a)
        if (hi == 0) then
            pos = 1
        elseif (new < a(lo)) then
            pos = 1
        else if (.not. (new < a(hi))) then
            pos = hi + 1
        else
            do while (lo /= hi)
                mid = (lo + hi) / 2
                if (new < a(mid)) then
                    hi = mid
                else
                    lo = mid + 1
                endif
            end do
            pos = hi
        endif
    end

    subroutine p_sorted_insert1(a, new, pos)
        VALUE_TYPE, intent(in) :: a(:), new(:)
        integer, intent(out) :: pos(:)

        integer :: i, j

        do i = 1, size(new)
            pos(i) = get_sorted_insert(a, new(i))
        enddo

        ! TODO this may be slower than first sorting and then adding one
        do i = 1, size(new)
            do j = i + 1, size(new)
                if (new(j) < new(i)) then
                    pos(i) = pos(i) + 1
                else
                    pos(j) = pos(j) + 1
                endif
            enddo
        enddo
    end

end module
