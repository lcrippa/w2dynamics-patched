#define INFO line=__LINE__

program test_permutation
    use iso_c_binding, only: c_double, c_long
    use MPermutation
    use MSorting
    use MRandomNumbers
    use Testing
    implicit none

    call test_cyclic_perm()
    call test_permsign()
    call test_grow_perm()
    call test_grow_perm_front()
    call test_inv_perm()
    call test_shrink_perm()
    call test_move_perm()

    call test_perm_remove()

    call test_sort_perm()
    call test_compound_sort()
    call test_sort_rand_large()
    call test_sorted_insert()

    call test_rand_perm()
    call test_rand_perm_swaps(4)
    call test_rand_perm_swaps(10)

    call test_is_perm()
    call test_group_perm()
contains
    subroutine test_permsign()
        call assert_equal(perm_parity((/ 1 /)), .false., INFO)
        call assert_equal(perm_parity((/ 1, 2 /)), .false., INFO)
        call assert_equal(perm_parity((/ 2, 1 /)), .true., INFO)
        call assert_equal(perm_parity((/ 6, 1, 5, 4, 3, 2 /)), .true., INFO)
    end

    subroutine test_cyclic_perm()
        call assert_equal(get_cyclic_perm(3, 0), &
                          (/ 1, 2, 3 /), INFO)
        call assert_equal(get_cyclic_perm(3, 1), &
                          (/ 2, 3, 1 /), INFO)
        call assert_equal(get_cyclic_perm(3, -1), &
                          (/ 3, 1, 2 /), INFO)
        call assert_equal(get_cyclic_perm(4, 3), &
                          (/ 4, 1, 2, 3 /), INFO)
    end

    subroutine test_grow_perm()
        call assert_equal(get_grow_perm(3, (/ 4 /)), &
                          (/ 1, 2, 3, 4 /), INFO)
        call assert_equal(get_grow_perm(3, (/ 2 /)), &
                          (/ 1, 4, 2, 3 /), INFO)
        call assert_equal(get_grow_perm(3, (/ 1, 3 /)), &
                          (/ 4, 1, 5, 2, 3 /), INFO)
        call assert_equal(get_grow_perm(4, (/ 3, 1 /)), &
                          (/ 6, 1, 5, 2, 3, 4 /), INFO)
        call assert_equal(get_grow_perm(3, (/ 5, 4 /)), &
                          (/ 1, 2, 3, 5, 4 /), INFO)
        call assert_equal(get_grow_perm(3, (/ 2, 3 /)), &
                          (/ 1, 4, 5, 2, 3 /), INFO)
        call assert_equal(get_grow_perm(3, (/ 3, 2 /)), &
                          (/ 1, 5, 4, 2, 3 /), INFO)
        call assert_equal(get_grow_perm(3, (/ 2, 5 /)), &
                          (/ 1, 4, 2, 3, 5 /), INFO)
        call assert_equal(get_grow_perm(3, (/ 5, 2 /)), &
                          (/ 1, 5, 2, 3, 4 /), INFO)
    end

    subroutine test_grow_perm_front()
        call assert_equal(get_grow_perm_front(3, (/ 1 /)), &
                          (/ 1, 2, 3, 4 /), INFO)
        call assert_equal(get_grow_perm_front(3, (/ 3 /)), &
                          (/ 2, 3, 1, 4 /), INFO)
        call assert_equal(get_grow_perm_front(3, (/ 4 /)), &
                          (/ 2, 3, 4, 1 /), INFO)
        call assert_equal(get_grow_perm_front(3, (/ 2, 4 /)), &
                          (/ 3, 1, 4, 2, 5 /), INFO)
        call assert_equal(get_grow_perm_front(3, (/ 4, 2 /)), &
                          (/ 3, 2, 4, 1, 5 /), INFO)
        call assert_equal(get_grow_perm_front(4, (/ 2, 1 /)), &
                          (/ 2, 1, 3, 4, 5, 6 /), INFO)
    end

    subroutine test_inv_perm()
        call assert_equal(get_inv_perm((/ 1 /)), &
                          (/ 1 /), INFO)
        call assert_equal(get_inv_perm((/ 1, 2 /)), &
                          (/ 1, 2 /), INFO)
        call assert_equal(get_inv_perm((/ 2, 1 /)), &
                          (/ 2, 1 /), INFO)
        call assert_equal(get_inv_perm((/ 6, 1, 5, 4, 3, 2 /)), &
                          (/ 2, 6, 5, 4, 3, 1 /), INFO)
    end

    subroutine test_shrink_perm()
        call assert_equal(get_shrink_perm(5, (/ 5 /)), &
                          (/ 1, 2, 3, 4, 5 /), INFO)
        call assert_equal(get_shrink_perm(5, (/ 3 /)), &
                          (/ 1, 2, 4, 5, 3 /), INFO)
        call assert_equal(get_shrink_perm(5, (/ 2, 3 /)), &
                          (/ 1, 4, 5, 2, 3 /), INFO)
        call assert_equal(get_shrink_perm(5, (/ 3, 2 /)), &
                          (/ 1, 4, 5, 3, 2 /), INFO)
        call assert_equal(get_shrink_perm(5, (/ 2, 4 /)), &
                          (/ 1, 3, 5, 2, 4 /), INFO)
        call assert_equal(get_shrink_perm(5, (/ 5, 2 /)), &
                          (/ 1, 3, 4, 5, 2 /), INFO)
    end

    subroutine test_move_perm()
        call assert_equal(get_move_perm(3, (/ 2 /), (/ 2 /)), &
                          (/ 1, 2, 3 /), INFO)
        call assert_equal(get_move_perm(3, (/ 2 /), (/ 3 /)), &
                          (/ 1, 3, 2 /), INFO)
        call assert_equal(get_move_perm(3, (/ 3 /), (/ 2 /)), &
                          (/ 1, 3, 2 /), INFO)
        call assert_equal(get_move_perm(7, (/ 2, 5 /), (/ 3, 4 /)), &
                          (/ 1, 3, 2, 5, 4, 6, 7 /), INFO)
        call assert_equal(get_move_perm(7, (/ 2, 5 /), (/ 4, 3 /)), &
                          (/ 1, 3, 5, 2, 4, 6, 7 /), INFO)
    end

    subroutine test_perm_remove()
        call assert_equal(get_perm_remove((/ 2, 1, 5, 3, 4/), (/ 3 /)), &
                          (/ 2, 1, 3, 4/))
        call assert_equal(get_perm_remove((/ 2, 1, 5, 3, 4/), (/ 4 /)), &
                          (/ 2, 1, 4, 3/))
        call assert_equal(get_perm_remove((/ 2, 1, 3 /), (/ 1, 2 /)), &
                          (/ 1 /))
        call assert_equal(get_perm_remove((/ 2, 1, 5, 3, 4 /), (/ 2, 3 /)), &
                          (/ 1, 2, 3 /))
    end subroutine

    subroutine test_sort_perm()
        real(c_double) :: x1(4) = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0 /)
        real(c_double) :: x2(3) = (/ 0.4d0, 0.1d0, 0.35d0 /)
        real(c_double) :: x3(5) = (/ 1.0d0, 5.0d0, 3.0d0, 2.0d0, 4.0d0 /)

        call assert_equal(get_sort_perm(x1), (/ 1, 2, 3, 4/), INFO)
        call assert_equal(get_sort_perm(x2), (/ 2, 3, 1 /), INFO)
        call assert_equal(get_sort_perm(x3), (/ 1, 4, 3, 5, 2 /), INFO)
    end

    subroutine test_compound_sort()
        type :: compound
            integer :: i
            real(c_double) :: x
        end type
        type(compound) :: arr(3)
        real(c_double) :: x(5) = (/ 1.0d0, 5.0d0, 4.0d0, 2.0d0, 3.0d0 /)

        arr(:)%i = (/ 1, 4, 5 /)
        arr(:)%x = (/ 0.4d0, 0.1d0, 0.35d0 /)

        ! Here, for some reason, gfortran wants to create a temporary :(
        call assert_equal(get_sort_perm(arr%x), (/ 2, 3, 1/), INFO)

        ! Here, no temporary is needed
        call assert_equal(get_sort_perm(x(1:5:2)), (/ 1, 3, 2 /), INFO)
    end

    subroutine test_sort_rand_large()
        class(RandomEngine), allocatable :: rng
        real(c_double), allocatable :: x(:), xsort(:)
        integer, allocatable :: p(:)
        integer :: i, n = 200000

        allocate(x(n))
        allocate(xsort(n))
        allocate(p(n))

        call init_mersenne_twister(rng, 4711)
        do i = 1, size(x)
            x(i) = random_real(rng)
        enddo

        xsort(:) = x(:)
        call sort_perm(xsort, p)
        call assert_equal(issorted(xsort), .true., INFO)
        call assert_close(x(p), xsort)
    end

    subroutine test_rand_perm()
        class(RandomEngine), allocatable :: rng
        integer, allocatable :: x(:)
        integer :: i, n = 20000

        allocate(x(n))
        call init_mersenne_twister(rng, 4711)
        call random_permutation(rng, n, x)
        call assert_equal(issorted(x), .false., INFO)
        call sort(x)
        call assert_equal(x, (/ (i, i=1,n) /), INFO)
    end

    subroutine test_is_perm()
        class(RandomEngine), allocatable :: rng
        integer :: x(200), i, j

        call init_xorshift_star(rng)
        do i = 1, 200
            call random_permutation(rng, i, x(:i))
            !write (0, "(200I4)") x(:i)
            call assert_equal(is_permutation(x(:i)), .true.)

            j = random_integer(rng, 1, i)
            x(j) = x(j) + 1
            call assert_equal(is_permutation(x(:i)), .false.)
        enddo
    end

    subroutine test_rand_perm_swaps(n)
        class(RandomEngine), allocatable :: rng
        integer, intent(in) :: n
        integer :: p(n), s(n), i, sprime(n)

        call init_mersenne_twister(rng, 4711)
        do i = 1, 1000
            call random_swaps(rng, n, s)
            call swaps_to_perm(s, p)
            call assert_equal(swap_parity(s), perm_parity(p))

            call perm_to_swaps(p, sprime)
            call assert_equal(s, sprime)
        enddo
    end

    subroutine test_sorted_insert()
        real(c_double) :: empty(0)
        integer :: emptyint(0)

        call assert_equal(get_sorted_insert(empty, 2.0d0), &
                          1, INFO)
        call assert_equal(get_sorted_insert((/ 3.0d0 /), 2.0d0), &
                          1, INFO)
        call assert_equal(get_sorted_insert((/ 3.0d0 /), 4.0d0), &
                          2, INFO)
        call assert_equal(get_sorted_insert((/ 3.0d0 /), 3.0d0), &
                          2, INFO)
        call assert_equal(get_sorted_insert((/ 2.0d0, 3.0d0, 3.1d0 /), 2.5d0), &
                          2, INFO)
        call assert_equal(get_sorted_insert((/ 2.0d0, 3.0d0, 3.1d0 /), 3.0d0), &
                          3, INFO)

        call assert_equal(get_sorted_insert(empty, &
                                            (/ 1.0d0, 0.5d0 /)), &
                                            (/ 2, 1 /))
        call assert_equal(get_sorted_insert((/ 1.0d0 /), &
                                            empty), &
                                            emptyint)
        call assert_equal(get_sorted_insert((/ 0.5d0 /), &
                                            (/ 1.0d0, 0.3d0 /)), &
                                            (/ 3, 1 /))
        call assert_equal(get_sorted_insert((/ 0.5d0, 0.8d0 /), &
                                            (/ 0.7d0, 0.7d0 /)), &
                                            (/ 2, 3 /))
        call assert_equal(get_sorted_insert((/ 0.2d0, 0.5d0, 1.5d0, 2.0d0 /), &
                                            (/ 0.4d0, 0.0d0, 2.3d0, 1.6d0 /)), &
                                            (/ 3, 1, 8, 6 /))
        call assert_equal(get_sorted_insert((/ 0.2d0, 0.5d0, 1.5d0, 2.0d0 /), &
                                            (/ 0.4d0, 0.0d0, 2.3d0, 0.4d0 /)), &
                                            (/ 3, 1, 8, 4 /))
    end

    subroutine test_group_perm()
        integer :: x(10), g(6), p(10)

        call assert_equal(get_group_perm((/ 2, 1, 1, 3, 1, 3 /)), &
                          (/ 2, 3, 5, 1, 4, 6 /))
        call assert_equal(get_group_perm((/ 1, 1 /)), &
                          (/ 1, 2 /))

        x(:) = (/ 2, 1, 5, 5, 1, 4, 4, 5, 1, 2 /)
        call group_perm(x, p, g)
        call assert_equal(x(p), (/ 1, 1, 1, 2, 2, 4, 4, 5, 5, 5 /))
        call assert_equal(g, (/ 3, 5, 5, 7, 10, 10 /))
    end
end program
