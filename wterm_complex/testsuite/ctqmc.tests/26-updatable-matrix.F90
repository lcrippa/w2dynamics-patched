#define INFO line=__LINE__

program test_updmatrix
    use iso_c_binding
    use MUpdatableMatrixD
    use MUpdatableMatrixZ
    use MExtRange
    use MBLAS
    use MLAPACK
    use MPrinting
    use MRandomNumbers
    use Testing
    implicit none

    call test_simple()
    call test_hilbert_grow()
    call test_hilbert_shrink()
    call test_hilbert_rank1()
    call test_hilbert_rows()
    call test_hilbert_cols()
    call test_hilbert_perm()

contains

    subroutine test_simple()
        type(DUpdatableMatrix) :: mat
        type(DMatrixReplace) :: repl

        call assert_close(limrange(get_det(mat)), 1.0d0, INFO)
        call verify(mat)

        call propose_replace_matrix(repl, mat, &
                reshape((/ 1.0d0, 1.0d0, 0.0d0, 1.0d0 /), (/ 2, 2 /)))
        call assert_close(limrange(get_det(mat)), 1.0d0, INFO)

        call accept_update(repl)
        call assert_close(getinv(mat), &
                reshape((/ 1.0d0, -1.0d0, 0.0d0, 1.0d0 /), (/ 2, 2 /)), INFO)
        call verify(mat)
    end

    subroutine test_hilbert_grow()
        type(DUpdatableMatrix) :: mat
        type(DMatrixGrow) :: grow
        real(c_double), parameter :: rtol = 1d-4   ! XXX

        real(c_double) :: h(10, 10)

        h = hilbert_matrix(10)
        call setmatrix(mat, h(1:8, 1:8))
        call verify(mat, rtol=rtol)

        call propose_grow_matrix(grow, mat, (/ 9, 10 /), (/ 9, 10 /), &
                                 h(9:, :8), h(:8, 9:), h(9:, 9:))
        call accept_update(grow)
        call verify(mat, rtol=rtol)
        call assert_close(getmatrix(mat), h)

        call setmatrix(mat, h(3:, 3:))
        call verify(mat, rtol=rtol)

        call propose_grow_matrix(grow, mat, (/ 1, 2 /), (/ 1, 2 /), &
                                 h(:2, 3:), h(3:, :2), h(:2, :2))
        call accept_update(grow)
        call verify(mat, rtol=rtol)  ! Hilbert is famously ill-conditioned

        call setmatrix(mat, h(2:10, 1:9))
        call verify(mat, rtol=rtol)
        call propose_grow_matrix(grow, mat, (/ 1 /), (/ 10 /), &
                                 h(:1, :9), h(2:, 10:), h(:1, 10:))
        call accept_update(grow)
        call verify(mat, rtol=rtol)  ! Hilbert is famously ill-conditioned

        call assert_close(getmatrix(mat), h)
    end

    subroutine test_hilbert_shrink()
        type(DUpdatableMatrix) :: mat
        type(DMatrixShrink) :: shrink
        real(c_double), parameter :: rtol = 1d-4

        real(c_double) :: h(10, 10)

        h = hilbert_matrix(10)
        call setmatrix(mat, h)
        call verify(mat, rtol=rtol)

        call propose_shrink_matrix(shrink, mat, (/ 1, 2 /), (/ 9, 10 /))
        call accept_update(shrink)
        call verify(mat, rtol=rtol)

        call assert_close(getmatrix(mat), h(3:, :8))
    end

    subroutine test_hilbert_rank1()
        type(DUpdatableMatrix) :: mat
        type(DLowRankUpdate) :: upd
        real(c_double), parameter :: rtol = 1d-4

        real(c_double) :: h(10, 10), u(10, 1), vdag(1, 10)
        integer :: i

        h = hilbert_matrix(10)
        call setmatrix(mat, h)
        call verify(mat, rtol=rtol)

        u(:, 1) = (/ (1.0 * i, i=1,10) /)
        vdag(1, :) = (/ (1.0 * i, i=1,10) /)

        call propose_low_rank_update(upd, mat, u, vdag, &
                                   (/ (i, i=1,10) /), (/ (i, i=1,10) /) )
        call accept_update(upd)
        call verify(mat, rtol=6d-2)
    end

    subroutine test_hilbert_rows()
        type(DUpdatableMatrix) :: mat
        type(DReplaceRows) :: upd
        real(c_double), parameter :: rtol = 1d-4

        real(c_double) :: h(10, 10), u(1, 10)
        integer :: i

        h = hilbert_matrix(10)
        call setmatrix(mat, h)
        call verify(mat, rtol=rtol)

        u(1, :) = (/ (1.0 * i, i=1,10) /)

        call propose_replace_rows(upd, mat, (/ 3 /), (/ 2 /), u)
        call accept_update(upd)
        call verify(mat, rtol=rtol)

        h(3, :) = h(2, :)
        h(2, :) = u(1, :)
        call assert_close(getmatrix(mat), h, rtol=6d-2)
    end

    subroutine test_hilbert_cols()
        type(DUpdatableMatrix) :: mat
        type(DReplaceCols) :: upd
        real(c_double), parameter :: rtol = 1d-4

        real(c_double) :: h(10, 10), u(10, 1)
        integer :: i

        h = hilbert_matrix(10)
        call setmatrix(mat, h)
        call verify(mat, rtol=rtol)

        u(:, 1) = (/ (1.0 * i, i=1,10) /)

        call propose_replace_cols(upd, mat, (/ 1 /), (/ 4 /), u)
        call accept_update(upd)
        call verify(mat, rtol=rtol)

        h(:, 1:3) = h(:, 2:4)
        h(:, 4) = u(:, 1)
        call assert_close(getmatrix(mat), h, rtol=6d-2)
    end

    function hilbert_matrix(n) result(h)
        integer, intent(in) :: n
        real(c_double), allocatable :: h(:, :)

        integer :: i, j
        allocate(h(n, n))

        do j = 1, n
            do i = 1, n
                h(i, j) = 1 / (i + j - 1.0d0)
            end do
        end do
    end

    subroutine test_hilbert_perm()
        type(DUpdatableMatrix) :: mat
        type(DMatrixPermute) :: perm
        real(c_double), parameter :: rtol = 1d-4

        real(c_double) :: h(5, 5), iscale(5), jscale(5)
        integer :: iperm(5), jperm(5)

        iscale = (/ 1.0, -1.0, -1.0, 1.0, -1.0 /)
        jscale = (/ 1.0, 1.0, -1.0, -1.0, 1.0 /)
        iperm = (/ 1, 3, 2, 5, 4 /)
        jperm = (/ 5, 1, 4, 3, 2 /)

        h = hilbert_matrix(10)
        call setmatrix(mat, h)
        call verify(mat, rtol=rtol)

        call propose_permute_matrix(perm, mat, iperm, jperm, iscale, jscale)
        call accept_update(perm)
        call verify(mat, rtol=rtol)
    end

end program
