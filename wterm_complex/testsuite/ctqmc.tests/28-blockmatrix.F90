#define INFO line=__LINE__

program test_block_matrix
    use iso_c_binding
    use MBlockMatrixD
    use MBlockMatrixZ
    use MExtRange
    use Testing
    implicit none

    call test_simple()
    call test_grow()
    call test_shrink()
    call test_rows()
    call test_cols()

contains
    subroutine test_simple()
        type(DBlockMatrix) :: mat
        type(DBlockReplace) :: repl

        call initialize(mat, 3)
        call assert_close(limrange(get_det(mat)), 1.0d0, INFO)
        call verify(mat, 1d-14)

        call propose_replace_matrix( &
                repl, mat, &
                (/ 1, 1 /), (/ 1, 2 /), &
                reshape((/ &
                    2d0, 0d0, &
                    0d0, 0d0 /), (/ 2, 2 /)) &
                )
        call assert_close(limrange(get_det(repl)), 0.0d0, INFO)
        call assert_equal(is_update_allowed(repl), .false., INFO)

        call propose_replace_matrix( &
                repl, mat, &
                (/ 1, 2, 2, 3 /), (/ 1, 2, 2, 3 /), &
                reshape((/ &
                    2d0, 0d0, 0d0, 0d0, &
                    0d0, 1d0, 1d0, 0d0, &
                    0d0, 0d0, 1d0, 0d0, &
                    0d0, 0d0, 0d0, 4d0 /), (/ 4, 4 /)) &
                )
        call assert_close(limrange(get_det(repl)), 8.0d0, INFO)
        call assert_equal(is_update_allowed(repl), .true., INFO)

        call accept_update(repl)
        call verify(mat, 1d-14)
        call assert_close( &
                getinv(mat), &
                reshape((/ &
                    0.5d0, 0d0,  0d0, 0d0, &
                    0d0,   1d0, -1d0, 0d0, &
                    0d0,   0d0,  1d0, 0d0, &
                    0d0,   0d0,  0d0, 0.25d0 /), (/ 4, 4 /)), &
                INFO)
    end

    subroutine test_grow()
        type(DBlockMatrix) :: mat
        type(DBlockGrow) :: move

        integer :: orig_blocks(4) = (/ 1, 2, 2, 3/)
        real(c_double) :: orig(4, 4) = reshape((/ &
                2d0, 0d0, 0d0, 0d0, &
                0d0, 1d0, 1d0, 0d0, &
                0d0, 5d0, 6d0, 0d0, &
                0d0, 0d0, 0d0, 4d0  &
                /), shape(orig), order=(/ 2, 1 /))
        real(c_double) :: grown(6, 6) = reshape((/ &
        !       vvv                 vvv
                4d0, 1d0, 0d0, 0d0, 0d0, 0d0, &   ! <<<
                1d0, 2d0, 0d0, 0d0, 0d0, 0d0, &
                0d0, 0d0, 1d0, 1d0, 3d0, 0d0, &
                0d0, 0d0, 0d0, 1d0, 4d0, 0d0, &   ! <<<
                0d0, 0d0, 5d0, 6d0, 7d0, 0d0, &
                0d0, 0d0, 0d0, 0d0, 0d0, 4d0  &
                /), shape(grown), order=(/ 2, 1 /))
        real(c_double) :: rows(2, 4) = reshape((/ &
                1d0, 0d0, 0d0, 0d0, &
                0d0, 0d0, 1d0, 0d0  &
                /), shape(rows), order=(/ 2, 1 /))
        real(c_double) :: cols(4, 2) = reshape((/ &
                1d0, 0d0, &
                0d0, 3d0, &
                0d0, 7d0, &
                0d0, 0d0  &
                /), shape(cols), order=(/ 2, 1 /))
        real(c_double) :: dot(2, 2) = reshape((/ &
                4d0, 0d0, &
                0d0, 4d0  &
                /), shape(dot), order=(/ 2, 1 /))
        integer :: rowind(2) = (/ 1, 4 /)
        integer :: colind(2) = (/ 1, 5 /)
        integer :: insblocks(2) = (/ 1, 2 /)
        real(c_double) :: grow_rec(6, 6)

        call initialize(mat, 3)
        call setmatrix(mat, orig_blocks, orig)
        call assert_close(limrange(get_det(mat)), 8.0d0, rtol=2d-14, INFO)
        call verify(mat, 2d-14)

        call propose_grow_matrix(move, mat, insblocks, insblocks, rowind, &
                                 colind, rows, cols, dot)
        call assert_close(limrange(get_det(move)), -336.0d0, rtol=2d-14, INFO)

        call accept_update(move)
        call verify(mat, 2d-14)
        grow_rec(:, :) = getmatrix(mat)
        call assert_close(grow_rec, grown, rtol=2d-14)
    end

    subroutine test_shrink()
        type(DBlockMatrix) :: mat
        type(DBlockShrink) :: move

        integer :: grown_blocks(6) = (/ 1, 1, 2, 2, 2, 3 /)
        real(c_double) :: orig(4, 4) = reshape((/ &
                2d0, 0d0, 0d0, 0d0, &
                0d0, 1d0, 1d0, 0d0, &
                0d0, 5d0, 6d0, 0d0, &
                0d0, 0d0, 0d0, 4d0  &
                /), shape(orig), order=(/ 2, 1 /))
        real(c_double) :: grown(6, 6) = reshape((/ &
        !       vvv                 vvv
                4d0, 1d0, 0d0, 0d0, 0d0, 0d0, &   ! <<<
                1d0, 2d0, 0d0, 0d0, 0d0, 0d0, &
                0d0, 0d0, 1d0, 1d0, 3d0, 0d0, &
                0d0, 0d0, 0d0, 1d0, 4d0, 0d0, &   ! <<<
                0d0, 0d0, 5d0, 6d0, 7d0, 0d0, &
                0d0, 0d0, 0d0, 0d0, 0d0, 4d0  &
                /), shape(grown), order=(/ 2, 1 /))
        integer :: rowind(2) = (/ 4, 1 /)
        integer :: colind(2) = (/ 1, 5 /)

        call initialize(mat, 3)
        call setmatrix(mat, grown_blocks, grown)
        call assert_close(limrange(get_det(mat)), -336.0d0, rtol=2d-14, INFO)
        call verify(mat, 2d-14)

        call propose_shrink_matrix(move, mat, (/ 1 /), (/ 3 /))
        call assert_close(limrange(get_det(move)), 0.0d0, INFO)

        call propose_shrink_matrix(move, mat, rowind, colind)
        call assert_close(limrange(get_det(move)), 8.0d0, rtol=2d-14, INFO)

        call accept_update(move)
        call verify(mat, 2d-14)

        call assert_close(getmatrix(mat), orig, rtol=2d-14)
    end

    subroutine test_rows()
        type(DBlockMatrix) :: mat
        type(DBlockReplaceRows) :: move

        integer :: grown_blocks(6) = (/ 1, 1, 2, 2, 2, 3 /)
        real(c_double) :: grown(6, 6) = reshape((/ &
                4d0, 1d0, 0d0, 0d0, 0d0, 0d0, &   ! <<<
                1d0, 2d0, 0d0, 0d0, 0d0, 0d0, &
                0d0, 0d0, 1d0, 1d0, 3d0, 0d0, &
                0d0, 0d0, 0d0, 1d0, 4d0, 0d0, &   ! <<<
                0d0, 0d0, 5d0, 6d0, 7d0, 0d0, &
                0d0, 0d0, 0d0, 0d0, 0d0, 4d0  &
                /), shape(grown), order=(/ 2, 1 /))
        integer :: rowind(2) = (/ 4, 1 /)

        real(c_double) :: rows(2, 6) = reshape((/ &
                0d0, 0d0, 0d0, 2d0, 9d0, 0d0, &   ! <<<
                5d0, 2d0, 0d0, 0d0, 0d0, 0d0 &   ! <<<
                /), shape(rows), order=(/ 2, 1 /))
        real(c_double) :: reshuffled(6, 6) = reshape((/ &
                1d0, 2d0, 0d0, 0d0, 0d0, 0d0, &
                5d0, 2d0, 0d0, 0d0, 0d0, 0d0, &   ! <<<
                0d0, 0d0, 0d0, 2d0, 9d0, 0d0, &   ! <<<
                0d0, 0d0, 1d0, 1d0, 3d0, 0d0, &
                0d0, 0d0, 5d0, 6d0, 7d0, 0d0, &
                0d0, 0d0, 0d0, 0d0, 0d0, 4d0  &
                /), shape(reshuffled), order=(/ 2, 1 /))
        real(c_double) :: reshuffled_rec(6, 6)

        call initialize(mat, 3)
        call setmatrix(mat, grown_blocks, grown)
        call assert_close(limrange(get_det(mat)), -336.0d0, rtol=1d-14, INFO)
        call verify(mat, 1d-14)

        call propose_replace_rows(move, mat, rowind, (/ 3, 2 /), rows)
        !call assert_close(limrange(get_det(move)), 8.0d0, rtol=1d-14, INFO)

        call accept_update(move)
        call verify(mat, 1d-14)

        reshuffled_rec(:, :) = getmatrix(mat)
        call assert_close(getmatrix(mat), reshuffled, rtol=1d-14)
    end

    subroutine test_cols()
        type(DBlockMatrix) :: mat
        type(DBlockReplaceCols) :: move

        integer :: grown_blocks(6) = (/ 1, 1, 2, 2, 2, 3 /)
        real(c_double) :: grown(6, 6) = reshape((/ &
                4d0, 1d0, 0d0, 0d0, 0d0, 0d0, &   ! <<<
                1d0, 2d0, 0d0, 0d0, 0d0, 0d0, &
                0d0, 0d0, 1d0, 1d0, 3d0, 0d0, &
                0d0, 0d0, 0d0, 1d0, 4d0, 0d0, &   ! <<<
                0d0, 0d0, 5d0, 6d0, 7d0, 0d0, &
                0d0, 0d0, 0d0, 0d0, 0d0, 4d0  &
                /), shape(grown), order=(/ 1, 2 /))
        integer :: colind(2) = (/ 4, 1 /)

        real(c_double) :: cols(6, 2) = reshape((/ &
                0d0, 0d0, 0d0, 2d0, 9d0, 0d0, &   ! <<<
                5d0, 2d0, 0d0, 0d0, 0d0, 0d0 &   ! <<<
                /), shape(cols), order=(/ 1, 2 /))
        real(c_double) :: reshuffled(6, 6) = reshape((/ &
                1d0, 2d0, 0d0, 0d0, 0d0, 0d0, &
                5d0, 2d0, 0d0, 0d0, 0d0, 0d0, &   ! <<<
                0d0, 0d0, 0d0, 2d0, 9d0, 0d0, &   ! <<<
                0d0, 0d0, 1d0, 1d0, 3d0, 0d0, &
                0d0, 0d0, 5d0, 6d0, 7d0, 0d0, &
                0d0, 0d0, 0d0, 0d0, 0d0, 4d0  &
                /), shape(reshuffled), order=(/ 1, 2 /))
        real(c_double) :: reshuffled_rec(6, 6)

        call initialize(mat, 3)
        call setmatrix(mat, grown_blocks, grown)
        call assert_close(limrange(get_det(mat)), -336.0d0, rtol=1d-14, INFO)
        call verify(mat, 1d-14)

        call propose_replace_cols(move, mat, colind, (/ 3, 2 /), cols)
        !call assert_close(limrange(get_det(move)), 8.0d0, rtol=1d-14, INFO)

        call accept_update(move)
        call verify(mat, 1d-14)

        reshuffled_rec(:, :) = getmatrix(mat)
        call assert_close(getmatrix(mat), reshuffled, rtol=1d-14)
    end

end program
