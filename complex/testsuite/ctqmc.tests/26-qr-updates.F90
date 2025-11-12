program test_qrupdates
    use iso_c_binding
    use MQRDecompositionD
    use MQRDecompositionZ
    use MUpdatableMatrixD
    use MUpdatableMatrixZ
    use MExtRange
    use MBLAS
    use MLAPACK
    use MPrinting
    use MRandomNumbers
    use Testing
    implicit none

    call test_qr_det()
    call test_qr_det_complex()
    call test_qr_hess()
    call test_qr_hess_complex()
    call test_qr_rot()
    call test_qr_rot_complex()
    call test_qr_update()
    call test_qr_update_complex()

contains

    subroutine test_qr_det()
        real(c_double) :: x(5, 5) = reshape((/ &
                70868.0d0, 39361.0d0, 71337.0d0, 17427.0d0, 84823.0d0, &
                43320.0d0, 24085.0d0,   101.0d0, 33262.0d0, 40862.0d0, &
                48744.0d0, 26480.0d0, 59551.0d0, 29495.0d0, 73639.0d0, &
                39509.0d0, 25141.0d0, 19686.0d0, 72501.0d0, 82889.0d0, &
                61481.0d0, 45640.0d0, 21339.0d0,   749.0d0, 33191.0d0  &
                /), (/ 5, 5 /))

        real(c_double) :: tau(5), work(25)
        type(DExtRange) :: det
        integer :: info

        call DGEQRF(5, 5, x(1,1), size(x, 1), tau(1), work(1), size(work), info)
        call check_info(info)

        det = qr_det(x, tau)
        call assert_close(limrange(det), 1.2249523060153992d+22, rtol=1d-14)
    end

    subroutine test_qr_hess()
        real(c_double) :: x(7, 7), r(7, 7), q(7, 7), x_rec(7, 7)
        type(DExtRange) :: det1, det2
        integer :: i, j
        integer, parameter :: k = 2
        type(DPlaneRotation) :: giv(k, 6)
        type(DUpdatableMatrix) :: mat
        class(RandomEngine), allocatable :: rng

        call init_xorshift_star(rng)
        x(:, :) = 0
        do j = 1, 7; do i = 1, min(j+k, 7)
            x(i, j) = random_real(rng, -1.0d0, 1.0d0)
        enddo; enddo
        call print_array(x, name='X')
        call setmatrix(mat, x)
        det1 = get_det(mat)

        r(:, :) = x(:, :)
        call qr_hessenberg(k, r, size(r, 1), giv)
        call print_array(r, name='R')
        call assert_equal(is_upper_triangular(r), .true.)
        call setmatrix(mat, r)
        det2 = get_det(mat)

        call assert_close(det1, det2, rtol=4 * epsilon(1.0d0))

        call identity(q)
        call rotate_q(q, size(q, 1), giv)
        call print_array(q, name='Q')
        !call matmul(transpose(q), q)

        x_rec = matmul(q, r)
        call print_array(x_rec, name='X_rec')
        call assert_close(x, x_rec, rtol=4 * epsilon(1.0d0))
    end subroutine

    subroutine test_qr_hess_complex()
        complex(c_double_complex) :: x(7, 7), r(7, 7), q(7, 7), x_rec(7, 7)
        type(ZExtRange) :: det1, det2
        integer :: i, j
        integer, parameter :: k = 2
        type(ZPlaneRotation) :: giv(k, 6)
        type(ZUpdatableMatrix) :: mat
        class(RandomEngine), allocatable :: rng

        call init_xorshift_star(rng)
        x(:, :) = 0
        do j = 1, 7; do i = 1, min(j+k, 7)
            x(i, j) = cmplx(random_real(rng, -1.0d0, 1.0d0), &
                            random_real(rng, -0.5d0, 0.5d0), c_double_complex)
        enddo; enddo
        call print_array(x, name='X')
        call setmatrix(mat, x)
        det1 = get_det(mat)

        r(:, :) = x(:, :)
        call qr_hessenberg(k, r, size(r, 1), giv)
        call print_array(r, name='R')
        call assert_equal(is_upper_triangular(r), .true.)
        call setmatrix(mat, r)
        det2 = get_det(mat)

        call assert_close(det1, det2, rtol=8 * epsilon(1.0d0))

        call identity(q)
        call rotate_q(q, size(q, 1), giv)
        call print_array(q, name='Q')
        !call matmul(transpose(q), q)

        x_rec = matmul(q, r)
        call print_array(x_rec, name='X_rec')
        call assert_close(x, x_rec, rtol=8 * epsilon(1.0d0))
    end subroutine

    subroutine test_qr_det_complex()
        complex(c_double_complex) :: x(4, 4) = reshape((/ &
            (564.0d0, +73.0d0), (586.0d0, 948.0d0), ( 24.0d0, 167.0d0), (660.0d0, 322.0d0), &
            (453.0d0, 536.0d0), (577.0d0, 631.0d0), (910.0d0, 301.0d0), (286.0d0, 319.0d0), &
            (478.0d0, 946.0d0), (510.0d0, 342.0d0), (316.0d0, 679.0d0), (902.0d0, 443.0d0), &
            (986.0d0, 791.0d0), (368.0d0, 546.0d0), (401.0d0, 820.0d0), ( 92.0d0, 525.0d0)  &
            /), (/ 4, 4 /))

        complex(c_double_complex) :: tau(4), work(16)
        type(ZExtRange) :: det
        integer :: info

        call ZGEQRF(4, 4, x(1,1), size(x, 1), tau(1), work(1), size(work), info)
        call check_info(info)

        det = qr_det(x, tau)
        call assert_close(limrange(det), (291942591691d0, -1034120419529d0), rtol=1d-12)
    end

    subroutine test_qr_rot()
        real(c_double) :: q(7, 7), qpart(3, 7), qref(3, 7)
        type(DPlaneRotation) :: giv(3, 6)
        class(RandomEngine), allocatable :: rng
        integer :: i, k

        write (0,*) 'test_qr_rot'
        call init_xorshift_star(rng)
        q(:, :) = random_orthogonal(rng, 7)
        call print_array(q, name='Q')

        qpart(1, :) = q(4, :)
        qpart(2, :) = q(6, :)
        qpart(3, :) = q(7, :)
        call q_diagonalize(qpart, size(qpart, 1), giv)
        call print_array(qpart, name="Q'")

        qref = merge(1, 0, reshape((/ ((i == k, k=1,3), i=1,7) /), (/ 3, 7 /)))
        call assert_close(abs(qpart), qref)
    end

    function random_orthogonal(rng, n) result(q)
        class(RandomEngine), intent(inout) :: rng
        integer, intent(in) :: n
        real(c_double), allocatable :: q(:, :)
        real(c_double) :: tau(n), work(n * n)
        integer :: i, j, info

        ! Random matrix
        allocate(q(n, n))
        do j = 1, n; do i = 1, n
            q(i, j) = random_real(rng)
        enddo; enddo

        call DGEQRF(n, n, q(1,1), size(q, 1), tau(1), work(1), size(work), info)
        call check_info(info)

        call DORGQR(n, n, n, q(1,1), size(q, 1), tau(1), work(1), size(work), info)
        call check_info(info)
    end

    subroutine test_qr_rot_complex()
        complex(c_double_complex) :: q(7, 7), qpart(3, 7), qref(3, 7)
        type(ZPlaneRotation) :: giv(3, 6)
        class(RandomEngine), allocatable :: rng
        integer :: i, k

        write (0,*) 'test_qr_rot_complex'
        call init_xorshift_star(rng)
        q(:, :) = random_unitary(rng, 7)
        call print_array(q, name='Q')

        qpart(1, :) = q(4, :)
        qpart(2, :) = q(6, :)
        qpart(3, :) = q(7, :)
        call q_diagonalize(qpart, size(qpart, 1), giv)
        call print_array(qpart, name="Q'")

        qref = merge(1, 0, reshape((/ ((i == k, k=1,3), i=1,7) /), (/ 3, 7 /)))
        call assert_close(abs(qpart), real(qref))
    end

    function random_unitary(rng, n) result(q)
        class(RandomEngine), intent(inout) :: rng
        integer, intent(in) :: n
        complex(c_double_complex), allocatable :: q(:, :)
        complex(c_double_complex) :: tau(n), work(n * n)
        integer :: i, j, info

        ! Random matrix
        allocate(q(n, n))
        do j = 1, n; do i = 1, n
            q(i, j) = cmplx(random_real(rng), random_real(rng), c_double_complex)
        enddo; enddo

        call ZGEQRF(n, n, q(1,1), size(q, 1), tau(1), work(1), size(work), info)
        call check_info(info)

        call ZUNGQR(n, n, n, q(1,1), size(q, 1), tau(1), work(1), size(work), info)
        call check_info(info)
    end

    subroutine test_qr_update()
        class(RandomEngine), allocatable :: rng
        real(c_double) :: x(7, 7), r(7, 7), q(7, 7), work(49), tau(7)
        real(c_double) :: qpart(3, 7), qsquare(3, 3), id(3, 3)
        type(DPlaneRotation) :: giv(3, 6)
        integer :: i, j, info

        call init_xorshift_star(rng)
        x(:, :) = 0
        do j = 1, 7; do i = 1, 7
            x(i, j) = random_real(rng, -1.0d0, 1.0d0)
        enddo; enddo
        call print_array(x, name='X')

        r(:, :) = x(:, :)
        call DGEQRF(7, 7, r(1,1), size(r, 1), tau(1), work(1), size(work), info)
        call check_info(info)

        q(:, :) = r(:, :)
        call DORGQR(7, 7, 7, q(1,1), size(q, 1), tau(1), work(1), size(work), info)
        call check_info(info)
        call print_array(q, name='Q')

        do j = 1, 7; do i = j+1, 7
            r(i, j) = 0
        enddo; enddo
        call print_array(r, name='R')
        call assert_close(matmul(q, r), x, rtol=4 * epsilon(1d0))

        qpart(:, :) = q((/3, 4, 7/), :)
        call q_diagonalize(qpart, size(qpart, 1), giv)

        call rotate_r(r, size(r, 1), giv)
        call print_array(r, name='H')

        call rotate_back_q(q, size(q, 1), giv)
        call assert_close(matmul(q, r), x, rtol=5 * epsilon(1d0))

        qsquare(:, :) = q((/3, 4, 7/), (/1, 2, 3/))
        call print_array(qsquare, name='Q part')
        call identity(id)
        call assert_close(abs(qsquare), id)
    end subroutine

    subroutine test_qr_update_complex()
        class(RandomEngine), allocatable :: rng
        complex(c_double_complex) :: x(7, 7), r(7, 7), q(7, 7), work(49), tau(7)
        complex(c_double_complex) :: qpart(3, 7), qsquare(3, 3), id(3, 3)
        type(ZPlaneRotation) :: giv(3, 6)
        integer :: i, j, info

        call init_xorshift_star(rng)
        x(:, :) = 0
        do j = 1, 7; do i = 1, 7
            x(i, j) = cmplx(random_real(rng, -1.0d0, 1.0d0), &
                            random_real(rng, -1.0d0, 1.0d0), c_double_complex)
        enddo; enddo
        call print_array(x, name='X')

        r(:, :) = x(:, :)
        call ZGEQRF(7, 7, r(1,1), size(r, 1), tau(1), work(1), size(work), info)
        call check_info(info)

        q(:, :) = r(:, :)
        call ZUNGQR(7, 7, 7, q(1,1), size(q, 1), tau(1), work(1), size(work), info)
        call check_info(info)
        call print_array(q, name='Q')

        do j = 1, 7; do i = j+1, 7
            r(i, j) = 0
        enddo; enddo
        call print_array(r, name='R')
        call assert_close(matmul(q, r), x, rtol=4 * epsilon(1d0))

        qpart(:, :) = q((/3, 4, 7/), :)
        call q_diagonalize(qpart, size(qpart, 1), giv)

        call rotate_r(r, size(r, 1), giv)
        call print_array(r, name='H')

        call rotate_back_q(q, size(q, 1), giv)
        call assert_close(matmul(q, r), x, rtol=4 * epsilon(1d0))

        qsquare(:, :) = q((/3, 4, 7/), (/1, 2, 3/))
        call print_array(qsquare, name='Q part')
        call identity(id)
        call assert_close(abs(qsquare), abs(id))
    end subroutine

end program
