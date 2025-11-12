module MLAPACK
    implicit none

    !> LU decomposition with partial pivoting
    interface
        subroutine DGETRF(M, N, A, LDA, IPIV, INFO)
            integer, intent(in)  :: LDA, M, N
            integer, intent(out) :: INFO, IPIV
            double precision, intent(inout) :: A
        end
        subroutine ZGETRF(M, N, A, LDA, IPIV, INFO)
            integer, intent(in)  :: LDA, M, N
            integer, intent(out) :: INFO, IPIV
            double complex, intent(inout) :: A
        end
    end interface

    !> Matrix from LU decomposition
    interface
        subroutine DGETRI(N, A, LDA, IPIV, WORK, LWORK, INFO)
            integer, intent(in) :: N, LDA, IPIV, LWORK
            integer, intent(out) :: INFO
            double precision, intent(inout) :: A, WORK(LWORK)
        end
        subroutine ZGETRI(N, A, LDA, IPIV, WORK, LWORK, INFO)
            integer, intent(in)  :: N, LDA, IPIV, LWORK
            integer, intent(out) :: INFO
            double complex, intent(inout) :: A, WORK(LWORK)
        end
    end interface

    !> LU decomposition with full pivoting
    interface
        subroutine DGETC2(n, a, lda, ipiv, jpiv, info)
            integer, intent(in) :: lda, n
            integer, intent(out) :: info, ipiv, jpiv
            double precision, intent(inout) ::  a
        end subroutine
        subroutine ZGETC2(n, a, lda, ipiv, jpiv, info)
            integer, intent(in) :: lda, n
            integer, intent(out) :: info, ipiv, jpiv
            double complex, intent(inout) ::  a
        end subroutine
    end interface

    !> Compute norm of matrix ('1', '2', 'I', 'M')
    interface
        function DLANGE(norm, m, n, a, lda, work)
            double precision :: dlange
            character(len=1), intent(in) :: norm
            integer, intent(in) :: lda, m, n
            double precision, intent(in) :: a
            double precision, intent(out) :: work
        end
        function ZLANGE(norm, m, n, a, lda, work)
            double precision :: zlange
            character(len=1), intent(in) :: norm
            integer, intent(in) :: lda, m, n
            double complex, intent(in) :: a
            double precision, intent(out) :: work
        end
    end interface

    !> Eigendecomposition of real symmetric matrix
    interface SYEV
        subroutine DSYEV(jobz, uplo, n, a, lda, w, work, lwork, info)
            character(len=1), intent(in) :: jobz, uplo
            integer, intent(in) :: lda, lwork, n
            integer, intent(out) :: info
            double precision, intent(inout) :: a
            double precision, intent(out) :: w
            double precision, intent(out) :: work
        end
    end interface

    !> Eigendecomposition of Hermitian matrix
    interface
        subroutine ZHEEV(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
            character(len=1), intent(in) :: jobz, uplo
            integer, intent(in) :: lda, lwork, n
            integer, intent(out) :: info
            double complex, intent(inout) :: a
            double precision, intent(out) :: w, rwork
            double complex, intent(out) :: work
        end
    end interface
    public :: DHEEV

    !> QR decomposition
    interface
        subroutine DGEQRF(m, n, a, lda, tau, work, lwork, info)
            integer, intent(in) :: lda, lwork, m, n
            integer, intent(out) :: info
            double precision, intent(inout) :: a
            double precision, intent(out) :: tau, work
        end subroutine
        subroutine ZGEQRF(m, n, a, lda, tau, work, lwork, info)
            integer, intent(in) :: lda, lwork, m, n
            integer, intent(out) :: info
            double complex, intent(inout) :: a
            double complex, intent(out) :: tau, work
        end subroutine
    end interface

    !> QR decomposition with column pivoting
    interface
        subroutine DGEQP3(m, n, a, lda, jpvt, tau, work, lwork, info)
            integer, intent(in) :: lda, lwork, m, n
            integer, intent(out) :: info
            integer, intent(inout) :: jpvt
            double precision, intent(inout) :: a
            double precision, intent(out) :: tau, work
        end subroutine
        subroutine ZGEQP3(m, n, a, lda, jpvt, tau, work, lwork, rwork, info)
            integer, intent(in) :: lda, lwork, m, n
            integer, intent(out) :: info
            integer, intent(inout) :: jpvt
            double complex, intent(inout) :: a
            double complex, intent(out) :: tau, work
            double precision, intent(out) :: rwork
        end subroutine
    end interface
    public DGEQP3X, ZGEQP3X

    !> Extract matrix Q from QR decomposition
    interface
        subroutine DORGQR(m, n, k, a, lda, tau, work, lwork, info)
            integer, intent(in) :: k, lda, lwork, m, n
            integer, intent(out) :: info
            double precision, intent(inout) :: a
            double precision, intent(in) :: tau
            double precision, intent(out) :: work
        end subroutine
        subroutine ZUNGQR(m, n, k, a, lda, tau, work, lwork, info)
            integer, intent(in) :: k, lda, lwork, m, n
            integer, intent(out) :: info
            double complex, intent(inout) :: a
            double complex, intent(in) :: tau
            double complex, intent(out) :: work
        end subroutine
    end interface

    !> Solve triangular system
    interface
        subroutine DTRTRS(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info)
            character(len=1), intent(in) :: diag, trans, uplo
            integer, intent(in) :: lda, ldb, n, nrhs
            integer, intent(out) :: info
            double precision, intent(in) :: a
            double precision, intent(inout) :: b
        end subroutine
        subroutine ZTRTRS(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info)
            character(len=1), intent(in) :: diag, trans, uplo
            integer, intent(in) :: lda, ldb, n, nrhs
            integer, intent(out) :: info
            double complex, intent(in) :: a
            double complex, intent(inout) :: b
        end subroutine
    end interface

    !> Generate plane rotation
    interface
        pure subroutine DLARTG(f, g, c, s, r)
            double precision, intent(in) :: f, g
            double precision, intent(out) :: c, s, r
        end subroutine
        pure subroutine ZLARTG(f, g, c, s, r)
            double complex, intent(in) :: f, g
            double precision, intent(out) :: c
            double complex, intent(out) :: s, r
        end subroutine
    end interface

    public :: check_info
contains
    subroutine DHEEV(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
        character(len=1), intent(in) :: jobz, uplo
        integer, intent(in) :: lda, lwork, n
        integer, intent(out) :: info
        double precision, intent(inout) :: a
        double precision, intent(out) :: w, rwork
        double precision, intent(out) :: work

        call DSYEV(jobz, uplo, n, a, lda, w, work, lwork, info)
    end

    subroutine DGEQP3X(m, n, a, lda, jpvt, tau, work, lwork, rwork, info)
        integer, intent(in) :: lda, lwork, m, n
        integer, intent(out) :: info
        integer, intent(inout) :: jpvt
        double precision, intent(inout) :: a
        double precision, intent(out) :: tau, work
        double precision, intent(out) :: rwork

        call DGEQP3(m, n, a, lda, jpvt, tau, work, lwork, info)
    end subroutine

    subroutine ZGEQP3X(m, n, a, lda, jpvt, tau, work, lwork, rwork, info)
        integer, intent(in) :: lda, lwork, m, n
        integer, intent(out) :: info
        integer, intent(inout) :: jpvt
        double complex, intent(inout) :: a
        double complex, intent(out) :: tau, work
        double precision, intent(out) :: rwork

        call ZGEQP3(m, n, a, lda, jpvt, tau, work, lwork, rwork, info)
    end subroutine

    subroutine check_info(info)
        integer, intent(in) :: info

        if (info == 0) &
            return
        write (0, *) 'LAPACK ERROR #', info
        error stop 'LAPACK ERROR'
    end
end module
