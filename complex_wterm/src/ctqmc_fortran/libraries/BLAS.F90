module MBLAS
    implicit none

    interface
        pure function DDOT(N, ZX, INCX, ZY, INCY) result(RES)
            double precision :: RES
            double precision, intent(in) :: ZX, ZY
            integer, intent(in) :: N, INCX, INCY
        end
        pure function ZDOTC(N, ZX, INCX, ZY, INCY) result(RES)
            double complex :: RES
            double complex, intent(in) :: ZX, ZY
            integer, intent(in) :: N, INCX, INCY
        end
    end interface

    interface
        subroutine DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
                         BETA, C, LDC)
            character, intent(in) :: TRANSA, TRANSB
            integer, intent(in) :: M, N, K, LDA, LDB, LDC
            double precision, intent(in) :: ALPHA, BETA, A, B
            double precision, intent(inout) :: C
        end
        subroutine ZGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
                         BETA, C, LDC)
            character, intent(in) :: TRANSA, TRANSB
            integer, intent(in) :: M, N, K, LDA, LDB, LDC
            double complex, intent(in) :: ALPHA, BETA, A, B
            double complex, intent(inout) :: C
        end
    end interface

    interface
        subroutine DGER(M, N, ALPHA, X, INCX, Y, INCY, A, LDA)
            integer, intent(in) :: M, N, INCX, INCY, LDA
            double precision, intent(in) :: ALPHA, X, Y
            double precision, intent(inout) :: A
        end
        subroutine ZGERU(M, N, ALPHA, X, INCX, Y, INCY, A, LDA)
            integer, intent(in) :: M, N, INCX, INCY, LDA
            double complex, intent(in) :: ALPHA, X, Y
            double complex, intent(inout) :: A
        end
    end interface

    interface
        subroutine DROT(n, dx, incx, dy, incy, c, s)
            integer, intent(in) :: n, incx, incy
            double precision, intent(inout) :: dx, dy
            double precision, intent(in) :: c, s
        end
        subroutine ZROT(n, dx, incx, dy, incy, c, s)
            integer, intent(in) :: n, incx, incy
            double complex, intent(inout) :: dx, dy
            double precision, intent(in) :: c
            double complex, intent(in) :: s
        end
    end interface

end module
