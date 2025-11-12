#define ASSERT_CLOSE(a, b) call assert_close(a, b)

program test_extrange
    use iso_c_binding, only: c_double, c_double_complex
    use MExtRange
    use Testing
    implicit none

    call test_large(1.0d300)
    call test_large(1.0d-300)
    call test_large(-1.0d-300)

    call test_largez((1.0d300, -0.5d100))
    call test_largez((-1.0d-100, 0.5d-300))
    call test_largez((-1.0d100, 0.5d-300))
    call test_largez((-1.0d-300, 0.5d200))

    call test_norm()

    call test_frexp(1.0d-27)

    call test_add(1.0d0, -1.0d0)
    call test_add(14.0d300, 1.0d0)

    call test_addz((1.0d0, 0.5d300), (-1.0d-100, 1.0d1))

    call test_scale(1d300)

    call test_exp_large(1d0)
    call test_exp_large(600d0)
    call test_exp_large(-650d0)
    call test_exp_large(642.23984023048d0)

contains
    subroutine test_large(num)
        real(c_double) :: num
        type(DExtRange) :: x

        x = extrange(num)
        call dump(x)
        x = x * (-num)
        call dump(x)
        x = x / (-extrange(num))
        call dump(x)
        ASSERT_CLOSE(limrange(x), num)

        call assert_equal(isclose(DExtRange(0.0d0), DExtRange(0.0d0)), .true.)
    end subroutine

    subroutine test_largez(num)
        complex(c_double_complex) :: num
        type(ZExtRange) :: x

        x = extrange(num)
        call dump(x)
        x = x * (-num)
        call dump(x)
        x = x / (-extrange(num))
        call dump(x)
        ASSERT_CLOSE(limrange(x), num)

        x = conjg(extrange(num))
        ASSERT_CLOSE(limrange(x), conjg(num))
    end subroutine

    subroutine test_norm()
        real(c_double), parameter :: annoying = 0.5 + (0.5)**(-10)
        type(DExtRange) :: x

        integer :: i

        x = extrange(annoying)
        do i = 1, 1000
            x = x * annoying
        enddo
        do i = 1, 1000
            x = x / annoying
        enddo
        call assert_close(limrange(x), annoying)
    end subroutine

    subroutine test_frexp(num)
        real(c_double) :: num
        type(DExtRange) :: x

        x = extrange(num)
        call assert_equal(int(exponent(x)), int(exponent(num)))
        call assert_close(fraction(x), fraction(num))
    end subroutine

    subroutine test_add(x, y)
        real(c_double) :: x, y

        call assert_close(x + y, limrange(extrange(x) + extrange(y)))
        call assert_close(x + y, limrange(extrange(y) + extrange(x)))

        call assert_close(x - y, limrange(extrange(x) - extrange(y)))
        call assert_close(y - x, limrange(extrange(y) - extrange(x)))
    end

    subroutine test_addz(x, y)
        complex(c_double_complex) :: x, y

        call assert_close(x + y, limrange(extrange(x) + extrange(y)))
        call assert_close(x + y, limrange(extrange(y) + extrange(x)))

        call assert_close(x - y, limrange(extrange(x) - extrange(y)))
        call assert_close(y - x, limrange(extrange(y) - extrange(x)))
    end

    subroutine test_scale(x)
        real(c_double) :: x

        call assert_close(scale(x, 10), &
                          limrange(scale(extrange(x), 10)))
        call assert_close(scale(x, -40), &
                          limrange(scale(scale(extrange(x), -100040), 100000)))
    end

    subroutine test_exp_large(x)
        real(c_double) :: x
        type(DExtRange) :: xx, exp_xx

        ! XXX tighten this tolerance
        real(c_double), parameter :: rtol = 1d-10

        xx = x
        exp_xx = exp(x)
        call assert_close(exp(xx), exp_xx)
        call assert_close(exp(xx * 2d0), exp_xx * exp_xx, rtol=rtol)

        ! XXX tighten bounds
        call assert_close( &
                exp(xx * 1000d0), exp(xx * 500d0) * exp(xx * 500d0), rtol=rtol)
    end

end program
