module MComplexNumbers
    use iso_c_binding, only: c_int64_t, c_long, c_double, c_double_complex
    implicit none
    private

    interface conjg
        module procedure :: ddconjg
    end interface
    public conjg

    interface sign
        module procedure :: zsign
    end interface
    public sign

    interface scale
       module procedure :: zscale
    end interface
    public scale

    interface nan
        module procedure :: dnan, znan
    end interface
    public nan

contains
    pure elemental function zsign(x, y) result(r)
        complex(c_double_complex), intent(in) :: x, y
        complex(c_double_complex) :: r

        r = abs(x)
        if (y /= 0) then
            r = r * (y / abs(y))
        endif
    end

    pure elemental function ddconjg(x) result(r)
        real(c_double), intent(in) :: x
        real(c_double) :: r

        r = x
    end

    pure elemental function zscale(z, exp) result(r)
        complex(c_double_complex), intent(in) :: z
        integer(c_long), intent(in) :: exp
        complex(c_double_complex) :: r

        r = cmplx(scale(real(z), exp), scale(aimag(z), exp), &
                  kind=c_double_complex)
    end

    pure function dnan(mold) result(r)
        real(c_double), intent(in) :: mold
        real(c_double) :: r

        ! this relies on IEEE 751, which we always have
        r = transfer(-2251799813685248_c_int64_t, 1.0_c_double)
    end

    pure function znan(mold) result(r)
        complex(c_double_complex), intent(in) :: mold
        complex(c_double_complex) :: r

        ! this relies on IEEE 751, which we always have
        r = cmplx(nan(1.0d0), nan(1.0d0), kind=c_double_complex)
    end
end module
