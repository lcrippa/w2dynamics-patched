module MExtRange
    use iso_c_binding, only: c_double, c_double_complex, c_long
    use MExtRangeD
    use MExtRangeZ
    implicit none

    interface ZExtRange
        module procedure zext_d
    end interface

    interface real
        module procedure zext_real
    end interface

    interface aimag
        module procedure zext_imag
    end interface

    interface assignment ( = )
        module procedure zassgn_d
    end interface
    public assignment ( = )

    interface operator ( * )
        module procedure dext_mul_z, dext_mul_complex
    end interface
    public operator ( * )

contains

    pure elemental function zext_d(x) result(r)
        type(DExtRange), intent(in) :: x
        type(ZExtRange) :: r

        r = ZExtRange(cmplx(fraction(x), kind=c_double_complex), exponent(x))
    end

    pure elemental function zext_dx(x) result(r)
        real(c_double), intent(in) :: x
        type(ZExtRange) :: r

        r = ZExtRange(cmplx(x, kind=c_double_complex))
    end

    pure elemental function zext_real(x) result(r)
        type(ZExtRange), intent(in) :: x
        type(DExtRange) :: r

        r = DExtRange(real(fraction(x)), exponent(x))
    end

    pure elemental function zext_imag(x) result(r)
        type(ZExtRange), intent(in) :: x
        type(DExtRange) :: r

        r = DExtRange(aimag(fraction(x)), exponent(x))
    end

    pure subroutine zassgn_d(self, x)
        type(ZExtRange), intent(out) :: self
        type(DExtRange), intent(in) :: x

        self = ZExtRange(x)
    end

    pure subroutine zassgn_dx(self, x)
        type(ZExtRange), intent(out) :: self
        real(c_double), intent(in) :: x

        self = ZExtRange(x)
    end

    pure elemental function dext_mul_complex(x, y) result(r)
        type(DExtRange), intent(in) :: x
        complex(c_double_complex), intent(in) :: y
        type(ZExtRange) :: r

        r = ZExtRange(fraction(x) * y, exponent(x))
    end

    pure elemental function dext_mul_z(x, y) result(r)
        type(DExtRange), intent(in) :: x
        type(ZExtRange), intent(in) :: y
        type(ZExtRange) :: r

        r = ZExtRange(fraction(x) * fraction(y), exponent(x) + exponent(y))
    end

end module
