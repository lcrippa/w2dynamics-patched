program test_accumulators
    use iso_c_binding
    use MAccumulator
    use testing
    implicit none

    type(TAccumulator) :: mean, corr

    call accumulator_init(mean, .false., ACC_MEAN, 2_c_int64_t)
    call test_arith(mean)

    call accumulator_init(corr, .false., ACC_BLOCKS, 2_c_int64_t, 3_c_int64_t)
    call test_arith(corr)

    call accumulator_init(mean, .true., ACC_MEAN, 2_c_int64_t)
    call test_arithz(mean)

    call accumulator_init(corr, .true., ACC_BLOCKS, 2_c_int64_t, 3_c_int64_t)
    call test_arithz(corr)

contains
    subroutine test_arith(acc)
        type(TAccumulator), intent(inout) :: acc
        real(c_double) :: curr(2), cmean(2) !, cvar(2,2)
        real(c_double), pointer :: buf(:)
        integer :: i

        call accumulator_buffer(acc, buf)

        ! arithmetic progression
        do i = 0, 10000
            curr = (/ i, 10000 - i /)

            buf = buf + curr
            call accumulator_add(acc)
        end do

        if (accumulator_count(acc) /= 10001) &
            call terminate('wrong number of items')

        ! check if Gauss was right
        if (accumulator_has_mean(acc)) then
            call accumulator_mean(acc, cmean)
            call assert_close(cmean, (/ 5000.0d0, 5000.0d0 /), &
                              line=__LINE__)
        end if
    end subroutine

    subroutine test_arithz(acc)
        type(TAccumulator), intent(inout) :: acc
        complex(c_double_complex) :: curr(2), cmean(2) !, cvar(2,2)
        complex(c_double_complex), pointer :: buf(:)
        integer :: i

        call accumulator_buffer(acc, buf)

        ! arithmetic progression
        do i = 0, 10000
            curr = (/ cmplx(1.0d0*i, 0.1d0*i, c_double_complex), &
                      cmplx(10000d0 - i, 0.2d0*i, c_double_complex) /)

            buf = buf + curr
            call accumulator_add(acc)
        end do

        if (accumulator_count(acc) /= 10001) &
            call terminate('wrong number of items')

        ! check if Gauss was right
        if (accumulator_has_mean(acc)) then
            call accumulator_mean(acc, cmean)
            call assert_close(cmean, (/ (5000d0, 500d0), (5000d0, 1000d0) /), &
                              rtol=1d-13, line=__LINE__)
        end if
    end subroutine

end program test_accumulators
