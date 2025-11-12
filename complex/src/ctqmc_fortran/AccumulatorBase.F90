module MAccumulatorBase
    use iso_c_binding, only: c_int64_t
    implicit none
    private

    !> Accumulator type tracking the number of observations
    integer, parameter, public :: ACC_NULL = 0

    !> Accumulator type tracking the sample mean
    integer, parameter, public :: ACC_MEAN = 1

    !> Accumulator type tracking a stack of sample variances
    integer, parameter, public :: ACC_BLOCKS = 3

end module
