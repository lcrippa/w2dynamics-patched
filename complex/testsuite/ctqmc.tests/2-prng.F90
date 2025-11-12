program test_prng
    use iso_c_binding
    use MRandomNumbers
    use testing
    implicit none

    call mt_ref
    call mt_std
    call mt_int
    call mt_int_large
    call mt_real

contains

    subroutine mt_std
        class(RandomEngine), allocatable :: rng
        integer(c_int32_t) :: i, num

        call init_mersenne_twister(rng)
        do i = 1, 10000
            num = rng_get_raw(rng)
        enddo

        ! Required by the C++ standard 28.5.6.3
        call assert_equal(num, -171307301)
    end subroutine

    subroutine mt_ref
        integer(c_int32_t), parameter :: ref_data(10) = (/ &
                -795755684,  581869302, -404620562, -708632711,  545404204, &
                -133711905, -372047867,  949333985,-1579004998, 1323567403 &
                /)
        class(RandomEngine), allocatable :: rng
        integer(c_int32_t) :: i, x(10)

        call init_mersenne_twister(rng)
        do i = 1, 10
            x(i) = rng_get_raw(rng)
        enddo
        call assert_equal(x, ref_data)
    end

    subroutine mt_int
        integer, parameter :: ref_data(20) = (/ &
                58, 9, 65, 60, 9, 69, 65, 15, 45, 22, &
                7, 39, 20, 13, 39, 71, 68, 71, 69, 69 /)
        class(RandomEngine), allocatable :: rng
        integer :: i, x(20)

        call init_mersenne_twister(rng)
        do i = 1, 20
            x(i) = random_integer(rng, 72);
        enddo
        call assert_equal(x, ref_data)
    end

    subroutine mt_int_large
        integer, parameter :: ref_data(40) = (/ &
                874802903,  145467325,  972586684,  1040313848, 237333496, &
                104733208,  587573641,  299035185,  202273606,  587209560, &
                1066098180, 1028115130, 1069942201, 1039054527, 169235752, &
                779363652,  1042166061, 1053458510, 1027750187, 117963156, &
                521168134,  859294615,  5136227,    452862732,  120757854, &
                983263532,  686940674,  1030247147, 540803682,  704096211, &
                856769327,  38345123,   387936480,  728786327,  428142225, &
                813617265,  795263923,  797932415,  509768251,  421150555 /)
        class(RandomEngine), allocatable :: rng
        integer :: i, x(40)

        call init_mersenne_twister(rng)
        do i = 1, 40
            x(i) = random_integer(rng, 1073741825)
        enddo
        call assert_equal(x, ref_data)
    end

    subroutine mt_real
        real(c_double), parameter :: ref_data(10) = (/ &
            0.1354770042967805d0, 0.8350085899945795d0, 0.9688677711242314d0, &
            0.2210340429827049d0, 0.3081670505070033d0, 0.5472205963678519d0, &
            0.1883819760471811d0, 0.9928813019178067d0, 0.9964613255480087d0, &
            0.9676949370105026d0 /)
        class(RandomEngine), allocatable :: rng
        integer :: i
        real(c_double) :: x(10)

        call init_mersenne_twister(rng)
        do i = 1, 10
            x(i) = random_real(rng)
        enddo
        call assert_close(x, ref_data)
    end

end program
