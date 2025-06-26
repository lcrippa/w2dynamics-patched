program trng_test
   use MRandomNumbers
   use iso_c_binding, only: c_double, c_int
   implicit none

   integer, parameter :: seed=12345, realiter=5000000, intiter=5000000, cmpiter=5000000, min=100, max=110
   integer            :: i
   class(RandomEngine), allocatable :: rng
   real(c_double)     :: dres
   integer(c_int)     :: ires

   call init_mersenne_twister(rng, seed)

   dres = -1.0_c_double
   do i = 1, realiter
      dres = random_real(rng)
!     write (*, *) "d ", dres
      if (dres < 0.0_c_double .or. dres >= 1.0_c_double) then
         write(*, *) "dres out of range"
         call abort()
      end if
   end do

   ires = -1_c_int
   do i = 1, intiter
      ires = random_integer(rng, min, max)
!     write (*, *) "i ", ires
      if (ires < min .or. ires > max) then
         write (*, *) "ires out of range"
         call abort()
      end if
   end do

end program trng_test
