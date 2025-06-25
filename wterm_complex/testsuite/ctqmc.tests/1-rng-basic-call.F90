PROGRAM MTRNGTEST1
  use iso_c_binding
  USE MRandomNumbers
  IMPLICIT NONE

  INTEGER, PARAMETER    :: seed = 42, realiter = 10000000, intiter = 10000000, min = 8, max = 16
  real(c_double) :: res
  INTEGER               :: ires, i
  class(RandomEngine), allocatable :: rng

  call init_mersenne_twister(rng, seed)
  do i = 1, realiter
     res = random_real(rng)
     if (res < 0.0d0 .or. res > 1.0d0) then
        write (*, *) "rand() out of range"
        call abort()
     end if
  end do
  do i = 1, intiter
     ires = random_integer(rng, min, max)
     if (ires < min .or. ires > max) then
        write (*, *) "randint() out of range"
        call abort()
     end if
  end do
  write (*,*) "success"

END PROGRAM MTRNGTEST1
