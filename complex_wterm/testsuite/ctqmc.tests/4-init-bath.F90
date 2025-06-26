! compile with
! gfortran -I /tmp/w2dynamics/lib/ -L /tmp/w2dynamics/lib/  4-init-bath.F90
! /tmp/w2dynamics/lib/tmp/w2dynamics/src/maxent/MaximumEntropy.o -llapack -lblas
! /tmp/w2dynamics/bin/maxent/CMakeFiles/MAXENTLIB.dir/

PROGRAM test4

   use MRandomNumbers
   use MBath
   IMPLICIT NONE

   type(TBath) :: bath

   integer, parameter :: Nftau = 250
   integer, parameter :: nbands = 3
   real(c_double), parameter :: beta = 10.0

   real(c_double) :: ftau_data(nbands * 2 * nbands * 2 * Nftau)
   complex(c_double_complex) :: ftau_full(nbands, 2, nbands, 2, Nftau)
   integer :: ib, is, i

   real(c_double) :: maxdev

   logical :: complexHyb
   class(ZHybridization), allocatable :: zhybr
   class(DHybridization), allocatable :: dhybr

   open (newunit=ib, file="ftau.dat", form="formatted", status="old")
   do is = 1, size(ftau_data)
   read (ib, *) ftau_data(is)
   end do
   close (ib)

   ftau_full = reshape(ftau_data, [nbands, 2, nbands, 2, Nftau])

   do i = 1,2
      if(i==1) complexHyb=.False.
      if(i==2) complexHyb=.True.

      !allocate(bath)

      !if(complexHyb)then
         !!!allocate(bath, source = zbath_create())
         !call zinit_Bath(bath, nbands, ftau_diag, Ftau_full, nftau, beta)
         !write(*,*) "zinit_Bath()"
      !else
         !!allocate(bath, source = zbath_create())
         !call dinit_Bath(bath, nbands, ftau_diag, Ftau_full, nftau, beta)
         !write(*,*) "dinit_Bath()"
      !endif

      if (complexHyb) then
        call init_linear_hybridization(zhybr, Ftau_full, beta, .true.)
        call init_bath(bath, zhybr)
      else
        call init_linear_hybridization(dhybr, real(Ftau_full), beta, .true.)
        call init_bath(bath, dhybr)
      endif

      !call init_Bath(bath, nbands, Ftau_full, nftau, beta, complexHyb)

      !write(*,*) "calling init bath..."
      !write(*,*) "nbands", nbands
      !write(*,*) "shape(ftau_diag) ", shape(ftau_diag)
      !write(*,*) "shape(ftau_full) ", shape(ftau_full)
      !write(*,*) "nftau", nftau
      !write(*,*) "beta", beta

      maxdev = abs(get_nbands(bath) - nbands)
      if(maxdev.gt.2d-14)then
         write(*,*) "maxdev3", maxdev
         stop 1
      endif

      maxdev = abs(get_beta(bath) - beta)
      if(maxdev.gt.2d-14)then
         write(*,*) "maxdev4", maxdev
         stop 1
      endif

      maxdev = abs(limrange(get_det(bath))) - 1.0d0
      if(maxdev.gt.2d-14)then
         write(*,*) "maxdev5", maxdev
         stop 1
      endif

      !maxdev = abs(bath%det%sign - 1.0d0)
      !if(maxdev.gt.2d-14)then
      !write(*,*) "maxdev6", maxdev
      !stop 1
      !endif

      maxdev = abs(size(bath))
      if(maxdev.gt.2d-14)then
         write(*,*) "maxdev7", maxdev
         stop 1
      endif

      !!!! TODO: this does not yet work...
      !!maxdev = abs(shape(bath%bathcreators) - 1000)
      !!if(maxdev.gt.2d-14)then
      !!write(*,*) "maxdev8", maxdev
      !!stop 1
      !!endif



   enddo

END PROGRAM test4
