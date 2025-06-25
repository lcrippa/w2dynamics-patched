! compile with
! gfortran -I /tmp/w2dynamics/lib/ -L /tmp/w2dynamics/lib/ 7-bathglobal-vs-bathshift.F90
! /tmp/w2dynamics/lib/tmp/w2dynamics/src/maxent/MaximumEntropy.o -llapack -lblas
! /tmp/w2dynamics/bin/maxent/CMakeFiles/MAXENTLIB.dir/

PROGRAM test8
   use MRandomNumbers
   use MBath
   use MCommon
   implicit none

   type(TBath) :: bath
   type(TBathMove) :: DBathMove
   type(TBathReplaceMove)  :: DBathReplaceMove

   integer, parameter :: Nftau = 250
   integer, parameter :: nbands = 3
   real(c_double), parameter :: beta = 20.0d0

   real(c_double) :: ftau_data(nbands * 2 * nbands * 2 * Nftau)
   complex(c_double_complex) :: ftau_full(nbands, 2, nbands, 2, Nftau)
   integer :: ib, is, i

   real(c_double) :: tmp, detrat

   integer, parameter :: noper=12
   type(TBathOperator), allocatable :: crea_orig(:), annh_orig(:)

   integer :: bath_repidx
   integer(kind(OpDummy)) :: bath_reptype
   real(c_double) :: wormtau
   integer :: wormorb, wormspin

   logical :: complexHyb
   class(ZHybridization), allocatable :: zhybr
   class(DHybridization), allocatable :: dhybr

   open (newunit=ib, file="ftau.dat", form="formatted", status="old")
   do is = 1, size(ftau_data)
   read (ib, *) ftau_data(is)
   end do
   close (ib)

   ftau_full = reshape(ftau_data, [nbands, 2, nbands, 2, Nftau])

   allocate(crea_orig(noper/2))
   allocate(annh_orig(noper/2))

   do i = 1,2
      if(i==1) complexHyb=.false.
      if(i==2) complexHyb=.true.

      !write(*,*) "ftau_full[1,1,1,1,1] ", ftau_full(1,1,1,1,1)
      !write(*,*) "ftau_full[1,2,1,2,1] ", ftau_full(1,2,1,2,1)
      !write(*,*) "ftau_full[1,1,1,1,1] ", ftau_full(1,1,1,1,1)
      !write(*,*) "ftau_full[2,1,2,1,1] ", ftau_full(2,1,2,1,1)
      !write(*,*) "ftau_full[1,1,1,1,1] ", ftau_full(1,1,1,1,1)
      !write(*,*) "ftau_full[3,1,3,1,1] ", ftau_full(3,1,3,1,1)

      crea_orig(:)%tau = (/ &
            2.3324616784743228d0, 2.3531238900464815d0, &
            4.5269482543169488d0, 5.4605337932302476d0, &
            5.6985151965150038d0, 9.6263303207052804d0 /)
      crea_orig(:)%orb = (/1, 1, 1, 3, 1, 1/)
      crea_orig(:)%sp = (/1, 2, 2, 1, 2, 2/)

      annh_orig(:)%tau = (/ &
            2.3486133841301595d0, 2.6734175374000624d0, &
            4.5473011655859281d0, 5.5261995668127710d0, &
            5.7353097185072279d0, 9.9615768336485182d0 /)
      annh_orig(:)%orb = (/1, 3, 1, 3, 1, 1/)
      annh_orig(:)%sp = (/1, 2, 2, 1, 2, 2/)

      wormtau = 4.5269482543169488d0
      wormorb = 1
      wormspin = 1

      !!!! init modules
      !allocate(bath)

      !if(complexHyb)then
         !write(*,*) "zinit_Bath(...):"
         !call zinit_Bath(bath, nbands, ftau_diag, Ftau_full, nftau, beta)
      !else
         !write(*,*) "dinit_Bath(...):"
         !call dinit_Bath(bath, nbands, ftau_diag, Ftau_full, nftau, beta)
      !endif

      if (complexHyb) then
        call init_linear_hybridization(zhybr, Ftau_full, beta, .true.)
        call init_bath(bath, zhybr)
      else
        call init_linear_hybridization(dhybr, real(Ftau_full), beta, .true.)
        call init_bath(bath, dhybr)
      endif

      !if(complexHyb)then
         !write(*,*) "zinit_bath_move(...):"
         !allocate(DBathMove, source = zinit_bath_move(bath))
         !allocate(DBathMove2, source = zinit_bath_move(bath))
      !else
         !write(*,*) "dinit_bath_move(...):"
         !allocate(DBathMove, source = dinit_bath_move(bath))
         !allocate(DBathMove2, source = dinit_bath_move(bath))
      !endif

      call init_bath_move(DBathMove, bath)

      !!! perform a global bath calculation
      call propose_bath_move(DBathMove, crea_orig, annh_orig)
      tmp = weight_ratio(DBathMove)

      call accept_update(DBathMove)
      call verify_bath(bath, 1d-12)

      bath_reptype = OpCrea
      bath_repidx = 3

      !!! now the fast update
      call init_bathreplacemove(DBathReplaceMove, bath)

      call propose_bathreplacemove(DBathReplaceMove, bath_repidx, bath_reptype,&
         TBathOperator(tau=wormtau, orb=wormorb, sp=wormspin))

      detrat = weight_ratio(DBathReplaceMove)
      write(*,*) "detrat", detrat

      !!! write this configuration as accepted configuration
      call accept_update(DBathReplaceMove)

      ! FIXME: improve stability
      call verify_bath(bath, 4d-2)
   enddo

END PROGRAM test8
