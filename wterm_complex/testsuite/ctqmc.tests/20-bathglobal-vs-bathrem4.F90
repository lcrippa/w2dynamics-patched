! compile with
! gfortran -I /tmp/w2dynamics/lib/ -L /tmp/w2dynamics/lib/ 7-bathglobal-vs-bathrem.F90
! /tmp/w2dynamics/lib/tmp/w2dynamics/src/maxent/MaximumEntropy.o -llapack -lblas
! /tmp/w2dynamics/bin/maxent/CMakeFiles/MAXENTLIB.dir/

PROGRAM test7

   use MRandomNumbers
   use testing
   use MBath
   IMPLICIT NONE

   type(TBath) :: bath
   type(TBathMove) :: DBathMove
   type(TBathRemMove)  :: DBathRemMove

   integer, parameter :: Nftau = 501
   integer, parameter :: norbitals = 2
   real(c_double), parameter :: beta = 30.0d0

   complex(c_double_complex), allocatable :: ftau_full(:, :, :, :, :)
   integer :: i

   real(c_double) :: detratfast

   real(c_double) :: rtemp

   integer :: Noper1
   type(TBathOperator), allocatable :: crea1(:), annh1(:)
   integer :: remidx(4)

   logical :: complexHyb
   class(ZHybridization), allocatable :: zhybr
   class(DHybridization), allocatable :: dhybr

   do i = 1,2

      if(i==1)then
         complexHyb=.False.
         call read_ftau("ftau2.dat", norbitals, beta, ftau_full, Nftau)
         call init_linear_hybridization(dhybr, real(Ftau_full), beta, .true.)
         call init_bath(bath, dhybr)
      endif
      if(i==2)then
         complexHyb=.True.
         call read_ftau("ftau3.dat", norbitals, beta, ftau_full, Nftau)
         call init_linear_hybridization(zhybr, Ftau_full, beta, .true.)
         call init_bath(bath, zhybr)
      endif

      call read_bathconfig("b14.dat", crea1, annh1, Noper1)

      !if(Nftau.ne.Nftau_check)then
         !write(*,*) "Nftau ", Nftau
         !write(*,*) "Nftau_check ", Nftau_check
         !write(*,*) "Nftau_check .ne. Nftau!"
         !stop 1
      !endif

      write(*,*) "ftau_full(1, 1, 1, 1, 1) ", ftau_full(1, 1, 1, 1, 1)
      write(*,*) "ftau_full(1, 1, 1, 2, 1) ", ftau_full(1, 1, 1, 2, 1)
      write(*,*) "ftau_full(1, 2, 1, 1, 1) ", ftau_full(1, 2, 1, 1, 1)
      write(*,*) "ftau_full(1, 2, 1, 2, 1) ", ftau_full(1, 2, 1, 2, 1)

      write(*,*) "ftau_full(1, 1, 1, 1, Nftau) ", ftau_full(1, 1, 1, 1, Nftau)
      write(*,*) "ftau_full(1, 1, 1, 2, Nftau) ", ftau_full(1, 1, 1, 2, Nftau)
      write(*,*) "ftau_full(1, 2, 1, 1, Nftau) ", ftau_full(1, 2, 1, 1, Nftau)
      write(*,*) "ftau_full(1, 2, 1, 2, Nftau) ", ftau_full(1, 2, 1, 2, Nftau)

      !!!!! init modules
      !allocate(bath)

      !if(complexHyb)then
         !write(*,*) "zinit_Bath(...):"
         !call zinit_Bath(bath, norbitals, ftau_diag, Ftau_full, nftau, beta)
      !else
         !write(*,*) "dinit_Bath(...):"
         !call dinit_Bath(bath, norbitals, ftau_diag, Ftau_full, nftau, beta)
      !endif

      !if(complexHyb)then
         !write(*,*) "zinit_bath_move(...):"
         !allocate(DBathMove, source = zinit_bath_move(bath))
      !else
         !write(*,*) "dinit_bath_move(...):"
         !allocate(DBathMove, source = dinit_bath_move(bath))
      !endif

      call init_bath_move(DBathMove, bath)

      !!! perform a global bath calculation
      !write(*,*) "shape(crea1) ", shape(crea1)
      !write(*,*) "shape(annh1) ", shape(annh1)
      call propose_bath_move(DBathMove, crea1, annh1)
      rtemp = weight_ratio(DBathMove)

      !!! write this configuration as accepted configuration
      call accept_update(DBathMove)
      call debug_print_bathconfig_highprec(bath)
      call verify_bath(bath, 1d-11)

      !!!! now the fast update
      !if(complexHyb)then
         !write(*,*) "zinit_bath_rem_move(...):"
         !allocate(DBathRemMove, source = zinit_bath_rem_move(bath))
      !else
         !write(*,*) "dinit_bath_rem_move(...):"
         !allocate(DBathRemMove, source = dinit_bath_rem_move(bath))
      !endif

      call init_bath_rem_move(DBathRemMove, bath)

      remidx(:) = (/ 19, 21, 21, 22 /)
      call propose_bath_rem_move(DBathRemMove, remidx)

      detratfast = weight_ratio(DBathRemMove)
      write(*,*) "detratfast", detratfast
      !det = DBathRemMove%det

      !!!! now generate and compare the fast matrices
      call accept_update(DBathRemMove)

      call verify_bath(bath, 1d-8)

      write(*,*) "destroy stuff:"
   enddo

END PROGRAM test7
