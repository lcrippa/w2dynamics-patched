! compile with
! gfortran -I /tmp/w2dynamics/lib/ -L /tmp/w2dynamics/lib/ 6-bathglobal-vs-bathadd.F90
! /tmp/w2dynamics/lib/tmp/w2dynamics/src/maxent/MaximumEntropy.o -llapack -lblas
! /tmp/w2dynamics/bin/maxent/CMakeFiles/MAXENTLIB.dir/

PROGRAM test6

   use MRandomNumbers
   use testing
   use MBath

   IMPLICIT NONE

   type(TBath) :: bath
   type(TBathMove)  :: DBathMove
   type(TBathAddMove)  :: DBathAddMove

   integer, parameter :: Nftau = 501
   integer, parameter :: norbitals = 2
   real(c_double), parameter :: beta = 30.0d0

   complex(c_double_complex), allocatable :: ftau_full(:, :, :, :, :)
   integer :: i

   real(c_double) :: detratfast

   real(c_double) :: rtemp

   integer :: Noper1, Noper2
   type(TBathOperator), allocatable :: crea1(:), annh1(:), crea2(:), annh2(:)
   type(TBathOperator), allocatable :: crea_tmp(:), annh_tmp(:)
   class(ZHybridization), allocatable :: zhybr
   class(DHybridization), allocatable :: dhybr
   logical :: complexHyb

   do i = 1,2

      if(i==1)then
         complexHyb=.False.
         call read_ftau("ftau2.dat", norbitals, beta, ftau_full, Nftau)
      endif
      if(i==2)then
         complexHyb=.True.
         call read_ftau("ftau3.dat", norbitals, beta, ftau_full, Nftau)
      endif

      call read_bathconfig("b5.dat", crea1, annh1, Noper1)
      call read_bathconfig("b6.dat", crea2, annh2, Noper2)

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

     !!!! init modules
      !allocate(bath)

      !if(complexHyb)then
         !write(*,*) "zinit_Bath(...):"
         !call zinit_Bath(bath, norbitals, ftau_diag, Ftau_full, nftau, beta)
      !else
         !write(*,*) "dinit_Bath(...):"
         !call dinit_Bath(bath, norbitals, ftau_diag, Ftau_full, nftau, beta)
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
      !else
         !write(*,*) "dinit_bath_move(...):"
         !allocate(DBathMove, source = dinit_bath_move(bath))
      !endif
      call init_bath_move(DBathMove, bath)

      !!! perform a global bath calculation
      !write(*,*) "shape(crea1) ", shape(crea1)
      !write(*,*) "shape(annh1) ", shape(annh1)
      call propose_bath_move(DBathMove, crea1, annh1)
      call compute_weight(DBathMove)
      rtemp = weight_ratio(DBathMove)

      !!! write this configuration as accepted configuration
      call accept_update(DBathMove)
      call debug_print_bathconfig_highprec(bath)
      call verify_bath(bath, 1d-12)

      !!!! now the fast update
      !if(complexHyb)then
         !write(*,*) "zinit_bath_add_move(...):"
         !allocate(DBathAddMove, source = zinit_bath_add_move(bath))
      !else
         !write(*,*) "dinit_bath_add_move(...):"
         !allocate(DBathAddMove, source = dinit_bath_add_move(bath))
      !endif

      call init_bath_add_move(DBathAddMove, bath)

      !!! these are the operators added
      crea_tmp = [ crea2(4) ]
      annh_tmp = [ annh2(4) ]
      call propose_bath_add_move(DBathAddMove, crea_tmp, annh_tmp)
      call compute_weight(DBathAddMove)
      detratfast = weight_ratio(DBathAddMove)

      !det = DBathAddMove%det

      !!! preform global update with full crea/annh arrays
      call propose_bath_move(DBathMove, crea2, annh2)
      call compute_weight(DBathMove)
      rtemp = weight_ratio(DBathMove)

      !!! call DBathMove%accept_update()
      !!! this must not be there, otherwise the fast update
      !!! will use the wrong hybridization matrix!
      call assert_close(detratfast, rtemp, rtol=1d-10)

      !!!! now generate and compare the fast matrices
      call accept_update(DBathAddMove)

      !!!! compare what was just written in bath,
      !!!! to what is still stored as correct/global result in DBathMove
      call verify_bath(bath, 1d-10)

      write(*,*) "destroy stuff:"

      deallocate(ftau_full)

      !!! TODO: it is maybe not good,
      !!! that the propose functions deallocate these arrays implicitly...
      write(*,*) "allocated(crea1) ", allocated(crea1)
      write(*,*) "allocated(crea2) ", allocated(crea2)
      write(*,*) "allocated(crea_tmp) ", allocated(crea_tmp)

   enddo

END PROGRAM test6
