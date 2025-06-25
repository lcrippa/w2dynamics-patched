! compile with
! gfortran -I /tmp/w2dynamics/lib/ -L /tmp/w2dynamics/lib/ 6-bathglobal-vs-bathadd.F90
! /tmp/w2dynamics/lib/tmp/w2dynamics/src/maxent/MaximumEntropy.o -llapack -lblas
! /tmp/w2dynamics/bin/maxent/CMakeFiles/MAXENTLIB.dir/

PROGRAM test22
   use MRandomNumbers
   use testing
   use MBath

   IMPLICIT NONE

   type(TBath) :: bath
   type(TBathMove) :: DBathMove
   type(TBathAddMove)  :: DBathAddMove

   integer, parameter :: Nftau = 501
   integer, parameter :: norbitals = 1
   real(c_double), parameter :: beta = 5.0d0

   complex(c_double_complex), allocatable :: ftau_full(:, :, :, :, :)

   real(c_double) :: detratfast

   real(c_double) :: rtemp

   integer :: Noper1, Noper2
   type(TBathOperator), allocatable :: crea1(:), annh1(:), crea2(:), annh2(:)
   type(TBathOperator), allocatable :: crea_tmp(:), annh_tmp(:)

   logical :: complexHyb
   class(ZHybridization), allocatable :: zhybr

   complex(c_double_complex) :: determ, sig, tr
   real(c_double) :: permsig

      complexHyb=.True.
      call read_ftau("ftau4.dat", norbitals, beta, ftau_full, Nftau)

      call read_bathconfig_highprec("b15.dat", crea1, annh1, noper1, tr, determ, sig, permsig)
      call read_bathconfig_highprec("b16.dat", crea2, annh2, noper2, tr, determ, sig, permsig)

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

      !!!! init modules
      call init_linear_hybridization(zhybr, Ftau_full, beta, .true.)
      call init_bath(bath, zhybr)

      call init_bath_move(DBathMove, bath)

      !!! perform a global bath calculation
      call propose_bath_move(DBathMove, crea1, annh1)
      rtemp = weight_ratio(DBathMove)

      !!! write this configuration as accepted configuration
      call accept_update(DBathMove)
      call debug_print_bathconfig_highprec(bath)
      call verify_bath(bath, 1d-12)

      !!! now the fast update
      call init_bath_add_move(DBathAddMove, bath)

      !!! these are the operators added
      crea_tmp = [ crea2(2) ]
      annh_tmp = [ annh2(1) ]
      call propose_bath_add_move(DBathAddMove, crea_tmp, annh_tmp)
      detratfast = weight_ratio(DBathAddMove)
      write(*,*) "detratfast", detratfast

      !!! preform global update with full crea/annh arrays
      call propose_bath_move(DBathMove, crea2, annh2)
      rtemp = weight_ratio(DBathMove)
      write(*,*) "rtemp", rtemp

      !!! call DBathMove%accept_update()
      !!! this must not be there, otherwise the fast update
      !!! will use the wrong hybridization matrix!
      call assert_close(detratfast, rtemp, rtol=1d-12)

      !!!! now generate and compare the fast matrices
      call accept_update(DBathAddMove)

      !!!! compare what was just written in bath,
      !!!! to what is still stored as correct/global result in DBathMove

      call verify_bath(bath, 1d-12)

      write(*,*) "destroy stuff:"

      deallocate(ftau_full)

      !!! TODO: it is maybe not good,
      !!! that the propose functions deallocate these arrays implicitly...
      write(*,*) "allocated(crea1) ", allocated(crea1)
      write(*,*) "allocated(crea2) ", allocated(crea2)
      write(*,*) "allocated(crea_tmp) ", allocated(crea_tmp)

END PROGRAM test22
