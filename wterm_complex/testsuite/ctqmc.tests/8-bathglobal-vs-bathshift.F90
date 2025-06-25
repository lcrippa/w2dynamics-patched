! compile with
! gfortran -I /tmp/w2dynamics/lib/ -L /tmp/w2dynamics/lib/ 7-bathglobal-vs-bathshift.F90
! /tmp/w2dynamics/lib/tmp/w2dynamics/src/maxent/MaximumEntropy.o -llapack -lblas
! /tmp/w2dynamics/bin/maxent/CMakeFiles/MAXENTLIB.dir/

PROGRAM test8
   use MRandomNumbers
   use testing
   use MBath
   use MCommon
   implicit none

   type(TBath) :: bath
   type(TBathMove) :: DBathMove, DBathMove2
   type(TBathShiftMove) :: DBathShiftMove

   integer, parameter :: Nftau = 501
   integer, parameter :: norbitals = 2
   real(c_double), parameter :: beta = 30.0d0

   complex(c_double_complex), allocatable :: ftau_full(:, :, :, :, :)
   integer :: i

   real(c_double) :: rtemp

   integer :: Noper1, Noper2

   !integer :: ref_wrapped_crea(nbands, 2), ref_wrapped_annh(nbands, 2)
   type(TBathOperator), allocatable :: crea_orig(:), annh_orig(:)
   type(TBathOperator), allocatable :: crea_shifted(:), annh_shifted(:)
   real(c_double) :: deltatau

   logical :: complexHyb
   class(ZHybridization), allocatable :: zhybr
   class(DHybridization), allocatable :: dhybr

   do i = 1,2

      if(i==1)then
         complexHyb=.False.
         call read_ftau("ftau2.dat", norbitals, beta, ftau_full, Nftau)
      endif
      if(i==2)then
         complexHyb=.True.
         call read_ftau("ftau3.dat", norbitals, beta, ftau_full, Nftau)
      endif

      call read_bathconfig("b11.dat", crea_orig, annh_orig, Noper1)
      call read_bathconfig("b12.dat", crea_shifted, annh_shifted, Noper2)

      deltatau = 4.8026006286193264

      !if(Nftau.ne.Nftau_check)then
         !write(*,*) "Nftau ", Nftau
         !write(*,*) "Nftau_check ", Nftau_check
         !write(*,*) "Nftau_check .ne. Nftau!"
         !stop 1
      !endif

      !!!!! init modules
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
         !allocate(DBathMove2, source = zinit_bath_move(bath))
      !else
         !write(*,*) "dinit_bath_move(...):"
         !allocate(DBathMove, source = dinit_bath_move(bath))
         !allocate(DBathMove2, source = dinit_bath_move(bath))
      !endif

      call init_bath_move(DBathMove, bath)
      call init_bath_move(DBathMove2, bath)

      !!! perform a global bath calculation
      call propose_bath_move(DBathMove, crea_orig, annh_orig)
      rtemp = weight_ratio(DBathMove)

      !!! write this configuration as accepted configuration
      call accept_update(DBathMove)
      call verify_bath(bath, 1d-12)

      !!!! now the fast update
      !if(complexHyb)then
         !write(*,*) "zinit_bath_shift_move(...):"
         !allocate(DBathShiftMove, source = zinit_bath_shift_move(bath))
      !else
         !write(*,*) "dinit_bath_shift_move(...):"
         !allocate(DBathShiftMove, source = dinit_bath_shift_move(bath))
      !endif

      call init_bath_shift_move(DBathShiftMove, bath)

      call propose_bath_shift_move(DBathShiftMove, deltatau)

      call assert_close(weight_ratio(DBathShiftMove), 1.0d0)

      call accept_update(DBathShiftMove)

      call verify_bath(bath, 1d-12)

      write(*,*) "destroy stuff:"

      write(*,*) "allocated(crea_orig) ", allocated(crea_orig)
      write(*,*) "allocated(crea_shifted) ", allocated(crea_shifted)
      write(*,*) "allocated(annh_orig) ", allocated(annh_orig)
      write(*,*) "allocated(annh_shifted) ", allocated(annh_shifted)

   enddo

END PROGRAM test8
