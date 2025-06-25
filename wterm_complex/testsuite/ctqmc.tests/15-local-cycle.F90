#define TEST_ERROR() call terminate(line=__LINE__)

program test
   use MBath
   use MDiagramMoves
   use MDiagram
   use MImpurity
   use MOperator
   use MPsi
   use MRandomNumbers
   use MStates
   use MSymmetryMoves
   use MWindow
   use MWormState
   use Testing
   use iso_c_binding
   implicit none

   logical, parameter :: b_offdiag = .true.

   integer, parameter :: Nftau = 501
   integer, parameter :: nbands = 2
   real(c_double), parameter :: beta = 30.0d0

   logical, parameter :: complexHyb=.false., complexLocTr=.false.

   complex(c_double_complex), allocatable :: ftau_full(:, :, :, :, :)
   real(c_double) :: u_matrix_data(16 * nbands**4)
   complex(c_double_complex) :: u_matrix(2 * nbands, 2 * nbands, 2 * nbands, 2 * nbands)
   real(c_double) :: muimp_data(4 * nbands**2)
   complex(c_double_complex) :: muimp_full(nbands, 2, nbands, 2)

   type(TQNSpec) :: qnspec
   class(RandomEngine), allocatable :: rng
   type(TStates) :: dstates
   type(TSymmetryMoves) :: sym_moves
   type(TDiagram) :: diagram
   type(TOperEigen) :: hloc_eb
   type(TPsis) :: psis_eb
   type(TImpurity) :: local
   type(TBath) :: bath
   type(TWindowState) :: window
   type(TWormState) :: worm
   integer :: i, ib, is, user_symmoves(0, 2 * nbands)
   integer :: TryAdd(NWormSectors), TryRem(NWormSectors)
   integer, parameter :: NSlide = 20, maxattempts = 100
   real(c_double) :: wormEta(NWormSectors)
   real(c_double), parameter :: Percentage4OperatorMove = 0.1, PercentageOuterMove = 0.1
   real(c_double), parameter :: PercentageTauShiftMove = 0.01, PercentageGlobalMove = 0.01
   real(c_double), parameter :: PercentageWormInsert = 0.0, PercentageWormReplace = 0.0

   real(c_double), parameter :: MINPROB = 5d-9, VERIFY_RTOL = 1e-7

   type(TAddMove) :: addmove, add4move
   type(TRemMove) :: remmove, rem4move
   type(TShiftMove) :: shiftmove
   type(TPermuteMove) :: permmove

   real(c_double) :: propprob1, propprob2, weightratio, rand, deltatau
   integer :: permcmp_outerbl, attempts
   integer, parameter :: RandPerm = 1, UserSym = 2, SpinFlip = 3, CAFlip = 4
   integer :: globtype, flavcaperm(0:2*nbands), outerbl
   logical :: allowed, do_spin_flip
   type(ZExtRange) :: trref
   type(TLocalOper), allocatable :: local_ops(:), newops(:)
   integer, allocatable :: local_remidx(:)
   class(ZHybridization), allocatable :: zhybr
   class(DHybridization), allocatable :: dhybr

   call init_mersenne_twister(rng)

   ! load U
   open (newunit=ib, file="u_matrix.dat", form="formatted", status="old")
   do is = 1, size(u_matrix_data)
      read (ib, *) u_matrix_data(is)
   end do
   close (ib)

   ! load muimp
   open (newunit=ib, file="muimp.dat", form="formatted", status="old")
   do is = 1, size(muimp_data)
      read (ib, *) muimp_data(is)
   end do
   close (ib)

   ! store U and muimp in correctly shaped arrays
   u_matrix = reshape(u_matrix_data, shape(u_matrix))
   muimp_full = reshape(muimp_data, shape(muimp_full))

   ! local initialization
   qnspec%Nt = .true.
   qnspec%Szt = .false.
   qnspec%Qzt = .false.
   call local_setup(.false.)

   outerbl = random_outer_state_by_weight(local, rng)
   call clear_impurity(local, outerbl)

   ! bath initialization
   if (complexHyb) then
      call read_ftau("ftau3.dat", nbands, beta, ftau_full, Nftau)
      call init_linear_hybridization(zhybr, Ftau_full, beta, .true.)
      call init_bath(bath, zhybr)
   else
      call read_ftau("ftau2.dat", nbands, beta, ftau_full, Nftau)
      call init_linear_hybridization(dhybr, real(Ftau_full), beta, .true.)
      call init_bath(bath, dhybr)
   endif

   ! worm initialization
   call enter_z_sector(worm)

   ! diagram initialization
   call init_diagram(diagram, local, bath, worm, b_offdiag)

   call init_window(window, diagram, beta/4.0)

   ! just sample around a bit first
   TryAdd(:) = 0
   TryRem(:) = 0
   do i = 1, 5000
      write(0, "('*')", advance='no')
      call test_ctqmc_sample_steps(1_c_int64_t, SecZ, 0, rng)
      call verify_diagram(diagram, VERIFY_RTOL)
   end do
   write (0,*)

   call init_addmove(addmove, diagram)
   call init_addmove(add4move, diagram)
   call init_remmove(remmove, diagram)
   call init_remmove(rem4move, diagram)
   call init_permutemove(permmove, diagram)
   call init_shiftmove(shiftmove, diagram)

   ! check equality of trace value after performing and reversing moves
   ! (FIXME: check probs if possible)
   do i = 1, 100
      write(0, "('*')", advance='no')
      ! sample around a bit to get different starting configurations for each iteration
      call test_ctqmc_sample_steps(100_c_int64_t, SecZ, 0, rng)
      call verify_diagram(diagram, VERIFY_RTOL)

      ! store local trace reference value
      trref = get_weight(diagram%local)

      ! test add single pair
      ! add move
      addloop: do
         if (allocated(newops)) deallocate(newops)
         allocate(newops(2))
         call generate_add(window, rng, 1, newops)
         propprob1 = 1.0d0 / proposal_probability_for_add(newops, window)
         if (propprob1 < MINPROB) cycle addloop
         call propose_addmove(addmove, newops)

         call choose_matching_outer_state(addmove%local, .false.)
         propprob2 = 1.0 !!!! FIXME: propprob_factor(addmove)
         if (propprob2 < MINPROB) cycle addloop
         weightratio = weight_ratio(addmove)
         if (weightratio < MINPROB) cycle addloop

         call accept_update(addmove)
         call verify_diagram(diagram, VERIFY_RTOL)
         exit addloop
      end do addloop

      ! FIXME: restore reverse move and numerical comparison if possible

      ! test add two pairs
      ! add 4-op move
      add4loop: do
         if (allocated(newops)) deallocate(newops)
         allocate(newops(4))
         call generate_add(window, rng, 2, newops)
         propprob1 = 1.0d0 / proposal_probability_for_add(newops, window)
         if (propprob1 < MINPROB) cycle add4loop
         call propose_addmove(add4move, newops)

         call choose_matching_outer_state(add4move%local, .false.)
         propprob2 = 1.0  !!!! FIXME: propprob_factor(add4move)
         if (propprob2 < MINPROB) cycle add4loop
         weightratio = weight_ratio(add4move)
         if (weightratio < MINPROB) cycle add4loop

         call accept_update(add4move)
         call verify_diagram(diagram, VERIFY_RTOL)
         exit add4loop
      end do add4loop

      ! ! revert add 4-op. move

      ! test rem one pair (FIXME: hoping there is at least one pair already...)
      ! rem move
      attempts = 0
      remloop: do
         if (allocated(local_remidx)) deallocate(local_remidx)
         allocate(local_remidx(2))
         local_ops = get_operators(diagram%local)
         call generate_remove( &
                        window, rng, local_ops, 1, &
                        local_remidx, allowed)
         if (.not. allowed) return
         propprob1 = proposal_probability_for_add( &
                        local_ops(local_remidx), window)

         if (propprob1 == 0) then
            ! there are no pairs
            write (*,*) &
                "WARNING: Cannot test one-pair removal and reversal ", &
                "(no pairs)"
            exit remloop
         end if
         attempts = attempts + 1
         if (attempts > maxattempts) then
            write (*,*) &
                "WARNING: Cannot test one-pair removal and reversal ", &
                "(no suitable pair found)"
            exit remloop
         end if
         call propose_remmove(remmove, local_remidx)

         call choose_matching_outer_state(remmove%local, .false.)
         propprob2 = 1.0 !!!! FIXME: propprob_factor(remmove)
         if (propprob2 < MINPROB) cycle remloop
         weightratio = weight_ratio(remmove)
         if (weightratio < MINPROB) cycle remloop

         call accept_update(remmove)
         call verify_diagram(diagram, VERIFY_RTOL)
         exit remloop
      end do remloop

      ! if (propprob1 /= 0 .and. attempts <= maxattempts) then
      !    ! revert rem move
      ! end if

      ! test rem one pair (FIXME: hoping there is at least one pair already...)
      ! rem move
      attempts = 0
      rem4loop: do
         if (allocated(local_remidx)) deallocate(local_remidx)
         allocate(local_remidx(4))
         local_ops = get_operators(diagram%local)
         call generate_remove( &
                    window, rng, local_ops, 2, &
                    local_remidx, allowed)
         if (.not. allowed) return
         propprob1 = proposal_probability_for_add( &
              local_ops(local_remidx), window)

         if (propprob1 < MINPROB) then
            ! there are not enough pairs
            write (*,*) &
                "WARNING: Cannot test two-pair removal and reversal ", &
                "(not enough pairs)"
            exit rem4loop
         end if
         attempts = attempts + 1
         if (attempts > maxattempts) then
            write (*,*) &
                "WARNING: Cannot test two-pair removal and reversal ", &
                "(no suitable pairs found)"
            exit rem4loop
         end if
         call propose_remmove(rem4move, local_remidx)

         call choose_matching_outer_state(rem4move%local, .false.)
         propprob2 = 1.0  !!!! FIXME: propprob_factor(rem4move)
         if (propprob2 < MINPROB) cycle rem4loop
         weightratio = weight_ratio(rem4move)
         if (weightratio < MINPROB) cycle rem4loop

         call accept_update(rem4move)
         call verify_diagram(diagram, VERIFY_RTOL)
         exit rem4loop
      end do rem4loop

      ! ! revert rem move

      call verify_diagram(diagram, VERIFY_RTOL)

      ! test shift
      ! shift move
      shiftloop: do
         deltatau = random_real(rng, beta)
         call propose_shiftmove(shiftmove, deltatau)
         propprob1 = 1

         weightratio = weight_ratio(shiftmove)
         if (weightratio < MINPROB) cycle shiftloop

         call accept_update(shiftmove)
         call verify_diagram(diagram, VERIFY_RTOL)
         exit shiftloop
      end do shiftloop

      ! reverse shift move
      deltatau = beta - deltatau
      call propose_shiftmove(shiftmove, deltatau)

      weightratio = weight_ratio(shiftmove)
      if (weightratio == 0) &
            TEST_ERROR()
      call accept_update(shiftmove)

      call verify_diagram(diagram, VERIFY_RTOL)

      ! test permute
      ! permute move
      permloop: do
         do_spin_flip = .false.
         rand = random_real(rng)
         if (rand < 0.3) then
            globtype = SpinFlip
            call generate_global_spinflip(get_nbands(dstates), flavcaperm)
            do_spin_flip = .true.
         else if (rand < 0.6) then
            globtype = CAFlip
            call generate_global_caflip(get_nbands(dstates), flavcaperm)
         else
            globtype = RandPerm
            call generate_global_flavca_permutation(get_nbands(dstates), flavcaperm, rng)
         end if

         call propose_permutemove(permmove, flavcaperm(0) /= 0, flavcaperm(1:))

         if (do_spin_flip) then
            outerbl = sst_for_flipped_spins(dstates, get_outerbl(diagram%local))
            if (is_allowed_outer_state(permmove%local, outerbl)) then
               propprob1 = 1
            else
               propprob1 = 0
            endif
         else
            ! one random number will be needed to choose the target outer
            ! sst/state, except for performing a spin-flip in
            ! superstatesampling only
            call make_random_outer_proposal( &
                   rng, permmove%local, diagram%local, outerbl, propprob1)
         endif
         if (propprob1 < MINPROB) cycle permloop

         call choose_outer_state(permmove%local, outerbl)

         weightratio = weight_ratio(permmove)
         if (weightratio < MINPROB) cycle permloop

         permcmp_outerbl = get_outerbl(diagram%local)
         call accept_update(permmove)
         call verify_diagram(diagram, 10*VERIFY_RTOL)
         exit permloop
      end do permloop

      ! ! reverse permute
   end do
   write (0,*)

   call dump(diagram%local)
   call verify_diagram(diagram, VERIFY_RTOL)

   ! Clear
   call clear_diagram(diagram)
   call verify_diagram(diagram, VERIFY_RTOL)

contains

    subroutine local_setup(b_complexLocTr)
        logical, intent(in)        :: b_complexLocTr

        type(TOperator)            :: dh
        real(c_double), parameter  :: EPSDEG = 1d-12

        call init_States(dstates, nbands, qnspec)
        call init_Psi(psis_eb, dstates, b_complexLocTr)
        call make_hamiltonian(muimp_full, u_matrix, dstates, psis_eb, dh, 0.0d0)

        call diag_Operator(dh, hloc_eb, EPSDEG)
        call transform_Psis(hloc_eb, psis_eb)

        call init_symmetry_moves(sym_moves, user_symmoves)
        call init_impurity(local, beta, hloc_eb, psis_eb, 0)
     end subroutine local_setup

    ! FIXME: maybe merge with ctqmc_sample_steps again when that one
    ! takes its configuration from arguments and not module variables
    subroutine test_ctqmc_sample_steps(Nsteps, allowed_worm_sector, iComponent, samprng)
        integer(c_int64_t), intent(in)          :: Nsteps
        integer, intent(in)                     :: allowed_worm_sector, iComponent
        class(RandomEngine), intent(inout)      :: samprng

        type(TDiagramMoveStats)                 :: movestat
        type(TAddMove) :: add_move
        type(TPermuteMove) :: permute_move
        type(TRemMove) :: rem_move
        type(TShiftMove) :: shift_move
        type(TWormAddMove) :: worm_add_move
        type(TWormRemMove) :: worm_rem_move
        type(TWormReplaceMove) :: worm_repl_move

        integer(c_int64_t) :: i
        integer :: move_type, current_sector

        call init_addmove(add_move, diagram)
        call init_permutemove(permute_move, diagram)
        call init_remmove(rem_move, diagram)
        call init_shiftmove(shift_move, diagram)
        call init_wormaddmove(worm_add_move, diagram)
        call init_wormremmove(worm_rem_move, diagram)
        call init_wormreplacemove(worm_repl_move, diagram)

        do i = 1, Nsteps
            current_sector = get_current_worm_sector(diagram%worm)
            call perform_move( &
                    samprng, movestat, allowed_worm_sector, iComponent, &
                    PercentageTauShiftMove, PercentageGlobalMove, &
                    PercentageWormInsert, PercentageWormReplace, &
                    PercentageOuterMove, Percentage4OperatorMove, &
                    wormEta(allowed_worm_sector), &
                    diagram, add_move, rem_move, permute_move, shift_move, &
                    worm_add_move, worm_rem_move, worm_repl_move, move_type, &
                    window, sym_moves, dstates)

            if (move_type == MOVE_TYPE_ADD) &
                TryAdd(current_sector) = TryAdd(current_sector) + 1
            if (move_type == MOVE_TYPE_REMOVE) &
                TryRem(current_sector) = TryRem(current_sector) + 1
            if (modulo(sum(TryAdd) + sum(TryRem), NSlide) == 0) &
                call shift_window(window)
        end do
    end subroutine test_ctqmc_sample_steps

end program test
