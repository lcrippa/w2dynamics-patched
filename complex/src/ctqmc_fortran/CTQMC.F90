!!! Module which performs a continuous time QMC simulation
! Working
!===============================================================================
module MCTQMC
!===============================================================================
use MParameters
!use MOperator
use MOperator
use MCompoundIndex
use MDiagram
use MDiagramMoves
use MBath
use MStates
use MRandomNumbers
use type_progress
use MSignal
!use iso_c_binding  ! GCC bug 103931 workaround
use MWormMeas
use MMeasurementInterface
use MMoveStatistics
use MAccumulator
use MUtilities
use MPsi
use MSamplerD
use MSamplerZ
use MSymmetryMoves
use MWindow
!use NFFT_Samplers
use testing
implicit none

   integer(C_INT64_T)         :: Nwarmups,Nmeas,NCorr
   integer                    :: NGtau,NGiw,NBands,NStates
   integer(c_int64_t)         :: NSlide
   real(c_double)                :: beta
   type(TProposalWeights)     :: proposal_weights

   !we make the umatrix global to qmc
   complex(c_double_complex),allocatable :: u_matrix(:,:,:,:)

   !we consider NWormSectors different spaces, which are defined in the sector variable
   !in the CTQMCsimulate subroutine
   !CntSampling: steps in each space taken
   !CntMeas: measurements in each space made
   integer(c_int64_t)                  :: CntSampling(NWormSectors), CntMeas(NWormSectors)

   ! state dimension index of the form (sst offset + in-sst index)
   integer(c_int32_t),allocatable, target        :: StatesSuperstates(:)          ! superstate of (eigen-)state
   real(c_double),allocatable, target :: EigenstatesEnergies(:)        ! energy of eigenstate
   integer(c_int32_t),allocatable, target        :: OccbasisMapping(:, :, :)      ! (dm-index, orb, spin) -> occ-eigval
   real(c_double),allocatable, target :: occ(:,:,:,:)  ! o1,s1,o2,s2
   real(c_double),allocatable, target :: rho2(:,:,:,:)  ! 4 contracted orbitals and spins
   real(c_double),allocatable    :: double_occ(:,:,:,:) ! o1, s1, o2, s2  ??? we have to exchange double_occ and single_occ in the old matrix functions
   real(c_double),allocatable    :: single_occ(:,:) ! o1, s1
   real(c_double),allocatable    :: Ssquare(:) ! # electrons
   real(c_double),allocatable    :: Ssquare_operator(:,:) ! size(Hilbert space), size(Hilbert space)
   complex(c_double_complex),allocatable, target    :: rho1(:,:)  ! 2 contracted orbitals and spins
   complex(c_double_complex),allocatable :: zDensityMatrix(:,:)
   real(c_double),allocatable    :: ExpResDensityMatrix(:,:,:,:,:)

   !we count the acceptance of hyb pairs in each sector including z sector
   real(c_double), target             :: AccRem(NWormSectors),AccAdd(NWormSectors),AccGlob,AccShift
   integer(c_int64_t)         :: TryRem(NWormSectors),TryAdd(NWormSectors)

   ! Separate acceptance counters for outer and inner moves (try =
   ! move attempts, acc = absolute number of accepted moves, accqn =
   ! absolute number of moves rejected due to quantum number checking)
   ! These counters are set but code for post-processing or writing
   ! them is not present by default.
   real(c_double)             :: AccRem4,AccAdd4
   real(c_double)             :: AccFlavc

   !worm accepts and tries
   !1 -> 1P Worm, 2-> 2P Worm, 3-> IE Sigma, 4->IE chi
   !(mind shift of 1 as we dont consider Z)
   real(c_double), target             :: AccWormRem(2:NWormSectors),AccWormAdd(2:NWormSectors),AccWormRep(2:NWormSectors)

   !FIXME: this needs to generalized for components
   real(c_double)                :: wormEta(1:NWormSectors)

   !sign measurement counter in z
   real(c_double)             :: time_calibration, time_warmup, time_sim, time_sim_steps

   logical                    :: fancy_prog = .false.

   integer                    :: FourPnt

   ! HACK: f2py is unable to handle derived types, but it is also unable to
   !       recognize Fortran 2003 polymorphic entities, `class(...)`, so we
   !       hide the structs from f2py by wrapping them like so.
   class(TSymmetryMoves), allocatable, target :: sym_moves
   class(TStates), allocatable, target :: DStates
   class(TOperEigen), allocatable, target :: Hloc
   class(TPsis), allocatable, target  :: psis_eb
   class(TDiagram), allocatable, target :: diagram
   class(TAddMove), allocatable, target :: add_move
   class(TPermuteMove), allocatable, target :: perm_move
   class(TRemMove), allocatable, target :: rem_move
   class(TShiftMove), allocatable, target :: shift_move
   class(TWormAddMove), allocatable, target :: worm_add_move
   class(TWormRemMove), allocatable, target :: worm_rem_move
   class(TWormReplaceMove), allocatable, target :: worm_repl_move
   class(TWindowState), allocatable, target :: window

   class(RandomEngine), allocatable :: samprng

   integer, parameter :: GF4_MATSUBARA_FULL=8

   logical :: b_Eigenbasis, b_Densitymatrix, b_meas_susz, b_Ssquare
   logical :: b_meas_susz_mat
   logical :: b_offdiag
   logical :: complexHyb, complexLocTr
   logical :: b_full_offdiag



   ! communicates to the caller if the QMC sampling was aborted by a signal.
   !  0 = no abort, 1 = aborted, some data, 2 = aborted, no data
   integer                    :: aborted = 0
   integer                    :: SimID = -1
   character(4)               :: simstr

contains

!===============================================================================
subroutine init_CTQMC()
!===============================================================================
   logical :: ZphConv, b_exch, GtauDetRat

   if (get_Integer_Parameter("segment") /= 0) then
      write (0,*) 'WARNING: segment is currently not supported'
   endif

   if (get_Integer_Parameter("statesampling") /= 0) then
      error stop 'State sampling not supported'
   end if

   Nwarmups=get_LongInteger_Parameter("Nwarmups")
   Nmeas=get_LongInteger_Parameter("Nmeas")
   NGtau=get_Integer_Parameter("Ntau")
   NBands=get_Integer_Parameter("Nd")
   NStates=4**NBands
   NCorr=get_Integer_Parameter("NCorr")
   NGiw=get_Integer_Parameter("Niw")
!   NBin=get_Integer_Parameter("NBin")

   beta=get_Real_Parameter("beta")
   FourPnt=get_Integer_Parameter("FourPnt")

   proposal_weights%GlobalMove = get_Real_Parameter("PercentageGlobalMove")
   proposal_weights%TauShiftMove = get_Real_Parameter("PercentageTauShiftMove")
   proposal_weights%OuterMove = get_Real_Parameter("PercentageOuterMove")
   proposal_weights%FourOperatorMove = get_Real_Parameter("Percentage4OperatorMove")
   proposal_weights%WormInsert = get_Real_Parameter("PercentageWormInsert")
   proposal_weights%WormReplace = get_Real_Parameter("PercentageWormReplace")

   wormEta = 0

   b_Eigenbasis = .true.
   b_Densitymatrix = get_Integer_Parameter("MeasDensityMatrix") /= 0

   ! choose offdiagonal or diagonal hybridisation
   if (get_Integer_Parameter("offdiag")==0) then
      b_offdiag = .false.
   else
      b_offdiag = .true.
   endif

   ! real or complex hybridisation function
   if (get_Integer_Parameter("complexHyb")==0) then
      complexHyb = .false.
      !call status_error("--> real F(tau).")
      write(*,*) "--> real F(tau)."
   else
      complexHyb = .true.
      !call status_error("--> complex F(tau).")
      write(*,*) "--> complex F(tau)."
   endif

   ! real or complex impurity
   if (get_Integer_Parameter("complexLocTr")==0) then
      complexLocTr = .false.
      !call status_error("--> real F(tau).")
      write(*,*) "--> real local impurity."
   else
      complexLocTr = .true.
      !call status_error("--> complex F(tau).")
      write(*,*) "--> complex local impurity."
   endif

   if(complexHyb.and..not.b_offdiag)then
      !call status_output( "complex hybridisation not possible with diagonal hybridisation!")
      write(*,*) "Flag complexHyb = 1 does not make sense with diagonal hybridisation function!"
      !stop
   endif

   if (get_Integer_Parameter("full_offdiag") /= 0) then
      write (0,*) "WARNING: Parameter full_offdiag has no effect"
   endif

   b_exch = get_Integer_Parameter("flavourchange_moves") /= 0
   GtauDetRat = get_integer_parameter("MeasGtauDetRat") /= 0

   b_meas_susz = .false.
   b_meas_susz_mat = get_Integer_Parameter("MeasSusz") /= 0

   ZphConv = get_Integer_Parameter("ZPHConvention") /= 0

   if(allocated(occ))then
     write(*,*)"occ already allocated"
   endif

   allocate(occ(NBands, 2,NBands, 2))
   allocate(double_occ(NBands, 2, NBands, 2))
   allocate(single_occ(NBands, 2))
   if (b_Ssquare) allocate(Ssquare(0:NBands*2))

   if (b_Densitymatrix) then
      allocate(zDensityMatrix(NStates,NStates))
      allocate(rho1(2*NBands,2*NBands))
      allocate(rho2(2*NBands,2*NBands,2*NBands,2*NBands))
   end if

   if (b_Ssquare) allocate(Ssquare_operator(NStates,NStates))

   if(get_Integer_Parameter("MeasExpResDensityMatrix").ne.0)&
     allocate(ExpResDensityMatrix(NBands,2,0:get_Integer_Parameter("MaxHisto"),NStates,NStates))

   allocate(StatesSuperstates(0:DStates%NStates-1))
   allocate(EigenstatesEnergies(0:DStates%NStates-1))
   allocate(OccbasisMapping(0:DStates%NStates-1, NBands, 2))

   occ=0d0
   double_occ=0d0
   single_occ=0d0
   if (allocated(Ssquare)) Ssquare=0d0

   call reset_qmc_counters()

   if(allocated(zDensityMatrix)) then
      zDensityMatrix=0d0
      rho1=dcmplx(0.0d0, 0.0d0)
      rho2=0d0
   end if
   if(allocated(Ssquare_operator)) Ssquare_operator=0d0
   if(allocated(ExpResDensityMatrix)) ExpResDensityMatrix=0d0


   !initialize sector to 1, i.e. start sampling in partition function space
   time_calibration = 0
   time_warmup = 0
   time_sim = 0
   time_sim_steps = 0

   aborted = 0
end subroutine init_CTQMC

!===============================================================================
subroutine reset_qmc_counters() bind(C, name='reset_qmc_counters')
   AccRem=0
   AccAdd=0
   AccGlob=0
   AccShift=0
   AccAdd4=0
   AccRem4=0
   AccFlavc=0
   AccWormAdd=0
   AccWormRem=0
   AccWormRep=0
   TryRem=0
   TryAdd=0

   CntSampling = 0
   CntMeas = 0
end subroutine reset_qmc_counters
!===============================================================================

!===============================================================================
subroutine dest_ctqmc() bind(C, name='dest_ctqmc')
   call dest_States(DStates)

   deallocate(DStates)
   deallocate(sym_moves)
   deallocate(psis_eb)
   deallocate(hloc)
   deallocate(diagram)
   deallocate(add_move)
   deallocate(perm_move)
   deallocate(rem_move)
   deallocate(shift_move)
   deallocate(worm_add_move)
   deallocate(worm_rem_move)
   deallocate(worm_repl_move)
   deallocate(window)

   call meastargets%decouple()

   if(allocated(StatesSuperstates))deallocate(StatesSuperstates)
   if(allocated(EigenstatesEnergies))deallocate(EigenstatesEnergies)
   if(allocated(OccbasisMapping))deallocate(OccbasisMapping)
   if(allocated(Occ))deallocate(Occ)
   if(allocated(double_occ))deallocate(double_occ)
   if(allocated(single_occ))deallocate(single_occ)
   if(allocated(Ssquare))deallocate(Ssquare)
   if(allocated(rho1))deallocate(rho1)
   if(allocated(rho2))deallocate(rho2)
   if(allocated(zDensityMatrix))deallocate(zDensityMatrix)
   if(allocated(Ssquare_operator))deallocate(Ssquare_operator)
   if(allocated(ExpResDensityMatrix))deallocate(ExpResDensityMatrix)

   if (allocated(u_matrix)) deallocate(u_matrix)
end subroutine dest_ctqmc

! Sample for Nwarmup + Nfix steps in the Z and iSector sectors and
! returns the number of the last Nfix steps in the Z and iSector
! sectors as stepsZ and stepsWorm, respectively
subroutine count_qmc_sector_steps(iSector, iComponent, Nfix, Nwarmup, &
     stepsWorm, stepsZ) bind(C, name='count_qmc_sector_steps')
  use type_progress

  integer(c_int32_t), value :: iSector, iComponent
  integer(c_int64_t), value :: Nfix, Nwarmup
  integer(c_int64_t), intent(out) :: stepsWorm, stepsZ

!local
  integer                    :: Sector

  type(TDiagramMoveStats)    :: movestat
  type(TLocalOper), allocatable :: ops(:)
  integer :: outer_sst

   !save starting configuration
   allocate(ops, source=get_operators(diagram%local))
   outer_sst = get_outerbl(diagram%local)
   call reset_qmc_counters()

   !starting in partition function space
   Sector = SecZ
   call wormstate_leave_wormsector(diagram%worm)

   if (Nwarmup > 0) then
      call ctqmc_worm_warmup(Nwarmup, iSector, iComponent)
      call reset_qmc_counters()
   end if

   call ctqmc_sample_steps(Nfix, movestat, iSector, iComponent)

   stepsWorm = CntSampling(iSector)
   stepsZ = CntSampling(SecZ)

   Sector = SecZ
   call reset_qmc_counters()
   !reset to initial configuraiton
   call clear_diagram(diagram)
   call set_window(window, 0)
   call set_diagram_from_local_ops(diagram, ops, outer_sst)
end subroutine count_qmc_sector_steps


!===============================================================================
subroutine init_qmc_solver(nbands, nftau, c_umat, c_ftau, c_muimp, c_screening) &
     bind(C, name='init_qmc_solver')
!===============================================================================
   use MSignal

!input
   integer(c_int32_t), value :: nbands, nftau
   type(c_ptr), value        :: c_umat, c_ftau, c_muimp, c_screening

   complex(c_double_complex), pointer :: u_matrix_in(:, :, :, :) ! 2*NBands,2*NBands,2*NBands,2*NBands
   complex(c_double_complex), pointer :: Ftau_full(:, :, :, :, :) ! NBands,2,NBands,2,Nftau
   complex(c_double_complex), pointer :: muimp_full(:, :, :, :) ! NBands,2,NBands,2
   real(c_double), pointer :: screening_function(:, :, :, :, :) ! NBands,2,NBands,2,Nftau
!local
   real(c_double)                :: qnthreshold, taudiff_max
   type(TPsis)                :: DPsis
   type(TOperator)            :: DH
   integer                    :: i, outer_sst
   !integer, allocatable       :: states2substates(:)
   class(RandomEngine), allocatable :: rng
   class(DHybridization), allocatable :: dhybr
   class(ZHybridization), allocatable :: zhybr
   type(TImpurity)          :: local
   type(TBath)                :: bath
   type(TWormState)           :: worm

   call c_f_pointer(c_umat, u_matrix_in, [2 * nbands, 2 * nbands, &
                                          2 * nbands, 2 * nbands])
   call c_f_pointer(c_ftau, Ftau_full, [nbands, 2, nbands, 2, nftau])
   call c_f_pointer(c_muimp, muimp_full, [nbands, 2, nbands, 2])
   call c_f_pointer(c_screening, screening_function, [nbands, 2, nbands, 2, &
                                                      nftau])

   allocate(u_matrix(2*NBands,2*NBands,2*NBands,2*NBands))
   u_matrix = u_matrix_in

   allocate(DStates)
   allocate(sym_moves)
   allocate(psis_eb)
   allocate(hloc)
   allocate(diagram)
   allocate(add_move)
   allocate(perm_move)
   allocate(rem_move)
   allocate(shift_move)
   allocate(worm_add_move)
   allocate(worm_rem_move)
   allocate(worm_repl_move)
   allocate(window)

   !!! this needs to be done before init_ctqmc
   if(get_Integer_Parameter("complexLocTr").eq.0)then
       complexLocTr = .false.
   else
       complexLocTr = .true.
   endif

   if (get_String_Parameter("EnableOuterTruncation") == "YES") then
      error stop 'Truncation is not supported'
   end if
   if (get_Integer_Parameter("Uw") == 1) then
      error stop 'Phonons are currently unsupported'
   endif

   beta = get_Real_Parameter("beta")

   ! make threshold relative to the magnitude of the largest element of the
   ! local Hamiltonian
   qnthreshold = get_Real_Parameter("EPSBLOCK")
   qnthreshold = qnthreshold * max(maxval(abs(muimp_full)),&
                                   maxval(abs(u_matrix_in)))

   call init_States(DStates, NBands, qnspec_from_parameters())
   call init_Psi(DPsis,DStates,complexLocTr)
   call make_hamiltonian(muimp_full, u_matrix, DStates, DPsis, DH, qnthreshold)

   call diag_Operator(DH,Hloc,get_Real_Parameter("EPSEVEC"))
   call init_Psi(psis_eb, DStates, complexLocTr)
   call transform_Psis(Hloc, psis_eb)
   call dest_Psis(DPsis)
   call dest_TOperator(DH)

   call init_CTQMC()
   call init_states_arrays(DStates)
   call init_symmetry_moves_from_parameters(sym_moves)
   call init_impurity(local, beta, hloc, psis_eb, 0)

   ! modularized: initialize local part of diagram with empty state
   !call new_trng_from_singleton(rng)
   if (.not. allocated(samprng)) then
      call system_clock(i)
      write (0,*) "WARNING: RNG not set up before init_solver."
      write (0,*) "         Will seed with system_clock:", i
      call init_qmc_rng(i)
   end if
   allocate(rng, source=samprng)
   outer_sst = random_outer_state_by_weight(local, rng)
   call clear_impurity(local, outer_sst)

   ! modularized: bath initialization (FIXME: this comment if needed)
   if (complexHyb) then
      call init_linear_hybridization(zhybr, Ftau_full, beta, .true.)
      call init_bath(bath, zhybr)
   else
      if (.not. all(aimag(Ftau_full) == 0)) &
         write (0,*) 'WARNING: Function has imaginary part, but should not'

      call init_linear_hybridization(dhybr, real(Ftau_full), beta, .true.)
      call init_bath(bath, dhybr)
   endif

   ! modularized: worm initialization
   ! (FIXME: possibly sector & component, probably not; this comment if needed)
   ! ^^^ TODO: figure out what on Earth this comment means
   call enter_z_sector(worm)
   call init_diagram(diagram, local, bath, worm, b_offdiag)

   ! Initialize window
   taudiff_max = get_Real_Parameter("TaudiffMax")
   if (taudiff_max <= 0) &
        taudiff_max = beta
   call init_window(window, diagram, taudiff_max)

   NSlide = max(2 * NBands, &
                floor(NCorr / (2.0d0 * num_window_positions(window))))

   ! Allocate moves
   call init_addmove(add_move, diagram)
   call init_permutemove(perm_move, diagram)
   call init_remmove(rem_move, diagram)
   call init_shiftmove(shift_move, diagram)
   call init_wormaddmove(worm_add_move, diagram)
   call init_wormremmove(worm_rem_move, diagram)
   call init_wormreplacemove(worm_repl_move, diagram)

   write(simstr,'(I4)') SimID

   call register_signal_handler(SIGNAL_TERMINATE, record_signal)
   call register_signal_handler(SIGNAL_INTERRUPT, record_signal)
   call register_signal_handler(SIGNAL_USER_1,    record_signal)

contains
   subroutine init_states_arrays(DStates)
       type(TStates), intent(in) :: DStates
       integer :: i, ib, is, start, finish

       do i = 0, DStates%NSStates - 1
           start = DStates%SubStates(i)%Offset
           finish = start + DStates%SubStates(i)%NStates - 1
           StatesSuperstates(start:finish) = i
           EigenstatesEnergies(start:finish) = get_eigenvalues(Hloc, i)
           do ib = 0, NBands - 1
           do is = 0, 1
               OccbasisMapping(start:finish, ib + 1, is + 1) = &
                   ibits(DStates%SubStates(i)%States(:), ib + is * NBands, 1)
           end do
           end do
       end do
   end subroutine
end subroutine init_qmc_solver

!===============================================================================
! accpairbuf: a zero-initialized buffer large enough to hold
! AccPairTau if list is true or AccPair if list is false.
subroutine ctqmc_calibrate(list, phase1pairnum, phase2taugrid, phase2pairnum, &
  nprogress, accpairbuf) &
     bind(C, name='ctqmc_calibrate')
!===============================================================================
   logical(c_bool), value :: list
   integer(c_int64_t), value :: phase1pairnum, phase2taugrid, phase2pairnum, nprogress
   type(c_ptr), value :: accpairbuf
   integer(c_int64_t), contiguous, pointer :: AccPair(:, :, :)
   real(c_double), contiguous, pointer :: AccPairTau(:)

   integer(c_int64_t)         :: i, last_progress
   integer                    :: j, k
   real(c_double)                :: time_start, time_end
   type(TDiagramMoveStats)    :: movestat

   if (list) then
      call init_tau_history_statistics(movestat, int((phase1pairnum + 1)/2))
      call c_f_pointer(accpairbuf, AccPairTau, [phase1pairnum])
      AccPair => null()
   else
      call init_tau_binning_statistics(movestat, &
                nbands, 2, int(phase2taugrid), max_pair_distance(window))
      call c_f_pointer(accpairbuf, AccPair, [NBands, 2, int(phase2taugrid)])
      AccPairTau => null()
   end if

   call cpu_time(time_start)

   if (list) then
      write (*, *) "Rank " // simstr // " starting window calibration phase 1"
   else
      write (*, *) "Rank " // simstr // " starting window calibration phase 2"
   end if

   i = 0
   last_progress = 0
   calibration_loop: do
      call ctqmc_sample_steps(int(min(1000, nprogress), c_int64_t), movestat, SecZ, 0)
      if (any_signal_fired()) exit
      i = i + min(1000, nprogress)

      if (modulo(i, min(1000, nprogress)) == 0) then
         if (modulo(i, nprogress) <= last_progress) then
            write (*, *) "Rank " // simstr // ": ", i, "steps"
         end if
         last_progress = modulo(i, nprogress)

         if (list) then
            if (movestat%addhybpair(1)%listindex + movestat%remhybpair(1)%listindex - 2 >= phase1pairnum) then
               write (*, *) "Rank " // simstr // " window calibration phase 1 complete"
               ! FIXME (very long-term): do not use module variable to pass res. to python
               j = int(min(phase1pairnum, movestat%remhybpair(1)%listindex - 1))  ! pair taudiffs taken from rem moves
               AccPairTau(1:j) = movestat%remhybpair(1)%acc_pair_list(1:j)
               AccPairTau(j + 1 : phase1pairnum) = movestat%addhybpair(1)%acc_pair_list(1 : phase1pairnum - j)
               exit calibration_loop
            end if
         else
            do j = 1, size(movestat%addhybpair(1)%acc_pair_taudiffbin, dim=1)
               do k = 1, size(movestat%addhybpair(1)%acc_pair_taudiffbin, dim=2)
                  if (sum(movestat%addhybpair(1)%acc_pair_taudiffbin(j, k, :))&
                       + sum(movestat%remhybpair(1)%acc_pair_taudiffbin(j, k, :)) < phase2pairnum)&
                     cycle calibration_loop
               end do
            end do
            write (*, *) "Rank " // simstr // " window calibration phase 2 complete"
            ! FIXME (very long-term): do not use module variable to pass res. to python

            do j = 1, size(movestat%addhybpair(1)%acc_pair_taudiffbin, dim=1)
               do k = 1, size(movestat%addhybpair(1)%acc_pair_taudiffbin, dim=2)
                  AccPair(j, k, :) = AccPair(j, k, :)&
                       + movestat%addhybpair(1)%acc_pair_taudiffbin(j, k, :)&
                       + movestat%remhybpair(1)%acc_pair_taudiffbin(j, k, :)
               end do
            end do
            exit calibration_loop
         end if
      end if
   end do calibration_loop

   call cpu_time(time_end)
   time_calibration = time_calibration + (time_end - time_start)
end subroutine ctqmc_calibrate

subroutine set_qmc_taudiffmax(taudiffmax) bind(C, name='set_qmc_taudiffmax')
   real(c_double), value :: taudiffmax

   if (.not. allocated(window)) &
      error stop 'Invalid window'

   call set_max_pair_distance(window, taudiffmax)
end subroutine set_qmc_taudiffmax

subroutine ctqmc_sample_steps(Nsteps, movestat, allowed_worm_sector, iComponent, p)
   use type_progress
   integer(c_int64_t), intent(in)          :: Nsteps
   type(TDiagramMoveStats), intent(inout)  :: movestat
   integer, intent(in)                     :: allowed_worm_sector, iComponent
   type(progress), intent(inout), optional :: p

   integer(c_int64_t) :: i, current_sector
   integer :: move_type

   do i = 1, Nsteps
      if (any_signal_fired()) return

      current_sector = get_current_worm_sector(diagram%worm)
      CntSampling(current_sector) = CntSampling(current_sector) + 1
      call perform_move( &
                samprng, movestat, allowed_worm_sector, iComponent, &
                proposal_weights, wormEta(allowed_worm_sector), &
                diagram, add_move, rem_move, perm_move, shift_move, &
                worm_add_move, worm_rem_move, worm_repl_move, move_type, &
                window, sym_moves, dstates)

      if (present(p)) call ptick(p)
      if (move_type == MOVE_TYPE_ADD) &
         TryAdd(current_sector) = TryAdd(current_sector) + 1
      if (move_type == MOVE_TYPE_REMOVE) &
         TryRem(current_sector) = TryRem(current_sector) + 1
      if (modulo(sum(TryAdd) + sum(TryRem), NSlide) == 0) &
         call shift_window(window)

      !call verify_bath(diagram%bath, 1.0d-5)
      !call verify_impurity(diagram%local, 1d-12)
   end do
end subroutine ctqmc_sample_steps


!===============================================================================
subroutine ctqmc_warmup(Nwarmups) bind(C, name='ctqmc_warmup')
!===============================================================================
   use type_progress

   integer(c_int64_t), value :: Nwarmups
!local
   type(progress)             :: p
   real(c_double)                :: time_start
   type(TDiagramMoveStats)    :: movestat


   call reset_qmc_counters()
   call cpu_time(time_start)

   call pstart(p, NWarmups, title='Warmup ' // simstr // ':', fancy=fancy_prog)

   call ctqmc_sample_steps(NWarmups, movestat, SecZ, 0, p)

   call cpu_time(time_warmup)
   time_warmup = time_warmup - time_start
end subroutine ctqmc_warmup

! Warmup including steps in worm sector
subroutine ctqmc_worm_warmup(Nwarmups, iSector, iComponent) &
     bind(C, name='ctqmc_worm_warmup')
   use type_progress
   integer(c_int64_t), value :: Nwarmups
   integer(c_int32_t), value :: iSector, iComponent
!local
   type(progress)             :: p
   type(TDiagramMoveStats)    :: movestat


   call reset_qmc_counters()

   call pstart(p, NWarmups, title='Wormup ' // simstr // ':', fancy=fancy_prog)

   call ctqmc_sample_steps(NWarmups, movestat, iSector, iComponent, p)

end subroutine ctqmc_worm_warmup

!ctqmc sampling with measurements
!Z-sampling/conventional CT-HYB: iSector=1,iComponent=1
!Worm-sampling/component samling: iSector set for worm
!estimator and iComponent set for flavor component
!===============================================================================
subroutine ctqmc_measure(iSector, iComponent, MaxMeasurementTime) &
     bind(C, name='ctqmc_measure')
!===============================================================================
   use type_progress
   use MSignal
   use DMOperator
   use ZMOperator

   integer(c_int32_t), value :: iSector, iComponent
   real(c_double), value :: MaxMeasurementTime
!local
   type(progress)             :: p
   integer(c_int64_t)         :: iNMeas
   integer(c_int64_t)         :: accucomps
   real(c_double)                :: time_start, time_temp
   real(c_double)                :: zScaling
   integer(kind(SecDummy))    :: Sector
   complex(c_double_complex), pointer    :: zbuffer1d(:)
   type(TDiagramMoveStats)    :: movestat
   complex(c_double_complex)             :: mean_sign, inv_mean_sign
   type(DTOperator)           :: rho_mb_d
   type(ZTOperator)           :: rho_mb_z

   if (iSector == SecZ) then
      if (SimID == 0 .and. proposal_weights%WormInsert /= 0.0d0) then
         write (*, *) &
            "WARNING: ignoring PercentageWormInsert != 0.0", &
            "because this is a pure Z-space run"
      end if
      proposal_weights%WormInsert = 0.0d0
   end if

   call reset_qmc_counters()
   call pstart(p, NMeas, title='Simulation ' // simstr, fancy=fancy_prog)

   call cpu_time(time_sim)

   Sector = SecDummy

   if (iSector /= SecZ .and. wormEta(iSector) == 0d0)&
        error stop "ctqmc_measure: eta is zero for the current worm sector"

   OUTER_LOOP: &
   do iNMeas=1,NMeas
      !correlation steps in z and worm space
      call cpu_time(time_start)

      call ctqmc_sample_steps(NCorr, movestat, iSector, iComponent)
      if (any_signal_fired()) exit OUTER_LOOP

      call cpu_time(time_temp)
      time_sim_steps = time_sim_steps + (time_temp - time_start)

      if(MaxMeasurementTime .gt. 0)then
      if(time_temp - time_sim .gt. MaxMeasurementTime)then
         write(0,*) "Exiting measurement due to MaxMeasurementTime!"
         exit OUTER_LOOP
      endif
      endif

      call ptick(p)

      Sector = get_current_worm_sector(diagram%worm)

      CntMeas(Sector) = CntMeas(Sector) + 1

      !partition function space measurements
      if(Sector == SecZ) then

         call meastargets%mSign%measure(diagram)
         call meastargets%mExpOrd%measure(diagram)
         call meastargets%mDM%measure(diagram)

         ! only measure if component samling not enabled
         ! (maybe reconsider this, cf. combined Z/worm G branch)
         if(iSector == SecZ) then
            call meastargets%mGreen%measure(diagram)
            !call meastargets%mG2iw%measure(diagram)
            if (IAND(FourPnt,GF4_MATSUBARA_FULL)/=0) then
               call meastargets%mG4iw_full%measure(diagram)
            end if
            !!! TODOCOMPL: reactivate ntau n0
            !if (b_meas_susz_mat) call meastargets%mNtauN0%measure(diagram)
         endif


      elseif (Sector == SecP2) then
         call meastargets%mWormP2ph%measure(diagram)
      elseif (Sector == SecP2pp) then
         call meastargets%mWormP2pp%measure(diagram)
      elseif (Sector == SecQQ) then
         call meastargets%mWormQQ%measure(diagram)
      elseif (Sector == SecRaman) then
         call meastargets%mWormRaman%measure(diagram)
      elseif (Sector == SecCustom) then
         call meastargets%mWormCustom%measure(diagram)
      endif

   enddo OUTER_LOOP

   call cpu_time(time_temp)
   time_sim = time_temp - time_sim

   if(any_signal_fired()) then
      write(0,"('Rank',I3,'> WARNING: Process caught signal')") SimID
      write(0,"('Rank',I3,'; measurement',I10,' of',I10)")&
         SimID, iNMeas, NMeas
      if (iNMeas == 1) then
         ! no data
         aborted = 2
      else
         ! some data
         aborted = 1
      endif
      ! DMFT.py stops if aborted; we don't want this if solver is
      ! stopped by MaxMeasurementTime
      if(MaxMeasurementTime .gt. 0)then
         aborted = 0
      endif
   endif

   ! ==================== POSTPROCESSING FROM HERE ON ========================

   if (SimID <= 0) &
      call dump(movestat)

   if (iSector /= SecZ) then
      write(0,*) "Sampling steps in Z and wormspaces", CntSampling(SecZ), CntSampling(iSector)
      write(0,*) "Measurements in Z and wormspaces  ", CntMeas(SecZ), CntMeas(iSector)
   else
      write(0,*) "Sampling steps in Z", CntSampling(SecZ)
      write(0,*) "Measurements in Z  ", CntMeas(SecZ)
   endif

   if (CntMeas(SecZ) == 0) then
      write (0,"('Rank', I3, '> ERROR: No measurements in Z space')") SimID
   else
      !CntMeas(1): measurements in partition function space
      !CntMeas(2): measurements in greens function space
      allocate(zbuffer1d(1))
      call accumulator_mean(meastargets%mSign%accumulator, zbuffer1d)
      mean_sign = zbuffer1d(1)
      if (mean_sign == 0) then
         write (*, *) "ERROR: mean sign was exactly zero on rank " &
                      // simstr // ", setting quantity scaling to zero"
         inv_mean_sign = 0
      else
         inv_mean_sign = 1 / mean_sign
      end if
      deallocate(zbuffer1d)

      !!! calculate static observables
      !!! FIXME: move this to python as discussed (pass eigenbasis psis)
      if (meastargets%mDM%has_mean()) then

         accucomps = int(accumulator_num_comp(meastargets%mDM%accumulator))
         allocate(zbuffer1d(accucomps))
         call accumulator_mean(meastargets%mDM%accumulator, zbuffer1d)

         zDensityMatrix = reshape(zbuffer1d, (/ nstates, nstates /))
         zDensityMatrix = zDensityMatrix * inv_mean_sign
         deallocate(zbuffer1d)

         if (complexLocTr) then
            call block_diagonal_oper(rho_mb_z, dstates)
            call oper_from_matrix(rho_mb_z, zDensityMatrix)
            rho1 = rho1_from_density_matrix( &
                        rho_mb_z, hloc%zeig, psis_eb%zpsis)
            occ = double_occ_from_density_matrix( &
                        rho_mb_z, hloc%zeig, psis_eb%zpsis)
        else
            call block_diagonal_oper(rho_mb_d, dstates)
            call oper_from_matrix(rho_mb_d, real(zDensityMatrix))

            rho1 = rho1_from_density_matrix( &
                        rho_mb_d, hloc%deig, psis_eb%dpsis)
            occ = double_occ_from_density_matrix( &
                        rho_mb_d, hloc%deig, psis_eb%dpsis)
        endif

      endif

      write(0,*) "Z Mean Sign", mean_sign

      !acceptance in partition function space
      AccGlob = acceptance_ratio(movestat%permute, 0)
      AccShift = acceptance_ratio(movestat%taushift, 0)
      AccAdd4 = acceptance_ratio(movestat%addhybpair(2), 0)
      AccRem4 = acceptance_ratio(movestat%remhybpair(2), 0)

      AccAdd(SecZ) = acceptance_ratio(movestat%addhybpair(1), SecZ)
      AccRem(SecZ) = acceptance_ratio(movestat%remhybpair(1), SecZ)
      if (iSector /= SecZ) then
         AccAdd(iSector) = acceptance_ratio(movestat%addhybpair(1), iSector)
         AccRem(iSector) = acceptance_ratio(movestat%remhybpair(1), iSector)
         AccWormAdd(iSector) = acceptance_ratio(movestat%addworm, iSector)
         AccWormRem(iSector) = acceptance_ratio(movestat%remworm, iSector)
         AccWormRep(iSector) = acceptance_ratio(movestat%repworm, iSector)
      end if

      if (abs(imag(mean_sign)) > 1e-10) &
          write (*,*) '=== WARNING: mean sign is not real:', mean_sign
      zScaling = -1./dble(CntMeas(SecZ)) * real(inv_mean_sign)

      if(allocated(expresdensitymatrix))then
        ! forall(ib = 1:NBands, is = 1:2, i=0:(size(Histo(1,1,:))-1), Histo(ib,is,i)/=0)
        !   expresdensitymatrix(ib,is,i,:,:)=expresdensitymatrix(ib,is,i,:,:)/Histo(ib,is,i) ! FIXME
        ! endforall
        write (0, *) &
            "WARNING: unimplemented: the measured expansion order resolved ", &
            "density matrix MUST be transformed from eigenbasis to occupation", &
            "number basis to get correct results!"
      endif

      write (0,"('Done post-processing, CTQMC ', A4, ' complete.')") simstr
   endif


   call clear_signals_fired()
   call clear_signal_handler(SIGNAL_TERMINATE)
   call clear_signal_handler(SIGNAL_INTERRUPT)
   call clear_signal_handler(SIGNAL_USER_1)
 contains
   pure integer function array_nonzero_elements(array)
     complex(c_double_complex), intent(in) :: array(:)
     integer                    :: in1

     array_nonzero_elements = 0

     do in1 = 1, size(array)
        if (array(in1) /= 0) array_nonzero_elements = array_nonzero_elements + 1
     end do
   end function array_nonzero_elements
end subroutine ctqmc_measure

integer(c_int32_t) function get_qmc_config_size() bind(C, name='get_qmc_config_size')
   get_qmc_config_size = size(diagram%local)
end function get_qmc_config_size

integer(c_int32_t) function get_qmc_config(N, taus, orbs, spins, cas, hashybs) &
     result(outer_sst) bind(C, name='get_qmc_config')
   integer(c_int32_t), value :: N
   real(c_double), dimension(N), intent(out)     :: taus
   integer(c_int32_t), dimension(N), intent(out) :: orbs, spins, cas, hashybs

   type(TLocalOper), allocatable :: ops(:)

   allocate(ops, source=get_operators(diagram%local))
   outer_sst = get_outerbl(diagram%local)
   if (size(ops) /= N) &
      error stop 'Size mismatch'
   call localops_to_arrays(ops, taus, orbs, spins, cas, hashybs)
end function get_qmc_config

subroutine set_qmc_config(N, taus, orbs, spins, cas, hashybs, outer_sst) &
     bind(C, name='set_qmc_config')
   integer(c_int32_t), value                    :: N
   real(c_double), dimension(N), intent(in)     :: taus
   integer(c_int32_t), dimension(N), intent(in) :: orbs, spins, cas, hashybs
   integer(c_int32_t), value                    :: outer_sst

   call set_window(window, 0)
   call set_diagram_from_local_ops(diagram, &
            local_opers_from_arrays(taus, orbs, spins, cas, hashybs), &
            outer_sst)
end subroutine set_qmc_config

subroutine set_empty_qmc_config(outer_sst) bind(C, name='set_empty_qmc_config')
   integer(c_int32_t), value              :: outer_sst
   type(TLocalOper)                       :: ops(0)

   call set_window(window, 0)
   call set_diagram_from_local_ops(diagram, ops, outer_sst)
end subroutine set_empty_qmc_config

!>\brief initialise the random number generator before starting the qmc for the
!! first time
subroutine init_qmc_rng(Nseed) bind(C, name='init_qmc_rng')
!input
  integer(c_int32_t), value :: Nseed

  if (.not. allocated(samprng)) then
    call init_mersenne_twister(samprng)
  endif
  call seed_rng(samprng, Nseed)
end subroutine init_qmc_rng

!> Set fancy progress bar (FIXME: doesn't do it properly anyway; change/remove progress?)
subroutine set_fancy_progress_bar(fancy) bind(C, name='set_fancy_progress_bar')
  logical(c_bool), value :: fancy

  fancy_prog = fancy
end subroutine set_fancy_progress_bar

!> Set QMC simulation ID
subroutine set_qmc_simid(id) bind(C, name='set_qmc_simid')
  integer(c_int32_t), value :: id

  SimID = id
end subroutine set_qmc_simid

!> Set sector worm eta, or all for negative sector
subroutine set_sector_eta(sector, eta) bind(C, name='set_sector_eta')
  integer(c_int32_t), value :: sector
  real(c_double), value :: eta

  if (sector < 0) then
     wormEta(:) = eta
  else if (sector < 1 .or. sector > NWormSectors) then
     error stop "set_sector_eta: sector out of range"
  else
     wormEta(sector) = eta
  end if
end subroutine set_sector_eta

!> Get sector worm eta
real(c_double) function get_sector_eta(sector) bind(C, name='get_sector_eta')
  integer(c_int32_t), value :: sector

  if (sector < 1 .or. sector > NWormSectors) then
     error stop "get_sector_eta: sector out of range"
  else
     get_sector_eta = wormEta(sector)
  end if
end function get_sector_eta

!> Get sector measurement count
integer(c_int64_t) function get_sector_meas_count(sector) &
     bind(C, name='get_sector_meas_count')
  integer(c_int32_t), value :: sector

  if (sector < 1 .or. sector > NWormSectors) then
     error stop "get_sector_meas_count: sector out of range"
  else
     get_sector_meas_count = CntMeas(sector)
  end if
end function get_sector_meas_count

!> Get sector sample count
integer(c_int64_t) function get_sector_sample_count(sector) &
     bind(C, name='get_sector_sample_count')
  integer(c_int32_t), value :: sector

  if (sector < 1 .or. sector > NWormSectors) then
     error stop "get_sector_sample_count: sector out of range"
  else
     get_sector_sample_count = CntSampling(sector)
  end if
end function get_sector_sample_count

!> Get integer indicating whether measurement aborted due to signal
integer(c_int32_t) function get_qmc_aborted() &
     bind(C, name='get_qmc_aborted')
  get_qmc_aborted = aborted
end function get_qmc_aborted

!>\brief initialise the parameters with values defined in c_paramstring
!! the format is the same as for the parameter file
subroutine init_qmc_parameters(length, c_paramstring) &
     bind(C, name='init_qmc_parameters')
   integer(c_int32_t), value :: length
   character(kind=c_char), intent(in) :: c_paramstring(*)
   character(len=:), allocatable :: paramstring
   integer :: i

   allocate(character(len=length) :: paramstring)
   do i = 1, length
      paramstring(i:i) = c_paramstring(i)
   end do
   call read_ParameterString(paramstring)
   deallocate(paramstring)
end subroutine init_qmc_parameters

!>\brief initialise the parameters with values defined in parastr
!! the format is the same as for the parameter file
subroutine dest_qmc_paras() bind(C, name='dest_qmc_paras')
!input
   call dest_Parameters()
end subroutine dest_qmc_paras

subroutine get_leftovers_into_ndarrays_FIXME(nsize, c_name, asize, array) &
     bind(C, name='get_leftovers_into_ndarrays_FIXME')
  integer(c_int32_t), value :: nsize, asize
  character(kind=c_char)    :: c_name(*)
  type(c_ptr), value        :: array

  real(c_double), pointer, contiguous :: f_ptr_r(:), p_ptr_r(:)
  complex(c_double_complex), pointer, contiguous :: f_ptr_c(:), p_ptr_c(:)
  integer(c_int32_t), pointer, contiguous :: f_ptr_i(:), p_ptr_i(:)

  character(len=:), allocatable :: name
  integer :: i

  allocate(character(len=nsize) :: name)
  do i = 1, nsize
     name(i:i) = c_name(i)
  end do
  select case (name)
  case ("StatesSuperstates")
     call checksize(size(StatesSuperstates))
     f_ptr_i => StatesSuperstates
     call c_f_pointer(array, p_ptr_i, [asize])
     p_ptr_i(:) = f_ptr_i
  case ("EigenstatesEnergies")
     call checksize(size(EigenstatesEnergies))
     f_ptr_r => EigenstatesEnergies
     call c_f_pointer(array, p_ptr_r, [asize])
     p_ptr_r(:) = f_ptr_r
  case ("OccbasisMapping")
     call checksize(size(OccbasisMapping))
     f_ptr_i(1:asize) => OccbasisMapping
     call c_f_pointer(array, p_ptr_i, [asize])
     p_ptr_i(:) = f_ptr_i
  case ("occ")
     call checksize(size(occ))
     f_ptr_r(1:asize) => occ
     call c_f_pointer(array, p_ptr_r, [asize])
     p_ptr_r(:) = f_ptr_r
  case ("AccAdd")
     call checksize(size(AccAdd))
     f_ptr_r => AccAdd
     call c_f_pointer(array, p_ptr_r, [asize])
     p_ptr_r(:) = f_ptr_r
  case ("AccAdd4")
     call checksize(1)
     call c_f_pointer(array, p_ptr_r, [asize])
     p_ptr_r(1) = AccAdd4
  case ("AccRem")
     call checksize(size(AccRem))
     f_ptr_r => AccRem
     call c_f_pointer(array, p_ptr_r, [asize])
     p_ptr_r(:) = f_ptr_r
  case ("AccRem4")
     call checksize(1)
     call c_f_pointer(array, p_ptr_r, [asize])
     p_ptr_r(1) = AccRem4
  case ("AccGlob")
     call checksize(1)
     call c_f_pointer(array, p_ptr_r, [asize])
     p_ptr_r(1) = AccGlob
  case ("AccShift")
     call checksize(1)
     call c_f_pointer(array, p_ptr_r, [asize])
     p_ptr_r(1) = AccShift
  case ("AccFlavc")
     call checksize(1)
     call c_f_pointer(array, p_ptr_r, [asize])
     p_ptr_r(1) = AccFlavc
  case ("time_warmup")
     call checksize(1)
     call c_f_pointer(array, p_ptr_r, [asize])
     p_ptr_r(1) = time_warmup
  case ("time_sim")
     call checksize(1)
     call c_f_pointer(array, p_ptr_r, [asize])
     p_ptr_r(1) = time_sim
  case ("time_sim_steps")
     call checksize(1)
     call c_f_pointer(array, p_ptr_r, [asize])
     p_ptr_r(1) = time_sim_steps
  case ("rho2")
     call checksize(size(rho2))
     f_ptr_r(1:asize) => rho2
     call c_f_pointer(array, p_ptr_r, [asize])
     p_ptr_r(:) = f_ptr_r
  case ("rho1")
     call checksize(size(rho1))
     f_ptr_c(1:asize) => rho1
     call c_f_pointer(array, p_ptr_c, [asize])
     p_ptr_c(:) = f_ptr_c
  case ("AccWormAdd")
     call checksize(size(AccWormAdd))
     f_ptr_r => AccWormAdd
     call c_f_pointer(array, p_ptr_r, [asize])
     p_ptr_r(:) = f_ptr_r
  case ("AccWormRem")
     call checksize(size(AccWormRem))
     f_ptr_r => AccWormRem
     call c_f_pointer(array, p_ptr_r, [asize])
     p_ptr_r(:) = f_ptr_r
  case ("AccWormRep")
     call checksize(size(AccWormRep))
     f_ptr_r => AccWormRep
     call c_f_pointer(array, p_ptr_r, [asize])
     p_ptr_r(:) = f_ptr_r
  case default
     error stop "get_leftovers_into_ndarrays_FIXME: no such quantity"
  end select

  deallocate(name)
contains
  subroutine checksize(truesize)
    integer(c_int32_t) :: truesize
    if (asize /= truesize) &
         error stop "get_leftovers_into_ndarrays_FIXME: size mismatch"
  end subroutine checksize
end subroutine get_leftovers_into_ndarrays_FIXME

type(c_ptr) function get_wormstate() bind (C, name='get_wormstate')
    get_wormstate = c_loc(diagram%worm)
end function get_wormstate

end module MCTQMC
