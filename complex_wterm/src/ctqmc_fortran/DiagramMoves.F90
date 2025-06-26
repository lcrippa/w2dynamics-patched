module MDiagramMoves
    use iso_c_binding, only: c_double, c_int64_t
    use MBath
    use MCommon
    use MDiagram
    use MImpurity
    use MMoveStatistics
    use MOperator
    use MPrinting, only: default
    use MCompoundIndex
    use MRandomNumbers
    use MStates
    use MSymmetryMoves
    use MWindow
    use MWormState
    implicit none
    private

    type, public :: TDiagramMoveStats
        type(TMoveStat) :: addhybpair(2), remhybpair(2)
        type(TMoveStat) :: addworm, remworm, repworm
        type(TMoveStat) :: taushift, changeouter, permute
    end type TDiagramMoveStats

    interface init_tau_binning_statistics
        module procedure diagmovestats_bin_mode
    end interface
    public init_tau_binning_statistics

    interface init_tau_history_statistics
        module procedure diagmovestats_list_mode
    end interface
    public init_tau_history_statistics

    interface reset_statistics
        module procedure diagmovestats_reset
    end interface
    public reset_statistics

    interface dump
        module procedure p_dump_statistics
    end interface
    public dump

    public :: perform_AddHybPair, perform_ChangeOuter, perform_Permute&
        &, perform_RemHybPair, perform_Shift, perform_WormAdd,&
        & perform_WormRem, perform_WormReplace, perform_move
    public :: local_opers_from_arrays, localops_to_arrays

    public :: generate_global_spinflip, generate_global_caflip, &
             generate_global_user_move, generate_global_flavca_permutation, &
             make_random_outer_proposal, random_outer_state_by_weight

    !> Collection of weights for proposing different move types
    type, public :: TProposalWeights
        real(c_double) :: GlobalMove
        real(c_double) :: TauShiftMove
        real(c_double) :: OuterMove
        real(c_double) :: FourOperatorMove
        real(c_double) :: WormInsert
        real(c_double) :: WormReplace
    end type

    ! Move types
    integer, parameter, public :: &
        MOVE_TYPE_ADD = 1,          MOVE_TYPE_ADD_4 = 2, &
        MOVE_TYPE_REMOVE = 3,       MOVE_TYPE_REMOVE_4 = 4, &
        MOVE_TYPE_PERMUTE = 5,      MOVE_TYPE_SHIFT = 6, &
        MOVE_TYPE_CHANGE_OUTER = 7, MOVE_TYPE_REPLACE = 8, &
        MOVE_TYPE_ADD_WORM = 9,     MOVE_TYPE_REMOVE_WORM = 10, &
        MOVE_TYPE_REPLACE_WORM = 11


contains

    !> Main routine for selecting and perfoming a Monte Carlo update
    subroutine perform_move( &
                    samprng, movestat, allowed_worm_sector, iComponent, &
                    prop_weights, wormEta, &
                    diagram, add_move, rem_move, perm_move, shift_move, &
                    worm_add_move, worm_rem_move, worm_repl_move, move_type, &
                    window, sym_moves, dstates)
        ! XXX the argument list is ridiculous ... clean it up!
        class(RandomEngine), intent(inout) :: samprng
        type(TDiagramMoveStats), intent(inout)  :: movestat
        integer, intent(in) :: allowed_worm_sector, iComponent
        type(TProposalWeights), intent(in) :: prop_weights
        real(c_double), intent(in) :: wormEta
        class(TDiagram), intent(inout) :: diagram
        class(TAddMove), intent(inout) :: add_move
        class(TPermuteMove), intent(inout) :: perm_move
        class(TRemMove), intent(inout) :: rem_move
        class(TShiftMove), intent(inout) :: shift_move
        class(TWormAddMove), intent(inout) :: worm_add_move
        class(TWormRemMove), intent(inout) :: worm_rem_move
        class(TWormReplaceMove), intent(inout) :: worm_repl_move
        integer, intent(out) :: move_type
        type(TWindowState), intent(in) :: window
        type(TSymmetryMoves), intent(in) :: sym_moves
        type(TStates), intent(in) :: dstates ! XXX remove this

        integer(c_int64_t) :: current_sector
        real(c_double) :: move_type_rnd, move_type_threshold, outer_rnd, fourop_rnd

        outer_rnd = 1
        current_sector = get_current_worm_sector(diagram%worm)
        move_type_rnd = random_real(samprng)
        move_type = 0

        ! currently in Z sector
        move_type_threshold = 1.0d0

        ! In Z sector, we also perform tau shift and global moves
        ! XXX figure out why only there
        if (current_sector == SecZ) then
            ! MOVE: tau shift (#operators > 0) / outer state change (#operators = 0)
            move_type_threshold = move_type_threshold - prop_weights%TauShiftMove
            if (move_type_rnd > move_type_threshold) then
                if (size(diagram%local) == 0) then
                    move_type = MOVE_TYPE_CHANGE_OUTER
                    call perform_ChangeOuter(diagram, samprng, movestat%changeouter)
                    return
                else
                    move_type = MOVE_TYPE_SHIFT
                    call perform_Shift(shift_move, diagram, samprng, movestat%taushift)
                    return
                end if
            end if

            ! MOVE: global (permute)
            move_type_threshold = move_type_threshold - prop_weights%GlobalMove
            if (move_type_rnd > move_type_threshold) then
                move_type = MOVE_TYPE_PERMUTE
                call perform_Permute(perm_move, diagram, samprng, sym_moves, &
                                     movestat%permute, dstates)
                return
            end if
        endif

        ! Worm moves: if we are in Z sector, propose insertion of worms, if we
        ! are in a worm sector, propose their removal or replacement
        if (current_sector == SecZ) then
            if (allowed_worm_sector > SecZ) then
                move_type_threshold = move_type_threshold - prop_weights%WormInsert
                if (move_type_rnd > move_type_threshold) then
                    move_type = MOVE_TYPE_ADD_WORM
                    call perform_WormAdd(worm_add_move, diagram, samprng, &
                            allowed_worm_sector, iComponent,&
                            wormEta, movestat%addworm)
                    return
                end if
            end if
        else
            ! XXX revisit this once we allow worm replacement moves to jump
            ! between sectors
            if (current_sector /= allowed_worm_sector) &
                error stop 'Not in an allowed worm sector'

            move_type_threshold = move_type_threshold - prop_weights%WormInsert
            if (move_type_rnd > move_type_threshold) then
                move_type = MOVE_TYPE_REMOVE_WORM
                call perform_WormRem(worm_rem_move, diagram, samprng,&
                        wormEta, movestat%remworm)
                return
            end if

            move_type_threshold = move_type_threshold - prop_weights%WormReplace
            if (move_type_rnd > move_type_threshold) then
                move_type = MOVE_TYPE_REPLACE_WORM
                call perform_WormReplace( &
                        worm_repl_move, diagram, samprng, movestat%repworm)
                return
            end if
        endif

        ! Otherwise, add/remove operators
        ! 1-/N-pair moves follow; need more random numbers
        fourop_rnd = random_real(samprng)
        if (fourop_rnd >= prop_weights%FourOperatorMove) &
            outer_rnd = random_real(samprng)

        move_type_threshold = move_type_threshold / 2.0d0
        if (move_type_rnd > move_type_threshold) then
            if (fourop_rnd >= prop_weights%FourOperatorMove) then
                ! 1-pair
                move_type = MOVE_TYPE_REMOVE
                call perform_RemHybPair(window, rem_move, diagram, 1, samprng,&
                        outer_rnd < prop_weights%OuterMove, movestat%remhybpair(1))
                return
            else
                ! 2-pair
                move_type = MOVE_TYPE_REMOVE_4
                call perform_RemHybPair(window, rem_move, diagram, 2, samprng,&
                        .false., movestat%remhybpair(2))
                return
            end if
        else
            if (fourop_rnd >= prop_weights%FourOperatorMove) then
                ! 1-pair
                move_type = MOVE_TYPE_ADD
                call perform_AddHybPair(window, add_move, diagram, 1, samprng,&
                        outer_rnd < prop_weights%OuterMove, movestat%addhybpair(1))
                return
            else
                ! 2-pair
                move_type = MOVE_TYPE_ADD_4
                call perform_AddHybPair(window, add_move, diagram, 2, samprng,&
                        .false., movestat%addhybpair(2))
                return
            end if
        end if

        error stop 'No move selected - move percentages likely add up to > 1'
    end subroutine

    !> Perform Monte Carlo configuration update removing (a) pair(s) of
    !> operators with connected hybridization event from diagram
    !> according to Metropolis update rule.
    !>
    !> @param diagram     Current diagram, will be modified if the move is accepted
    !> @param rng         Random number generator
    !> @param outer       Whether to perform outer or inner move
    ! FIXME: name when cleaning up
    ! FIXME: pull outer decision inside
    subroutine perform_RemHybPair(window, move, diagram, npairs, rng, outer, stat)
        type(TWindowState), intent(in) :: window
        type(TDiagram), intent(inout), target :: diagram
        integer, intent(in)           :: npairs
        class(RandomEngine), intent(inout)     :: rng
        logical, intent(in)           :: outer
        type(TMoveStat), intent(inout) :: stat

        type(TRemMove), intent(inout) :: move
        integer                       :: outcome
        real(c_double)                   :: propprob1, propprob3
        real(c_double)                   :: weightratio, rand
        logical                       :: allowed
        type(TLocalOper), allocatable :: ops(:)
        type(TLocalOper)              :: rem_ops(2*npairs)
        integer                       :: local_remidx(2*npairs)

        allocate(ops, source=get_operators(diagram%local))
        outcome = MOVE_FAILED_TO_GENERATE

        call generate_remove(window, rng, ops, npairs, local_remidx, allowed)
        if (.not. allowed) goto 99

        outcome = MOVE_NOT_ALLOWED
        rem_ops = ops(local_remidx)

        call propose_remmove(move, local_remidx)
        call choose_matching_outer_state(move%local, outer)

        if (.not. is_update_allowed(move%local)) &
            goto 99
        if (.not. is_update_allowed(move%bath)) &
            goto 99

        propprob3 = proposal_probability_for_remove(local_remidx, ops, window)
        if (propprob3 == 0) goto 99

        propprob3 = 1.0 / propprob3

        propprob1 = proposal_probability_for_add(rem_ops, window)

        call compute_weight(move%local)
        call compute_weight(move%bath)
        weightratio = weight_ratio(move)
        if (weightratio == 0) goto 99

        rand = random_real(rng)
        if (rand < abs(propprob1 * propprob3 * weightratio)) then
            call accept_update(move)
            outcome = MOVE_ACCEPTED
        else
            outcome = MOVE_REJECTED
        endif

    99 call record_move( &
                stat, get_current_worm_sector(diagram%worm), outcome, rem_ops)
    end subroutine perform_RemHybPair

    !> Perform Monte Carlo configuration update adding (a) pair(s) of
    !> operators with connected hybridization event to the diagram
    !> according to Metropolis update rule.
    !>
    !> @param diagram     Current diagram, will be modified if the move is accepted
    !> @param rng         Random number generator
    !> @param outer       Whether to perform outer or inner move
    ! FIXME: name when cleaning up
    ! FIXME: pull outer decision inside
    subroutine perform_AddHybPair(window, move, diagram, npairs, rng, outer, stat)
        type(TWindowState), intent(in) :: window
        type(TAddMove), intent(inout) :: move
        type(TDiagram), intent(inout), target :: diagram
        integer, intent(in)           :: npairs
        class(RandomEngine), intent(inout)     :: rng
        logical, intent(in)           :: outer
        type(TMoveStat), intent(inout) :: stat

        type(TLocalOper)              :: newops(2*npairs)
        real(c_double)                :: propprob1, propprob3, weightratio, rand
        integer                       :: outcome
        integer, allocatable          :: newidx(:)

        ! Subroutine to perform an addition move of an operator pair
        ! connected to hybridization events
        outcome = MOVE_NOT_ALLOWED

        call generate_add(window, rng, npairs, newops)

        call propose_addmove(move, newops)
        call choose_matching_outer_state(move%local, outer)

        if (.not. is_update_allowed(move%local)) &
            goto 99
        if (.not. is_update_allowed(move%bath)) &
            goto 99

        propprob1 = 1.0d0 / proposal_probability_for_add(newops, window)

        call get_newidx_addorder(move%local, newidx)
        propprob3 = proposal_probability_for_remove( &
                            newidx, get_operators(move%local), window)

        call compute_weight(move%local)
        call compute_weight(move%bath)
        weightratio = weight_ratio(move)
        if (weightratio == 0) goto 99

        rand = random_real(rng)
        if (rand < abs(propprob1 * propprob3 * weightratio)) then
            call accept_update(move)
            outcome = MOVE_ACCEPTED
        else
            outcome = MOVE_REJECTED
        end if

    99 call record_move(stat, get_current_worm_sector(diagram%worm), outcome, &
                       newops)
    end subroutine perform_AddHybPair

    !> Function returning the proposal probability factor for addition
    !> of worm operators
    !>
    !> @param sector     Worm sector
    !> @param component  Worm component
    !> @param beta       Inverse temperature
    !> @param eta        Sector change proposal probability factor eta
    !> @param nops       Number of worm operators  (FIXME: how counted for composite when/if added?)
    !> @param nbands     Number of orbitals
    real(c_double) function propprob_factor_wormsector( &
                    sector, component, beta, eta, ntaus, nbands) result(factor)
        integer(kind(SecDummy)), intent(in) :: sector
        integer, intent(in) :: component
        real(c_double), intent(in) :: beta, eta
        integer, intent(in) :: ntaus, nbands

        ! note: the factor nbands is only present to ensure
        ! comparability of eta values with old code versions
        ! FIXME: keep it?
        select case (sector)
        case (SecP2, SecP2pp)
            factor = eta * beta**2 * nbands**4  ! 2 times, 4 flavors
        case (SecQQ)
            factor = eta * beta**2 * nbands**6  ! 2 times, 6 flavors
        case (SecRaman)
            factor = eta * beta**3 * nbands**6  ! 3 times, 6 flavors
        case (SecCustom)
            factor = eta * beta**ntaus
        case default
            error stop "propprob_factor_wormsector: Sector unimplemented"
        end select
    end function propprob_factor_wormsector


    !> Perform Monte Carlo configuration update adding worm operators
    !> to the diagram according to Metropolis update rule.
    !>
    !> @param diagram           Current diagram, will be modified if the move is accepted
    !> @param rng               Random number generator
    !> @param target_sector     Worm sector to propose
    !> @param target_component  Worm component to propose
    !> @param eta               Sector change proposal probability factor eta
    ! FIXME: name when cleaning up
    subroutine perform_WormAdd(move, diagram, rng, target_sector, target_component, eta, stat)
        type(TDiagram), intent(inout), target :: diagram
        class(RandomEngine), intent(inout)     :: rng
        integer(kind(SecDummy)), intent(in) :: target_sector
        integer, intent(in) :: target_component
        real(c_double), intent(in) :: eta
        type(TMoveStat), intent(inout) :: stat

        type(TWormAddMoveSpec)        :: spec
        type(TWormAddMove), intent(inout) :: move
        real(c_double)                   :: propprob1, weightratio, rand
        integer                       :: nwormops, nwormtaus, nbands, outcome

        nbands = get_nbands(diagram%local)
        outcome = MOVE_FAILED_TO_GENERATE

        spec%target_sector = target_sector
        spec%target_component = target_component

        nwormops = get_sector_nworm(diagram%worm, target_sector)
        nwormtaus = get_sector_ntaus(diagram%worm, target_sector)

        allocate(spec%newops(nwormops))
        call set_optypes(spec%newops, nbands, target_sector, target_component)
        call set_optimes(spec%newops, diagram%beta, target_sector, rng)
        propprob1 =  propprob_factor_wormsector(target_sector,&
                                              target_component,&
                                              diagram%beta,&
                                              eta,&
                                              nwormtaus,&
                                              nbands)
        if (propprob1 == 0) goto 99

        outcome = MOVE_NOT_ALLOWED
        call propose_wormaddmove(move, spec)

        ! FIXME: in principle we could also do outer moves
        call choose_matching_outer_state(move%local, .false.)

        if (.not. is_update_allowed(move%local)) &
            goto 99

        call compute_weight(move%local)
        weightratio = weight_ratio(move)
        if (weightratio == 0) goto 99

        rand = random_real(rng)
        if (rand < abs(propprob1 * weightratio)) then
            call accept_update(move)
            outcome = MOVE_ACCEPTED
        else
            outcome = MOVE_REJECTED
        end if

    99 call record_move(stat, target_sector, outcome)
    contains
        ! Subroutine setting proper types (type = crea/annh/composite
        ! and flavor(s)) of new worm operators according to sector and
        ! component.
        subroutine set_optypes(newops, nbands, sector, component)
            type(TLocalOper), intent(inout)     :: newops(:)
            integer, intent(in)                 :: nbands
            integer(kind(SecDummy)), intent(in) :: sector
            integer, intent(in)                 :: component
            integer                             :: i
            integer :: bs(size(newops)), b(size(newops)), s(size(newops))

            ! convert component index to band spin pattern of length size(newops)
            call index2component_general(nbands, size(newops), component, bs, b, s)

            select case (sector)
            case (SecP2)
                newops%type = (/ OpCreaW, OpAnnhW, OpCreaW, OpAnnhW /)
                do i = 1, size(newops)
                    newops(i)%orb = b(size(newops) + 1 - i)
                    newops(i)%sp = s(size(newops) + 1 - i)
                end do
            case (SecP2pp)
                newops%type = (/ OpAnnhW, OpAnnhW, OpCreaW, OpCreaW /)
                do i = 1, size(newops)
                    newops(i)%orb = b(size(newops) + 1 - i)
                    newops(i)%sp = s(size(newops) + 1 - i)
                end do
            case (SecQQ)
                ! order of ops: from smallest to largest time,
                ! i.e. inverse order than in the formulas on paper.
                ! The left three operators are q^\dag (again, reversed order),
                ! the right three operators are q.
                newops%type = (/ OpAnnhW, OpCreaW, OpCreaW, &
                             OpAnnhW, OpAnnhW, OpCreaW /)
                do i = 1, size(newops)
                    newops(i)%orb = b(size(newops) + 1 - i)
                    newops(i)%sp = s(size(newops) + 1 - i)
                end do
            case (SecRaman)
                newops%type = (/ OpAnnhW, OpCreaW, OpAnnhW, &
                             OpCreaW, OpAnnhW, OpCreaW /)
                do i = 1, size(newops)
                    newops(i)%orb = b(size(newops) + 1 - i)
                    newops(i)%sp = s(size(newops) + 1 - i)
                end do
            case (SecCustom)
                do i = 1, size(newops)
                    newops(i)%type = get_custom_worm_optype(diagram%worm, sector, i)
                    newops(i)%orb = b(i)
                    newops(i)%sp = s(i)
                end do
            case default
                error stop "perform_WormAdd: set_optypes: Sector unimplemented"
            end select
        end subroutine set_optypes

        ! Subroutine generating as many random times as needed for a
        ! given sector and assigning them to new worm operators.
        subroutine set_optimes(newops, beta, sector, rng)
            type(TLocalOper), intent(inout) :: newops(:)
            real(c_double), intent(in) :: beta
            integer(kind(SecDummy)), intent(in) :: sector
            class(RandomEngine), intent(inout) :: rng
            integer :: i, j
            real(c_double) :: eqtime_tau

            ! Sectors Z, 1P-GF, 2P-GF:
            ! N operators with N random times
            select case (sector)
            case (SecZ)
                do i = 1, size(newops)
                    newops(i)%tau = random_real(rng, beta)
                end do
            case (SecP2, SecP2pp)  ! 2 random times for 4 operators
                newops(1)%tau = random_real(rng, beta)
                newops(2)%tau = newops(1)%tau

                newops(3)%tau = random_real(rng, beta)
                newops(4)%tau = newops(3)%tau
            case (SecQQ)  ! 2 random times for 6 operators

                ! time argument corresponding to q^\dag -> tau2
                newops(1)%tau = random_real(rng, beta)
                newops(2)%tau = newops(1)%tau
                newops(3)%tau = newops(1)%tau

                newops(4)%tau = random_real(rng, beta)
                newops(5)%tau = newops(4)%tau
                newops(6)%tau = newops(4)%tau
            case (SecRaman)  ! 3 random times for 6 operators
                newops(1)%tau = random_real(rng, beta)
                newops(2)%tau = newops(1)%tau

                newops(3)%tau = random_real(rng, beta)
                newops(4)%tau = newops(3)%tau

                newops(5)%tau = random_real(rng, beta)
                newops(6)%tau = newops(5)%tau
            case (SecCustom)

                j = 1

                do i = 1, get_custom_worm_ntaus(diagram%worm, sector)

                    eqtime_tau = random_real(rng, beta)

                    do j = j, j + get_custom_worm_ntauops(diagram%worm, sector, i) - 1
                        newops(j)%tau = eqtime_tau
                    end do

                end do

            case default
                error stop "perform_WormAdd: set_optimes: Sector unimplemented"
            end select
        end subroutine set_optimes
    end subroutine perform_WormAdd

    !> Perform Monte Carlo configuration update removing the worm operators from the diagram.
    !>
    !> @param diagram     Current diagram, will be modified if the move is accepted
    !> @param rng         Random number generator
    !> @param eta         Sector change proposal probability factor eta
    ! FIXME: name when cleaning up
    subroutine perform_WormRem(move, diagram, rng, eta, stat)
        type(TDiagram), intent(inout), target :: diagram
        class(RandomEngine), intent(inout)     :: rng
        real(c_double), intent(in) :: eta
        type(TMoveStat), intent(inout) :: stat

        type(TWormRemMove), intent(inout):: move
        real(c_double)                   :: propprob1, weightratio, rand
        integer                       :: sector, nwormtaus, nbands, outcome

        nbands = get_nbands(diagram%local)
        sector = get_current_worm_sector(diagram%worm)
        outcome = MOVE_FAILED_TO_GENERATE

        nwormtaus = get_sector_ntaus(diagram%worm, sector)

        ! Proposal probability is inverse of worm insertion proposal probability
        propprob1 =  1.0/propprob_factor_wormsector(sector,&
                                                  get_current_worm_comp(diagram%worm),&
                                                  diagram%beta,&
                                                  eta,&
                                                  nwormtaus,&
                                                  nbands)

        ! FIXME: can this happen and if so, maybe print warning?
        if (propprob1 == 0) goto 99

        outcome = MOVE_NOT_ALLOWED
        call propose_wormremmove(move)

        ! FIXME: in principle, we could also do outer proposal
        call choose_matching_outer_state(move%local, .false.)

        if (.not. is_update_allowed(move%local)) &
            goto 99

        call compute_weight(move%local)
        weightratio = weight_ratio(move)
        if (weightratio == 0) goto 99

        outcome = MOVE_REJECTED
        rand = random_real(rng)
        if (rand < abs(propprob1 * weightratio)) then
            call accept_update(move)
            outcome = MOVE_ACCEPTED
        end if

    99 call record_move(stat, sector, outcome)
    end subroutine perform_WormRem

    !> Perform Monte Carlo configuration update switching operator
    !> types between a worm operator and an operator with hybridization
    !> according to Metropolis update rule.
    !>
    !> @param diagram     Current diagram, will be modified if the move is accepted
    !> @param rng         Random number generator
    subroutine perform_WormReplace(move, diagram, rng, stat)
        type(TDiagram), intent(inout), target :: diagram
        class(RandomEngine), intent(inout)     :: rng
        type(TMoveStat), intent(inout) :: stat

        type(TWormReplaceMoveSpec)    :: spec
        type(TWormReplaceMove), intent(inout) :: move
        real(c_double)                   :: propprob, weightratio, rand
        integer                       :: outcome

        outcome = MOVE_FAILED_TO_GENERATE
        call generate_worm_replace(diagram,&
                                    spec,&
                                    propprob,&
                                    rng)
        if (propprob == 0) goto 99

        outcome = MOVE_NOT_ALLOWED
        call propose_wormreplacemove(move, spec)
        if (.not. is_update_allowed(move%bath)) &
            goto 99

        ! Only in this case we have a local move to make
        if (move%neqtauops > 1) then
            ! FIXME: in principle we could also do outer moves
            call choose_matching_outer_state(move%local, .false.)

            if (.not. is_update_allowed(move%local)) &
                goto 99

            call compute_weight(move%local)
        endif

        call compute_weight(move%bath)
        weightratio = weight_ratio(move)
        if (weightratio == 0) goto 99

        outcome = MOVE_REJECTED
        rand = random_real(rng)
        if (rand < abs(weightratio)) then
            call accept_update(move)
            outcome = MOVE_ACCEPTED
        end if

    99 call record_move(stat, get_current_worm_sector(diagram%worm), outcome)
    end subroutine perform_WormReplace

    !> Generates a worm replacement move.
    !>
    !> @param base Base diagram from which to move away from
    !> @param newspec Specification of the generated replacement
    !> @param rng Random number generator
    subroutine generate_worm_replace(base, newspec, propprob, rng)
        type(TDiagram), intent(in)      :: base
        type(TWormReplaceMoveSpec), intent(out) :: newspec
        real(c_double), intent(out)        :: propprob
        class(RandomEngine), intent(inout)       :: rng

        type(TLocalOper) :: wormop
        integer :: nhybops, nhybop, hybop_idx

        type(TLocalOper), allocatable :: local_ops(:)

        allocate(local_ops, source=get_operators(base%local))

        propprob = 1.0d0
        newspec%worm_reptypepos = random_integer(rng, 1, get_nworm(base%worm))
        newspec%local_wormidx = &
                get_worm_oper_position(base%worm, newspec%worm_reptypepos)
        wormop = local_ops(newspec%local_wormidx)

        nhybops = get_nosoper(base%bath, wormop%orb, wormop%sp)/2  ! FIXME: offdiag
        if (nhybops == 0) then
                propprob = 0
                return
        end if
        nhybop = random_integer(rng, 1, nhybops)

        select case (wormop%type)
        case (OpCreaW)
                do hybop_idx = 1, size(local_ops)
                if (local_ops(hybop_idx)%type == OpCrea&
                              .and. local_ops(hybop_idx)%orb == wormop%orb&
                              .and. local_ops(hybop_idx)%sp == wormop%sp) then
                    if (nhybop <= 1) exit
                    nhybop = nhybop - 1
                end if
                end do
        case (OpAnnhW)
                do hybop_idx = 1, size(local_ops)
                if (local_ops(hybop_idx)%type == OpAnnh&
                              .and. local_ops(hybop_idx)%orb == wormop%orb&
                              .and. local_ops(hybop_idx)%sp == wormop%sp) then
                    if (nhybop <= 1) exit
                    nhybop = nhybop - 1
                end if
                end do
        case default
                error stop 'Invalid worm operator type'
        end select
        newspec%local_repidx = hybop_idx
    end subroutine generate_worm_replace

    !> Perform Monte Carlo configuration update cyclically shifting all
    !> operators by the same Delta tau according to Metropolis update
    !> rule.
    !>
    !> @param diagram     Current diagram, will be modified if the move is accepted
    !> @param rng         Random number generator
    ! FIXME: name when cleaning up
    subroutine perform_Shift(move, diagram, rng, stat)
        type(TDiagram), intent(inout), target :: diagram
        class(RandomEngine), intent(inout)     :: rng
        type(TMoveStat), intent(inout) :: stat

        type(TShiftMove), intent(inout) :: move
        real(c_double)                   :: propprob, weightratio, rand, delta_tau
        integer                       :: outcome

        outcome = MOVE_NOT_ALLOWED

        ! Subroutine to perform a tau-shift move
        delta_tau = random_real(rng, diagram%beta)

        call propose_shiftmove(move, delta_tau)
        propprob = 1

        if (.not. is_update_allowed(move%local)) &
            goto 99
        if (.not. is_update_allowed(move%bath)) &
            goto 99

        call compute_weight(move%local)
        call compute_weight(move%bath)
        weightratio = weight_ratio(move)
        if (weightratio == 0) goto 99

        outcome = MOVE_REJECTED
        rand = random_real(rng)
        if (rand < abs(propprob * weightratio)) then
            call accept_update(move)
            outcome = MOVE_ACCEPTED
        end if

    99 call record_move(stat, get_current_worm_sector(diagram%worm), outcome)
    end subroutine perform_Shift

    !> Change the outer superstate (and state if applicabe) of a
    !> configuration without operators randomly.
    !>
    !> @param diagram     Current diagram, will be modified if the move is accepted
    !> @param rng         Random number generator
    ! FIXME: name when cleaning up
    subroutine perform_ChangeOuter(diagram, rng, stat)
        type(TDiagram), intent(inout), target :: diagram
        class(RandomEngine), intent(inout)     :: rng
        type(TMoveStat), intent(inout) :: stat
        integer :: outerbl

        if (size(diagram%local) /= 0) &
            error stop 'Diagram must be empty'

        outerbl = random_outer_state_by_weight(diagram%local, rng)
        call clear_diagram(diagram, outerbl)

        call record_move(stat, get_current_worm_sector(diagram%worm), MOVE_ACCEPTED)
    end subroutine perform_ChangeOuter

    !> Generate a local state localstate with empty configuration and
    !> random outer superstate. Useful for (re-)initializing a
    !> superstate-sampling simulation.
    function random_outer_state_by_weight(localstate, rng) result(outer)
        type(TImpurity), intent(in)    :: localstate
        class(RandomEngine), intent(inout) :: rng
        real(c_double), allocatable      :: weights(:)
        real(c_double)                   :: accumulator, rand
        integer                          :: i, outer

        allocate(weights, source=thermal_weights(localstate))

        accumulator = 0
        rand = random_real(rng)
        do i = 1, size(weights)
          accumulator = accumulator + weights(i)
          if (accumulator >= rand .and. weights(i) /= 0) then
              outer = i - 1
              return
          end if
        end do

        ! Fall-back
        outer = 0
    end function

    !> Perform Monte Carlo configuration update permuting orb/sp/ca
    !> according to Metropolis update rule.
    !>
    !> @param diagram     Current diagram, will be modified if the move is accepted
    !> @param rng         Random number generator
    ! FIXME: name when cleaning up
    subroutine perform_Permute(move, diagram, rng, sym_moves, stat, dstates)
        type(TDiagram), intent(inout), target :: diagram
        class(RandomEngine), intent(inout)       :: rng
        type(TMoveStat), intent(inout) :: stat
        type(TSymmetryMoves), intent(in) :: sym_moves
        type(TStates), intent(in) :: dstates

        integer, parameter              :: RandPerm = 1, UserSym = 2, SpinFlip = 3, CAFlip = 4
        integer                         :: globtype, outer_state
        integer, allocatable            :: flavcaperm(:)
        logical                         :: do_spin_flip
        type(TPermuteMove), intent(inout) :: move
        real(c_double)                     :: propprob, weightratio, rand
        integer                         :: nbands, outcome

        nbands = get_nbands(diagram%local)
        outcome = MOVE_FAILED_TO_GENERATE
        if (size(diagram%bath) <= 0) goto 99

        outcome = MOVE_NOT_ALLOWED
        do_spin_flip = .false.
        allocate(flavcaperm(0 : 2 * nbands))
        rand = random_real(rng)
        if (rand < 0.3) then
            globtype = SpinFlip
            call generate_global_spinflip(nbands, flavcaperm)
            do_spin_flip = dstates%qnspec%Szt
        else if (rand < 0.6) then
            globtype = CAFlip
            call generate_global_caflip(nbands, flavcaperm)
        else if (rand < 0.8 .and. sym_moves%NSymMove > 0) then
            globtype = UserSym
            call generate_global_user_move(nbands, flavcaperm, sym_moves, rng)
        else
            globtype = RandPerm
            call generate_global_flavca_permutation(nbands, flavcaperm, rng)
        end if

        ! XXX use two out arguments in generate_...
        call propose_permutemove(move, flavcaperm(0) /= 0, flavcaperm(1:))

        if (do_spin_flip) then
            outer_state = sst_for_flipped_spins(dstates, get_outerbl(diagram%local))
            if (is_allowed_outer_state(move%local, outer_state)) then
                propprob = 1
            else
                propprob = 0
            endif
        else
            call make_random_outer_proposal( &
                rng, move%local, diagram%local, outer_state, propprob)
        endif
        if (propprob == 0) goto 99

        call choose_outer_state(move%local, outer_state)
        if (.not. is_update_allowed(move%local)) &
            goto 99
        if (.not. is_update_allowed(move%bath)) &
            goto 99

        call compute_weight(move%local)
        call compute_weight(move%bath)
        weightratio = weight_ratio(move)
        if (weightratio == 0) goto 99

        outcome = MOVE_REJECTED
        rand = random_real(rng)
        if (rand < abs(propprob * weightratio)) then
            call accept_update(move)
            outcome = MOVE_ACCEPTED
        endif

    99 call record_move(stat, get_current_worm_sector(diagram%worm), outcome)
    end subroutine perform_Permute

    !> Generate a uniformly random flavour permutation possibly randomly
    !> including a type exchange.
    !>
    !> @param nbands Number of orbitals in the calculation
    !> @param flavcaperm Array representing the generated flavor permutation
    !> @param rng Random number generator
    subroutine generate_global_flavca_permutation(nbands, flavcaperm, rng)
        integer, intent(in)       :: nbands
        integer, intent(out)      :: flavcaperm(0 : 2 * nbands)
        class(RandomEngine), intent(inout) :: rng

        ! permute flavors randomly; note one-based indexing of arguments
        call random_permutation(rng, 2*nbands, flavcaperm(1:2*nbands))

        ! choose randomly whether to flip creators and annihilators
        flavcaperm(0) = random_integer(rng, 0, 1)
    end subroutine generate_global_flavca_permutation

    !> Generate the flavor permutation array for a spin-flip
    !>
    !> @param nbands Number of orbitals in the calculation
    !> @param flavcaperm Array representing the spin-flip flavor permutation
    pure subroutine generate_global_spinflip(nbands, flavcaperm)
        integer, intent(in)       :: nbands
        integer, intent(out)      :: flavcaperm(0 : 2 * nbands)
        integer                   :: i

        flavcaperm(0) = 0
        do i = 1, nbands
            flavcaperm(i) = i + nbands
        end do
        do i = nbands + 1, 2 * nbands
            flavcaperm(i) = i - nbands
        end do
    end subroutine generate_global_spinflip

    !> Generate the flavor permutation matrix for a creator/annihilator exchange
    !>
    !> @param nbands Number of orbitals in the calculation
    !> @param flavcaperm Array representing the generated type permutation without flavor change
    pure subroutine generate_global_caflip(nbands, flavcaperm)
        integer, intent(in)       :: nbands
        integer, intent(out)      :: flavcaperm(0 : 2 * nbands)
        integer                   :: i

        ! exchange creators and annihilators
        flavcaperm(0) = 1
        ! just set flavor exchange entries to flavor numbers
        do i = 1, 2 * nbands
            flavcaperm(i) = i
        end do
    end subroutine generate_global_caflip

    !> Generate the flavor permutation matrix for a random user-specified flavor permutation.
    !>
    !> @param nbands Number of orbitals in the calculation
    !> @param flavcaperm Array representing the generated flavor permutation
    !> @param sym_moves TSymmetryMoves object containing the arrays of the user-specified permutations
    !> @param rng Random number generator
    subroutine generate_global_user_move(nbands, flavcaperm, sym_moves, rng)
        integer, intent(in)       :: nbands
        integer, intent(out)      :: flavcaperm(0 : 2 * nbands)
        type(TSymmetryMoves), intent(in)  :: sym_moves
        class(RandomEngine), intent(inout) :: rng
        integer                   :: isym

        ! do not exchange creators and annihilators
        flavcaperm(0) = 0
        ! set flavor permutation array to one of the user specified ones
        isym = random_integer(rng, 1, sym_moves%NSymMove)
        flavcaperm(1 : 2 * nbands) = sym_moves%SymMoves(isym, :)
    end subroutine generate_global_user_move

    !> Turns a proposal cfgprop for a new operator list into a complete
    !> local configuration proposal cfgprop with one superstate allowed
    !> by quantum numbers randomly chosen, typically used for global
    !> moves.
    subroutine make_random_outer_proposal(rng, move, lcfg, outer_state, prob)
     class(RandomEngine), intent(inout)     :: rng
     type(TImpurityReplace), intent(in)    :: move
     type(TImpurity), intent(in)          :: lcfg
     integer, intent(out)                   :: outer_state
     real(c_double), intent(out)            :: prob

     integer, allocatable                   :: ssts_before(:), ssts(:)
     integer                                :: npssts_before, npssts, k
     real(c_double)                         :: rand

     allocate(ssts_before, source=allowed_outer_states(lcfg))
     npssts_before = size(ssts_before)

     allocate(ssts, source=allowed_outer_states(move))
     npssts = size(ssts)

     ! XXX we always make random number even if not needed
     rand = random_real(rng)

     if (npssts <= 0) then
            prob = 0
            return
     end if

     prob = real(npssts, c_double) / real(npssts_before, c_double)

     ! XXX this is an inefficient and biased way of getting a random int
     k = floor(real(npssts, c_double) * rand) + 1
     outer_state = ssts(k)
    end subroutine make_random_outer_proposal



    function local_opers_from_arrays(taus, orbs, spins, cas, hashybs) result(ops)
     real(c_double), dimension(:), intent(in) :: taus
     integer, dimension(size(taus)), intent(in) :: orbs, spins, cas, hashybs
     type(TLocalOper), allocatable :: ops(:)
     integer :: i

     allocate(ops(size(taus)))
     do i = 1, size(taus)
        ops(i)%tau = taus(i)
        ops(i)%orb = orbs(i)
        ops(i)%sp = spins(i)
        if (hashybs(i) == 0) then
           select case (cas(i))
           case (1)
              ops(i)%type = OpCreaW
           case (2)
              ops(i)%type = OpAnnhW
           case default
              error stop 'Invalid worm type'
           end select
        else
           select case (cas(i))
           case (1)
              ops(i)%type = OpCrea
           case (2)
              ops(i)%type = OpAnnh
           case default
              error stop 'Invalid Z type'
           end select
        endif
     enddo
    end

    !> Dump sufficient information about the operator sequence into arrays to
    !  recreate Z-space configurations.
    subroutine localops_to_arrays(ops, taus, orbs, spins, cas, hashybs)
     type(TLocalOper), intent(in) :: ops(:)
     real(c_double), dimension(size(ops)), intent(out) :: taus
     integer, dimension(size(ops)), intent(out) :: orbs, spins, cas, hashybs
     integer :: i

     do i = 1, size(ops)
        taus(i) = ops(i)%tau
        orbs(i) = ops(i)%orb
        spins(i) = ops(i)%sp
        select case (ops(i)%type)
        case (OpCrea)
           cas(i) = 1
           hashybs(i) = -1
        case (OpAnnh)
           cas(i) = 2
           hashybs(i) = -1
        case (OpCreaW)
           cas(i) = 1
           hashybs(i) = 0
        case (OpAnnhW)
           cas(i) = 2
           hashybs(i) = 0
        case default
           error stop 'Invalid type of operator'
        end select
     enddo
    end

    pure subroutine diagmovestats_list_mode(coll, nlist)
        class(TDiagramMoveStats), intent(inout) :: coll
        integer, intent(in) :: nlist
        integer :: i

        do i = 1, 1  ! size(coll%addhybpair)
            call init_tau_history_statistics(coll%addhybpair(i), nlist)
        end do

        do i = 1, 1  ! size(coll%remhybpair)
            call init_tau_history_statistics(coll%remhybpair(i), nlist)
        end do

        ! call init_tau_history_statistics(coll%addworm, nlist)
        ! call init_tau_history_statistics(coll%remworm, nlist)
        ! call init_tau_history_statistics(coll%repworm, nlist)
        ! call init_tau_history_statistics(coll%taushift, nlist)
        ! call init_tau_history_statistics(coll%changeouter, nlist)
        ! call init_tau_history_statistics(coll%permute, nlist)
    end subroutine diagmovestats_list_mode

    pure subroutine diagmovestats_bin_mode(coll, orbs, spins, nbins, limit)
        class(TDiagramMoveStats), intent(inout) :: coll
        integer, intent(in) :: orbs, spins, nbins
        real(c_double), intent(in) :: limit
        integer :: i

        do i = 1, 1  ! size(coll%addhybpair)
            call init_tau_binning_statistics(coll%addhybpair(i), orbs, spins, nbins, limit)
        end do

        do i = 1, 1  ! size(coll%remhybpair)
            call init_tau_binning_statistics(coll%remhybpair(i), orbs, spins, nbins, limit)
        end do

        ! call init_tau_binning_statistics(coll%addworm, orbs, spins, nbins, limit)
        ! call init_tau_binning_statistics(coll%remworm, orbs, spins, nbins, limit)
        ! call init_tau_binning_statistics(coll%repworm, orbs, spins, nbins, limit)
        ! call init_tau_binning_statistics(coll%taushift, orbs, spins, nbins, limit)
        ! call init_tau_binning_statistics(coll%changeouter, orbs, spins, nbins, limit)
        ! call init_tau_binning_statistics(coll%permute, orbs, spins, nbins, limit)
    end subroutine diagmovestats_bin_mode

    pure subroutine diagmovestats_reset(coll)
        class(TDiagramMoveStats), intent(inout) :: coll
        integer :: i

        do i = 1, size(coll%addhybpair)
            call reset_statistics(coll%addhybpair(i))
        end do

        do i = 1, size(coll%remhybpair)
        call reset_statistics(coll%remhybpair(i))
        end do

        call reset_statistics(coll%addworm)
        call reset_statistics(coll%remworm)
        call reset_statistics(coll%repworm)
        call reset_statistics(coll%taushift)
        call reset_statistics(coll%changeouter)
        call reset_statistics(coll%permute)
    end subroutine diagmovestats_reset

    subroutine p_dump_statistics(self, unit)
        type(TDiagramMoveStats), intent(in) :: self
        integer, intent(in), optional :: unit
        integer(c_int64_t) :: total
        integer :: u, i

        u = default(unit, 0)

        total = sum(self%addworm%pretry) + sum(self%remworm%pretry) &
              + sum(self%repworm%pretry) + sum(self%taushift%pretry) &
              + sum(self%changeouter%pretry) + sum(self%permute%pretry)
        do i = 1, size(self%addhybpair)
            total = total + sum(self%addhybpair(i)%pretry)
        enddo
        do i = 1, size(self%remhybpair)
            total = total + sum(self%remhybpair(i)%pretry)
        enddo

        write(u, "(60('='))")
        write(u, "(A, 23X, A, 5X, A, 4X, A)") &
            "MOVE", "PROPOSALS", "ACCEPTED", "ILLEGAL"
        !write(u, "(61('-'))")
        call dump_array_stat("Add pairs", self%addhybpair, total, u)
        call dump_array_stat("Remove pairs", self%remhybpair, total, u)
        call dump_single_stat("Add worm ops", self%addworm, total, u)
        call dump_single_stat("Remove worm", self%remworm, total, u)
        call dump_single_stat("Replace worm", self%repworm, total, u)
        call dump_single_stat("Shift tau", self%taushift, total, u)
        call dump_single_stat("Switch outer in empty", self%changeouter, total, u)
        call dump_single_stat("Permute flavours", self%permute, total, u)
        write(u, "(60('='))")
    end subroutine

    subroutine dump_array_stat(name, arr, total, u)
        character(len=*), intent(in) :: name
        type(TMoveStat), intent(in) :: arr(:)
        integer, intent(in) :: u
        integer(c_int64_t), intent(in) :: total
        character(len=len(name)+2) :: fname
        integer :: i

        do i = 1, size(arr)
            write(fname, "(A, ' ', I0)") name, i
            call dump_single_stat(fname, arr(i), total, u)
        enddo
    end subroutine

    subroutine dump_single_stat(name, self, xtotal, u)
        character(len=*), intent(in) :: name
        type(TMoveStat), intent(in) :: self
        integer, intent(in) :: u
        integer(c_int64_t), intent(in) :: xtotal
        integer(c_int64_t) :: total, valid, accepted
        real(c_double) :: factor

        total = sum(self%pretry)
        if (total == 0) &
            return

        valid = sum(self%try)
        accepted = sum(self%acc)
        factor = 100.0d0 / max(total, 1)
        write(u, 99) name, repeat('.', 25-len(name)), &
                     100.0d0 * total / xtotal, &
                     factor * accepted, factor * (total - valid)
    99  format (A, X, A, X, F7.3, ' %', 2X, 3(3X, F6.2, ' %'))
    end subroutine

end module MDiagramMoves
