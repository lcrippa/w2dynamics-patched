program test
    use iso_c_binding
    use MDiagram
    use MDiagramMoves
    use MImpurity
    use MMoveStatistics
    use MOperator
    use MPrinting
    use MRandomNumbers
    use MStates
    use MSymmetryMoves
    use MWindow
    use MWormState, only: SecZ
    use Testing
    implicit none

    real(c_double), parameter :: beta = 3, u = 1.3, mu = u/2
    real(c_double), parameter :: epsk(5) = (/ &
        -3.631194257431713d-01, -9.049035188745458d-02, &
        -3.291563979729537d-10, 9.049035017337467d-02, &
        3.631194218303095d-01 /)
    real(c_double), parameter :: vk(5) = (/ &
        0.4225631253617651d0, 0.2538903070113485d0, 0.15496346803614264d0, &
        0.25389030555412373d0, 0.4225631244300106d0 /)
    integer, parameter :: nmeas = 200000, nwarmup = 100, ncorr = 10
    integer, parameter :: user_sym_moves(1, 2) = reshape((/ 2, 1 /), (/ 1, 2 /))

    real(c_double), parameter :: perc_outer = 0.21d0

    type(TStates) :: states
    type(TOperEigen) :: hloc
    type(TDiagram) :: diagram
    type(TAddMove) :: add_move
    type(TRemMove) :: rem_move
    type(TShiftMove) :: shift_move
    type(TPermuteMove) :: perm_move
    type(TWormAddMove) :: worm_add_move
    type(TWormRemMove) :: worm_rem_move
    type(TWormReplaceMove) :: worm_repl_move
    class(RandomEngine), allocatable :: rng
    type(TDiagramMoveStats) :: movestat
    type(TSymmetryMoves) :: sym_moves
    type(TWindowState) :: window
    integer :: imeas, istep, move_type

    type(TOperator) :: mean_rho
    complex(c_double_complex) :: mean_sign, weight
    real(c_double), allocatable :: mean_rho_mat(:, :), mean_n, mean_d
    real(c_double), parameter :: mean_d_exact = 0.16672365210061169d0, mean_n_exact = 1

    write (0,*) 'Constructing diagram ...'
    call single_orbital(beta, u, mu, vk, epsk, states, hloc, diagram)

    write (0,*) 'Initialize sampling ...'
    call init_addmove(add_move, diagram)
    call init_remmove(rem_move, diagram)
    call init_shiftmove(shift_move, diagram)
    call init_permutemove(perm_move, diagram)
    call init_wormaddmove(worm_add_move, diagram)
    call init_wormremmove(worm_rem_move, diagram)
    call init_wormreplacemove(worm_repl_move, diagram)

    call init_mersenne_twister(rng)
    call init_symmetry_moves(sym_moves, user_sym_moves)
    call init_window(window, diagram, beta/4)

    call init_operator_like(mean_rho, hloc)
    mean_sign = 0

    write (0,*) 'Simulation ...'
    do imeas = -nwarmup, nmeas
        if (mod(80 * imeas, nmeas) == 0) &
            write (0, "('*')", advance='no')
        do istep = 1, ncorr
            call perform_move( &
                    rng, movestat, SecZ, 0, &
                    0.1d0, 0.01d0, 0.0d0, 0.0d0, perc_outer, 0.2d0, 0.0d0, &
                    diagram, add_move, rem_move, perm_move, shift_move, &
                    worm_add_move, worm_rem_move, worm_repl_move, move_type, &
                    window, sym_moves, states)
        enddo
        if (mod(imeas, 100) == 0) &
            call shift_window(window)

        if (imeas > 0) then
            weight = get_diagram_sign(diagram) / nmeas
            mean_sign = mean_sign + weight
            call compute_density_matrix(diagram%local, weight, mean_rho)
        endif
    enddo
    write (0,*)
    call dump(movestat)
    write (0,*)

    call print_array(real(mean_sign), name='sign')
    mean_rho_mat = getmatrix(mean_rho%dop) / real(mean_sign)
    call print_array(mean_rho_mat, name='rho')

    mean_n = mean_rho_mat(2, 2) + mean_rho_mat(3, 3) + 2 * mean_rho_mat(4, 4)
    mean_d = (mean_rho_mat(1, 1) + mean_rho_mat(4, 4)) / 2

    write (0, 98) mean_n, mean_n_exact
    write (0, 99) mean_d, mean_d_exact
98  format ('<n> =          ', F10.6, ' (exact =', F10.6, ')')
99  format ('<(n - <n>)²> = ', F10.6, ' (exact =', F10.6, ')')

    call assert_close(mean_n, mean_n_exact, rtol=1d-2)
    call assert_close(mean_d, mean_d_exact, rtol=1d-2)

contains
    subroutine single_orbital(beta, u, mu, vk, epsk, states, hloc_eb, diagram)
        use MBath
        use MHybridizationD
        use MWormState
        use MPsi
        real(c_double), intent(in) :: beta, u, mu, vk(:), epsk(:)
        type(TStates), intent(out) :: states
        type(TOperEigen), intent(out) :: hloc_eb
        type(TDiagram), intent(out) :: diagram

        class(DHybridization), allocatable :: hybr
        type(TBath) :: bath
        type(TImpurity) :: imp
        type(TPsis) :: psi
        type(TOperator) :: hloc
        type(TWormState) :: worm
        real(c_double) :: vkfull(1, 2, 1, 2, size(vk))
        !real(c_double) :: ftau(1, 2, 1, 2, 0:100)
        complex(c_double_complex) :: u_matrix(2, 2, 2, 2), muimp(1, 2, 1, 2)
        !integer :: i

        vkfull(:,:,:,:,:) = 0
        vkfull(1,1,1,1,:) = -vk**2
        vkfull(1,2,1,2,:) = -vk**2
        call init_log_hybridization(hybr, vkfull, -epsk, beta, .true.)
        ! call discretize_propagator(vkfull, -epsk, beta, .true., ftau)
        !do i = 0, ubound(ftau,5)
        !    write (99, '(2F16.6)') i * beta / ubound(ftau,5), ftau(1, 1, 1, 1, i)
        !enddo
        !call init_linear_hybridization(hybr, ftau, beta, .true.)
        call init_bath(bath, hybr)

        call init_states(states, 1, TQNSpec(Nt=.true., Szt=.true., All=.true.))
        call init_psi(psi, states, .false.)

        u_matrix(:,:,:,:) = 0
        u_matrix(1,2,1,2) = u
        u_matrix(2,1,2,1) = u
        muimp(:,:,:,:) = 0
        muimp(1,1,1,1) = -mu
        muimp(1,2,1,2) = -mu
        call make_hamiltonian(muimp, u_matrix, states, psi, hloc, 0.d0)
        call diag_operator(hloc, hloc_eb, 0.0d0)
        call transform_psis(hloc_eb, psi)
        call init_impurity(imp, beta, hloc_eb, psi, 0)

        call enter_z_sector(worm)
        call init_diagram(diagram, imp, bath, worm, .false.)
    end subroutine

end program
