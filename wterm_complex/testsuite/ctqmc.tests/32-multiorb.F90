program test
    use iso_c_binding
    use MDiagram
    use MDiagramMoves
    use MImpurity
    use MMoveStatistics
    use MOperator
    use MPrinting
    use MPsi
    use MRandomNumbers
    use MStates
    use MSymmetryMoves
    use MWindow
    use MWormState, only: SecZ
    use Testing
    implicit none

    integer, parameter :: norb = 3
    real(c_double), parameter :: beta = 25
    real(c_double), parameter :: u = 3.4, j = 0.47, v = 2.37, mu = -0.73798412
    real(c_double), parameter :: ftau(100) = (/ &
            0.28658684, 0.21109435, 0.15918748, 0.12290588, 0.0971054 ,&
            0.0784298 , 0.06466584, 0.05433736, 0.04644801, 0.04031689,&
            0.03547261, 0.03158452, 0.02841759, 0.02580246, 0.0236155 ,&
            0.02176525, 0.02018317, 0.01881726, 0.01762761, 0.01658318,&
            0.01565965, 0.01483767, 0.01410177, 0.01343941, 0.01284038,&
            0.01229626, 0.01180006, 0.01134594, 0.01092901, 0.0105451 ,&
            0.01019065, 0.00986261, 0.00955838, 0.00927565, 0.00901245,&
            0.00876704, 0.00853789, 0.00832363, 0.00812308, 0.00793516,&
            0.00775891, 0.00759349, 0.00743811, 0.00729209, 0.0071548 ,&
            0.00702566, 0.00690417, 0.00678983, 0.00668223, 0.00658097,&
            0.00648569, 0.00639605, 0.00631175, 0.00623251, 0.00615807,&
            0.00608818, 0.00602264, 0.00596122, 0.00590375, 0.00585004,&
            0.00579994, 0.00575328, 0.00570993, 0.00566976, 0.00563264,&
            0.00559847, 0.00556713, 0.00553852, 0.00551256, 0.00548916,&
            0.00546823, 0.00544971, 0.00543353, 0.00541961, 0.00540791,&
            0.00539837, 0.00539093, 0.00538554, 0.00538217, 0.00538078,&
            0.00538131, 0.00538375, 0.00538805, 0.0053942 , 0.00540216,&
            0.0054119 , 0.00542342, 0.00543668, 0.00545167, 0.00546838,&
            0.00548679, 0.00550689, 0.00552867, 0.00555212, 0.00557725,&
            0.00560403, 0.00563248, 0.00566258, 0.00569435, 0.00572782 /)
    integer, parameter :: nmeas = 25000, nwarmup = 1000, ncorr = 20
    integer, parameter :: user_sym_moves(1, 6) = reshape( &
                (/ 2, 1, 3, 5, 6, 4 /), &
                (/ 1, 6 /))

    real(c_double), parameter :: perc_outer = 0.21d0

    type(TStates) :: states
    type(TPsis) :: psi
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
    real(c_double), allocatable :: rho1(:, :)

    write (0,*) 'Constructing diagram ...'
    call diag_orbitals(beta, norb, u, v, j, mu, ftau, states, psi, hloc, diagram)

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
    rho1 = rho1_from_density_matrix(mean_rho%dop, hloc%deig, psi%dpsis) / real(mean_sign)
    call print_array(rho1, name='rho1')

contains

    function umatrix_kanamori(n, u, v, j) result(umat)
        integer, intent(in) :: n
        real(c_double), intent(in) :: u, v, j
        real(c_double), allocatable :: umat(:,:,:,:)
        integer :: b1, b2
        integer, parameter :: up = 1, dn = 2

        allocate(umat(2*n, 2*n, 2*n, 2*n))
        do b1 = 0, 2*n-1, 2
            do b2 = 0, 2*n-1, 2
                if (b1 == b2) then
                    umat(b1+up, b1+dn, b1+up, b1+dn) = u
                    umat(b1+dn, b1+up, b1+dn, b1+up) = u
                else
                    ! density-density
                    umat(b1+up, b1+up, b2+up, b2+up) = v
                    umat(b1+up, b1+dn, b2+up, b2+dn) = v
                    umat(b1+dn, b1+up, b2+dn, b2+up) = v
                    umat(b1+dn, b1+dn, b2+dn, b2+dn) = v

                    umat(b1+up, b2+up, b2+up, b1+up) = j
                    umat(b1+dn, b2+dn, b2+dn, b1+dn) = j

                    ! pair hopping
                    umat(b1+up, b1+dn, b2+up, b2+dn) = j
                    umat(b1+dn, b1+up, b2+dn, b2+up) = j

                    ! spin flip
                    umat(b1+dn, b2+up, b2+dn, b1+up) = j
                    umat(b1+up, b2+dn, b2+up, b1+dn) = j
                endif
            enddo
        enddo
    end function

    subroutine diag_orbitals( &
                    beta, norb, u, v, j, mu, ftau, states, psi, hloc_eb, diagram)
        use MBath
        use MHybridizationD
        use MWormState
        integer, intent(in) :: norb
        real(c_double), intent(in) :: beta, u, v, j, mu, ftau(:)
        type(TStates), intent(out) :: states
        type(TPsis), intent(out) :: psi
        type(TOperEigen), intent(out) :: hloc_eb
        type(TDiagram), intent(out) :: diagram

        class(DHybridization), allocatable :: hybr
        type(TBath) :: bath
        type(TImpurity) :: imp
        type(TOperator) :: hloc
        type(TWormState) :: worm
        real(c_double) :: ftau_full(norb, 2, norb, 2, size(ftau))
        complex(c_double_complex) :: u_matrix(2*norb, 2*norb, 2*norb, 2*norb)
        complex(c_double_complex) :: muimp(norb, 2, norb, 2)

        integer :: i, s

        do i = 1, norb
            do s = 1, 2
                ftau_full(i, s, i, s, :) = ftau
            enddo
        enddo
        call init_linear_hybridization(hybr, ftau_full, beta, .true.)
        call init_bath(bath, hybr)

        call init_states(states, norb, TQNSpec(Nt=.true., Szt=.true., Qzt=.true.))
        call init_psi(psi, states, .false.)

        u_matrix = umatrix_kanamori(norb, u, v, j)
        do i = 1, norb
            do s = 1, 2
                muimp(i, s, i, s) = -mu
            enddo
        enddo

        call make_hamiltonian(muimp, u_matrix, states, psi, hloc, 0.d0)
        call diag_operator(hloc, hloc_eb, 0.0d0)
        call transform_psis(hloc_eb, psi)
        call init_impurity(imp, beta, hloc_eb, psi, 0)

        call enter_z_sector(worm)
        call init_diagram(diagram, imp, bath, worm, .false.)
    end subroutine

end program
