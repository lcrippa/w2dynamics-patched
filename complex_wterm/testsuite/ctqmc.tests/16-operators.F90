program test_operators
    use iso_c_binding, only: c_double, c_double_complex
    use MExtRange
    use MImpurityD
    use MImpurityZ
    use MOperator
    use MPrinting
    use MPsi
    use MRandomNumbers
    use MStates
    use Testing
    implicit none

    call check_op_lifecycle()
    call check_partition_fn(.false.)
    call check_partition_fn(.true.)
    call check_connect_list()
    call check_adjoint()
    call check_integral()

contains

    subroutine init_problem(is_complex, psi, ham_eb, states)
        logical, intent(in) :: is_complex
        type(TPsis), intent(out) :: psi
        type(TOperEigen), intent(out) :: ham_eb
        type(TStates), intent(out) :: states

        class(RandomEngine), allocatable :: rng
        type(TQNSpec) :: qn_spec
        type(TOperator) :: ham

        complex(c_double_complex) :: muimp(3,2,3,2), u_matrix(6,6,6,6)
        integer :: i, j, s

        call init_xorshift_star(rng)
        muimp = 0
        do i = 1, 3; do j = 1, 3; do s = 1, 2
            muimp(i, s, j, s) = random_real(rng)
        enddo; enddo; enddo
        u_matrix = 0

        qn_spec%Nt = .true.
        qn_spec%Szt = .true.
        call init_states(states, 3, qn_spec)
        call init_psi(psi, states, is_complex)
        call make_hamiltonian(muimp, u_matrix, states, psi, ham, 0.0d0)
        call diag_operator(ham, ham_eb)
        call transform_psis(ham_eb, psi)
    end

    subroutine check_op_lifecycle()
        type(TOperator) :: op
        type(TPsis) :: psi
        type(TOperEigen) :: ham_eb
        type(TStates) :: states

        call init_problem(.false., psi, ham_eb, states)
        call assert_equal(is_initialized(op), .false.)
        call init_operator_like(op, ham_eb)
        call assert_equal(is_initialized(op), .true.)
    end

    subroutine check_partition_fn(is_complex)
        logical, intent(in) :: is_complex
        type(TPsis) :: psi
        type(TOperEigen) :: ham_eb
        type(TStates) :: states

        call init_problem(is_complex, psi, ham_eb, states)
        write (0, "('E0 =', ES20.10)") get_egs(ham_eb)

        if (is_complex) then
            call check_partition_fn_complex(psi%zpsis, ham_eb%zeig)
        else
            call check_partition_fn_real(psi%dpsis, ham_eb%deig)
        endif
    end

    subroutine check_partition_fn_complex(psi, ham_eb)
        type(ZTPsis), intent(in) :: psi
        type(ZOperEigen), intent(in) :: ham_eb

        type(ZExtRange) :: tr_naive, tr_naive_full, tr_expl, tr_expl_full, tr_imp, tr_fn
        type(ZTKet) :: ket
        type(ZTBra) :: bra
        type(ZImpurity) :: imp
        integer :: sst
        real(c_double), allocatable :: weights(:)
        real(c_double), parameter :: beta = 14.5

        tr_naive_full = cmplx(0.0, kind=c_double_complex)
        do sst = 0, get_num_superstates(ham_eb)-1
            call set_to_identity(ket, ham_eb, sst)
            call set_to_identity(bra, ham_eb, sst)
            call time_evolve(ham_eb, beta, ket)
            tr_naive = braket_dotprod(bra, ket)
            tr_naive_full = tr_naive_full + tr_naive
            tr_expl = partition_function(ham_eb, beta, sst)
            call assert_close(tr_naive, tr_expl, rtol=1d-14)

            call init_impurity(imp, beta, ham_eb, psi, sst)
            tr_imp = get_weight(imp)
            call assert_close(tr_imp, tr_expl, rtol=1d-14)

            call set_to_identity(ket, ham_eb, sst)
            call time_evolve(ham_eb, 0.75*beta, ket)
            tr_fn = trace_boltzmann(ham_eb, 0.25*beta, ket)
            call assert_close(tr_fn, tr_expl, rtol=1d-14)
        enddo

        tr_expl_full = partition_function(ham_eb, beta)
        call assert_close(tr_naive_full, tr_expl_full, rtol=1d-14)
        call assert_close(sum(thermal_weights(ham_eb, beta)), 1.0d0)

        weights = thermal_weights(ham_eb, beta)
        call assert_close(sum(weights), 1.0d0)
        do sst = 0, get_num_superstates(ham_eb)-1
            tr_expl = partition_function(ham_eb, beta, sst)
            call assert_close( &
                    weights(sst+1), real(limrange(tr_expl / tr_expl_full)), atol=1d-15)
        enddo
    end subroutine

    subroutine check_partition_fn_real(psi, ham_eb)
        type(DTPsis), intent(in) :: psi
        type(DOperEigen), intent(in) :: ham_eb

        type(DExtRange) :: tr_naive, tr_naive_full, tr_expl, tr_expl_full, tr_imp, tr_fn
        type(DTKet) :: ket
        type(DTBra) :: bra
        type(DImpurity) :: imp
        integer :: sst
        real(c_double), allocatable :: weights(:)
        real(c_double), parameter :: beta = 14.5

        tr_naive_full = 0.0d0
        do sst = 0, get_num_superstates(ham_eb)-1
            call set_to_identity(ket, ham_eb, sst)
            call set_to_identity(bra, ham_eb, sst)
            call time_evolve(ham_eb, beta, ket)
            tr_naive = braket_dotprod(bra, ket)
            tr_naive_full = tr_naive_full + tr_naive
            tr_expl = partition_function(ham_eb, beta, sst)
            call assert_close(tr_naive, tr_expl, rtol=1d-14)

            call init_impurity(imp, beta, ham_eb, psi, sst)
            tr_imp = get_weight(imp)
            call assert_close(tr_imp, tr_expl, rtol=1d-14)

            call set_to_identity(ket, ham_eb, sst)
            tr_fn = trace_boltzmann(ham_eb, beta, ket)
            call assert_close(tr_fn, tr_expl, rtol=1d-14)
        enddo

        tr_expl_full = partition_function(ham_eb, beta)
        call assert_close(tr_naive_full, tr_expl_full, rtol=1d-14)

        weights = thermal_weights(ham_eb, beta)
        call assert_close(sum(weights), 1.0d0)
        do sst = 0, get_num_superstates(ham_eb)-1
            tr_expl = partition_function(ham_eb, beta, sst)
            call assert_close( &
                    weights(sst+1), limrange(tr_expl / tr_expl_full), atol=1d-15)
        enddo
    end subroutine

    subroutine check_integral()
        type(TPsis), target :: psi_any
        type(TOperEigen), target :: ham_eb_any
        type(TStates) :: states
        type(DTPsis), pointer :: psi
        type(DOperEigen), pointer :: ham_eb
        type(DTOperator) :: op, rho, rho_ref, work
        type(DTKet) :: ket, ket_t
        type(DTBra) :: bra, bra_t
        type(DExtRange) :: alpha

        real(c_double), parameter :: tau = 0.4
        integer, parameter :: n = 400
        integer :: i

        call init_problem(.false., psi_any, ham_eb_any, states)
        psi => psi_any%dpsis
        ham_eb => ham_eb_any%deig

        op = matmul(psi%psis(1, 1, 1), psi%psis(1, 1, 2))
        call set_to_operator(ket, op, 2)
        call dump(ket)

        op = matmul(psi%psis(2, 1, 1), psi%psis(2, 1, 2))
        call set_to_operator(bra, op, 2)
        call dump(bra)

        alpha = 3.0d0 * tau / n
        call block_diagonal_oper(rho_ref, states)
        do i = 0, n - 1
            call copy(ket, ket_t)
            call copy(bra, bra_t)
            call time_evolve(ham_eb, 1.0d0 * i / n * tau, ket_t)
            call time_evolve(ham_eb, 1.0d0 * (n - i) / n * tau, bra_t)
            call outer_product(rho_ref, bra_t, ket_t, alpha, 1.0d0)
        enddo

        alpha = 3.0d0
        call block_diagonal_oper(rho, states)
        call block_diagonal_oper(work, states)
        call outer_integral(rho, bra, ket, tau, ham_eb, alpha, work)

        call assert_close(getmatrix(rho), getmatrix(rho_ref), rtol=1d-3)
    end

    subroutine check_connect_list()
        type(TPsis), target :: psi_any
        type(TOperEigen), target :: ham_eb_any
        type(TStates) :: states
        type(DTPsis), pointer :: psi
        type(DOperEigen), pointer :: ham_eb

        integer :: conn(20), nconn

        call init_problem(.false., psi_any, ham_eb_any, states)
        psi => psi_any%dpsis
        ham_eb => ham_eb_any%deig

        call start_sst_list(ham_eb, conn, nconn)
        call assert_equal(nconn, 16)
        call connecting_sst_list(psi%psis(1,1,1), conn, nconn)
        call connecting_sst_list(psi%psis(1,2,1), conn, nconn)
        call connecting_sst_list(psi%psis(2,1,1), conn, nconn)
        call connecting_sst_list(psi%psis(2,2,1), conn, nconn)
        call connecting_sst_list(psi%psis(3,1,1), conn, nconn)
        call connecting_sst_list(psi%psis(3,2,1), conn, nconn)
        call connecting_sst_list(psi%psis(1,1,1), conn, nconn)
        call finish_sst_list(ham_eb, conn, nconn)
        call assert_equal(nconn, 0)
    end subroutine

    subroutine check_adjoint()
        type(TPsis), target :: psi_any
        type(TOperEigen), target :: ham_eb_any
        type(TStates) :: states
        type(DTPsis), pointer :: psi
        type(DOperEigen), pointer :: ham_eb
        type(DTOperator) :: op
        type(DTBra) :: bra, bra2, bra3
        type(DTKet) :: ket, ket2, ket3
        type(DExtRange) :: r1, r2

        call init_problem(.false., psi_any, ham_eb_any, states)
        psi => psi_any%dpsis
        ham_eb => ham_eb_any%deig

        op = matmul(psi%psis(1, 1, 1), psi%psis(3, 1, 2))
        call set_to_identity(bra, ham_eb, 2)
        call time_evolve(ham_eb, 1.25d0, bra)
        call set_to_identity(ket, ham_eb, 2)
        call apply_operator(ket2, op, ket)
        r1 = limrange(braket_dotprod(bra, ket2))

        call set_to_identity(ket, ham_eb, 2)
        call time_evolve(ham_eb, 1.25d0, ket)
        call set_to_identity(bra, ham_eb, 2)
        call apply_adjoint(bra2, adjoint(adjoint(op)), bra)
        r2 = limrange(braket_dotprod(bra2, ket))

        call assert_close(r1, r2)

        ! -- check set_to_operator --

        call set_to_operator(ket, psi%psis(1, 1, 1), 2)
        call relocate_state(ham_eb, 2, 4.0d0, ket)

        call set_to_identity(ket2, ham_eb, 2)
        call time_evolve(ham_eb, 4.0d0, ket2)
        call apply_operator(ket3, psi%psis(1, 1, 1), ket2)
        call assert_close(ket3, ket, 1d-15)

        ! -- for bras

        call set_to_operator(bra, psi%psis(1, 1, 1), 2)
        call relocate_state(ham_eb, 2, 4.0d0, bra)

        call set_to_identity(bra2, ham_eb, 2)
        call time_evolve(ham_eb, 4.0d0, bra2)
        call apply_adjoint(bra3, psi%psis(1, 1, 1), bra2)
        call assert_close(bra3, bra, 1d-15)

    end subroutine

end program