program test_operators
    use iso_c_binding, only: c_double, c_double_complex
    use MExtRange
    use MImpurity
    use MImpurityD
    use MImpurityZ
    use MOperator
    use MPsi
    use MPrinting
    use MRandomNumbers
    use MStates
    use Testing
    implicit none

    class(RandomEngine), allocatable :: rng
    type(TImpurity) :: imp

    call init_xorshift_star(rng)
    call init_test_impurity(imp, 2, .true., .false.)
    call test_add_remove(imp)

    call init_xorshift_star(rng)
    call init_test_impurity(imp, 2, .true., .true.)
    call test_add_remove(imp)

    call init_xorshift_star(rng)
    call init_test_impurity(imp, 3, .false., .false.)
    call test_add_remove(imp)

contains

    subroutine init_test_impurity(imp, norb, dd, is_complex)
        type(TImpurity), intent(out) :: imp
        integer, intent(in) :: norb
        logical, intent(in) :: dd, is_complex

        type(TQNSpec) :: qn_spec
        type(TStates) :: states
        type(TOperator) :: ham
        type(TPsis) :: psi
        type(TOperEigen) :: ham_eb

        complex(c_double_complex) :: muimp(3,2,3,2)
        complex(c_double_complex) :: uval(6,6,6,6)
        integer :: i, j, s

        muimp = 0
        uval = 0
        do i = 0, norb-1
            uval(2*i+1, 2*i+2, 2*i+1, 2*i+2) = random_real(rng)
        enddo
        do i = 1, norb; do s = 1, 2
            muimp(i, s, i, s) = random_real(rng)
        enddo; enddo
        if (.not. dd) then
            do i = 1, norb; do j = i+1, norb; do s = 1, 2
                if (is_complex) then
                    muimp(i, s, j, s) = &
                        cmplx(random_real(rng), random_real(rng), c_double_complex)
                else
                    muimp(i, s, j, s) = random_real(rng)
                endif
                muimp(j, s, i, s) = conjg(muimp(i, s, j, s))
            enddo; enddo; enddo
        endif

        qn_spec%Nt = .true.
        qn_spec%Szt = .true.
        qn_spec%Azt = dd
        call init_states(states, norb, qn_spec)
        call init_psi(psi, states, is_complex)
        call make_hamiltonian(muimp, uval, states, psi, ham, 0.0d0)
        call diag_operator(ham, ham_eb)
        call transform_psis(ham_eb, psi)
        call init_impurity(imp, 14.0d0, ham_eb, psi, 0)
        write (0, "('E0 =', ES23.15)") get_egs(ham_eb)
    end

    subroutine test_add_remove(imp)
        type(TImpurity), intent(inout) :: imp
        class(TImpurityUpdateBase), pointer :: move
        type(TImpurityGrow), target :: grow_move
        type(TImpurityShrink), target :: shrink_move
        type(TImpurityShift), target :: shift_move
        type(TImpurityReplace), target :: permute_move
        type(TImpurityUpdate), target :: update_move

        ! XXX understand where the loss of precision comes from :/
        real(c_double), parameter :: rtol = 1e5 * epsilon(1.0d0)

        integer, parameter :: NMOVE_TYPES = 5
        complex(c_double_complex), parameter :: ONE = 1
        integer :: istep, move_type, rank, rempos(4), sst
        integer :: accepts(0:NMOVE_TYPES-1)
        integer, allocatable :: perm(:)
        integer, pointer :: outer(:)
        real(c_double) :: delta_tau
        logical :: caflip
        type(TLocalOper) :: newops(4)
        type(ZExtRange) :: wr
        type(TOperator) :: rho_op

        accepts(:) = 0
        allocate(perm(2 * get_nbands(imp)))
        call verify_impurity(imp, rtol)
        call init_grow_impurity(grow_move, imp)
        call init_shrink_impurity(shrink_move, imp)
        call init_shift_impurity(shift_move, imp)
        call init_replace_impurity(permute_move, imp)
        call init_update_impurity(update_move, imp)

        do istep = 1, 100000
            move_type = random_integer(rng, NMOVE_TYPES)
            rank = 2 * random_integer(rng, 1, 2)
            if (move_type == 0) then
                !write (0, "(/, I0, '+', I0, ': ')", advance="no") size(imp), rank
                call random_op(get_nbands(imp), get_beta(imp), newops(:rank))
                call propose_grow_impurity(grow_move, newops(:rank))
                move => grow_move
            else if (move_type == 1) then
                !write (0, "(/, I0, '-', I0, ': ')", advance="no") size(imp), rank
                if (rank > size(imp)) cycle
                call random_sample(rng, size(imp), rank, rempos(:rank))
                !write (0, "(10(I3, ', '))", advance="no") rempos(:rank)
                call propose_shrink_impurity(shrink_move, rempos(:rank))
                move => shrink_move
            else if (move_type == 2) then
                !write (0, "(/, I0, '>: ')", advance="no") size(imp)
                delta_tau = random_real(rng, get_beta(imp))
                call propose_shift_impurity(shift_move, delta_tau)
                !call assert_close(weight_ratio(shift_move), 1.0d0, rtol)
                move => shift_move
            else if (move_type == 3) then
                call random_permutation(rng, size(perm), perm)
                caflip = random_logical(rng)
                call propose_permute_impurity_flavours(permute_move, caflip, perm)
                move => permute_move
            else if (move_type == 4) then
                if (rank > size(imp)) cycle
                call random_op(get_nbands(imp), get_beta(imp), newops(:rank))
                call random_sample(rng, size(imp), rank, rempos(:rank))
                call propose_update_impurity(update_move, rempos(:rank), newops(:rank))
                move => update_move
            else
                error stop 'XXX'
            endif

            if (move_type == 3) then
                ! Only for permuter moves
                outer => allowed_outer_states(permute_move)
                if (size(outer) == 0) &
                    cycle

                sst = outer(random_integer(rng, 1, size(outer)))
                call choose_outer_state(permute_move, sst)
            else if (move_type == 2) then
                continue
            else
                if (random_real(rng) <= 0.2 .and. move_type /= 4) then
                    call choose_matching_outer_state(move, .true.)
                else
                    call choose_matching_outer_state(move, .false.)
                endif
            endif

            if (.not. is_update_allowed(move)) &
                cycle

            call compute_weight(move)
            wr = weight_ratio(move)
            !write (0, "('wr =')", advance="no")
            !call dump(wr)

            if (limrange(wr) == 0) &
                cycle

            if (random_integer(rng, 2) < 1) &
                cycle

            call accept_update(move)
            call verify_impurity(imp, rtol)
            accepts(move_type) = accepts(move_type) + 1

            call compute_density_matrix(imp, ONE, rho_op)
        enddo
        call dump(imp)
        call print_array(accepts, name='accepts')
        call assert_equal(all(accepts > 15), .true.)
    end

    subroutine random_op(norb, beta, spec)
        integer, intent(in) :: norb
        real(c_double), intent(in) :: beta
        type(TLocalOper), intent(out) :: spec(:)
        integer :: i

        do i = 1, size(spec)
            spec(i)%orb = random_integer(rng, 1, norb)
            spec(i)%sp = random_integer(rng, 1, 2)
            spec(i)%tau = random_real(rng, beta)
            spec(i)%type = 2*random_integer(rng, 0, 1) - 1
        enddo
    end

end program
