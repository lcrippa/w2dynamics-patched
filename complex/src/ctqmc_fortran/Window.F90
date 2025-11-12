!> The idea between sliding windows is to have a small window in imaginary
!! time, which slides back and forth every so many steps:
!!
!!     pos   |.....:.....:.....:.....:.....|     window size
!!       1   |...........|                      [       0, 0.4*beta]
!!       2         |...........|                [0.2*beta, 0.6*beta]
!!       3               |...........|          [0.4*beta, 0.8*beta]
!!       4                     |...........|    [0.6*beta,     beta]
!!       5               |...........|          [0.4*beta, 0.8*beta]
!!       6         |...........|                [0.2*beta, 0.6*beta]
!!     ...    ...
!!
!! (Non-global) updates are then confined only to within the window. The
!! idea is that most updates in CT-HYB are short-lived excursions into
!! high-lying states, which are well-covered by sliding windows [1].
!!
!! [1] H. Shinaoka et. al, J. Stat. Mech. P06012 (2014).
!!
module MWindow
    use iso_c_binding, only: c_double
    use MCommon
    use MDiagram
    use MLocalBase, only: isworm, TLocalOper
    use MPermutation
    use MRandomNumbers
    use MSortingD
    implicit none
    private

    public :: TLocalOper

    !> A type containing the configuration and state of the sliding
    !> window.
    type, public :: TWindowState
        private
        !> The number of positions for the sliding window (which is
        !> _not_ the divisor NWin for the window width as in width =
        !> beta/NWin).
        integer     :: num_windows = 0

        !> The current window position expressed as number of half
        !> window widths from tau = 0 to the beginning of the window,
        !> which is always at least 0 and less than num_windows.
        integer     :: window_position = -HUGE(1)

        !> The maximum tau difference between operators for which a pair
        !> insertion or removal may be proposed.
        real(c_double) :: taudiff_max = -HUGE(1.0d0)

        !> Inverse temperature
        real(c_double) :: beta = -HUGE(1.0d0)

        !> The lower and upper edge of the sliding window in the current
        !> position.
        real(c_double) :: tauwin_min = -HUGE(1.0d0), tauwin_max = -HUGE(1.0d0)

        !> Number of bands
        integer :: nbands

        !> Consider only pairs which are diagonal in spin/orbital
        logical :: diag_only
    end type TWindowState

    public :: init_window, set_max_pair_distance, set_window, shift_window
    public :: num_window_positions, max_pair_distance
    public :: tau_range, window_range
    public :: proposal_probability_for_add, proposal_probability_for_remove
    public :: generate_add, generate_remove, forms_pair

contains

    !> (Re-)initialize the sliding window.
    subroutine init_window(ws, diagram, taudiff_max)
        type(TWindowState), intent(out) :: ws
        type(TDiagram), intent(in) :: diagram
        real(c_double), intent(in)  :: taudiff_max

        ws%beta = diagram%beta
        ws%nbands = diagram%nbands
        ws%diag_only = .not. diagram%b_offdiag
        call set_max_pair_distance(ws, taudiff_max)
    end subroutine init_window

    !> Set maximum pair distance to taudiff_max and reset window
    subroutine set_max_pair_distance(self, taudiff_max)
        type(TWindowState), intent(inout) :: self
        real(c_double), intent(in)  :: taudiff_max

        if (.not. (taudiff_max > 0)) &
            error stop 'taudiff_max must be positive'

        self%taudiff_max = min(taudiff_max, self%beta)
        self%num_windows = max(floor(self%beta / self%taudiff_max) - 1, 1)

        ! Reset window to apply new pair distance
        call set_window(self, 0)
    end subroutine

    !> Adjust related window state to a newly set position.
    !  window_position must be a valid position
    pure subroutine set_window(ws, window_position)
        type(TWindowState), intent(inout) :: ws
        integer, intent(in)           :: window_position

        ws%window_position = max(window_position, 0)
        if (ws%window_position >= ws%num_windows) then
            ws%window_position = 0
        end if

        ws%tauwin_min = (dble(ws%window_position)/ &
                (dble(ws%num_windows + 1))) * ws%beta
        ws%tauwin_max = (dble(ws%window_position + 2) / &
                (dble(ws%num_windows + 1))) * ws%beta
    end subroutine set_window

    !> Shift the sliding window to the next position.
    pure subroutine shift_window(ws)
        type(TWindowState), intent(inout) :: ws

        call set_window(ws, ws%window_position + 1)
    end subroutine shift_window

    !> Number of distinct positions of the sliding window
    pure integer function num_window_positions(self)
        type(TWindowState), intent(in) :: self

        num_window_positions = self%num_windows
    end function

    !> Maximum tau distance between pairs considered for insertion/removal
    pure function max_pair_distance(self) result(delta_tau)
        type(TWindowState), intent(in) :: self
        real(c_double) :: delta_tau

        delta_tau = self%taudiff_max
    end function

    ! -------------------------------------------------------------------------
    ! PAIR COUNTING FUNCTIONS

    !> Check whether two operators form a pair matching the window
    !!
    !! Two operators form a pair if they both reside inside the current window,
    !! if they are close enough together, if both of them have hybridization
    !! lines, and one of them is a creator and the othe one is an annihilator.
    !! If diag_only is set, then their flavours must also match.
    pure logical function forms_pair(op1, op2, ws) result(p)
        type(TLocalOper), intent(in) :: op1, op2
        type(TWindowState), intent(in) :: ws

        p = .false.

        ! Check times
        if (op1%tau > ws%tauwin_max) return
        if (op2%tau > ws%tauwin_max) return

        if (op1%tau < ws%tauwin_min) return
        if (op2%tau < ws%tauwin_min) return

        if (op2%tau > op1%tau + ws%taudiff_max) return
        if (op1%tau > op2%tau + ws%taudiff_max) return

        ! Check worm
        if (isworm(op1)) return

        ! Check rest
        p = forms_operator_pair(op1, op2, ws%diag_only)
    end function

    !> Choose specific pair of creators/annihilators in window except in list.
    !!
    !! Builds the set of pairs of operators (one creator and one annihilator,
    !! each with a hybridization line). To form a pair, both operators must
    !! reside inside the current window, be close enough together, and neither
    !! must have an index in the exclude list. If diag_only is set, then they
    !! also must have the same flavour.
    !!
    !! This function chooses the ipair'th pair of operators (A1(tau1), A2(tau2))
    !! where tau1 < tau2 and the pairs are ordered lexicographically.
    subroutine choose_pair_except(ipair, ops, ws, exclude, pair)
        integer, intent(in)            :: ipair
        type(TLocalOper), intent(in)   :: ops(:)
        type(TWindowState), intent(in) :: ws
        integer, intent(in)            :: exclude(:)
        integer, intent(out)           :: pair(2)
        integer                        :: i, j, imin, imax, jmax, icur

        icur = 0
        call indices_in_tau_range(ops, ws%tauwin_min, ws%tauwin_max, imin, imax)

        do i = imin, imax
            if (isworm(ops(i))) &
                cycle
            if (is_element_of(i, exclude)) &
                cycle

            jmax = end_of_range(ops(i:imax), ws%taudiff_max) + i - 1
            do j = i + 1, jmax
                if (is_element_of(j, exclude)) &
                    cycle
                if (.not. forms_operator_pair(ops(i), ops(j), ws%diag_only)) &
                    cycle

                icur = icur + 1
                if (icur == ipair) then
                    pair(:) = (/ i, j /)
                    return
                endif
            enddo
        enddo
        error stop 'No pair with that index exists'
    end

    !> Return set of indices {imin, ..., imax} of operators inside window
    subroutine indices_in_tau_range(ops, tau_min, tau_max, imin, imax)
        type(TLocalOper), intent(in) :: ops(:)
        real(c_double), intent(in) :: tau_min, tau_max
        integer, intent(out) :: imin, imax
        real(c_double) :: tau(size(ops))

        tau(:) = ops(:)%tau
        imin = get_sorted_insert(tau, tau_min)
        imax = get_sorted_insert(tau(imin:), tau_max) + imin - 1
        if (imax > size(ops)) then
           imax = imax - 1
        else if (ops(imax)%tau /= tau_max) then
           imax = imax - 1
        endif
    end subroutine

    ! XXX: work with Move objects once they've been sanitized

    !> Generate N creation and annihilation operators at random times
    !! inside of the window with a maximum distance of taudiffmax between
    !! each pair (1&2, 3&4, ...). Operators are stored in (c, cdag) order.
    !!
    !! @param ws Window state
    !! @param rng Random number generator
    !! @param N Number of operators to generate
    !! @param ops Array where the operators are stored
    subroutine generate_add(ws, rng, npairs, ops)
        type(TWindowState), intent(in)      :: ws
        class(RandomEngine), intent(inout)  :: rng
        integer, intent(in)                 :: npairs
        type(TLocalOper), intent(out)       :: ops(2*npairs)

        real(c_double)                         :: tau_min, tau_max
        integer                             :: i
        type(TLocalOper)                    :: op

        do i = 1, npairs
            ! first of pair, random in window and used as reference
            call window_range(ws, tau_min, tau_max)
            op%type = OpAnnh
            op%tau = random_real(rng, tau_min, tau_max)
            op%orb = random_integer(rng, 1, ws%NBands)
            op%sp = random_integer(rng, 1, 2)
            ops(2*i-1) = op

            ! second of pair, random in window within taudiffmax of first
            call tau_range(ops(2*i - 1)%tau, ws, tau_min, tau_max)
            op%type = OpCrea
            op%tau = random_real(rng, tau_min, tau_max)
            if (.not. ws%diag_only) then
                ! Override annihilator with new flavours
                op%orb = random_integer(rng, 1, ws%NBands)
                op%sp = random_integer(rng, 1, 2)
            endif
            ops(2*i) = op
        end do
    end subroutine

    function proposal_probability_for_add(ops, ws) result(p)
        type(TLocalOper), intent(in) :: ops(:)
        type(TWindowState), intent(in) :: ws
        real(c_double) :: p, inv_p
        integer :: i, npairs, comb_fact

        npairs = get_npairs(ops, ws)

        ! The proposal probability is a bit tricky: we first choose the time
        ! t1 of the annihilation operator uniformly inside the [tmin, tmax]
        ! and then t2 uniformly inside tau_range(t1). Thus, if x1 and x2 are
        ! our random numbers [0,1), we have:
        !
        !     dx1 dx2 = abs(det(dx/dt)) dt1 dt2
        !             = dt1 dt2 / ((tmax - tmin) * tau_range(t1))
        !
        ! The proposal probabilities for the orbitals simply depend on whether
        ! we choose diagonal or general.
        inv_p = 1.0
        do i = 1, npairs
            inv_p = inv_p * tau_range_size(ops(2*i-1)%tau, ws)
            inv_p = inv_p * window_size(ws)
            inv_p = inv_p * (2 * ws%nbands)
            if (.not. ws%diag_only) &
                inv_p = inv_p * (2 * ws%nbands)
        enddo

        ! There is a labelling ambiguity: we can permute the order of choosing
        ! in the pairs, but since the resulting diagram is always understood
        ! to be time-ordered, this results in the same diagram. We correct for
        ! this here
        comb_fact = factorial(npairs)

        ! There is another potential ambiguity for multiple pairs: there may
        ! be more than one way to form legal pairs from the inserted
        ! operators. We have to correct for that.
        if (npairs == 2) then
            if (forms_pair(ops(1), ops(4), ws) .and. &
                forms_pair(ops(2), ops(3), ws) &
            ) then
                comb_fact = 2 * comb_fact
            endif
        elseif (npairs >= 3) then
            error stop 'Not supported yet'
        endif

        ! Multiply the weight of the current proposal by the number of
        ! equivalent proposals to arrive at the full probability
        p = comb_fact / inv_p
    end function

    !> Picks one randomly chosen removable operator pair of an
    !! annihilator and creator inside the window if possible and
    !! returns their indices in remidx.
    !!
    !! @param ws Window state
    !! @param ops Base diagram from which to move away from
    !! @param remidx Indices of the operators to be removes
    !! @param rng Random number generator
    subroutine generate_remove(ws, rng, ops, nrempairs, remidx, allowed)
        type(TWindowState), intent(in)  :: ws
        class(RandomEngine), intent(inout) :: rng
        type(TLocalOper), intent(in) :: ops(:)
        integer, intent(in) :: nrempairs
        integer, intent(out) :: remidx(2*nrempairs)
        logical, intent(out) :: allowed

        integer :: pairs(2, (size(ops)/2)**2), npairs(num_classes(ws))
        integer :: cl(size(ops)), npairs_total
        integer :: npairs_remaining, i, offset12, cl12

        call compute_pairs(ops, ws, cl, pairs, npairs, npairs_total)

        allowed = .false.

        ! Choose first pair
        if (nrempairs >= 1) then
            if (npairs_total == 0) &
                return

            i = random_integer(rng, 1, npairs_total)

            ! XXX reuse computed pairs
            ! Choose pair, making sure that annihilator is the first in each
            ! pair, consistent with generate_add.
            call choose_pair_except(i, ops, ws, remidx(1:0), remidx(1:2))
            if (ops(remidx(1))%type /= OpAnnh) &
                call swap(remidx(1), remidx(2))
        endif

        ! Choose second pair
        if (nrempairs >= 2) then
            cl12 = cl(remidx(1))
            offset12 = sum(npairs(:cl12-1))
            npairs_remaining = npairs_total - count_pairs_containing( &
                remidx(1:2), pairs(:, offset12+1:offset12+npairs(cl12)))
            if (npairs_remaining == 0) &
                return

            i = random_integer(rng, 1, npairs_remaining)
            call choose_pair_except(i, ops, ws, remidx(1:2), remidx(3:4))
            if (ops(remidx(3))%type /= OpAnnh) &
                call swap(remidx(3), remidx(4))
        endif

        if (nrempairs >= 3) &
            error stop 'Only up to 2 pairs supported for now'

        allowed = .true.
    end subroutine

    !> Computes the probability of choosing a remove move.
    function proposal_probability_for_remove(idx, ops, ws) result(p)
        integer, intent(in) :: idx(:)
        type(TLocalOper), intent(in) :: ops(:)
        type(TWindowState), intent(in) :: ws
        real(c_double) :: p

        integer :: pairs(2, (size(ops)/2)**2), npairs(num_classes(ws))
        integer :: cl(size(ops))
        integer :: npairs_total, cl12, cl34, offset12, offset34, k
        integer :: nexcl12, nexcl34, nexcl14, nexcl32

        k = get_npairs(ops(idx), ws)
        call compute_pairs(ops, ws, cl, pairs, npairs, npairs_total)

        ! The way we choose pairs is a bit peculiar: we choose the first pair,
        ! then for each of these choices, we choose the next pairs, and so on.
        ! Let `i, j, k ...` be the pairs we are choosing, then we have:
        !
        !     p(i,j,k,...) = p(i) * p(j | i) * p(k | i,j) ...
        !
        if (k == 1) then
            p = 1.0d0 / npairs_total
            !write (0, "('NEW',5I4)") npairs_total
        else if (k == 2) then
            cl12 = cl(idx(1))
            cl34 = cl(idx(3))

            ! Operators with indices (1,2) and (3,4) must each form a pair.
            ! Once we have chosen one of them with some probability, we must
            ! determine how likely it is to choose the other (thes probabilites
            ! may be different because of taudiff_max)
            offset12 = sum(npairs(:cl12-1))
            nexcl12 = npairs_total - count_pairs_containing( &
                        (/ idx(1), idx(2) /), &
                        pairs(:, offset12+1:offset12+npairs(cl12)))

            offset34 = sum(npairs(:cl34-1))
            nexcl34 = npairs_total - count_pairs_containing( &
                        (/ idx(3), idx(4) /), &
                        pairs(:, offset34+1:offset34+npairs(cl34)))

            ! It may be that (1,4) and (2,3) also form a pair. In this case
            ! this provides another possibility for ending up in the same
            ! configuration.
            if (cl12 == cl34 &
                .and. abs(ops(idx(1))%tau - ops(idx(4))%tau) < ws%taudiff_max &
                .and. abs(ops(idx(2))%tau - ops(idx(3))%tau) < ws%taudiff_max &
            ) then
                nexcl14 = npairs_total - count_pairs_containing( &
                                (/ idx(1), idx(4) /), &
                                pairs(:, offset12+1:offset12+npairs(cl12)))
                nexcl32 = npairs_total - count_pairs_containing( &
                                (/ idx(3), idx(2) /), &
                                pairs(:, offset12+1:offset12+npairs(cl12)))
            else
                nexcl14 = 0
                nexcl32 = 0
            endif

            p = 0
            if (nexcl12 /= 0) &
                p = 1.0d0 / nexcl12
            if (nexcl34 /= 0) &
                p = p + 1.0d0 / nexcl34
            if (nexcl14 /= 0) &
                p = p + 1.0d0 / nexcl14
            if (nexcl32 /= 0) &
                p = p + 1.0d0 / nexcl32

            p = p / npairs_total
        else
            error stop 'Not implemented (new)'
        endif
    end

    function count_pairs_containing(idx, pairs) result(npairs)
        integer, intent(in) :: idx(2), pairs(:,:)
        integer :: i, npairs

        !write (0, "(2I4, X, 100(I0, ',', I0, X))") idx, pairs
        npairs = 0
        do i = 1, size(pairs,2)
            if (pairs(1,i) == idx(1) .or. pairs(2,i) == idx(2)) &
                npairs = npairs + 1
        enddo
    end

    !> Verify that the operators form proper pairs and return their count
    function get_npairs(ops, ws) result(npairs)
        type(TLocalOper), intent(in) :: ops(:)
        type(TWindowState), intent(in) :: ws
        integer :: i, npairs

        npairs = size(ops) / 2
        if (size(ops) /= 2 * npairs) &
            error stop 'Invalid number of operators'

        do i = 1, npairs
            if (ops(2*i - 1)%type /= OpAnnh) &
                error stop 'Expecting annihilator'
            if (.not. forms_pair(ops(2*i - 1), ops(2*i), ws)) &
                error stop 'Expecting pair'
        enddo
    end function

    !> Get index i of last operator such that ops(1:i) spans at most tau_span
    integer function end_of_range(ops, tau_span) result(i)
        type(TLocalOper), intent(in) :: ops(:)
        real(c_double), intent(in) :: tau_span
        real(c_double) :: tau_end

        if (size(ops) == 0) &
            error stop 'Empty range'

        tau_end = ops(1)%tau + tau_span
        do i = 2, size(ops)
            if (ops(i)%tau > tau_end) &
                exit
        enddo
        i = i - 1
    end function

    !> Return if operators form a valid pair
    pure logical function forms_operator_pair(lhs, rhs, diag_only) result(t)
        type(TLocalOper), intent(in) :: lhs, rhs
        logical, intent(in) :: diag_only

        t = .true.
        if (lhs%type /= -rhs%type) &
            t = .false.
        if (diag_only) then
            if (lhs%orb /= rhs%orb .or. lhs%sp /= rhs%sp) &
                t = .false.
        endif
    end function

    !> Returns true if and only if x is equal to an element in list
    pure logical function is_element_of(x, list) result(t)
        integer, intent(in) :: x, list(:)
        integer :: i

        t = .false.
        do i = 1, size(list)
            if (x == list(i)) then
                t = .true.
                return
            end if
        end do
    end function

    ! -------------------------------------------------------------------------
    ! PAIR COMPUTATION

    subroutine compute_pairs(ops, ws, cl, pairs, npairs, npairs_total)
        type(TLocalOper), intent(in) :: ops(:)
        type(TWindowState), intent(in) :: ws
        integer, intent(out) :: cl(:), pairs(:,:), npairs(:), npairs_total

        integer :: i, cstart, cstop, astart, astop, nclasses
        integer, dimension(size(ops)) :: clperm
        integer :: clstop(0:2*num_classes(ws))

        ! First, we group operators into classes, where only operators with
        ! matching class can form pairs, and all operators which cannot form
        ! pairs are put into class 0. There is some linear effort here, which
        ! will pay off later.
        do i = 1, size(ops)
            cl(i) = operator_class(ops(i), ws) + 1
        enddo
        call group_perm(cl, clperm, clstop)
        cl(:) = cl(:) - 1

        ! Now, instead of checking every operator against every other, we
        ! only check operators in matching classes if they are close enough
        ! together.
        nclasses = num_classes(ws)
        astart = clstop(0) + 1
        cstart = clstop(nclasses) + 1
        npairs_total = 0
        do i = 1, nclasses
            astop = clstop(i)
            cstop = clstop(i + nclasses)
            call compute_pairs_class( &
                    ops, clperm(astart:astop), clperm(cstart:cstop), &
                    ws%taudiff_max, pairs(:, npairs_total+1:), npairs(i))

            astart = astop + 1
            cstart = cstop + 1
            npairs_total = npairs_total + npairs(i)
        enddo
    end

    subroutine compute_pairs_class(ops, aidx, cidx, taudiff_max, pairs, npairs)
        type(TLocalOper), intent(in) :: ops(:)
        real(c_double), intent(in) :: taudiff_max
        integer, intent(in) :: cidx(:), aidx(:)
        integer, intent(out) :: pairs(:,:), npairs

        integer :: c, a
        real(c_double) :: taudiff

        npairs = 0
        do c = 1, size(cidx)
            do a = 1, size(aidx)
                taudiff = ops(aidx(a))%tau - ops(cidx(c))%tau
                if (abs(taudiff) > taudiff_max) &
                    cycle

                npairs = npairs + 1
                pairs(1, npairs) = aidx(a)
                pairs(2, npairs) = cidx(c)
            enddo
        enddo
    end subroutine

    !> Returns the class of an operator for the purpose of pairing
    pure integer function operator_class(op, ws) result(cl)
        type(TLocalOper), intent(in) :: op
        type(TWindowState), intent(in) :: ws

        cl = 0
        if (.not. (op%tau >= ws%tauwin_min .and. op%tau <= ws%tauwin_max)) &
            return
        if (op%type /= OpCrea .and. op%type /= OpAnnh) &
            return

        if (op%type == OpCrea) then
            cl = 1
        endif
        if (ws%diag_only) then
            cl = ws%nbands * cl + (op%orb - 1)
            cl = 2 * cl + (op%sp - 1)
        endif
        cl = cl + 1
    end function

    pure integer function num_classes(ws)
        type(TWindowState), intent(in) :: ws

        if (ws%diag_only) then
            num_classes = 2 * ws%nbands
        else
            num_classes = 1
        endif
    end

    ! -------------------------------------------------------------------------
    ! WINDOW INQUIRY

    !> Returns the valid tau range around some point in the window.
    !!
    !! Given a time `tau`, which must lie inside the current window `ws`,
    !! compute the largest interval `[tau_min, tau_max]` which is still
    !! contained inside the window and not further than the maximum distance
    !! away from `tau`.
    subroutine tau_range(tau, ws, tau_min, tau_max)
        real(c_double), intent(in) :: tau
        type(TWindowState), intent(in) :: ws
        real(c_double), intent(out) :: tau_min, tau_max

        tau_min = tau - ws%taudiff_max
        if (tau_min < ws%tauwin_min) then
            if (tau < ws%tauwin_min) &
                error stop 'tau is too small for the window'
            tau_min = ws%tauwin_min
        endif

        tau_max = tau + ws%taudiff_max
        if (tau_max > ws%tauwin_max) then
            if (tau > ws%tauwin_max) &
                error stop 'tau is too big for the window'
            tau_max = ws%tauwin_max
        endif
    end subroutine

    !> Return size of tau range around some point in the window.
    function tau_range_size(tau, ws) result(r)
        real(c_double), intent(in) :: tau
        type(TWindowState), intent(in) :: ws
        real(c_double) :: r, tau_min, tau_max

        call tau_range(tau, ws, tau_min, tau_max)
        r = tau_max - tau_min
    end function

    !> Fill variables with the range [tau_min, tau_max] of current window
    pure subroutine window_range(ws, tau_min, tau_max)
        type(TWindowState), intent(in) :: ws
        real(c_double), intent(out) :: tau_min, tau_max

        tau_min = ws%tauwin_min
        tau_max = ws%tauwin_max
    end

    !> Return size of the current window
    pure function window_size(ws) result(r)
        type(TWindowState), intent(in) :: ws
        real(c_double) :: r

        r = ws%tauwin_max - ws%tauwin_min
    end function

    !> Factorial of non-negative integer
    pure elemental integer function factorial(n)
        integer, intent(in) :: n
        integer :: i

        factorial = max(n, 1)
        do i = n - 1, 2, -1
            factorial = factorial * i
        enddo
    end function
end module
