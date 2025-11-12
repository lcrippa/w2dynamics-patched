module MWormState
    use iso_c_binding, only: c_double
    use MCommon
    use MSorting
    use MLocalBase
    implicit none
    private

    enum, bind(c)
        enumerator :: SecDummy, SecZ=1, SecG, SecGSigma, SecG4, SecH4, SecP2
        enumerator :: SecP2pp, SecP3, SecP3pp, SecQQ, SecQ4, SecNQQdag, SecQQdd
        enumerator :: SecUcaca, SecUccaa, SecQUDdag, SecRaman, SecCustom
    end enum

    public :: SecDummy, SecZ, SecG, SecGSigma, SecG4, SecH4, SecP2
    public :: SecP2pp, SecP3, SecP3pp, SecQQ, SecQ4, SecNQQdag, SecQQdd
    public :: SecUcaca, SecUccaa, SecQUDdag, SecRaman, SecCustom

    !> number of worm sectors
    integer, parameter, public :: NWormSectors = SecCustom - SecZ + 1

    !> number of worm operators in specific sector
    integer, parameter, public :: NOperWorm(NWormSectors) = &
            (/ 0, 2, 4, 4, 6, 4, 4, 4, 4, 6, 12, 8, 8, 4, 4, 4, 6, -1/)

    type, public :: TCustomWormConfig
        private
        integer :: ntaus=-1, ntotalops=-1
        integer, allocatable :: ntauops(:)
        integer(kind(OpDummy)), allocatable :: optypes(:)
        integer, allocatable :: offsets(:)
    end type TCustomWormConfig

    type, public :: TWormState
        private
        !> Sector the current diagram belongs to
        integer(kind(SecDummy)) :: current_sector = SecDummy
        !> Component the current diagram belongs to
        integer :: current_component = -1
        type(TCustomWormConfig) :: custom_config(size(NOperWorm))

        !> Indices of the worm operators in the time-ordered local operator list,
        !! where index `i` refers to the `i`-th worm operator in the "natural"
        !! order (how they appear in the perturbation series expansion).
        !! E.g., for the worm `c(t) c+(t')`, locpos_typeorder(1) would be the
        !! position of the `c(t)` operator in the time-ordered local operator
        !! list.
        integer, allocatable :: locpos_typeorder(:)

        !> Indices of the worm operators in the time-ordered local operator list,
        !! where index `i` refers to the `i`-th worm operator in order of tau
        !! order, i.e., this is a sorted version of locpos_typeorder.
        integer, allocatable :: locpos_ownorder(:)

        !> 'indices' in the sector-dependent fixed order, in ascending local
        !! operator list index order, such that
        !!
        !!      locpos_ownorder(i) == locpos_typeorder(typepos_locposorder(i))
        !!
        !! in other words, this is the permutation that sorts the array
        !! `locpos_typeorder`.  Suppose, e.g., that  `typeorder = [5, 6, 7, 2]`,
        !! then `locposorder = [4, 1, 2, 3]` such that `ownorder = [2, 5, 6, 7]`.
        integer, allocatable :: typepos_locposorder(:)
    end type TWormState

    ! FIXME: remove / replace everything that is unnecessary, here and
    ! the implementations
    public :: enter_z_sector
    public :: get_current_worm_sector, get_current_worm_comp
    public :: get_worm_oper_positions_time_ordered, get_worm_oper_position
    public :: get_tau_oper_position
    public :: get_nworm, get_ntaus
    public :: verify_worm
    public :: get_custom_worm_ntaus, get_custom_worm_ntotalops
    public :: get_custom_worm_ntauops, get_custom_worm_optype
    public :: get_sector_nworm, get_sector_ntaus

    public :: wormstate_update_pos
    public :: assign_pos_from_locpos_typeorder
    public :: wormstate_taushift, wormstate_remove_hybops
    public :: update_wormstate_replace, update_wormstate_add_opers
    public :: update_wormstate_replace_eqtauop
    public :: wormstate_enter_wormsector, wormstate_leave_wormsector
    public :: natural_order_worm_perm

contains

    subroutine enter_z_sector(self)
        type(TWormState), intent(inout) :: self

        if (allocated(self%typepos_locposorder)) then
            deallocate(self%locpos_ownorder)
            deallocate(self%locpos_typeorder)
            deallocate(self%typepos_locposorder)
        endif

        self%current_sector = SecZ
        self%current_component = 1
    end subroutine

    function get_current_worm_sector(self) result(res)
        type(TWormState), intent(in) :: self
        integer :: res

        ! FIXME: enable eventually
        ! if (self%current_sector /= SecZ .and. self%current_sector /= SecCustom) &
        !      error stop 'get_current_worm_sector: Invalid sector'
        res = self%current_sector
    end function

    function get_current_worm_comp(self) result(res)
        type(TWormState), intent(in) :: self
        integer :: res

        res = self%current_component
    end function

    function get_worm_oper_positions_time_ordered(self) result(res)
        type(TWormState), intent(in) :: self
        integer, allocatable :: res(:)

        if (self%current_sector == SecDummy) &
          error stop 'Invalid worm state'
        if (self%current_sector == SecZ) then
        allocate(res(0))
        else
        allocate(res, source=self%locpos_ownorder)
        endif
    end

    function get_worm_oper_position(self, i) result(res)
        type(TWormState), intent(in) :: self
        integer, intent(in) :: i
        integer :: res

        if (self%current_sector == SecDummy) &
          error stop 'Invalid worm state'
        if (self%current_sector == SecZ) &
          error stop 'Not in a worm sector'
        if (i < 1 .or. i > size(self%locpos_typeorder)) &
          error stop 'Invalid position'
        res = self%locpos_typeorder(i)
    end

    ! get the (local operator list) index of an operator at (worm time)
    ! tau index ntau
    pure integer function get_tau_oper_position(self, ntau) result(locind)
        type(TWormState), intent(in) :: self
        integer, intent(in) :: ntau
        integer :: offsets(4), typeind

        select case (self%current_sector)
        case (SecGSigma, SecQQ, SecQUDdag)
            offsets(1:2) = [1, 4]
            typeind = offsets(ntau)
        case (SecH4)
            offsets(1:4) = [1, 4, 5, 6]
            typeind = offsets(ntau)
        case (SecP2, SecP2pp, SecUccaa, SecUcaca)
            offsets(1:2) = [1, 3]
            typeind = offsets(ntau)
        case (SecNQQdag, SecQQdd)
            offsets(1:3) = [1, 4, 7]
            typeind = offsets(ntau)
        case (SecQ4)
            offsets(1:4) = [1, 4, 7, 10]
            typeind = offsets(ntau)
        case (SecRaman)
            offsets(1:3) = [1, 3, 5]
            typeind = offsets(ntau)
        case (SecCustom)
            typeind = self%custom_config(self%current_sector)%offsets(ntau)
        case default
            error stop 'get_tau_oper_position: invalid sector'
        end select

        locind = self%locpos_typeorder(typeind)
    end function get_tau_oper_position

    pure integer function get_nworm(worm_state) result(nworm)
        type(TWormState), intent(in) :: worm_state

        ! TODO: figure out whether size(...) is also a good measure in
        !       the Z sectors
        if (worm_state%current_sector == SecZ) then
                nworm = 0
        else
                nworm = size(worm_state%locpos_typeorder)
        endif
    end

    pure integer function get_ntaus(worm_state) result(ntaus)
        type(TWormState), intent(in) :: worm_state

        ntaus = get_sector_ntaus(worm_state, worm_state%current_sector)
    end function get_ntaus

    subroutine natural_order_perm(ops, perm)
        type(TLocalOper), intent(in) :: ops(:)
        integer, intent(out) :: perm(size(ops))

        integer :: i, nannh, ncrea, nworm

        ! First pass: figure out worm operators
        nworm = 0
        do i = 1, size(ops)
                if (ops(i)%type /= OpCrea .and. ops(i)%type /= OpAnnh) then
                nworm = nworm + 1
                perm(i) = nworm
                endif
        enddo

        ! Second pass: figure out non-worm operators
        ncrea = 0
        nannh = 0
        do i = 1, size(ops)
                if (ops(i)%type == OpCrea) then
                ncrea = ncrea + 1
                perm(i) = nworm + 2 * ncrea - 1
                elseif (ops(i)%type == OpAnnh) then
                nannh = nannh + 1
                perm(i) = nworm + 2 * nannh
                endif
        enddo

        if (nannh /= ncrea) &
                error stop 'Creators and annihilators are not balanced'
    end

    subroutine natural_order_worm_perm(ops, worm_state, perm)
        type(TLocalOper), intent(in) :: ops(:)
        type(TWormState), intent(in) :: worm_state
        integer, intent(out) :: perm(size(ops))

        integer :: i, nworm

        call natural_order_perm(ops, perm)
        nworm = get_nworm(worm_state)
        do i = 1, nworm
                perm(worm_state%locpos_typeorder(i)) = i
        enddo
    end

    subroutine verify_worm(self, ops)
        type(TWormState), intent(in) :: self
        type(TLocalOper), intent(in) :: ops(:)

        integer :: i, iworm, nworm
        real(c_double) :: tau

        if (self%current_sector == SecDummy) &
            error stop 'verify_worm: Dummy sector'
        if (self%current_sector > size(NOperWorm) &
           .or. self%current_sector < SecDummy) &
            error stop 'verify_worm: Nonexistent sector'

        nworm = get_sector_nworm(self, self%current_sector)

        if (self%current_sector /= SecZ) then
            if (.not. allocated(self%locpos_ownorder)) &
                error stop 'verify_worm: No ownorder'
            if (.not. allocated(self%locpos_typeorder)) &
                error stop 'verify_worm: No typeorder'
            if (.not. allocated(self%typepos_locposorder)) &
                error stop 'verify_worm: No permutation'

            if (size(self%locpos_ownorder) /= nworm) &
                error stop 'verify_worm: Invalid size ownorder'
            if (size(self%locpos_typeorder) /= nworm) &
                error stop 'verify_worm: Invalid size typeorder'
            if (size(self%typepos_locposorder) /= nworm) &
                error stop 'verify_worm: Invalid size locposorder'
        else
            if (allocated(self%locpos_ownorder)) &
                error stop 'verify_worm: ownorder in sector Z'
            if (allocated(self%locpos_typeorder)) &
                error stop 'verify_worm: typeorder in sector Z'
            if (allocated(self%typepos_locposorder)) &
                error stop 'verify_worm: permutation in sector Z'
        endif

        do i = 1, nworm
            if (self%locpos_typeorder(i) < 1 .or. self%locpos_typeorder(i) > size(ops)) &
                error stop 'verify_worm: Invalid worm loc'
            if (self%typepos_locposorder(i) < 1 .or. self%typepos_locposorder(i) > nworm) &
                error stop 'verify_worm: Invalid worm permutation'
            if (self%locpos_ownorder(i) &
                        /= self%locpos_typeorder(self%typepos_locposorder(i))) &
                error stop 'verify_worm: Mismatch'
            if (ops(self%locpos_typeorder(i))%type == OpAnnh &
              .or. ops(self%locpos_typeorder(i))%type == OpCrea) &
                error stop 'verify_worm: non-worm operator in worm array'
        enddo
        if (nworm /= 0) then
            if (.not. issorted(self%locpos_ownorder)) &
                error stop 'verify_worm: locpos ownorder must be sorted'
        endif

        iworm = 0
        do i = 1, size(ops)
            if (ops(i)%type == OpCreaW .or. ops(i)%type == OpAnnhW) then
                iworm = iworm + 1
                if (iworm > nworm) &
                    error stop 'verify_worm: Too many worms'
                if (self%locpos_ownorder(iworm) /= i) &
                    error stop 'verify_worm: Worm location mismatch'
            endif
        enddo
        if (iworm /= nworm) &
            error stop 'verify_worm: Not enough worms'

        ! check products at equal time (for adjacent positions and equal times)
        select case (self%current_sector)
        ! FIXME: check this after implementing the missing sectors
        case (SecGSigma, SecH4, SecQUDdag)
            if (self%locpos_typeorder(1) + 1 /= self%locpos_typeorder(2) &
              .or. self%locpos_typeorder(2) + 1 /= self%locpos_typeorder(3)) &
              error stop 'verify_worm: equal time product not contiguous (1)'

            tau = ops(self%locpos_typeorder(1))%tau
            if (ops(self%locpos_typeorder(2))%tau /= tau &
              .or. ops(self%locpos_typeorder(3))%tau /= tau) &
              error stop 'verify_worm: operators not at equal time (1)'
        case (SecP2, SecP2pp, SecUccaa, SecUcaca, SecRaman)
            do i = 0, nworm / 2 - 1
                if (self%locpos_typeorder(2 * i + 1) + 1 &
                 /= self%locpos_typeorder(2 * i + 2)) &
                 error stop 'verify_worm: equal time product not contiguous (2)'

                tau = ops(self%locpos_typeorder(2 * i + 1))%tau
                if (ops(self%locpos_typeorder(2 * i + 2))%tau /= tau) &
                 error stop 'verify_worm: operators not at equal time (2)'
            end do
        case (SecQQ, SecQ4)
            do i = 0, nworm / 3 - 1
                if (self%locpos_typeorder(3 * i + 1) + 1 &
                 /= self%locpos_typeorder(3 * i + 2) &
                 .or. self%locpos_typeorder(3 * i + 2) + 1 /= &
                 self%locpos_typeorder(3 * i + 3)) &
                 error stop 'verify_worm: equal time product not contiguous (4)'

                tau = ops(self%locpos_typeorder(3 * i + 1))%tau
                if (ops(self%locpos_typeorder(3 * i + 2))%tau /= tau &
                 .or. ops(self%locpos_typeorder(3 * i + 3))%tau /= tau) &
                 error stop 'verify_worm: operators not at equal time (4)'
            end do
        case (SecNQQdag, SecQQdd)
            do i = 0, 1
                if (self%locpos_typeorder(3 * i + 1) + 1 &
                 /= self%locpos_typeorder(3 * i + 2) &
                 .or. self%locpos_typeorder(3 * i + 2) + 1 /= &
                 self%locpos_typeorder(3 * i + 3)) &
                 error stop 'verify_worm: equal time product not contiguous (5/1)'

                tau = ops(self%locpos_typeorder(3 * i + 1))%tau
                if (ops(self%locpos_typeorder(3 * i + 2))%tau /= tau &
                 .or. ops(self%locpos_typeorder(3 * i + 3))%tau /= tau) &
                 error stop 'verify_worm: operators not at equal time (5/1)'
            end do

            if (self%locpos_typeorder(7) + 1 /= self%locpos_typeorder(8)) &
              error stop 'verify_worm: equal time product not contiguous (5/2)'

            tau = ops(self%locpos_typeorder(7))%tau
            if (ops(self%locpos_typeorder(8))%tau /= tau) &
              error stop 'verify_worm: operators not at equal time (5/2)'
        case (SecCustom)

            i = 1

            do iworm = 1, get_custom_worm_ntaus(self, self%current_sector)

                tau = ops(self%locpos_typeorder(i))%tau

                do i = i + 1, i - 1 + get_custom_worm_ntauops(self, &
                                                          self%current_sector, &
                                                          iworm)
                    if (self%locpos_typeorder(i - 1) + 1 /= self%locpos_typeorder(i)) &
                    error stop 'verify_worm: equal time product not contiguous (6)'
                    if (ops(self%locpos_typeorder(i))%tau /= tau) &
                    error stop 'verify_worm: operators not at equal time (6)'
                end do

            end do

        end select
    end

    ! Assigns all position arrays in wormst using the explicitly given
    ! array of local indices in 'type order' (sector-dependent) and
    ! mergesort
    subroutine assign_pos_from_locpos_typeorder(wormst, locpos_typeorder)
        type(TWormState), intent(inout)     :: wormst
        integer, allocatable, intent(inout) :: locpos_typeorder(:)  ! allocation is MOVED OUT

        if (allocated(wormst%locpos_typeorder)) deallocate(wormst%locpos_typeorder)
        if (allocated(wormst%locpos_ownorder)) deallocate(wormst%locpos_ownorder)
        if (allocated(wormst%typepos_locposorder)) deallocate(wormst%typepos_locposorder)

        call move_alloc(locpos_typeorder, wormst%locpos_typeorder)
        allocate(wormst%locpos_ownorder(size(wormst%locpos_typeorder)))
        allocate(wormst%typepos_locposorder(size(wormst%locpos_typeorder)))

        ! assign local positions in type order as initial value
        ! to local positions in local position order, then sort the
        ! array and perform the same operations on another array
        ! containing integers in ascending order to get the positions of
        ! the worm operator 'types' in that array
        wormst%locpos_ownorder(:) = wormst%locpos_typeorder(:)
        call sort_perm(wormst%locpos_ownorder, wormst%typepos_locposorder)
    end subroutine assign_pos_from_locpos_typeorder

    subroutine update_wormstate_replace(self, index_in_typeorder, new_value)
        type(TWormState), intent(inout) :: self
        integer, intent(in) :: index_in_typeorder, new_value

        integer, allocatable :: worm_new_locpos_typeorder(:)

        ! adjust worm operator position information
        call move_alloc(self%locpos_typeorder, worm_new_locpos_typeorder)
        worm_new_locpos_typeorder(index_in_typeorder) = new_value
        call assign_pos_from_locpos_typeorder(self, worm_new_locpos_typeorder)
    end subroutine

    subroutine update_wormstate_replace_eqtauop(self, &
                                               eqtau_typeidx, eqtau_locidx, &
                                               old_wormidx, old_repidx)
        type(TWormState), intent(inout) :: self
        ! type indices and new local indices of the moved operators
        integer, intent(in) :: eqtau_typeidx(:), eqtau_locidx(:)
        ! old indices of the operators selected for replacement
        integer, intent(in) :: old_wormidx, old_repidx

        integer, allocatable :: worm_new_locpos_typeorder(:)
        integer :: i

        ! adjust worm operator position information
        call move_alloc(self%locpos_typeorder, worm_new_locpos_typeorder)

        do i = 1, size(worm_new_locpos_typeorder)
            ! fix other worm indices by taking the move of neqtauops - 1
            ! operators from around wormidx to around repidx into account
            if (old_wormidx < worm_new_locpos_typeorder(i) &
              .and. worm_new_locpos_typeorder(i) < old_repidx) then
                worm_new_locpos_typeorder(i) = &
                 worm_new_locpos_typeorder(i) - size(eqtau_typeidx) + 1
            else if (old_repidx < worm_new_locpos_typeorder(i) &
              .and. worm_new_locpos_typeorder(i) < old_wormidx) then
                worm_new_locpos_typeorder(i) = &
                 worm_new_locpos_typeorder(i) + size(eqtau_typeidx) - 1
            end if
        end do

        do i = 1, size(eqtau_typeidx)
            ! fix moved operators indices which are passed as arguments
            worm_new_locpos_typeorder(eqtau_typeidx(i)) = eqtau_locidx(i)
        end do

        call assign_pos_from_locpos_typeorder(self, worm_new_locpos_typeorder)
    end subroutine update_wormstate_replace_eqtauop

    ! Convenience routine for entering a different worm sector
    subroutine wormstate_enter_wormsector(wormst, sector, component, locpos_typeorder)
        type(TWormState), intent(inout)     :: wormst
        integer(kind(SecDummy)), intent(in) :: sector
        integer, intent(in)                 :: component
        integer, allocatable, intent(inout) :: locpos_typeorder(:)  ! allocation is MOVED OUT

        wormst%current_sector = sector
        wormst%current_component = component
        call assign_pos_from_locpos_typeorder(wormst, locpos_typeorder)
    end subroutine wormstate_enter_wormsector

    ! Convenience routine for leaving a worm sector for Z
    pure subroutine wormstate_leave_wormsector(wormst)
        type(TWormState), intent(inout)     :: wormst

        wormst%current_sector = SecZ
        wormst%current_component = 1
        if (allocated(wormst%locpos_typeorder)) deallocate(wormst%locpos_typeorder)
        if (allocated(wormst%locpos_ownorder)) deallocate(wormst%locpos_ownorder)
        if (allocated(wormst%typepos_locposorder)) deallocate(wormst%typepos_locposorder)
    end subroutine wormstate_leave_wormsector

    ! Update worm operator positions stored in a worm state with given
    ! index differences, e.g. due to hyb operator insertion/removal.
    subroutine wormstate_update_pos(wormst, posdiff)
        type(TWormState), intent(inout) :: wormst
        integer, intent(in)             :: posdiff(1:size(wormst%locpos_ownorder))

        if (.not. (allocated(wormst%locpos_typeorder)&
                 .and. allocated(wormst%locpos_typeorder)&
                 .and. allocated(wormst%locpos_typeorder))) &
          error stop 'Worm state is not initialized'

        ! Update the current index in both arrays. For the array in type
        ! order, the correct update is simply the one that keeps the
        ! invariant noted in the type definition valid.
        wormst%locpos_ownorder(:) = wormst%locpos_ownorder(:) + posdiff(:)
        wormst%locpos_typeorder(wormst%typepos_locposorder(:))&
           = wormst%locpos_typeorder(wormst%typepos_locposorder(:)) + posdiff(:)
    end subroutine wormstate_update_pos

    ! Get a difference array 'diff' to update an ascending array of
    ! indices 'idxarray' into the local operator array to apply to the
    ! operator list after acceptance of the move 'move'.
    pure subroutine updatediff_local_index_array(new_idx, idxarray, diff)
        integer, intent(in)              :: new_idx(:)
        integer, intent(in)              :: idxarray(:)
        integer, intent(out)             :: diff(:)

        integer :: inspos(size(new_idx))
        integer :: i, j

        if (size(idxarray) /= size(diff)) &
            error stop 'idxarry and diff must be of same size'

        ! Compute insertion positions from new indices by shifting the new
        ! positions by the number of operators inserted before the current
        ! one
        do i = 1, size(new_idx)
            inspos(i) = new_idx(i)
            do j = 1, size(new_idx)
                if (new_idx(j) < new_idx(i)) &
                    inspos(i) = inspos(i) - 1
            enddo
        enddo

        ! Compute the shifts by counting number of insertions before the
        ! current worm position
        do i = 1, size(idxarray)
            diff(i) = count(inspos <= idxarray(i))
        enddo
    end subroutine updatediff_local_index_array

    subroutine update_wormstate_add_opers(self, new_idx)
     type(TWormState), intent(inout) :: self
     integer, intent(in) :: new_idx(:)

     integer, allocatable :: wormposdiff(:)

     if (self%current_sector == SecDummy) &
        error stop 'Invalid worm state'
     if (self%current_sector == SecZ) &
        return

     allocate(wormposdiff(size(self%locpos_ownorder)))
     call updatediff_local_index_array(&
                new_idx, self%locpos_ownorder, wormposdiff)
     call wormstate_update_pos(self, wormposdiff)
    end subroutine

    ! Update worm operator positions stored in a worm state on removal
    ! of hybridization operators
    pure subroutine wormstate_remove_hybops(wormst, rem_idx)
        type(TWormState), intent(inout) :: wormst
        integer, intent(in) :: rem_idx(:)
        integer :: i, j, subpos

        if (.not. (allocated(wormst%locpos_typeorder)&
                 .and. allocated(wormst%locpos_typeorder)&
                 .and. allocated(wormst%locpos_typeorder))) return


        ! Go through the worm local index list (indexed by j) and update
        ! entries due to the removal of operators at old indices given
        ! in rem_idx (indexed by i). subpos keeps track of the shift to
        ! be subtracted from the current (at j) worm operator local
        ! index due to removal of operators before it.
        do j = 1, size(wormst%locpos_ownorder)
            subpos = 0
            do i = 1, size(rem_idx)
                ! While the next removed index is not greater than the
                ! current worm local index, increment the shift and
                ! proceed through the list of removed indices.
                if (rem_idx(i) < wormst%locpos_ownorder(j)) &
                subpos = subpos + 1
            end do

            ! Update the current index in both arrays. For the array in
            ! type order, the correct update is simply the one that keeps
            ! the invariant noted in the type definition valid.
            wormst%locpos_ownorder(j) = wormst%locpos_ownorder(j) - subpos
            wormst%locpos_typeorder(wormst%typepos_locposorder(j))&
                = wormst%locpos_typeorder(wormst%typepos_locposorder(j)) - subpos
        end do
    end subroutine wormstate_remove_hybops

    ! Update worm operator positions stored in a worm state on tau
    ! shift
    pure subroutine wormstate_taushift(wormst, localops_total, localops_wrapped)
        type(TWormState), intent(inout) :: wormst
        integer, intent(in) :: localops_total, localops_wrapped
        integer :: i, j, last_wrapped
        integer, allocatable :: locpos_ownorder_temp(:), typepos_locposorder_temp(:)

        if (.not. (allocated(wormst%locpos_typeorder)&
                 .and. allocated(wormst%locpos_typeorder)&
                 .and. allocated(wormst%locpos_typeorder))) return

        ! Go through the worm local index list (indexed by j) and update
        ! entries due to a tau shift that wrapped localops_wrapped local
        ! operators around 0=beta out of a total of localops_total.
        allocate(locpos_ownorder_temp(size(wormst%locpos_ownorder)))
        allocate(typepos_locposorder_temp(size(wormst%typepos_locposorder)))
        last_wrapped = 0
        do i = 1, size(wormst%locpos_ownorder)
            ! Write updated local indices into a temporary array and keep
            ! track of how many entries need to be wrapped around 0
            ! (using last_wrapped).
            if (wormst%locpos_ownorder(i) <= localops_wrapped) then
                locpos_ownorder_temp(i) = wormst%locpos_ownorder(i) + localops_total - localops_wrapped
                last_wrapped = i
            else
                locpos_ownorder_temp(i) = wormst%locpos_ownorder(i) - localops_wrapped
            end if
            ! Just copy, as these are not local indices, but in local
            ! index order and so need to be wrapped around 0 as well.
            typepos_locposorder_temp(i) = wormst%typepos_locposorder(i)
        end do

        ! Copy unwrapped part back from the temporary arrays to the
        ! beginning of the state arrays.
        j = 1
        if (last_wrapped < size(wormst%locpos_ownorder)) then
            do i = last_wrapped + 1, size(wormst%locpos_ownorder)
                wormst%locpos_ownorder(j) = locpos_ownorder_temp(i)
                wormst%typepos_locposorder(j) = typepos_locposorder_temp(i)
                j = j + 1
            end do
        end if

        ! Copy back wrapped part after that.
        if (last_wrapped > 0) then
            do i = 1, last_wrapped
                wormst%locpos_ownorder(j) = locpos_ownorder_temp(i)
                wormst%typepos_locposorder(j) = typepos_locposorder_temp(i)
                j = j + 1
            end do
        end if

        ! Update entries in type ordered list using invariant noted in
        ! type definition.
        do i = 1, size(wormst%locpos_ownorder)
            wormst%locpos_typeorder(wormst%typepos_locposorder(i)) = wormst%locpos_ownorder(i)
        end do
    end subroutine wormstate_taushift

    integer function get_custom_worm_ntaus(ws, sector) result(ntaus)
        type(TWormState), intent(in) :: ws
        integer, intent(in) :: sector

        if (.not. allocated(ws%custom_config(sector)%ntauops)) &
           error stop "ntaus: No custom worm config for sector"

        ntaus = ws%custom_config(sector)%ntaus
    end function get_custom_worm_ntaus

    integer function get_custom_worm_ntotalops(ws, sector) result(ntotalops)
        type(TWormState), intent(in) :: ws
        integer, intent(in) :: sector

        if (.not. allocated(ws%custom_config(sector)%ntauops)) &
           error stop "ntotalops: No custom worm config for sector"

        ntotalops = ws%custom_config(sector)%ntotalops
    end function get_custom_worm_ntotalops

    integer function get_custom_worm_ntauops(ws, sector, ntau) result(ntauops)
        type(TWormState), intent(in) :: ws
        integer, intent(in) :: sector, ntau

        if (.not. allocated(ws%custom_config(sector)%ntauops)) &
           error stop "ntauops: No custom worm config for sector"

        ntauops = ws%custom_config(sector)%ntauops(ntau)
    end function get_custom_worm_ntauops

    integer(kind(OpDummy)) function get_custom_worm_optype(ws, sector, nop) &
        result(type)
        type(TWormState), intent(in) :: ws
        integer, intent(in) :: sector, nop

        if (.not. allocated(ws%custom_config(sector)%optypes)) &
           error stop "optype: No custom worm config for sector"

        type = ws%custom_config(sector)%optypes(nop)
    end function get_custom_worm_optype

    integer function get_sector_nworm(ws, sector) result(nworm)
        type(TWormState), intent(in) :: ws
        integer, intent(in) :: sector

        select case (sector)
        case (SecCustom)

            if (.not. allocated(ws%custom_config(sector)%ntauops)) &
              error stop "No custom worm config for sector"

            nworm = ws%custom_config(sector)%ntotalops

        case default

            nworm = NOperWorm(sector)

        end select

    end function get_sector_nworm

    pure integer function get_sector_ntaus(ws, sector) result(ntaus)
        type(TWormState), intent(in) :: ws
        integer, intent(in) :: sector

        select case (sector)
        case (SecP2, SecP2pp, SecGSigma, SecQQ, &
                SecUcaca, SecUccaa, SecQUDdag)
            ntaus = 2
        case (SecNQQdag, SecQQdd, SecRaman)
            ntaus = 3
        case (SecH4, SecQ4)
            ntaus = 4
        case (SecCustom)

            if (.not. allocated(ws%custom_config(sector)%ntauops)) &
              error stop "get_sector_ntaus: No custom worm config for sector"

            ntaus = ws%custom_config(sector)%ntaus

        case default
            error stop "get_sector_ntaus: Sector not implemented"
        end select

    end function get_sector_ntaus

    subroutine set_custom_worm(c_ws, ntaus, c_ntauops, c_optypes) &
        bind (C, name='set_custom_worm')
     use iso_c_binding
     integer(c_int), value :: ntaus
     type(c_ptr), value :: c_ws, c_ntauops, c_optypes

     type(TWormState), pointer :: ws
     integer(c_int), pointer :: ntauops(:), optypes(:)
     integer :: ntotalops, i


     call c_f_pointer(c_ws, ws)
     call c_f_pointer(c_ntauops, ntauops, [ntaus])

     ws%custom_config(SecCustom)%ntaus = ntaus

     if (allocated(ws%custom_config(SecCustom)%ntauops)) &
          deallocate(ws%custom_config(SecCustom)%ntauops)
     allocate(ws%custom_config(SecCustom)%ntauops(ntaus))

     ws%custom_config(SecCustom)%ntauops(:) = ntauops(:)

     if (allocated(ws%custom_config(SecCustom)%offsets)) &
          deallocate(ws%custom_config(SecCustom)%offsets)
     allocate(ws%custom_config(SecCustom)%offsets(ntaus))

     ntotalops = 1
     do i = 1, ntaus
        ws%custom_config(SecCustom)%offsets(i) = ntotalops
        ntotalops = ntotalops + ntauops(i)
     end do
     ntotalops = ntotalops - 1

     ws%custom_config(SecCustom)%ntotalops = ntotalops

     call c_f_pointer(c_optypes, optypes, [ntotalops])

     if (allocated(ws%custom_config(SecCustom)%optypes)) &
          deallocate(ws%custom_config(SecCustom)%optypes)
     allocate(ws%custom_config(SecCustom)%optypes(ntotalops))

     ws%custom_config(SecCustom)%optypes(:) = optypes(:)
    end subroutine set_custom_worm

end module MWormState
