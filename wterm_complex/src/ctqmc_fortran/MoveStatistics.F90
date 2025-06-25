module MMoveStatistics
    use iso_c_binding
    use MLocalBase
    use MWormState
    implicit none
    private

    integer, parameter, public :: MOVE_FAILED_TO_GENERATE = 1
    integer, parameter, public :: MOVE_NOT_ALLOWED = 2
    integer, parameter, public :: MOVE_REJECTED = 3
    integer, parameter, public :: MOVE_ACCEPTED = 4

    type, public :: TMoveStat
        ! trial and acceptance statistics; pretry counts all attempts
        ! due to move type selection, even if no sensible move spec can
        ! then be generated
        integer(c_int64_t) :: pretry(NWormSectors) = 0, try(NWormSectors) = 0, acc(NWormSectors) = 0

        ! cache to avoid passing info multiple times per move
        integer :: cache_sector = -1, cache_orb = -1, cache_sp = -1
        real(c_double) :: cache_taudiff = -huge(1.0d0)

        ! taudiff pair trial and acceptance statistics
        integer :: ntaudiffbins = 0, listindex = 0
        real(c_double) :: taudifflimit = huge(0)
        real(c_double), allocatable :: acc_pair_list(:)
        integer, allocatable :: try_pair_taudiffbin(:, :, :), acc_pair_taudiffbin(:, :, :)
    end type TMoveStat

    public :: record_move
    public :: acceptance_ratio

    interface init_tau_binning_statistics
        module procedure movestat_taudiff_bin_mode
    end interface
    public init_tau_binning_statistics

    interface init_tau_history_statistics
        module procedure movestat_taudiff_list_mode
    end interface
    public init_tau_history_statistics

    interface reset_statistics
        module procedure movestat_reset_statistics
    end interface
    public reset_statistics

contains

    !> Update statistics with the outcome of a move.
    !!
    !! @param self     Move statics
    !! @param sector   Worm sector
    !! @param outcome  Outcome of the move, one of MOVE_*
    !! @param ops      Operators touched by this move (optional)
    subroutine record_move(self, sector, outcome, ops)
        type(TMoveStat), intent(inout) :: self
        integer, intent(in) :: sector, outcome
        type(TLocalOper), intent(in), optional :: ops(:)

        if (sector < 1 .or. sector > NWormSectors) &
                error stop 'Invalid sector'

        call movestat_count_pretry(self, sector)
        if (outcome /= MOVE_FAILED_TO_GENERATE) then
                if (present(ops)) then
                if (size(ops) == 2) then
                    call movestat_count_try( &
                                self, abs(ops(2)%tau - ops(1)%tau), &
                                ops(1)%orb, ops(1)%sp)
                else
                    call movestat_count_try(self)
                endif
                else
                call movestat_count_try(self)
                endif
                if (outcome /= MOVE_NOT_ALLOWED) then
                if (outcome == MOVE_REJECTED) then
                    call movestat_count_rej(self)
                else
                    call movestat_count_acc(self)
                endif
                endif
        endif
    end subroutine

    pure subroutine movestat_count_pretry(stat, sector)
        type(TMoveStat), intent(inout) :: stat
        integer, intent(in) :: sector

        stat%pretry(sector) = stat%pretry(sector) + 1
        stat%cache_sector = sector
    end subroutine movestat_count_pretry

    pure subroutine movestat_count_try(stat, taudiff, orb, sp, sector)
        type(TMoveStat), intent(inout)    :: stat
        real(c_double), intent(in), optional :: taudiff
        integer, intent(in), optional     :: orb, sp
        integer, intent(in), optional     :: sector
        integer                           :: bin

        if (present(sector)) stat%cache_sector = sector

        stat%try(stat%cache_sector) = stat%try(stat%cache_sector) + 1

        if (present(taudiff) .and. present(orb) .and. present(sp)) then
            if (stat%listindex == 0 .and. stat%ntaudiffbins >= 1) then  ! tau bin mode
                bin = min(stat%ntaudiffbins, floor(taudiff/stat%taudifflimit * real(stat%ntaudiffbins, c_double)) + 1)
                stat%try_pair_taudiffbin(orb, sp, bin) = stat%try_pair_taudiffbin(orb, sp, bin) + 1
            end if

            stat%cache_taudiff = taudiff
            stat%cache_orb = orb
            stat%cache_sp = sp
        else
            stat%cache_taudiff = -huge(1.0d0)
            stat%cache_orb = -1
            stat%cache_sp = -1
        end if
    end subroutine movestat_count_try

    pure subroutine movestat_count_acc(stat, sector, taudiff, orb, sp)
        type(TMoveStat), intent(inout)    :: stat
        integer, intent(in), optional     :: sector
        real(c_double), intent(in), optional :: taudiff
        integer, intent(in), optional     :: orb, sp
        real(c_double)                       :: td
        integer                           :: o, s, bin

        if (present(sector)) stat%cache_sector = sector

        stat%acc(stat%cache_sector) = stat%acc(stat%cache_sector) + 1

        if (present(taudiff) .and. present(orb) .and. present(sp)) then
            td = taudiff
            o = orb
            s = sp
        else
            td = stat%cache_taudiff
            o = stat%cache_orb
            s = stat%cache_sp
        end if

        if (td /= -huge(1.0d0)) then
            if (stat%listindex >= 1) then  ! list mode
                if (stat%listindex <= size(stat%acc_pair_list)) then
                    stat%acc_pair_list(stat%listindex) = td
                    stat%listindex = stat%listindex + 1
                end if
            else  ! tau bin mode
                if (stat%ntaudiffbins >= 1) then
                    bin = min(stat%ntaudiffbins, floor(td/stat%taudifflimit * real(stat%ntaudiffbins, c_double)) + 1)
                    stat%acc_pair_taudiffbin(o, s, bin) = stat%acc_pair_taudiffbin(o, s, bin) + 1
                end if
            end if
        end if

        stat%cache_sector = -1
        stat%cache_taudiff = -huge(1.0d0)
        stat%cache_orb = -1
        stat%cache_sp = -1
    end subroutine movestat_count_acc

    pure subroutine movestat_count_rej(stat)
        type(TMoveStat), intent(inout)    :: stat

        stat%cache_sector = -1
        stat%cache_taudiff = -huge(1.0d0)
        stat%cache_orb = -1
        stat%cache_sp = -1
    end subroutine movestat_count_rej

    ! Acceptance ratio (number of acceptance / number of pretries) for
    ! a statistics collector stat, either for one specific sector (Z /
    ! worm) if sector is a valid sector index or summed over all
    ! sectors otherwise. Returns zero in case of zero pretries.
    pure real(c_double) function acceptance_ratio(stat, sector) result(res)
     class(TMoveStat), intent(in)    :: stat
     integer, intent(in)             :: sector
     integer(c_int64_t)              :: pretry, acc

     res = 0
     if (sector >= 1 .and. sector <= NWormSectors) then
        pretry = stat%pretry(sector)
        acc = stat%acc(sector)
     else
        pretry = sum(stat%pretry)
        acc = sum(stat%acc)
     end if
     if (pretry > 0) res = real(acc, c_double) / real(pretry, c_double)
    end function acceptance_ratio

    pure subroutine movestat_taudiff_bin_mode(stat, orbs, spins, nbins, limit)
        class(TMoveStat), intent(inout) :: stat
        integer, intent(in) :: orbs, spins, nbins
        real(c_double), intent(in) :: limit

        call tmovestat_scalar_dest(stat)
        stat%listindex = 0

        allocate(stat%try_pair_taudiffbin(orbs, spins, nbins))
        stat%try_pair_taudiffbin = 0
        allocate(stat%acc_pair_taudiffbin(orbs, spins, nbins))
        stat%acc_pair_taudiffbin = 0
        stat%ntaudiffbins = nbins
        stat%taudifflimit = limit
    end subroutine movestat_taudiff_bin_mode

    pure subroutine movestat_taudiff_list_mode(stat, nlist)
        class(TMoveStat), intent(inout) :: stat
        integer, intent(in) :: nlist

        stat%ntaudiffbins = 0
        stat%taudifflimit = huge(stat%taudifflimit)
        call tmovestat_scalar_dest(stat)

        allocate(stat%acc_pair_list(nlist))
        stat%listindex = 1
    end subroutine movestat_taudiff_list_mode

    pure subroutine movestat_reset_statistics(stat)
        class(TMoveStat), intent(inout) :: stat

        call tmovestat_scalar_dest(stat)
        stat%pretry = 0
        stat%try = 0
        stat%acc = 0

        stat%ntaudiffbins = 0
        stat%listindex = 0
        stat%taudifflimit = huge(stat%taudifflimit)
    end subroutine

    pure subroutine tmovestat_scalar_dest(stat)
        type(TMoveStat), intent(inout) :: stat

        if (allocated(stat%try_pair_taudiffbin)) deallocate(stat%try_pair_taudiffbin)
        if (allocated(stat%acc_pair_taudiffbin)) deallocate(stat%acc_pair_taudiffbin)
        if (allocated(stat%acc_pair_list)) deallocate(stat%acc_pair_list)
    end subroutine tmovestat_scalar_dest

end module MMoveStatistics
