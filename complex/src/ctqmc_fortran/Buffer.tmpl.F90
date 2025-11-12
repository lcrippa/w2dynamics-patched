module MBUFFER
    use iso_c_binding, only: c_double, c_double_complex
    implicit none
    private

    interface reserve
        module procedure p_reserve_1, p_reserve_2
#if VALUE_IS_INTEGER
        module procedure p_reserve_l1
#endif
    end interface
    public reserve

    interface copy
        module procedure p_copy_1
    end interface
    public copy

#if VALUE_IS_INTEGER
    VALUE_TYPE, parameter :: default_value = -4711
#else
    REAL_TYPE, parameter :: ONE = 1
    VALUE_TYPE, parameter :: default_value = huge(ONE)
#endif

contains

    subroutine p_reserve_1(buf, rowcap, preserve)
        VALUE_TYPE, allocatable, intent(inout) :: buf(:)
        integer, intent(in) :: rowcap
        logical, intent(in), optional :: preserve
        VALUE_TYPE, allocatable :: tmp(:)

        if (allocated(buf)) then
            if (rowcap <= size(buf, 1)) &
                return
            if (value_or_false(preserve)) then
                call move_alloc(buf, tmp)
            else
                deallocate(buf)
            endif
        endif
        allocate(buf(rowcap), source=default_value)
        if (allocated(tmp)) then
            buf(:size(tmp, 1)) = tmp(:)
        endif
    end subroutine

    subroutine p_reserve_2(buf, rowcap, colcap, preserve)
        VALUE_TYPE, allocatable, intent(inout) :: buf(:,:)
        integer, intent(in) :: rowcap, colcap
        logical, intent(in), optional :: preserve
        VALUE_TYPE, allocatable :: tmp(:,:)

        if (allocated(buf)) then
            if (rowcap <= size(buf, 1) .and. colcap <= size(buf, 2)) &
                return
            if (value_or_false(preserve)) then
                call move_alloc(buf, tmp)
            else
                deallocate(buf)
            endif
        endif
        allocate(buf(rowcap, colcap), source=default_value)
        if (allocated(tmp)) then
            buf(:size(tmp, 1), :size(tmp, 2)) = tmp(:, :)
        endif
    end subroutine

#if VALUE_IS_INTEGER
    subroutine p_reserve_l1(buf, rowcap, preserve)
        logical, allocatable, intent(inout) :: buf(:)
        integer, intent(in) :: rowcap
        logical, intent(in), optional :: preserve
        logical, allocatable :: tmp(:)

        if (allocated(buf)) then
            if (rowcap <= size(buf, 1)) &
                return
            if (value_or_false(preserve)) then
                call move_alloc(buf, tmp)
            else
                deallocate(buf)
            endif
        endif
        allocate(buf(rowcap))
        if (allocated(tmp)) then
            buf(:size(tmp, 1)) = tmp(:)
        endif
    end subroutine
#endif

    subroutine p_copy_1(val, buf)
        VALUE_TYPE, intent(in) :: val(:)
        VALUE_TYPE, allocatable, intent(inout) :: buf(:)

        if (allocated(buf)) then
            if (size(buf) /= size(val)) &
                deallocate(buf)
        endif
        if (allocated(buf)) then
            buf(:) = val(:)
        else
            allocate(buf, source=val)
        endif
    end subroutine

    pure logical function value_or_false(preserve)
        logical, intent(in), optional :: preserve

        if (present(preserve)) then
            value_or_false = preserve
        else
            value_or_false = .false.
        endif
    end function

end module
