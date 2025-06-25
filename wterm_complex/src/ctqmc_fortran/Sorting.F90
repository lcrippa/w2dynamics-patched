module MSorting
    ! Type-specific sortings
    use MSortingI
    use MSortingD

    !>      issorted(arr, isless)
    !!
    !! Return if array of values are sorted from smallest to largest. NaNs,
    !! if present, are ignored.
    interface issorted
        module procedure p_issorted
    end interface
    public issorted

contains

    pure logical function p_issorted(x, isless) result(answer)
        character, intent(in) :: x(:, :)
        interface
            pure logical function isless(a, b)
                character, intent(in) :: a(*), b(*)
            end function
        end interface
        integer :: i

        do i = 1, size(x) - 1
            if (isless(x(:, i+1), x(:, i))) then
                answer = .false.
                return
            endif
        enddo
        answer = .true.
    end

end module
