module MUtilities
    use iso_c_binding, only: c_int64_t
    implicit none
    private

    public :: connected_components
    public :: ilog2

contains

    !> Extract connected components from a graph given as adjacency matrix
    !!
    !! Expects an `n x n` adjacency matrix `adjmat`. This corresponds to a
    !! directed graph with vertices `{1, ..., n}`, where `adjmat(i,j)` is true
    !! if and only if there is an edge leading from vertex `i` to `j`.
    !!
    !! Let us assume the graph has `k` connected components `{1, ..., k}`. The
    !! function returns a vector `whichblock` of size `n`, where
    !! `whichblock(i) == j` if the vertex `i` belongs to the connected
    !! component `j`.
    function connected_components(adjmat) result(whichblock)
        logical, intent(in) :: adjmat(:, :)
        integer, allocatable :: whichblock(:)

        integer :: i, j, bl

        if (size(adjmat, 1) /= size(adjmat, 2)) &
            stop 'Matrix must be square'

        allocate(whichblock(1:size(adjmat, 1)))
        whichblock(:) = 0
        bl = 0
        do j = 1, size(whichblock)
            if (whichblock(j) == 0) then
                bl = bl + 1
                whichblock(j) = bl
            endif
            do i = j+1, size(whichblock)
                if (whichblock(i) == 0 .and. adjmat(i, j)) &
                    whichblock(i) = whichblock(j)
            enddo
            do i = j+1, size(whichblock)
                if (whichblock(i) == 0 .and. adjmat(j, i)) &
                    whichblock(i) = whichblock(j)
            enddo
        enddo
    end

    pure function ilog2(val) result(res)
        integer(c_int64_t), intent(in) :: val
        integer(c_int64_t) :: res, tmp

        res = -1
        if (val < 1) return

        tmp = val
        do while (tmp /= 0)
            res = res + 1
            tmp = tmp / 2
        enddo
    end function

end
