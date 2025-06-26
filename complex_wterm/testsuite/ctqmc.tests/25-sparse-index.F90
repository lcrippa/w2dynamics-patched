#define INFO line=__LINE__

program test_sparse_index
    use MUtilities
    use Testing
    implicit none

    call test_conn_comp()
contains

    subroutine test_conn_comp()
        integer :: x(4, 4) = reshape((/ &
            1, 0, 0, 0, &
            0, 1, 0, 0, &
            1, 0, 0, 1, &
            0, 0, 1, 1  &
            /), (/ 4, 4 /))

        call assert_equal(connected_components(x /= 0), &
                          (/ 1, 2, 1, 1/), INFO)
    end
end program
