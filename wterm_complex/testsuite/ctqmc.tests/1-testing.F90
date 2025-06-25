program test_extrange
    use Testing
    use MPrinting
    implicit none

    call test_close()
    call test_fmt
contains
    subroutine test_close()
        if (isclose(1.0d0, 0.99999d0)) &
            error stop 'Should not be close'
        if (isclose(0.0d0, 1d-200)) &
            error stop 'Should not be close'
        if (.not. isclose(1.0d0, 1.0d0 + epsilon(0.0d0))) &
            error stop 'Should be close'
    end subroutine

    subroutine test_fmt()
        write (*,*) format_matrix(3, 2, 'ES11.3')
        write (*,*) format_matrix(0, 1, 'ES11.3')

        write (*, format_matrix(3, 2, 'I3')) 'X', (/ 1, 2, 3, 4, 5, 6 /)

        call print_array((/ 1.0d0, 2.0d0 /))
        call print_array(reshape((/ 1, 2, 3, 4, 5, 6/), (/ 2, 3 /)))
        call print_array(reshape((/ 1, 2, 3, 4, 5, 6/), (/ 6, 1 /)))
    end subroutine
end program
