module MBathBase
    use iso_c_binding, only: c_double, c_double_complex
    use MCommon
    implicit none
    private

    interface operator ( == )
        module procedure equalBathOperator
    end interface
    public operator ( == )

    !!> A helper type to store all operators with hybridization lines in one array.
    type, public, extends(TBaseOper) :: TBathOperator
    end type

    ! XXX: these are also somewhat misplaced here
    public read_bathconfig
    public read_bathconfig_highprec

contains

    pure elemental function equalBathOperator(op1, op2) result(eq)
        type(TBathOperator), intent(in)  :: op1, op2
        logical                          :: eq
        eq =  (op1%tau == op2%tau&
            .and. op1%orb == op2%orb&
            .and. op1%sp == op2%sp)
    end function equalBathOperator

    !> reads bath config from file with format as created in subroutine above
    subroutine read_bathconfig(filename, bc, ba, Noper)
        character(len=*), intent(in) :: filename
        type(TBathOperator), allocatable, intent(out) :: bc(:), ba(:)
        integer :: i, j
        integer, intent(out) :: Noper
        character(len=10) :: c1
        character(len=4) :: c2

        !write(*,*) "filename:", filename
        open (66, file=filename, form="formatted", status="old")
        read(66, "(A10, I5, A4)") c1, Noper, c2
        !write(*,*) "Noper:", Noper
        allocate(bc(Noper/2))
        allocate(ba(Noper/2))
        read(66, *)
        read(66, *)
        do i = 1, Noper/2
            read(66, "(I5, F12.7, I2, I2)") j, bc(i)%tau, bc(i)%orb, bc(i)%sp
        enddo
        read(66, *)
        read(66, *)
        do i = 1, Noper/2
            read(66, "(I5, F12.7, I2, I2)") j, ba(i)%tau, ba(i)%orb, ba(i)%sp
        enddo
        close(66)
    end subroutine read_bathconfig

    !> reads bath config from file with format as created in subroutine above
    subroutine read_bathconfig_highprec(filename, bc, ba, Noper, tr, det, sig, permsig)
        character(len=*), intent(in) :: filename
        type(TBathOperator), allocatable, intent(out) :: bc(:), ba(:)
        integer :: i, j, orb, sp
        integer, intent(out) :: Noper
        character(len=10) :: c10
        character(len=4) :: c4
        character(len=3) :: c3
        real(c_double) :: tau
        real(c_double) :: retr, imtr, redet, imdet, resig, imsig, permsig
        complex(c_double_complex) :: tr, det, sig

        !write(*,*) "filename:", filename
        open (66, file=filename, form="formatted", status="old")
        read(66, "(A10, I5, A4)") c10, Noper, c4
        !read(66, "(A9)") c9
        !write(*,*) "Noper:", Noper
        !write(*,*) "noper/2", noper/2
        allocate(bc(Noper/2))
        allocate(ba(Noper/2))
        !write(*,*) "size(ba)", size(ba)
        read(66, *)
        read(66, *)
        do i = 1, Noper/2
            read(66, "(I5, A3, F19.16, I3, I3)") j, c3, tau, orb, sp
            !write(*,*) "tau, orb, sp", tau, orb, sp
            bc(i)%tau = tau
            bc(i)%orb = orb
            bc(i)%sp = sp
        enddo
        read(66, *)
        read(66, *)
        do i = 1, Noper/2
            read(66, "(I5, A3, F19.16, I3, I3)") j, c3, tau, orb, sp
            !write(*,*) "tau, orb, sp", tau, orb, sp
            ba(i)%tau = tau
            ba(i)%orb = orb
            ba(i)%sp = sp
        enddo

        read(66, "(E26.16)") redet
        read(66, "(E26.16)") imdet
        read(66, "(E26.16)") retr
        read(66, "(E26.16)") imtr
        read(66, "(E26.16)") resig
        read(66, "(E26.16)") imsig
        read(66, *)
        read(66, "(E26.16)") permsig

        tr = dcmplx(retr, imtr)
        det = dcmplx(redet, imdet)
        sig = dcmplx(resig, imsig)
        !write(*,*) "det", det
        !write(*,*) "sig", sig

        close(66)
    end subroutine read_bathconfig_highprec

end module