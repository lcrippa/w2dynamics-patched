module MPsi
    use MOperator
    use MPsiD
    use MPsiZ
    use MStates
    implicit none

    public :: rho1_from_density_matrix
    public :: double_occ_from_density_matrix

    type, public :: TPsis
        !private
        type(DTPsis), allocatable :: dpsis
        type(ZTPsis), allocatable :: zpsis

        ! This duplicates the structure of dpsis/zpsis with toperator
        ! XXX this is fairly hacky
        type(TOperator), allocatable :: psis(:,:,:)
    end type

    interface block_diagonal_oper
        module procedure p_block_diag
    end interface
    public block_diagonal_oper

    interface make_hamiltonian
        module procedure p_make_hamiltonian
    end interface
    public make_hamiltonian

    interface transform_Psis
       module procedure p_transform_Psis
    end interface
    public transform_Psis

    interface init_psi
        module procedure p_init_psi
    end interface
    public init_psi

    public is_complex_psis
    public dest_Psis

contains

    subroutine p_block_diag(this, DStates, iscomplex)
        type(TOperator), intent(out) :: this
        type(TStates), intent(in)    :: DStates
        logical, intent(in)          :: iscomplex

        if (iscomplex) then
            allocate(this%zop)
            call block_diagonal_oper(this%zop, DStates)
        else
            allocate(this%dop)
            call block_diagonal_oper(this%dop, DStates)
        endif
    end subroutine

    logical function is_complex_psis(this)
        type(TPsis), intent(in) :: this

        if (allocated(this%dpsis) .eqv. allocated(this%zpsis)) &
            error stop 'Psis not allocated'
        is_complex_psis = allocated(this%zpsis)
    end


    subroutine p_init_Psi(this,DStates,iscomplex)
        type(TPsis), intent(out)  :: this
        type(TStates), intent(in) :: DStates
        logical, intent(in)       :: iscomplex

        if (iscomplex) then
            allocate(this%zpsis)
            call init_Psi(this%zpsis, DStates)
        else
            allocate(this%dpsis)
            call init_Psi(this%dpsis, DStates)
        endif

        call copy_inner_psis(this)
    end subroutine

    subroutine copy_inner_psis(this)
        type(TPsis), intent(inout) :: this
        integer :: norb, iorb, isp, ca

        if (is_complex_psis(this)) then
            norb = size(this%zpsis%psis, 1)
        else
            norb = size(this%dpsis%psis, 1)
        endif
        allocate(this%psis(norb, 2, 2))

        do ca = 1, 2
            do isp = 1, 2
                do iorb = 1, norb
                    if (is_complex_psis(this)) then
                        allocate(this%psis(iorb, isp, ca)%zop, &
                                 source=this%zpsis%psis(iorb, isp, ca))
                    else
                        allocate(this%psis(iorb, isp, ca)%dop, &
                                 source=this%dpsis%psis(iorb, isp, ca))
                    endif
                enddo
            enddo
        enddo
    end subroutine

    subroutine p_make_hamiltonian(muimp, u_matrix, states, psi, hloc, qnthr)
        complex(c_double_complex), intent(in) :: muimp(:, :, :, :)
        complex(c_double_complex), intent(in) :: u_matrix(:, :, :, :)
        real(c_double), intent(in) :: qnthr
        type(TStates), intent(inout) :: states
        type(TPsis), intent(inout) :: psi
        type(TOperator), intent(out) :: hloc

        if (is_complex_psis(psi)) then
            allocate(hloc%zop)
            call make_hamiltonian( &
                        muimp, u_matrix, states, psi%zpsis, hloc%zop, qnthr)
        else
            if (any(imag(muimp) /= 0)) &
                error stop 'muimp must be real'
            if (any(imag(u_matrix) /= 0)) &
                error stop 'U matrix must be real'

            allocate(hloc%dop)
            call make_hamiltonian( &
                        real(muimp), real(u_matrix), states, &
                        psi%dpsis, hloc%dop, qnthr)
        endif
    end subroutine

    subroutine p_transform_Psis(Eigen, transformed_Psis)
        type(TOperEigen), intent(in) :: Eigen
        type(TPsis), intent(inout)   :: transformed_Psis

        if (is_complex_opeigen(Eigen)) then
            if (.not. is_complex_psis(transformed_Psis)) &
                error stop 'Invalid psis'
            call transform_Psis(Eigen%zeig, transformed_Psis%zpsis)
        else
            if (is_complex_psis(transformed_Psis)) &
                error stop 'Invalid psis'
            call transform_Psis(Eigen%deig, transformed_Psis%dpsis)
        endif
    end

    ! XXX remove these routines
    subroutine dest_Psis(this)
        type(TPsis)                               :: this

        if (allocated(this%dpsis)) deallocate(this%dpsis)
        if (allocated(this%zpsis)) deallocate(this%zpsis)
    end subroutine

end module
