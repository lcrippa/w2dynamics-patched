
module MOperator
    use iso_c_binding, only: c_double, c_double_complex
    use DMOperator
    use ZMOperator
    use MComplexNumbers
    implicit none
    !private

    type, public :: TOperator
        !private
        type(DTOperator), allocatable :: dop
        type(ZTOperator), allocatable :: zop
    end type

    type, public :: TOperEigen
        !private
        type(DOperEigen), allocatable :: deig
        type(ZOperEigen), allocatable :: zeig
    end type

    interface is_initialized
        module procedure p_is_initialized
    end interface
    public is_initialized

    interface init_operator_like
        module procedure p_init_operator_like
    end interface
    public init_operator_like

    interface force_hermitian
        module procedure p_force_hermitian
    end interface
    public force_hermitian

    interface diag_operator
        module procedure p_diag_operator
    end interface
    public diag_operator

    interface connecting_sst
       module procedure p_connecting_sst
    end interface
    public connecting_sst

    interface connecting_sst_backwards
       module procedure p_connecting_sst_backwards
    end interface
    public connecting_sst_backwards

    interface get_eigenvalues
       module procedure p_get_eigenvalues
    end interface
    public get_eigenvalues

    interface get_num_superstates
       module procedure p_get_num_superstates, p_get_num_superstates_eig
    end interface
    public get_num_superstates

    interface get_egs
       module procedure p_get_egs
    end interface
    public get_egs

    interface thermal_weights
       module procedure t_thermal_weights
    end interface
    public thermal_weights

    interface max_sst_size
        module procedure p_max_sst_size
    end interface
    public max_sst_size

    interface dest_toperator
        module procedure p_dest_toperator
    end interface
    public dest_toperator

    public is_complex_opeigen
    public is_complex_op

contains
    logical function is_complex_op(this)
        type(TOperator), intent(in) :: this

        if (allocated(this%dop) .eqv. allocated(this%zop)) &
            error stop 'Invalid state of operator'
        is_complex_op = allocated(this%zop)
    end

    pure logical function p_is_initialized(self)
        type(TOperator), intent(in) :: self

        if (allocated(self%dop) .eqv. allocated(self%zop)) then
            p_is_initialized = .false.
        else if (allocated(self%zop)) then
            p_is_initialized = is_initialized(self%zop)
        else
            p_is_initialized = is_initialized(self%dop)
        endif
    end

    subroutine p_init_operator_like(self, ref)
        type(TOperator), intent(out) :: self
        type(TOperEigen), intent(in) :: ref

        if (is_complex_opeigen(ref)) then
            allocate(self%zop)
            call init_operator_like(self%zop, ref%zeig)
        else
            allocate(self%dop)
            call init_operator_like(self%dop, ref%deig)
        endif
    end subroutine

    subroutine p_force_hermitian(H)
        type(TOperator), intent(inout) :: H

        if (is_complex_op(h)) then
             call force_hermitian(H%zop)
        else
             call force_hermitian(H%dop)
        endif
    end subroutine

    subroutine p_diag_Operator(this,Eigen,epsevec)
        type(TOperator), intent(in)            :: this
        type(TOperEigen), intent(out)          :: Eigen
        real(c_double), optional, intent(in)   :: epsevec

        if (is_complex_op(this)) then
            allocate(Eigen%zeig)
            call diag_Operator(this%zop,Eigen%zeig,epsevec)
        else
            allocate(Eigen%deig)
            call diag_Operator(this%dop,Eigen%deig,epsevec)
        endif
    end subroutine

    logical function is_complex_opeigen(this)
        type(TOperEigen), intent(in) :: this

        if (allocated(this%deig) .eqv. allocated(this%zeig)) &
            error stop 'Invalid state of operator-eigen'
        is_complex_opeigen = allocated(this%zeig)
    end

    function p_connecting_sst(self, sst) result(next_sst)
        type(TOperator), intent(in) :: self
        integer, intent(in) :: sst
        integer :: next_sst

        if (is_complex_op(self)) then
            next_sst = connecting_sst(self%zop, sst)
        else
            next_sst = connecting_sst(self%dop, sst)
        endif
    end function

    function p_connecting_sst_backwards(self, sst) result(next_sst)
        type(TOperator), intent(in) :: self
        integer, intent(in) :: sst
        integer :: next_sst

        if (is_complex_op(self)) then
            next_sst = connecting_sst_backwards(self%zop, sst)
        else
            next_sst = connecting_sst_backwards(self%dop, sst)
        endif
    end function

    function p_get_eigenvalues(self, sst) result(evals)
        type(TOperEigen), intent(in) :: self
        integer, intent(in) :: sst
        real(c_double), allocatable :: evals(:)

        if (is_complex_opeigen(self)) then
            evals = get_eigenvalues(self%zeig, sst)
        else
            evals = get_eigenvalues(self%deig, sst)
        endif
    end function

    function p_get_egs(self) result(egs)
        type(TOperEigen), intent(in) :: self
        real(c_double), allocatable :: egs

        if (is_complex_opeigen(self)) then
            egs = get_egs(self%zeig)
        else
            egs = get_egs(self%deig)
        endif
    end function

    integer function p_get_num_superstates_eig(self) result(n)
        type(TOperEigen), intent(in) :: self

        if (is_complex_opeigen(self)) then
            n = get_num_superstates(self%zeig)
        else
            n = get_num_superstates(self%deig)
        endif
    end function

    integer function p_get_num_superstates(self) result(n)
        type(TOperator), intent(in) :: self

        if (is_complex_op(self)) then
            n = get_num_superstates(self%zop)
        else
            n = get_num_superstates(self%dop)
        endif
    end function

    function p_max_sst_size(self) result(n)
        type(TOperEigen), intent(in) :: self
        integer :: n

        if (is_complex_opeigen(self)) then
            n = max_sst_size(self%zeig)
        else
            n = max_sst_size(self%deig)
        endif
    end function

    function t_thermal_weights(H_eigen, beta) result(weights)
        type(TOperEigen), intent(in) :: H_eigen
        real(c_double), intent(in) :: beta
        real(c_double), allocatable :: weights(:)

        if (is_complex_opeigen(H_eigen)) then
            weights = thermal_weights(H_eigen%zeig, beta)
        else
            weights = thermal_weights(H_eigen%deig, beta)
        endif
    end function

    subroutine p_dest_TOperator(this)
        type(TOperator)                        :: this

        if (allocated(this%dop)) deallocate(this%dop)
        if (allocated(this%zop)) deallocate(this%zop)
    end subroutine

end module