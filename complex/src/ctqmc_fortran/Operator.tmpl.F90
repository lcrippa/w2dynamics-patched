#if 0
! This is a templated Fortran file.
!
! By invoking `make`, several modules are generated in sub-directories
! double/ etc, for each valid value type (real, complex, etc.).
! In each generated module, the value VALUE_TYPE is replaced by the
! corresponding value type, and each uppercase derived type is replaced by
! a typed version.
#endif

!> This module generates the Hamiltonian, creation and annihilation Operators.
!! They can then be stored in the type specified in MStates.
module MOPERATOR
    use iso_c_binding, only: c_double, c_double_complex
    use MBLAS
    use MLAPACK
    use MLibC
    use MComplexNumbers, only: nan
    use MExtRange
    use MPermutation
    use MPrinting
    implicit none
    private

    type :: TSubOperator
        VALUE_TYPE, allocatable :: op(:,:)
    end type TSubOperator

    !> Many-body operator (matrix) in block form
    !!
    !! This stores a block-diagonal (Hamiltonian etc.) or block-permutation
    !! (creators/annihilators) matrix. One can use `apply_oper_left()` to
    !! apply the operator to a state vector.
    type, public :: OPERATOR
        private
        type(TSubOperator), allocatable :: SubOps(:)
        integer, allocatable :: connect(:)
        integer, allocatable :: offset(:)
    end type OPERATOR

    type :: TSubEigen
        ! evec(:,i) is the i-th eigenvector,
        VALUE_TYPE, allocatable :: evec(:,:)
        REAL_TYPE, allocatable :: eval(:)
    end type

    !> Eigenvalues and eigenvectors of a block matrix
    !!
    !! Stores the eigenvalues and eigenvectors for each block of a given
    !! matrix (Hamiltonian).  Once constructed, one can perform time evolution
    !! using `time_evolve()`.
    type, public :: OPER_EIGEN
        private
        type(TSubEigen), allocatable :: SubEigen(:)
        integer, allocatable :: offset(:)
        REAL_TYPE :: Egs
    end type

    !> A type containing the local state vector (ket or bra) in one
    !! superstate as a dynamically allocated array and a common exponent
    !! of all components as well as the index of that superstate and its
    !! size.
    type, public :: STATE_VECTOR
        private
        integer :: sst = -999   ! current superstate
        integer :: sstsize      ! size of current superstate
        integer :: nvec         ! number of vectors
        REAL_TYPE :: eta        ! additional damping `exp(eta)`
        VALUE_TYPE, allocatable :: vec(:,:)
                                ! buffer (potentially larger)
    end type STATE_VECTOR

    type, public, extends(STATE_VECTOR) :: KET
    end type KET

    type, public, extends(STATE_VECTOR) :: BRA
    end type BRA

    ! OPERATORS

    interface is_initialized
        module procedure p_is_initialized
    end interface
    public is_initialized

    interface init_operator
       module procedure p_init_operator_sizes
    end interface
    public init_operator

    interface init_operator_like
        module procedure p_init_operator_like
    end interface
    public init_operator_like

    interface dest_toperator
        module procedure p_dest_operator
    end interface
    public dest_toperator

    interface matmul
        module procedure p_mul_opers
    end interface
    public matmul

    interface add_to_oper
       module procedure p_add_to_oper
    end interface
    public add_to_oper

    interface oper_set_block
       module procedure p_oper_set_block
    end interface
    public oper_set_block

    interface oper_get_block
       module procedure p_oper_get_block
    end interface
    public oper_get_block

    interface adjoint
        module procedure p_adjoint
    end interface
    public adjoint

    !> Hermiticize an operator
    interface force_hermitian
        module procedure p_force_hermitian
    end interface
    public force_hermitian

    !> Diagonalize an operator
    interface diag_operator
        module procedure p_diag_operator
    end interface
    public diag_operator

    !> Get left/outgoing superstate (block) for given right/incoming one.
    interface connecting_sst
       module procedure p_connecting_sst
    end interface
    public connecting_sst

    interface start_sst_list
        module procedure p_start_sst_list
    end interface
    public start_sst_list

    interface connecting_sst_list
        module procedure p_connect_sst_list
    end interface
    public connecting_sst_list

    interface finish_sst_list
        module procedure p_finish_sst_list
    end interface
    public finish_sst_list

    !> Get right/incoming superstate (block) for given left/outgoing one.
    interface connecting_sst_backwards
       module procedure p_connecting_sst_backwards
    end interface
    public connecting_sst_backwards

    interface transform_inv_eigenbasis
       module procedure p_transform_inv_eigenbasis
    end interface
    public transform_inv_eigenbasis

    !> Get eigenvalues for a given block
    interface get_eigenvalues
       module procedure p_get_eigenvalues
    end interface
    public get_eigenvalues

    !> Get ground state energy
    interface get_egs
       module procedure p_get_egs
    end interface
    public get_egs

    !> Get size of largest superstate (block)
    interface max_sst_size
        module procedure p_max_sst_size, p_max_sst_size_op
    end interface
    public max_sst_size

    interface get_num_superstates
       module procedure p_get_num_superstates, p_get_num_superstates_op
    end interface
    public get_num_superstates

    interface get_sst_size
       module procedure p_get_sst_size, p_get_sst_size_op
    end interface
    public get_sst_size

    !> Return operator as matrix
    interface getmatrix
        module procedure p_get_matrix
    end interface
    public getmatrix

    !> Fill operator from given matrix
    interface oper_from_matrix
        module procedure p_from_matrix
    end interface
    public oper_from_matrix

    ! STATES

    !> Swap contents of two states
    interface swap
        module procedure :: t_swap_bra, t_swap_ket
    end interface
    public :: swap

    !> Set state to given contents
    interface set_state
        module procedure :: t_set
#if VALUE_IS_COMPLEX
        module procedure :: t_set_real
#endif
    end interface
    public :: set_state

    interface copy
        module procedure :: t_copy_ket, t_copy_bra
    end interface
    public :: copy

    interface get_sst
       module procedure :: state_get_sst
    end interface
    public get_sst

    !> Return value of inner product <bra|ket>
    interface braket_dotprod
        module procedure :: t_braket_dotprod
    end interface braket_dotprod
    public braket_dotprod

    !> Fill a state with an eigenstate in the eigenbasis (indicator)
    interface eigenstate_in_eigenbasis
        module procedure :: t_eigenstate_in_eigenbasis
    end interface
    public eigenstate_in_eigenbasis

    interface set_to_identity
        module procedure :: t_state_set_identity
    end interface
    public set_to_identity

    interface set_to_operator
        module procedure :: t_ket_set_operator, t_bra_set_operator
    end interface
    public set_to_operator


    !> Perform time evolution exp(-H*tau) |ket>
    interface time_evolve
        module procedure :: t_time_evolve
    end interface time_evolve
    public time_evolve

    !> Perform time evolution on the "right side", |ket_i> exp(-e_i*tau)
    interface relocate_state
        module procedure :: t_relocate_state
    end interface relocate_state
    public relocate_state

    !> Perform trace tr(exp(-H*beta) |ket>)
    interface trace_boltzmann
       module procedure :: ket_trace_boltzmann, bra_trace_boltzmann
    end interface
    public trace_boltzmann

    !> Compute partition function tr(exp(-beta * H))
    interface partition_function
       module procedure :: t_partition_function
    end interface
    public partition_function

    !> Return array of thermal weights for each block
    interface thermal_weights
       module procedure t_thermal_weights
    end interface
    public thermal_weights

    !> Override output vector with operator applied to given |ket>
    interface apply_operator
        module procedure :: t_apply_oper
    end interface
    public apply_operator

    !> Override output vector with adjoint operator applied to given <bra|
    interface apply_adjoint
        module procedure :: t_apply_adjoint
    end interface
    public apply_adjoint

    !> Perform outer product |bra><ket| and add it to an operator.
    interface outer_product
        module procedure t_outer_product
    end interface
    public outer_product

    !> Perform outer product |bra><ket|, integrated over some time interval.
    interface outer_integral
        module procedure t_outer_integral
    end interface
    public outer_integral

    interface dump
       module procedure t_dump_state
    end interface
    public dump

    interface assert_close
        module procedure ket_assert_close, bra_assert_close
    end interface
    public assert_close

#ifndef VALUE_IS_COMPLEX
    interface conjg
        module procedure conjg_real
    end interface
#endif

    VALUE_TYPE, parameter :: ONE = 1, ZERO = 0
    REAL_TYPE, parameter :: ONE_REAL = 1

    ! Smallest x such that (1 + exp(x)) is different from x. In double
    ! precision, this comes out to be about -36.7, for single, about -16.6.
    REAL_TYPE, parameter :: EXP_MINARG = log(0.5 * epsilon(ONE_REAL))

contains

    pure function p_is_initialized(self) result(test)
        type(OPERATOR), intent(in) :: self
        logical :: test

        test = allocated(self%offset)
    end

    subroutine p_init_operator_sizes(self, block_size, connect)
        type(OPERATOR), intent(out) :: self
        integer, intent(in) :: block_size(0:)
        integer, intent(in), optional :: connect(0:)
        integer :: i, n, lsize, rsize

        n = size(block_size)

        allocate(self%offset(0:n))
        self%offset(0) = 0
        do i = 0, n-1
            if (block_size(i) < 0) &
                error stop 'Invalid block size'
            self%offset(i+1) = self%offset(i) + block_size(i)
        enddo

        allocate(self%connect(0:n-1))
        if (present(connect)) then
            if (size(connect) /= n) &
                error stop 'Invalid size of connect'

            do i = 0, n-1
                if (connect(i) < -1 .or. connect(i) >= n) &
                    error stop 'Invalid value of connect'

                self%connect(i) = connect(i)
            enddo
        else
            self%connect(:) = -1
        endif

        allocate(self%subops(0:n-1))
        do i = 0, n-1
            rsize = block_size(i)
            if (self%connect(i) == -1) then
                lsize = 0
            else
                lsize = block_size(self%connect(i))
            endif

            allocate(self%subops(i)%op(0:lsize-1, 0:rsize-1))
            self%subops(i)%op(:, :) = 0
        enddo
    end subroutine

    subroutine p_init_operator_like(self, ref)
        type(OPERATOR), intent(out) :: self
        type(OPER_EIGEN), intent(in) :: ref
        integer :: i, n, b

        n = size(ref%offset) - 1
        allocate(self%offset, source=ref%offset)

        allocate(self%connect(0:n-1))
        allocate(self%SubOps(0:n-1))
        do i = 0, n-1
            self%connect(i) = i

            b = self%offset(i+1) - self%offset(i)
            allocate(self%subops(i)%op(0:b-1, 0:b-1))
            self%SubOps(i)%op(:, :) = 0
        enddo
    end subroutine

    subroutine p_dest_operator(self)
        type(OPERATOR), intent(inout) :: self

        if (is_initialized(self)) then
            deallocate(self%connect, self%offset, self%SubOps)
        endif
    end subroutine

    subroutine p_oper_set_block(self, i, j, val)
        type(OPERATOR), intent(inout) :: self
        integer, intent(in) :: i, j
        VALUE_TYPE, intent(in) :: val(0:, 0:)
        integer :: isize, jsize, jprev

        if (i < 0 .or. i >= size(self%connect)) &
            error stop 'Invalid left block'
        if (j < 0 .or. j >= size(self%connect)) &
            error stop 'Invalid right block'

        isize = self%offset(i+1) - self%offset(i)
        jsize = self%offset(j+1) - self%offset(j)
        if (size(val, 1) /= isize .or. size(val, 2) /= jsize) &
            error stop 'Size mismatch between values and block'

        ! Check if anything already connects to i
        jprev = connecting_sst_backwards(self, i)
        if (jprev == -1) then
            self%connect(j) = i
            deallocate(self%subops(j)%op)
            allocate(self%subops(j)%op(0:size(val,1)-1, 0:size(val,2)-1))
            self%subops(j)%op =val
        elseif (jprev == j) then
            self%subops(j)%op(:,:) = val
        else
            error stop 'Block collision'
        endif
    end subroutine

    subroutine p_oper_get_block(self, i, j, op)
        type(OPERATOR), intent(in) :: self
        integer, intent(in) :: i, j
        VALUE_TYPE, allocatable, intent(inout) :: op(:,:)
        integer :: nsstates

        nsstates = size(self%connect)
        if (i < 0 .or. i >= nsstates) &
            error stop 'Invalid left superstate'
        if (j < 0 .or. j >= nsstates) &
            error stop 'Invalid right superstate'

        if (allocated(op)) &
            deallocate(op)
        allocate(op(0:self%offset(i+1) - self%offset(i) - 1, &
                    0:self%offset(j+1) - self%offset(j) - 1))
        if (i == self%connect(j)) then
            op(:,:) = self%SubOps(j)%op(:,:)
        else
            op(:,:) = 0
        endif
    end subroutine

    function p_mul_opers(x, y) result(r)
        type(OPERATOR), target, intent(in) :: x, y
        type(OPERATOR) :: r
        integer :: i, nsstates, rconn, isize
        VALUE_TYPE, pointer :: xop(:, :), yop(:, :)

        if (.not. all(x%offset == y%offset)) &
            error stop 'Inconsistent sizes for operator'

        allocate(r%offset, source=x%offset)

        nsstates = size(x%connect)
        allocate(r%connect(0:nsstates-1))
        allocate(r%SubOps(0:nsstates-1))
        do i = 0, nsstates - 1
            isize = x%offset(i+1) - x%offset(i)
            rconn = y%connect(i)
            if (rconn /= -1) then
                rconn = x%connect(rconn)
            endif
            r%connect(i) = rconn

            if (rconn == -1) then
                allocate(r%SubOps(i)%op(0:-1, 0:isize-1))
            else
                xop => x%SubOps(y%connect(i))%op
                yop => y%SubOps(i)%op
                allocate(r%SubOps(i)%op(0:size(xop,1)-1, 0:size(yop,2)-1))
                r%SubOps(i)%op = matmul(xop, yop)
            endif
        enddo
    end function

    subroutine p_add_to_oper(r, term, scale)
        type(OPERATOR), intent(inout) :: r
        type(OPERATOR), intent(in) :: term
        VALUE_TYPE :: scale
        integer :: i, nsstates

        if (.not. all(r%offset == term%offset)) &
            error stop 'Inconsistent sizes for operator'

        nsstates = size(r%connect)
        do i = 0, nsstates - 1
            if (term%connect(i) == -1) then
                cycle
            elseif (r%connect(i) == -1) then
                r%connect(i) = term%connect(i)
                deallocate(r%SubOps(i)%op)
                allocate(r%SubOps(i)%op(0:size(term%SubOps(i)%op,1)-1, &
                                        0:size(term%SubOps(i)%op,2)-1))
                r%SubOps(i)%op(:, :) = scale * term%SubOps(i)%op
            elseif (r%connect(i) == term%connect(i)) then
                r%SubOps(i)%op = r%SubOps(i)%op + scale * term%SubOps(i)%op
            else
                error stop 'Block structure is violated'
            endif
        enddo
    end subroutine

    subroutine p_diag_Operator(this, Eigen, epsevec)
        type(OPERATOR), intent(in) :: this
        type(OPER_EIGEN), intent(out) :: Eigen
        REAL_TYPE, intent(in), optional :: epsevec
    !output
    !local
        VALUE_TYPE, allocatable :: c_WORK(:)
        REAL_TYPE, allocatable :: RWORK(:)
        integer :: iSSt, INFO, NElements, NStatesMax, NSStates


        NStatesMax = max_sst_size(this)
        allocate(c_WORK(max(1, 3*NStatesMax)))
        allocate(RWORK(max(1, 3*NStatesMax-2)))

        NSStates = size(this%connect)
        allocate(Eigen%SubEigen(0:NSStates-1))

        do iSSt = 0, NSStates - 1
            if (this%connect(iSSt) /= iSSt) &
                error stop 'Operator is not block-diagonal'
        enddo

        !!! here comes the complex version
        do iSSt=0,NSStates-1
            NElements=size(this%SubOps(iSSt)%op, 2)
            allocate(Eigen%SubEigen(iSSt)%Eval(0:NElements-1))
            allocate(Eigen%SubEigen(iSSt)%Evec, source=this%SubOps(iSSt)%op)

            ! TODO: higher precision can be reached by using ?syevr?
            call HEEV('V', 'U', NElements, &
                    Eigen%SubEigen(iSSt)%Evec(0, 0), NElements,&
                    Eigen%SubEigen(iSSt)%Eval(0), &
                    c_WORK(1), 3*NElements, RWORK(1), INFO)

            if(INFO.ne.0)then
                write(0, *) "HEEV failed with INFO =", INFO
                error stop 'Failed to diagonalize'
            endif

            if (present(epsevec)) &
                call process_eigenvectors(Eigen%SubEigen(iSSt)%Evec, epsevec)
        enddo

        call recompute_egs(Eigen)
        allocate(Eigen%offset, source=this%offset)
    end subroutine

    function p_adjoint(self) result(adj)
        type(OPERATOR), intent(in) :: self
        type(OPERATOR) :: adj
        integer :: i, j, n

        n = get_num_superstates(self)
        call init_operator(adj, &
                (/ (self%offset(i+1) - self%offset(i), i = 0, n-1) /), &
                (/ (connecting_sst_backwards(self, i), i = 0, n-1) /) &
                )

        do i = 0, n - 1
            j = self%connect(i)
            if (j == -1) &
                cycle

            adj%SubOps(j)%op(:, :) = conjg(transpose(self%SubOps(i)%op))
        enddo
    end function

    subroutine p_force_hermitian(H)
        type(OPERATOR), intent(inout) :: H
        integer :: sst

        do sst = 0, size(H%SubOps) - 1
            call hermiticize(H%SubOps(sst)%op)
        end do
    end subroutine

    subroutine hermiticize(mat)
        VALUE_TYPE, intent(inout) :: mat(:, :)
        VALUE_TYPE :: elem
        integer :: i, j

        if (size(mat, 1) /= size(mat, 2)) &
            error stop 'Matrix is not square'

        do j = 1, size(mat, 2)
            do i = j, size(mat, 1)
                elem = 0.5 * (mat(i, j) + conjg(mat(j, i)))
                mat(i, j) = elem
                mat(j, i) = conjg(elem)
            enddo
        enddo
    end subroutine

    subroutine process_eigenvectors(vecs, epsevec)
        VALUE_TYPE, intent(inout) :: vecs(0:, 0:)
        REAL_TYPE, intent(in) :: epsevec
        integer :: i, j, n

        if (epsevec < 0) &
            error stop 'Invalid value of epsevec'
        if (size(vecs, 1) /= size(vecs, 2)) &
            error stop 'Matrix must be square'

        n = size(vecs, 1)
        do i=0, n-1
            do j=0, n-1
                if (abs(vecs(j, i)) < epsevec) then
                    vecs(j,i)=0d0
                endif
            enddo
            vecs(:,i) = vecs(:,i) / sum(abs(vecs(:,i))**2)
        enddo
    end subroutine

    function p_connecting_sst(self, sst) result(next_sst)
        type(OPERATOR), intent(in) :: self
        integer, intent(in) :: sst
        integer :: next_sst

        if (.not. allocated(self%connect)) &
            error stop 'Operator is not initialized'

        if (sst == -1) then
            ! A super state of -1 is mapped to itself
            next_sst = -1
        else if (sst < 0 .or. sst >= size(self%connect)) then
            error stop 'Invalid incoming sst'
        else
            next_sst = self%connect(sst)
        endif
    end function

    subroutine p_start_sst_list(self, valid_out, nvalid_out)
        type(OPER_EIGEN), intent(in) :: self
        integer, intent(out) :: valid_out(*), nvalid_out
        integer :: sst_in

        if (.not. allocated(self%offset)) &
            error stop 'Operator is not initialized'

        nvalid_out = size(self%offset)-1
        do sst_in = 0, nvalid_out-1
            valid_out(sst_in+1) = sst_in
        enddo
    end subroutine

    subroutine p_connect_sst_list(self, valid, nvalid)
        type(OPERATOR), intent(in) :: self
        integer, intent(inout) :: valid(nvalid)
        integer, intent(in) :: nvalid
        integer :: i, sst_in, sst_out

        if (.not. allocated(self%connect)) &
            error stop 'Operator is not initialized'

        !write (0,'(A,100I4)') 'CONNECT', self%connect
        do i = 1, nvalid
            sst_in = valid(i)
            if (sst_in < 0) &
                cycle
            sst_out = self%connect(sst_in)
            valid(i) = sst_out
        enddo
    end

    subroutine p_finish_sst_list(self, valid, nvalid)
        type(OPER_EIGEN), intent(in) :: self
        integer, intent(inout) :: valid(nvalid), nvalid
        integer :: i, sst_in, nvalid_out

        !write (0,'(A,100I4)') 'CONNECT', self%connect
        nvalid_out = 0
        do i = 1, nvalid
            sst_in = valid(i)
            if (sst_in == i - 1) then
                nvalid_out = nvalid_out + 1
                valid(nvalid_out) = sst_in
            endif
        enddo
        nvalid = nvalid_out
    end

    function p_connecting_sst_backwards(self, sst) result(next_sst)
        type(OPERATOR), intent(in) :: self
        integer, intent(in) :: sst
        integer :: next_sst

        if (.not. allocated(self%connect)) &
            error stop 'Operator is not initialized'

        if (sst == -1) then
            ! A super state of -1 is mapped to itself
            next_sst = -1
        else if (sst < 0 .or. sst >= size(self%connect)) then
            error stop 'Invalid incoming sst'
        else
            do next_sst = 0, size(self%connect) - 1
                if (self%connect(next_sst) == sst) &
                    return
            enddo
            next_sst = -1
        endif
    end function

    subroutine p_transform_inv_eigenbasis(eigen, op)
        type(OPER_EIGEN), intent(in) :: Eigen
        type(OPERATOR), intent(inout) :: op
        integer :: i, j, nmax
        VALUE_TYPE, allocatable :: work(:, :)

        nmax = max_sst_size(eigen)
        allocate(work(nmax, nmax))

        do j = 0, size(op%connect) - 1
            i = op%connect(j)
            if (i == -1) &
                cycle

            call p_transform_inv_eigenbasis_mat(&
                    op%SubOps(j)%op, &
                    Eigen%SubEigen(i)%evec, &
                    Eigen%SubEigen(j)%evec, &
                    work)
        enddo
    end subroutine

    subroutine p_transform_inv_eigenbasis_mat(op, lbasis, rbasis, work)
        VALUE_TYPE, intent(inout) :: op(:,:), work(:,:)
        VALUE_TYPE, intent(in) :: lbasis(:,:), rbasis(:,:)
        integer :: m, n

        if (size(lbasis, 1) /= size(lbasis, 2)) &
            error stop 'lbasis must be square'
        if (size(rbasis, 1) /= size(rbasis, 2)) &
            error stop 'rbasis must be square'
        if (size(op, 1) /= size(lbasis, 1) .or. size(op, 2) /= size(rbasis, 1)) &
            error stop 'mismatch between op and basis'

        m = size(op, 1)
        n = size(op, 2)
        work(:,:) = 0
        call GEMM('N', 'N', m, n, n, &
                  ONE, op(1, 1), size(op, 1), rbasis(1, 1), size(rbasis, 1), &
                  ZERO, work(1, 1), size(work, 1))
        call GEMM('C', 'N', m, n, m, &
                  ONE, lbasis(1, 1), size(lbasis, 1), work(1, 1), size(work, 1), &
                  ZERO, op(1, 1), size(op, 1))
    end subroutine

    function p_get_eigenvalues(self, sst) result(evals)
        type(OPER_EIGEN), intent(in) :: self
        integer, intent(in) :: sst
        REAL_TYPE, allocatable :: evals(:)

        allocate(evals, source=self%SubEigen(sst)%eval)
    end function

    function p_get_egs(self) result(egs)
        type(OPER_EIGEN), intent(in) :: self
        REAL_TYPE, allocatable :: egs

        egs = self%egs
    end function

    subroutine recompute_egs(self)
        type(OPER_EIGEN), intent(inout) :: self
        integer :: isst

        self%egs = self%SubEigen(0)%EVal(0)
        do iSSt = 0, size(self%SubEigen) - 1
           if (self%SubEigen(iSSt)%EVal(0) < self%Egs)&
              self%Egs = self%SubEigen(iSSt)%EVal(0)
        end do
    end

    function p_max_sst_size(self) result(n)
        type(OPER_EIGEN), intent(in) :: self
        integer :: n

        n = p_max_sst_size_offset(self%offset)
    end

    function p_max_sst_size_op(self) result(n)
        type(OPERATOR), intent(in) :: self
        integer :: n

        n = p_max_sst_size_offset(self%offset)
    end

    function p_max_sst_size_offset(offset) result(n)
        integer, intent(in) :: offset(0:)
        integer :: n, isst

        n = offset(0)
        do isst = 1, ubound(offset,1) - 1
            n = max(n, offset(isst+1) - offset(isst))
        enddo
    end

    function p_get_num_superstates(self) result(n)
        type(OPER_EIGEN), intent(in) :: self
        integer :: n

        n = size(self%SubEigen)
    end

    function p_get_num_superstates_op(self) result(n)
        type(OPERATOR), intent(in) :: self
        integer :: n

        n = size(self%SubOps)
    end

    function p_get_sst_size(self, sst) result(n)
        type(OPER_EIGEN), intent(in) :: self
        integer, intent(in) :: sst
        integer :: n

        n = p_get_sst_size_offset(self%offset, sst)
    end

    function p_get_sst_size_op(self, sst) result(n)
        type(OPERATOR), intent(in) :: self
        integer, intent(in) :: sst
        integer :: n

        n = p_get_sst_size_offset(self%offset, sst)
    end

    function p_get_sst_size_offset(offset, sst) result(n)
        integer, intent(in) :: offset(0:)
        integer, intent(in) :: sst
        integer :: n

        if (sst < 0 .or. sst >= ubound(offset, 1)) &
            error stop 'Invalid superstate'

        n = offset(sst+1) - offset(sst)
    end

    function p_get_matrix(self) result(mat)
        type(OPERATOR), intent(in) :: self
        VALUE_TYPE, allocatable :: mat(:, :)
        integer :: nn, is, js

        nn = self%offset(size(self%offset) - 1)
        allocate(mat(nn, nn))
        mat(:, :) = 0

        do js = 0, size(self%SubOps) - 1
            is = self%connect(js)
            mat(self%offset(is)+1:self%offset(is+1), &
                self%offset(js)+1:self%offset(js+1)) = self%SubOps(js)%op
        enddo
    end

    subroutine p_from_matrix(self, mat)
        type(OPERATOR), intent(inout) :: self
        VALUE_TYPE, intent(in) :: mat(:,:)
        integer :: nn, is, js, imin, imax, jmin, jmax

        nn = self%offset(size(self%offset) - 1)
        if (size(mat, 1) /= nn .or. size(mat, 2) /= nn) &
            error stop 'Size mismatch'

        do js = 0, size(self%SubOps) - 1
            is = self%connect(js)
            imin = self%offset(is) + 1
            imax = self%offset(is + 1)
            jmin = self%offset(js) + 1
            jmax = self%offset(js + 1)

            if (any(mat(:imin-1, jmin:jmax) /= 0)) &
                error stop 'Values outside of block'
            if (any(mat(imax+1:, jmin:jmax) /= 0)) &
                error stop 'Values outside of block'

            self%SubOps(js)%op(:,:) = self%SubOps(js)%op(:,:) &
                    + mat(imin:imax, jmin:jmax)
        enddo
    end

    ! -------------------------------------------------------------------------
    ! STATES

    pure subroutine t_swap(x, y)
        class(STATE_VECTOR), intent(inout) :: x, y

        call swap(x%sst, y%sst)
        call swap(x%sstsize, y%sstsize)
        call swap(x%nvec, y%nvec)
        call swap_real(x%eta, y%eta)
        call swap_vec(x%vec, y%vec)
    end subroutine

    pure subroutine t_swap_ket(x, y)
        type(KET), intent(inout) :: x, y

        call t_swap(x, y)
    end subroutine

    pure subroutine t_swap_bra(x, y)
        type(BRA), intent(inout) :: x, y

        call t_swap(x, y)
    end subroutine

    pure subroutine swap_vec(x, y)
        VALUE_TYPE, allocatable, intent(inout) :: x(:,:), y(:,:)
        VALUE_TYPE, allocatable :: tmp(:,:)

        call move_alloc(x, tmp)
        call move_alloc(y, x)
        call move_alloc(tmp, y)
    end

    pure subroutine swap_real(a, b)
        REAL_TYPE, intent(inout) :: a, b
        REAL_TYPE :: tmp

        tmp = a
        a = b
        b = tmp
    end subroutine

    subroutine t_set(self, sst, val, exp)
        class(STATE_VECTOR), intent(inout) :: self
        integer, intent(in) :: sst
        integer, intent(in), optional :: exp
        VALUE_TYPE, intent(in) :: val(:)

        call p_resize_state(self, size(val), 1)
        self%sst = sst
        self%vec(:size(val), 1) = val
        if (present(exp)) then
            self%eta = log(2.0d0) * exp
        else
            self%eta = 0
        endif
    end subroutine

#if VALUE_IS_COMPLEX
    subroutine t_set_real(self, sst, val, exp)
        class(STATE_VECTOR), intent(inout) :: self
        integer, intent(in) :: sst
        integer, intent(in), optional :: exp
        REAL_TYPE, intent(in) :: val(:)
        VALUE_TYPE :: val_complex(size(val))

        val_complex = val
        call set_state(self, sst, val_complex, exp)
    end subroutine
#endif

    subroutine t_copy_ket(state_in, state_out)
        class(KET), intent(in) :: state_in
        class(KET), intent(inout) :: state_out

        call t_copy(state_in, state_out)
    end

    subroutine t_copy_bra(state_in, state_out)
        class(BRA), intent(in) :: state_in
        class(BRA), intent(inout) :: state_out

        call t_copy(state_in, state_out)
    end

    subroutine t_copy(state_in, state_out)
        class(STATE_VECTOR), intent(in) :: state_in
        class(STATE_VECTOR), intent(inout) :: state_out

        call p_resize_state(state_out, state_in%sstsize, state_in%nvec)
        state_out%sst = state_in%sst
        state_out%sstsize = state_in%sstsize
        state_out%nvec = state_in%nvec
        state_out%eta = state_in%eta
        state_out%vec(:state_in%sstsize, :state_in%nvec) = &
            state_in%vec(:state_in%sstsize, :state_in%nvec)
    end subroutine

    subroutine t_add_canary(self)
        class(STATE_VECTOR), intent(inout) :: self

        if (self%sstsize < size(self%vec, 1)) then
            self%vec(self%sstsize+1, 1) = nan(1.0d0)
        endif
        if (self%nvec < size(self%vec, 2)) then
            self%vec(1, self%nvec+1) = nan(1.0d0)
        endif
    end

    subroutine t_remove_canary(self)
        class(STATE_VECTOR), intent(inout) :: self
        VALUE_TYPE :: v

        if (self%sstsize < size(self%vec, 1)) then
            v = self%vec(self%sstsize+1, 1)
            if (v == v) &
                error stop 'Canary missing along row dimension'
            self%vec(self%sstsize+1, 1) = ZERO
        endif
        if (self%nvec < size(self%vec, 2)) then
            v = self%vec(1, self%nvec+1)
            if (v == v) &
                error stop 'Canary missing along column dimension'
            self%vec(1, self%nvec+1) = ZERO
        endif
    end

    subroutine p_resize_state(self, sstsize, nvec)
        class(STATE_VECTOR), intent(inout) :: self
        integer, intent(in) :: sstsize, nvec

        if (sstsize < 0 .or. nvec < 1) &
            error stop 'Invalid size of sst or nvec'

        if (.not. allocated(self%vec)) then
            allocate(self%vec(max(sstsize, 1), nvec))
            self%vec(:, :) = 0
        else
            call t_remove_canary(self)
            if (sstsize > size(self%vec, 1) .or. nvec > size(self%vec, 2)) then
                deallocate(self%vec)
                allocate(self%vec(max(sstsize, 1), nvec))
                self%vec(:, :) = 0
            endif
        endif

        self%sstsize = sstsize
        self%nvec = nvec
        call t_add_canary(self)
    end subroutine

    pure elemental integer function state_get_sst(state) result(sst)
        class(STATE_VECTOR), intent(in) :: state

        sst = state%sst
    end

    ! Computes the scalar product of a bra state b and a ket state k
    ! with separately stored exponents resulting in a number with
    ! extended range of type EXT_RANGE.
    !
    ! The sst and sstsize of the bra and ket must be the same,
    ! but this is unchecked.
    function t_braket_dotprod(b, k) result(bk)
        type(BRA), intent(in) :: b
        type(KET), intent(in) :: k
        type(EXT_RANGE) :: bk
        VALUE_TYPE :: braket
        integer :: i

        if (b%sst /= k%sst) then
            bk = ZERO
        else
            if (b%sstsize /= k%sstsize .or. b%nvec /= k%nvec) &
                error stop 'size of vector differs'
            braket = ZERO
            do i = 1, k%nvec
                braket = braket + DOTC(b%sstsize, b%vec(1, i), 1, k%vec(1, i), 1)
            enddo
            bk = extended_exp(b%eta + k%eta) * braket
        endif
    end

    ! Computes the time evolution of a state with energies taken relative to E0.
    ! For a bra state, this computes
    !    <state| := <state| exp(-(H - E0) * dtau)
    ! and for a ket state, this computes
    !    |state> := exp(-(H - E0) * dtau) |state>
    ! (Note that the operator is the same in both cases, so to evolve a bra from
    ! higher imaginary time to lower imaginary time dtau must be positive.)
    pure subroutine t_time_evolve(H_eigen, dtau, state)
        type(OPER_EIGEN), intent(in) :: H_eigen
        REAL_TYPE, intent(in) :: dtau
        class(STATE_VECTOR), intent(inout) :: state
        REAL_TYPE :: sstE0, damping
        integer :: i

        if (.not. (dtau >= 0)) &
            error stop 'Time evolution can only move forward'

        ! If sst is -1, we are in the QN violated block - time evolution is
        ! trivial in that case
        if (state%sst == -1) &
            return

        ! use sst's lowest eigenvalue as intermediate 0 reference value. the
        ! difference between this and the absolute ground state is
        ! directly added to logarithmic accounting. (potential
        ! improvement: ideally, the value used here would keep the norm
        ! of vecout as close as possible to that of vecin, possibly also
        ! making the renormalizig unnecessary, but that would also
        ! depend on vecin and likely not be worth it in total)
        sstE0 = H_eigen%SubEigen(state%sst)%Eval(0)  ! E_int

        ! logarithmic accounting for difference between lowest energy in
        ! sst and ground state energy
        state%eta = state%eta - (sstE0 - H_eigen%Egs) * dtau
        do i = 1, state%sstsize - 1
            damping = exp(-dtau * (H_eigen%SubEigen(state%sst)%Eval(i) - sstE0))
            state%vec(i+1, :) = damping * state%vec(i+1, :)
        enddo
    end

    subroutine t_relocate_state(H, sst_right, dtau, state)
        type(OPER_EIGEN), intent(in) :: H
        integer, intent(in) :: sst_right
        REAL_TYPE, intent(in) :: dtau
        class(STATE_VECTOR), intent(inout) :: state

        REAL_TYPE :: E0, damping
        integer :: i

        if (.not. (dtau >= 0)) &
            error stop 'Time evolution can only move forward'
        if (sst_right < 0) &
            error stop 'right sst is invalid'
        if (state%nvec /= size(H%SubEigen(sst_right)%Eval)) &
            error stop 'size mismatch on right side'

        E0 = H%SubEigen(sst_right)%Eval(0)

        state%eta = state%eta - (E0 - H%Egs) * dtau
        do i = 1, state%nvec - 1
            damping = exp(-dtau * (H%SubEigen(sst_right)%Eval(i) - E0))
            state%vec(:, i+1) = damping * state%vec(:, i+1)
        enddo
    end

    function ket_trace_boltzmann(H_eigen, beta, state) result(tr)
        type(OPER_EIGEN), target, intent(in) :: H_eigen
        REAL_TYPE, intent(in) :: beta
        type(KET), intent(in) :: state
        type(EXT_RANGE) :: tr

        tr = state_trace_boltzmann(H_eigen, beta, state)
    end

    function bra_trace_boltzmann(H_eigen, beta, state) result(tr)
        type(OPER_EIGEN), target, intent(in) :: H_eigen
        REAL_TYPE, intent(in) :: beta
        type(BRA), intent(in) :: state
        type(EXT_RANGE) :: tr

        tr = conjg(state_trace_boltzmann(H_eigen, beta, state))
    end

    function state_trace_boltzmann(H_eigen, beta, state) result(tr)
        type(OPER_EIGEN), target, intent(in) :: H_eigen
        REAL_TYPE, intent(in) :: beta
        class(STATE_VECTOR), intent(in) :: state
        type(EXT_RANGE) :: tr
        REAL_TYPE :: E0
        VALUE_TYPE :: tr_diff

        ! XXX revisit once we can "drop" states
        if (state%sstsize /= state%nvec) &
            error stop 'state must be square'
        if (state%sst == -1) then
            tr = ZERO
            return
        endif

        tr_diff = trace_boltzmann_shifted( &
                        H_eigen%SubEigen(state%sst), beta, state%vec, E0)
        tr = extended_exp(-beta * (E0 - H_eigen%Egs) + state%eta) * tr_diff
    end

    function trace_boltzmann_shifted(H_sub, beta, vec, E0) result(tr_diff)
        type(TSubEigen), intent(in) :: H_sub
        REAL_TYPE, intent(in) :: beta
        REAL_TYPE, intent(out) :: E0
        VALUE_TYPE, intent(in) :: vec(:, :)
        VALUE_TYPE :: tr_diff
        REAL_TYPE :: arg
        integer :: i

        tr_diff = ZERO
        E0 = H_sub%eval(0)
        do i = 0, size(H_sub%eval)-1
            arg = -beta * (H_sub%eval(i) - E0)
            if (arg < EXP_MINARG) &
                exit
            tr_diff = tr_diff + exp(arg) * vec(i+1, i+1)
        enddo
    end

    function t_partition_function(H_eigen, beta, outer_sst) result(Z)
        type(OPER_EIGEN), target, intent(in) :: H_eigen
        REAL_TYPE, intent(in) :: beta
        integer, optional, intent(in) :: outer_sst
        type(DExtRange) :: Z
        REAL_TYPE :: Zdiff, E0
        integer :: sst

        if (.not. (beta >= 0 .and. beta <= huge(beta))) &
            error stop 'Invalid inverse temperature'

        ! We have that:
        !
        !   Z = tr exp(-beta * (H - E0)))
        !     = exp(-beta * (E1 - E0)) * sum(exp(-beta * (E - E1))),
        !
        ! The sum is then 1 plus smaller terms, which is guaranteed to not
        ! overflow.
        if (present(outer_sst)) then
            if (outer_sst < 0 .or. outer_sst >= size(H_eigen%SubEigen)) &
                error stop 'Invalid sst'

            E0 = H_eigen%SubEigen(outer_sst)%eval(0)
            Zdiff = p_partition_function_part( &
                        H_eigen%SubEigen(outer_sst), beta, E0)
            Z = exp(DExtRange(-beta * (E0 - H_eigen%Egs))) * Zdiff
        else
            E0 = H_eigen%Egs
            Zdiff = 0
            do sst = 0, size(H_eigen%SubEigen)-1
                Zdiff = Zdiff + p_partition_function_part( &
                                    H_eigen%SubEigen(sst), beta, H_eigen%Egs)
            enddo
            Z = Zdiff
        endif
    end

    function t_thermal_weights(H_eigen, beta) result(weights)
        type(OPER_EIGEN), intent(in) :: H_eigen
        REAL_TYPE, intent(in) :: beta
        REAL_TYPE, allocatable :: weights(:)
        REAL_TYPE :: Z
        integer :: sst

        if (.not. (beta >= 0 .and. beta <= huge(beta))) &
            error stop 'Invalid inverse temperature'

        allocate(weights(0:size(H_eigen%SubEigen)-1))

        Z = 0
        do sst = 0, size(weights) - 1
            weights(sst) = p_partition_function_part( &
                                H_eigen%SubEigen(sst), beta, H_eigen%Egs)
            Z = Z + weights(sst)
        enddo
        weights(:) = weights(:) / Z
    end

    pure function p_partition_function_part(Hsub, beta, E0) result(Zdiff)
        type(TSubEigen), intent(in) :: Hsub
        REAL_TYPE, intent(in) :: beta, E0
        REAL_TYPE :: Zdiff, arg
        integer :: i

        Zdiff = 0
        do i = lbound(Hsub%eval,1), ubound(Hsub%eval,1)
            arg = -beta * (Hsub%eval(i) - E0)
            if (arg < EXP_MINARG) &
               exit
            Zdiff = Zdiff + exp(arg)
        enddo
    end

    ! Apply operator op to state state,
    !    |state> := op |state>
    subroutine t_apply_oper(vout, op, vin)
        type(KET), intent(inout) :: vout
        type(OPERATOR), intent(in), target :: op
        type(KET), intent(in) :: vin
        VALUE_TYPE, pointer :: block_op(:,:)

        if (.not. allocated(vin%vec)) &
            error stop 'vin not initialized'

        vout%sst = connecting_sst(op, vin%sst)
        if (vout%sst == -1) then
            call p_resize_state(vout, 0, vin%nvec)
            vout%eta = 0
        else
            block_op => op%SubOps(vin%sst)%op
            call p_resize_state(vout, size(block_op, 1), vin%nvec)
            vout%eta = vin%eta

            call GEMM('N', 'N', vout%sstsize, vout%nvec, vin%sstsize, ONE, &
                      op%SubOps(vin%sst)%op(0, 0), size(block_op, 1), &
                      vin%vec(1, 1), size(vin%vec, 1), &
                      ZERO, vout%vec(1, 1), size(vout%vec, 1))
        endif
    end

    ! Apply operator op to state state,
    !    <state| := <state| op = ( op^\dag |state> )^\dag
    subroutine t_apply_adjoint(vout, op, vin)
        type(BRA), intent(inout) :: vout
        type(OPERATOR), intent(in), target :: op
        type(BRA), intent(in) :: vin
        VALUE_TYPE, pointer :: block_op(:,:)

        if (.not. allocated(vin%vec)) &
            error stop 'vin not initialized'

        vout%sst = connecting_sst_backwards(op, vin%sst)
        if (vout%sst == -1) then
            call p_resize_state(vout, 0, vin%nvec)
            vout%eta = 0
        else
            block_op => op%SubOps(vout%sst)%op
            call p_resize_state(vout, size(block_op, 2), vin%nvec)
            vout%eta = vin%eta

            call GEMM('C', 'N', vout%sstsize, vout%nvec, vin%sstsize, ONE, &
                      op%SubOps(vout%sst)%op(0, 0), size(block_op, 1), &
                      vin%vec(1, 1), size(vin%vec, 1), &
                      ZERO, vout%vec(1, 1), size(vout%vec, 1))
        endif
    end

    ! XXX remove
    subroutine t_eigenstate_in_eigenbasis(stvec, eigen, isst, ist)
        class(STATE_VECTOR), intent(inout) :: stvec
        type(OPER_EIGEN), intent(in) :: eigen
        integer, intent(in) :: isst, ist

        if (isst < 0 .or. isst >= size(eigen%SubEigen)) &
            error stop 'Invalid superstate number'
        if (ist < 1 .or. ist > size(eigen%SubEigen(isst)%eval)) &
            error stop 'Invalid state number'

        stvec%sst = isst
        call p_resize_state(stvec, size(eigen%SubEigen(isst)%eval), 1)
        stvec%eta = 0
        stvec%vec(:stvec%sstsize, :stvec%nvec) = ZERO
        stvec%vec(ist, 1) = ONE
    end

    subroutine t_state_set_identity(self, ham, sst)
        class(STATE_VECTOR), intent(inout) :: self
        class(OPER_EIGEN) :: ham
        integer, intent(in) :: sst
        integer :: i, sstsize

        if (sst < -1) &
            error stop 'Invalid sst'

        sstsize = get_sst_size(ham, sst)
        call p_resize_state(self, sstsize, sstsize)
        self%sst = sst
        self%eta = 0
        self%vec(:sstsize, :sstsize) = ZERO
        do i = 1, sstsize
            self%vec(i, i) = ONE
        enddo
    end subroutine

    subroutine t_ket_set_operator(self, op, sst_in)
        type(KET), intent(inout) :: self
        type(OPERATOR), intent(in) :: op
        integer, intent(in) :: sst_in
        integer :: rows, cols

        if (sst_in < 0 .or. sst_in >= size(op%connect)) &
            error stop 'Invalid SST'

        rows = size(op%SubOps(sst_in)%op, 1)
        cols = size(op%SubOps(sst_in)%op, 2)
        call p_resize_state(self, rows, cols)
        self%sst = op%connect(sst_in)
        self%eta = 0
        self%vec(:rows, :cols) = op%SubOps(sst_in)%op(:, :)
    end

    subroutine t_bra_set_operator(self, op, sst_in)
        type(BRA), intent(inout) :: self
        type(OPERATOR), intent(in) :: op
        integer, intent(in) :: sst_in
        integer :: rows, cols, sst_adj

        if (sst_in < 0 .or. sst_in >= size(op%connect)) &
            error stop 'Invalid SST'

        sst_adj = connecting_sst_backwards(op, sst_in)
        rows = size(op%SubOps(sst_adj)%op, 2)
        cols = size(op%SubOps(sst_adj)%op, 1)
        call p_resize_state(self, rows, cols)
        self%sst = sst_adj
        self%eta = 0
        self%vec(:rows, :cols) = conjg(transpose(op%SubOps(sst_adj)%op))
    end subroutine

    subroutine t_outer_integral(op, b, k, tau, H, alpha, work)
        type(OPERATOR), intent(inout) :: op, work
        type(KET), intent(in) :: k
        type(BRA), intent(in) :: b
        REAL_TYPE, intent(in) :: tau
        type(OPER_EIGEN), intent(in) :: H
        type(EXT_RANGE), intent(in) :: alpha

        if (k%sst /= b%sst) &
            error stop 'Ket and bra must belong to the same block for now'

        !alpha_x = alpha * EXT_RANGE(exp(DExtRange(-tau * H%egs)))
        call outer_product(work, b, k, alpha, ZERO)
        call weigh_with_kernel( &
                H%SubEigen(k%sst)%eval(:) - H%Egs, tau, work%SubOps(k%sst)%op)
        op%SubOps(k%sst)%op(:, :) = op%SubOps(k%sst)%op(:, :) &
                + work%SubOps(k%sst)%op(:, :)
    end

    subroutine t_outer_product(op, b, k, alpha, gamma)
        type(OPERATOR), intent(inout), target :: op
        type(KET), intent(in) :: k
        type(BRA), intent(in) :: b
        type(EXT_RANGE), intent(in) :: alpha
        VALUE_TYPE, intent(in) :: gamma
        VALUE_TYPE :: alpha_d
        VALUE_TYPE, pointer :: subop(:,:)

        if (connecting_sst(op, b%sst) /= k%sst) &
            error stop 'Not compatible with operator'
        if (k%nvec /= b%nvec) &
            error stop 'Number of components mismatch'

        subop => op%SubOps(b%sst)%op
        alpha_d = limrange(extended_exp(k%eta + b%eta) * alpha)

        ! Perform op = gamma * op + alpha * |k> * <b|
        call GEMM('N', 'C', k%sstsize, b%sstsize, k%nvec, alpha_d, &
                  k%vec(1, 1), size(k%vec, 1), b%vec(1, 1), size(b%vec, 1), &
                  gamma, subop(0, 0), size(subop, 1))
    end

    subroutine t_dump_state(x, unit)
        class(STATE_VECTOR), intent(in) :: x
        integer, intent(in), optional :: unit
        integer :: runit

        runit = default(unit, 0)
        write (runit, 100) x%sst, x%eta, size(x%vec, 1), size(x%vec, 2)
        call print_array(x%vec(:x%sstsize, :x%nvec), file=runit)

    100 format ('sst=', I0, ', eta=', F0.6, ' (cap=', I0, ' x ', I0, '):')
    end

    subroutine ket_assert_close(a, b, rtol)
        type(KET), intent(in) :: a, b
        REAL_TYPE, intent(in) :: rtol

        call state_assert_close(a, b, rtol)
    end subroutine

    subroutine bra_assert_close(a, b, rtol)
        type(BRA), intent(in) :: a, b
        REAL_TYPE, intent(in) :: rtol

        call state_assert_close(a, b, rtol)
    end subroutine

    subroutine state_assert_close(a, b, rtol)
        class(STATE_VECTOR), intent(in) :: a, b
        REAL_TYPE, intent(in) :: rtol
        character(len=70) :: reason

        REAL_TYPE :: alpha
        VALUE_TYPE :: diff(a%sstsize, a%nvec)

        if (a%sst /= b%sst) then
            reason = 'sst is different'
        else if (a%sstsize /= b%sstsize) then
            reason = 'sstsize is different'
        else if (a%nvec /= b%nvec) then
            reason = 'nvec is different'
        else
            alpha = exp(a%eta - b%eta)
            diff = alpha * a%vec(:a%sstsize,:a%nvec) - b%vec(:a%sstsize,:a%nvec)
            if (all(abs(diff) <= rtol)) &
                return
            reason = 'elements are different'
        endif

        call dump(a)
        call dump(b)
        write (0,*) reason
        error stop 'Assertion failed - states are different'
    end subroutine

    !> Computes the kernel for a bosonic zero-frequency object.
    !!
    !! Computes the following integral exactly:
    !! $$
    !!   K_{ij}(\tau) :=
    !!          \int_0^\tau d\tau' \exp(-\tau' E_i) \exp(-(\tau - \tau') E_j),
    !! $$
    !! where `E` gives the energies and `tau` is some non-negative imaginary
    !! time distance.
    subroutine weigh_with_kernel(e, tau, out)
        REAL_TYPE, intent(in) :: e(:), tau
        VALUE_TYPE, intent(inout) :: out(:, :)

        integer :: i, j
        REAL_TYPE :: arg, z(size(e))

        !write (0,*) "E=", E

        if (size(out, 1) /= size(e) .or. size(out, 2) /= size(e)) &
            error stop 'Invalid size of out array'
        if (.not. (tau >= 0 .and. tau <= huge(tau))) &
            error stop 'Time distance must be positive'

        ! The value of the kernel is
        !   K(i,j) = (exp(-tau*E(i)) - exp(-tau*E(j))) / (E(j) - E(i))
        ! Precompute the exponentials, because they are expensive
        z(:) = exp(-tau * e(:))
        do j = 1, size(e)
            do i = 1, size(e)
                arg = tau * (e(j) - e(i))
                if (arg == 0) then
                    out(i, j) = out(i, j) * tau * z(j)
                elseif (abs(arg) < 0.25) then
                    out(i, j) = out(i, j) * z(j) * libc_expm1(arg) / (e(j) - e(i))
                else
                    out(i, j) = out(i, j) * (z(i) - z(j)) / (e(j) - e(i))
                endif
            enddo
        enddo
    end subroutine

#ifndef VALUE_IS_COMPLEX
    pure elemental function conjg_real(x) result(r)
        REAL_TYPE, value :: x
        REAL_TYPE :: r
        r = x
    end
#endif

end module
