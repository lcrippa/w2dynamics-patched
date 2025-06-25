module MAccumulator
    use iso_c_binding, only: c_int64_t, c_double, c_double_complex
    use MAccumulatorBase
    use MAccumulatorD
    use MAccumulatorZ
    private

    ! Re-exports from the typed modules
    public DAccumulator, ZAccumulator

    ! Re-exports from the base module
    public ACC_NULL, ACC_MEAN, ACC_BLOCKS

    type, public :: TAccumulator
        private
        type(ZAccumulator), allocatable :: zacc
        type(DAccumulator), allocatable :: dacc
    end type

    interface accumulator_init
        module procedure p_accumulator_init
    end interface
    public accumulator_init

    interface accumulator_add
        module procedure p_accumulator_add
    end interface
    public accumulator_add

    interface accumulator_buffer
       module procedure p_accumulator_dbuffer, p_accumulator_zbuffer
    end interface
    public accumulator_buffer

    interface accumulator_reset
        module procedure p_accumulator_reset
    end interface
    public accumulator_reset

    interface accumulator_delete
        module procedure p_accumulator_delete
    end interface
    public accumulator_delete

    interface accumulator_num_comp
        module procedure p_accumulator_num_comp
    end interface
    public accumulator_num_comp

    interface accumulator_num_blocks
        module procedure p_accumulator_num_blocks
    end interface
    public accumulator_num_blocks

    interface accumulator_has_mean
        module procedure p_accumulator_has_mean
    end interface
    public accumulator_has_mean

    interface accumulator_has_blocks
        module procedure p_accumulator_has_blocks
    end interface
    public accumulator_has_blocks

    interface accumulator_count
       module procedure p_accumulator_count
    end interface
    public accumulator_count

    interface accumulator_mean
        module procedure p_accumulator_dmean, p_accumulator_zmean
    end interface
    public accumulator_mean

    interface accumulator_blocks
        module procedure p_accumulator_blocks
    end interface
    public accumulator_blocks

    interface accumulator_is_complex
       module procedure p_accumulator_is_complex
    end interface
    public accumulator_is_complex

contains

    subroutine p_accumulator_init(self, iscomplex, typ, ncomp, nlevels)
        type(TAccumulator), intent(inout) :: self
        logical, intent(in) :: iscomplex
        integer, intent(in) :: typ
        integer(c_int64_t), intent(in) :: ncomp
        integer(c_int64_t), intent(in), optional :: nlevels

        call accumulator_delete(self)
        if (iscomplex) then
            allocate(self%zacc)
            call accumulator_init(self%zacc, typ, ncomp, nlevels)
        else
            allocate(self%dacc)
            call accumulator_init(self%dacc, typ, ncomp, nlevels)
        endif
    end

    subroutine p_accumulator_reset(self)
        type(TAccumulator), intent(inout) :: self

        if (accumulator_is_complex(self)) then
            call accumulator_reset(self%zacc)
        else
            call accumulator_reset(self%dacc)
        endif
    end subroutine

    subroutine p_accumulator_delete(self)
        type(TAccumulator), intent(inout) :: self

        if (allocated(self%zacc)) &
            deallocate(self%zacc)
        if (allocated(self%dacc)) &
            deallocate(self%dacc)
    end subroutine

    subroutine p_accumulator_add(self)
        type(TAccumulator), intent(inout) :: self

        if (accumulator_is_complex(self)) then
            call accumulator_add(self%zacc)
        else
            call accumulator_add(self%dacc)
        endif
    end subroutine

    subroutine p_accumulator_dbuffer(this, buf)
        type(TAccumulator), intent(inout), target :: this
        real(c_double), intent(out), pointer :: buf(:)

        if (accumulator_is_complex(this)) then
            call backtrace
            error stop 'Want real, but accumulator is complex'
        else
            call accumulator_buffer(this%dacc, buf)
        endif
    end subroutine

    subroutine p_accumulator_zbuffer(this, buf)
        type(TAccumulator), intent(inout), target :: this
        complex(c_double_complex), intent(out), pointer :: buf(:)

        if (accumulator_is_complex(this)) then
            call accumulator_buffer(this%zacc, buf)
        else
            call backtrace
            error stop 'Want complex, but accumulator is real'
        endif
    end subroutine

    function p_accumulator_num_comp(self) result(answer)
        type(TAccumulator), intent(in) :: self
        integer(c_int64_t) :: answer

        if (accumulator_is_complex(self)) then
            answer = accumulator_num_comp(self%zacc)
        else
            answer = accumulator_num_comp(self%dacc)
        endif
    end function

    function p_accumulator_num_blocks(self) result(answer)
        type(TAccumulator), intent(in) :: self
        integer(c_int64_t) :: answer

        if (accumulator_is_complex(self)) then
            answer = accumulator_num_blocks(self%zacc)
        else
            answer = accumulator_num_blocks(self%dacc)
        endif
    end function

    function p_accumulator_has_mean(self) result(answer)
        type(TAccumulator), intent(in) :: self
        logical :: answer

        if (accumulator_is_complex(self)) then
            answer = accumulator_has_mean(self%zacc)
        else
            answer = accumulator_has_mean(self%dacc)
        endif
    end function

    function p_accumulator_has_blocks(self) result(answer)
        type(TAccumulator), intent(in) :: self
        logical :: answer

        if (accumulator_is_complex(self)) then
            answer = accumulator_has_blocks(self%zacc)
        else
            answer = accumulator_has_blocks(self%dacc)
        endif
    end function

    function p_accumulator_count(self) result(n)
        type(TAccumulator), intent(in) :: self
        integer(c_int64_t) :: n

        if (accumulator_is_complex(self)) then
            n = accumulator_count(self%zacc)
        else
            n = accumulator_count(self%dacc)
        endif
    end function

    subroutine p_accumulator_dmean(self, out_)
        type(TAccumulator), intent(in) :: self
        real(c_double), intent(out) :: out_(:)

        if (accumulator_is_complex(self)) then
            error stop 'Expecting double mean, but accumulator is complex'
        else
            call accumulator_mean(self%dacc, out_)
        endif
    end subroutine

    subroutine p_accumulator_zmean(self, out_)
        type(TAccumulator), intent(in) :: self
        complex(c_double_complex), intent(out) :: out_(:)

        if (accumulator_is_complex(self)) then
            call accumulator_mean(self%zacc, out_)
        else
            error stop 'Expecting complex mean, but accumulator is real'
        endif
    end subroutine

    subroutine p_accumulator_blocks(self, out_)
        type(TAccumulator), intent(in) :: self
        real(c_double), intent(out) :: out_(:,:)

        if (accumulator_is_complex(self)) then
            call accumulator_blocks(self%zacc, out_)
        else
            call accumulator_blocks(self%dacc, out_)
        endif
    end subroutine

    logical function p_accumulator_is_complex(self) result(test)
        type(TAccumulator), intent(in) :: self

        if (allocated(self%zacc) .eqv. allocated(self%dacc)) &
            stop 'Invalid allocator'

        test = allocated(self%zacc)
    end function

end module
