module MSamplerRegister
    use iso_c_binding, only: c_double_complex
    use MAccumulator
    use MSamplerD
    use MSamplerZ
    private

    type, public :: sampler_register_entry
        class(ZSampler), pointer :: p => null()
        type(TAccumulator), pointer :: acc => null()
    end type sampler_register_entry

    type, public :: sampler_register
        type(sampler_register_entry), allocatable :: ps(:)
        logical                                   :: empty = .true.
    contains
        procedure, pass(reg) :: register => sampler_register_register
        procedure, pass(reg) :: add => sampler_register_add
        procedure, pass(reg) :: clear => sampler_register_clear
    end type sampler_register

    type, public :: sampler2d_register_entry
        class(ZSampler2D), pointer :: p => null()
        type(TAccumulator), pointer :: acc => null()
    end type sampler2d_register_entry

    type, public :: sampler2d_register
        type(sampler2d_register_entry), allocatable :: ps(:)
        logical                                   :: empty = .true.
    contains
        procedure, pass(reg) :: register => sampler2d_register_register
        procedure, pass(reg) :: add => sampler2d_register_add
        procedure, pass(reg) :: clear => sampler2d_register_clear
    end type sampler2d_register

    type, public :: sampler3d_register_entry
        class(ZSampler3D), pointer :: p => null()
        type(TAccumulator), pointer :: acc => null()
     end type sampler3d_register_entry

    type, public :: sampler3d_register
        type(sampler3d_register_entry), allocatable :: ps(:)
        logical                                   :: empty = .true.
    contains
        procedure, pass(reg) :: register => sampler3d_register_register
        procedure, pass(reg) :: add => sampler3d_register_add
        procedure, pass(reg) :: clear => sampler3d_register_clear
    end type sampler3d_register

contains

    subroutine sampler_register_register(reg, samp, acc)
        class(sampler_register), intent(inout) :: reg
        class(ZSampler), target, intent(in)     :: samp
        type(TAccumulator), target, intent(in) :: acc
        type(sampler_register_entry), allocatable      :: temp(:)

        if (allocated(reg%ps)) then
            allocate(temp(size(reg%ps) + 1))
            temp(1 : size(temp) - 1) = reg%ps(:)
            deallocate(reg%ps)
        else
            allocate(temp(1))
        end if
        temp(size(temp))%p => samp
        temp(size(temp))%acc => acc
        call move_alloc(temp, reg%ps)
        reg%empty = .false.
    end subroutine sampler_register_register

    subroutine sampler_register_add(reg, sample)
        class(sampler_register), intent(inout) :: reg
        type(ZSamplePoint), intent(inout)       :: sample(:)
        complex(c_double_complex), pointer                :: ibuf(:), buf(:, :)
        integer :: i

        do i = 1, size(reg%ps)
            call accumulator_buffer(reg%ps(i)%acc, ibuf)
            buf(1:reg%ps(i)%p%nacc, 1:reg%ps(i)%p%ncomp) => ibuf
            call sampler_add(reg%ps(i)%p, sample, buf)
            call accumulator_add(reg%ps(i)%acc)
        end do
    end subroutine sampler_register_add

    subroutine sampler_register_clear(reg)
        class(sampler_register), intent(inout) :: reg

        if (allocated(reg%ps)) deallocate(reg%ps)
        reg%empty = .true.
    end subroutine sampler_register_clear

    subroutine sampler2d_register_register(reg, samp, acc)
        class(sampler2d_register), intent(inout) :: reg
        class(ZSampler2D), target, intent(in)     :: samp
        type(TAccumulator), target, intent(in) :: acc
        type(sampler2d_register_entry), allocatable      :: temp(:)

        if (allocated(reg%ps)) then
            allocate(temp(size(reg%ps) + 1))
            temp(1 : size(temp) - 1) = reg%ps(:)
            deallocate(reg%ps)
        else
            allocate(temp(1))
        end if
        temp(size(temp))%p => samp
        temp(size(temp))%acc => acc
        call move_alloc(temp, reg%ps)
        reg%empty = .false.
    end subroutine sampler2d_register_register

    subroutine sampler2d_register_add(reg, sample)
        class(sampler2d_register), intent(inout) :: reg
        type(ZSamplePoint2D), intent(inout)       :: sample(:)
        complex(c_double_complex), pointer                :: ibuf(:), buf(:, :)
        integer :: i

        do i = 1, size(reg%ps)
            call accumulator_buffer(reg%ps(i)%acc, ibuf)
            buf(1:reg%ps(i)%p%nacc, 1:reg%ps(i)%p%ncomp) => ibuf
            call sampler_add(reg%ps(i)%p, sample, buf)
            call accumulator_add(reg%ps(i)%acc)
        end do
    end subroutine sampler2d_register_add

    subroutine sampler2d_register_clear(reg)
        class(sampler2d_register), intent(inout) :: reg

        if (allocated(reg%ps)) deallocate(reg%ps)
        reg%empty = .true.
    end subroutine sampler2d_register_clear

    subroutine sampler3d_register_register(reg, samp, acc)
        class(sampler3d_register), intent(inout) :: reg
        class(ZSampler3D), target, intent(in)     :: samp
        type(TAccumulator), target, intent(in) :: acc
        type(sampler3d_register_entry), allocatable      :: temp(:)

        if (allocated(reg%ps)) then
            allocate(temp(size(reg%ps) + 1))
            temp(1 : size(temp) - 1) = reg%ps(:)
            deallocate(reg%ps)
        else
            allocate(temp(1))
        end if
        temp(size(temp))%p => samp
        temp(size(temp))%acc => acc
        call move_alloc(temp, reg%ps)
        reg%empty = .false.
    end subroutine sampler3d_register_register

    subroutine sampler3d_register_add(reg, sample)
        class(sampler3d_register), intent(inout) :: reg
        type(ZSamplePoint3D), intent(inout)       :: sample(:)
        complex(c_double_complex), pointer                :: ibuf(:), buf(:, :)
        integer :: i

        do i = 1, size(reg%ps)
            call accumulator_buffer(reg%ps(i)%acc, ibuf)
            buf(1:reg%ps(i)%p%nacc, 1:reg%ps(i)%p%ncomp) => ibuf
            call sampler_add(reg%ps(i)%p, sample, buf)
            call accumulator_add(reg%ps(i)%acc)
        end do
    end subroutine sampler3d_register_add

    subroutine sampler3d_register_clear(reg)
        class(sampler3d_register), intent(inout) :: reg

        if (allocated(reg%ps)) deallocate(reg%ps)
        reg%empty = .true.
    end subroutine sampler3d_register_clear
end module
