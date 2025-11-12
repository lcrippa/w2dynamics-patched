
!> This is the
!! module for generating all possible configurations for a given number of
!! bands and electrons - it is possible to sort the states by N_t (total
!! number of electrons) and Sz_t (total spin in z-direction)
!! this module is intended to be used with the module Hloc to create a
!! Hamiltonian and creation, annihilation matrices suitable for a
!! continuous time QMC
module MStates
    use iso_c_binding, only: c_double
    use MParameters
    implicit none
    private

    type, public :: TQNSpec
        logical :: Nt = .false.
        logical :: Szt = .false.
        logical :: Qzt = .false.
        logical :: Azt = .false.
        logical :: Lzt = .false.
        logical :: Jzt = .false.
        logical :: All = .false.
    end type TQNSpec

    type :: TSubStates
        !> Number of states in this superstate (block)
        integer                          :: NStates,Offset

        !> List of states in the occupation number basis in this superstate
        integer, allocatable             :: States(:)
    end type TSubStates

    !> \brief Structure holding the system´s occupancy number basis
    !!
    !! If there are good quantum numbers in this basis (the Hamiltonian takes a
    !! block-diagonal form), then the basis states of the disconnected subspaces
    !! are grouped together using \c TSubStates.
    type, public :: TStates
        !> the number of orbitals
        integer                          :: NBands
        !> the total number of states in the fock space
        integer                          :: NStates
        !> the number of blocks in the hamiltonian
        integer                          :: NSStates
        !> the number of states in the biggest block
        integer                          :: NStatesMax
        !> conserved quantum numbers
        type(TQNSpec)                    :: qnspec

        !> this lookup array has three elements:
        !! (:, 1) tell you which index belongs to which state
        !! (:, 2) it tells you which state (index) belongs to which block value
        !! (:, 3) gives a integer-representation of a ket ( 0 0 ... 0 0 1 0 0 .. 0 0)
        integer, allocatable             :: StatesSStates(:,:)
        Type(TSubStates), allocatable    :: SubStates(:)
    end type TStates

    public :: init_states, dest_States
    public :: qnspec_from_parameters
    public :: sst_for_flipped_spins

    interface get_nbands
       module procedure p_getnbands
    end interface
    public :: get_nbands

    interface get_num_superstates
       module procedure p_get_num_superstates
    end interface
    public get_num_superstates

contains
!===============================================================================
subroutine init_States(this, nbands, spec, states2substates)
!===============================================================================
    type(TStates), intent(out)          :: this
    integer, intent(in)                 :: nbands
    type(TQNSpec), intent(in)           :: spec
    integer, intent(in), optional       :: states2substates(:)

    !local
    integer                             :: iSt, nsst
    integer, allocatable                :: states2substates_from_qn(:)

    this%NBands = nbands
    this%NStates = 2**(2*this%NBands)
    this%qnspec = spec

    allocate(this%StatesSStates(0:this%NStates-1, 3))
    this%StatesSStates(:, :) = -1
    do iSt = 0, this%NStates-1
        this%StatesSStates(iSt, 1) = iSt
    enddo

    if (present(states2substates)) then
        nsst = maxval(states2substates) + 1
        call init_substates(this, nsst, states2substates)
    else
        call qns2substates(this, nsst, states2substates_from_qn)
        call init_substates(this, nsst, states2substates_from_qn)
    endif
end subroutine init_States

pure integer function p_getnbands(this) result(nbands)
    type(TStates), intent(in) :: this

    nbands = this%NBands
end

pure integer function p_get_num_superstates(this) result(nsst)
    type(TStates), intent(in) :: this

    nsst = this%NSStates
end

!> this subroutine initializes the substates via an array which tells
!! the subroutine which state (index) belongs to which substate (value)
subroutine init_substates(this, nsubstates, states2substates)
    type(TStates)                       :: this
    integer, intent(in)                 :: nsubstates, states2substates(:)
    !> loop variables
    integer                             :: st1, st2, isst, n
    !> index of the state in a block
    integer                             :: blockidx
    !> saves the states belonging to one substate
    integer                             :: temp(0:this%NStates)

    this%NSStates = nsubstates
    this%StatesSStates(:,2) = states2substates(:)

    allocate(this%SubStates(0:this%NSStates-1))
    this%SubStates(0)%Offset=0
    do st1 = 0, this%NSStates-1
        blockidx = 0
        do st2 = 0, this%NStates-1
            if(this%StatesSStates(this%statessstates(st2,1),2).eq.st1)then
                temp(blockidx) = this%StatesSStates(st2,1)
                blockidx = blockidx + 1
            endif
        enddo
        this%SubStates(st1)%NStates = blockidx
        if(st1.gt.0)&
            this%SubStates(st1)%Offset=this%SubStates(st1-1)%Offset+this%SubStates(st1-1)%NStates
        allocate(this%SubStates(st1)%States(0:blockidx-1))
        this%SubStates(st1)%States(0:blockidx-1)=temp(0:blockidx-1)
    enddo

        !StatesSStates(:,3) gives a integer-representation of a ket ( 0 0 ... 0 0 1 0 0 .. 0 0)
    n=0
    do isst=0,this%NSStates-1
!          write(unit=*,fmt=*) "--> isst", isst
        do st1=0,this%substates(isst)%NStates-1
!        write(unit=*,fmt=*) "n", n
!        write(unit=*,fmt=*) "st1", st1
!        write(unit=*,fmt=*) "zustand", this%substates(isst)%States(st1)
            this%StatesSStates(this%substates(isst)%States(st1),3)=n
            n=n+1
        enddo
    enddo

end subroutine

!> Calculates the number of superstates and a lookup array telling you
!> which state (index) belongs to which block from the specification
!> of conserved quantum numbers.
pure subroutine qns2substates(this, nsubstates, states2substates)
    type(TStates), intent(in)      :: this
    !> the number of blocks in the hamiltonian = number of substates
    integer, intent(out)              :: nsubstates
    !> saves the states belonging to one substate
    integer, allocatable, intent(out) :: states2substates(:)
    !> loop variables
    integer                           :: st1, st2
    integer                           :: hst1, hst2

    allocate(states2substates(0 : this%nstates-1))
    states2substates(:) = -1
    nsubstates = 0

    outer: do st1 = 0, this%nstates-1
        hst1 = this%StatesSStates(st1,1)
        if(states2substates(hst1).ne.-1)cycle outer
        nsubstates = nsubstates + 1
        states2substates(hst1) = nsubstates - 1
        inner: do st2 = st1, this%nstates-1
            hst2 = this%StatesSStates(st2,1)
            if(states2substates(hst2).ne.-1)cycle inner
            if(this%qnspec%Nt .and.(get_Nt(this,hst1) .ne. get_Nt(this,hst2)))cycle inner
            if(this%qnspec%Szt .and..not.eq(get_Szt(this,hst1),get_Szt(this,hst2)))cycle inner
            if(this%qnspec%Qzt .and..not.eq(get_Qzt(this,hst1),get_Qzt(this,hst2)))cycle inner
            if(this%qnspec%Azt .and..not.eq(get_Azt(this,hst1),get_Azt(this,hst2)))cycle inner
            if(this%qnspec%Lzt .and..not.eq(get_Lzt(this,hst1),get_Lzt(this,hst2)))cycle inner
            if(this%qnspec%Jzt .and..not.eq(get_Jzt(this,hst1),get_Jzt(this,hst2)))cycle inner
            states2substates(hst2) = nsubstates-1
        enddo inner
    enddo outer
end subroutine qns2substates

type(TQNSpec) function qnspec_from_parameters() result(spec)
    character(vallength)                :: QNs

    QNs = get_String_Parameter("QuantumNumbers")
    if (index(QNs, "Nt") /= 0) spec%Nt = .true.
    if (index(QNs, "Szt") /= 0) spec%Szt = .true.
    if (index(QNs, "Qzt") /= 0) spec%Qzt = .true.
    if (index(QNs, "Azt") /= 0) spec%Azt = .true.
    if (index(QNs, "Lzt") /= 0) spec%Lzt = .true.
    if (index(QNs, "Jzt") /= 0) spec%Jzt = .true.
    if (index(QNs, "All") /= 0) spec%All = .true.
end function qnspec_from_parameters

logical elemental function eq(real1,real2)
    real(c_double), intent(in) :: real1,real2

    ! XXX arbitrary cutoffs, absolute no less
    eq = abs(real1 - real2) .lt. 1d-12
end function eq

!===============================================================================
pure integer function get_Nt(this, State)
!===============================================================================
    type(TStates), intent(in)           :: this
!input
    integer,intent(in)                  :: State
!local
    integer                             :: iB

    get_Nt = 0
    do iB = 0, this%NBands-1
        if(btest(State,iB              )) get_Nt=get_Nt+1
        if(btest(State,iB + this%NBands)) get_Nt=get_Nt+1
    enddo
end function get_Nt

!===============================================================================
pure real(c_double) function get_Szt(this,State)
!===============================================================================
    type(TStates), intent(in)           :: this
!input
    integer,intent(in)                  :: State
!local
    integer                             :: iB
    get_Szt=0d0
    do iB=0,this%NBands-1
        if(btest(State,iB))get_Szt=get_Szt-0.5d0
        if(btest(State,iB+this%NBands))get_Szt=get_Szt+0.5d0
    enddo
end function get_Szt

!===============================================================================
pure real(c_double) function get_Qzt(this,State)
!===============================================================================
    type(TStates), intent(in)           :: this
!input
    integer,intent(in)                  :: State
!local
    integer                             :: iB
    get_Qzt=0d0
    do iB=0,this%NBands-1
        if((btest(State,iB).and..not.btest(State,iB+this%NBands))&
            .or.(.not.btest(State,iB).and.btest(State,iB+this%NBands)))get_Qzt=get_Qzt+10**iB
    enddo
end function get_Qzt

!===============================================================================
pure real(c_double) function get_Azt(this,State)
!===============================================================================
    type(TStates), intent(in)           :: this
!input
    integer,intent(in)                  :: State
!local
    integer                             :: iB
    get_Azt=0d0
    do iB=0,this%NBands-1
        if(btest(State,iB))get_Azt=get_Azt+10d0**(iB)
        if(btest(State,iB+this%NBands))get_Azt=get_Azt+10d0**(iB+this%NBands)
    enddo
end function get_Azt

!===============================================================================
pure integer function get_eg_up(State)
!===============================================================================
!input
    integer,intent(in)                  :: State
!local
    integer                             :: iB
    get_eg_up=0
    do iB=3,4
        if(btest(State,iB))get_eg_up=get_eg_up+1
    enddo
end function get_eg_up

!===============================================================================
pure integer function get_eg_do(this,State)
!===============================================================================
    type(TStates), intent(in)           :: this
!input
    integer,intent(in)                  :: State
!local
    integer                             :: iB
    get_eg_do=0
    do iB=3,4
        if(btest(State,iB+this%NBands))get_eg_do=get_eg_do+1
    enddo
end function get_eg_do

!===============================================================================
pure integer function get_t2g_up(State)
!===============================================================================
!input
    integer,intent(in)                  :: State
!local
    integer                             :: iB
    get_t2g_up=0
    do iB=0,2
        if(btest(State,iB))get_t2g_up=get_t2g_up+1
    enddo
end function get_t2g_up

!===============================================================================
pure integer function get_t2g_do(this,State)
!===============================================================================
    type(TStates), intent(in)           :: this
!input
    integer,intent(in)                  :: State
!local
    integer                             :: iB
    get_t2g_do=0
    do iB=0,2
        if(btest(State,iB+this%NBands))get_t2g_do=get_t2g_do+1
    enddo
end function get_t2g_do

!===============================================================================
pure real(c_double) function get_Lzt(this,State)
!===============================================================================
    type(TStates), intent(in)           :: this
!input
    integer,intent(in)                  :: State
!local
    integer                             :: iB,iS
    get_Lzt=0d0
    do iS=1,2
        do iB=1,this%NBands
            if(btest(State,(iB-1)+this%NBands*(iS-1)))get_Lzt=get_Lzt+dble(iB)-dble(this%NBands+1)/2d0
        enddo
    enddo
end function get_Lzt

!===============================================================================
pure real(c_double) function get_Jzt(this,State)
!===============================================================================
    type(TStates), intent(in)           :: this
!input
    integer,intent(in)                  :: State
    get_Jzt=get_Lzt(this,State)+get_Szt(this,State)
end function get_Jzt

!===============================================================================
integer pure function get_occ(this,state,b,sp)
!===============================================================================
    type(TStates),intent(in)            :: this
!input
    integer,intent(in)                  :: state,b,sp
    get_occ=ibits(state,(b-1)+this%nbands*(sp-1),1)
end function get_occ

!===============================================================================
integer function sst_for_flipped_spins(this, sst) result(sst_flipped)
    type(TStates), intent(in) :: this
    integer, intent(in) :: sst
    logical, save :: warn = .true.

    integer :: old_state, spin_up, spin_down, new_state

    if (sst < 0 .or. sst >= this%NSStates) &
        error stop 'Invalid superstate'
    if (.not. this%qnspec%Szt .and. warn) then
        write (0,*) 'WARNING: spin is not conserved; unsure of spin flip'
        warn = .false.
    endif

    ! get first state in the current outer superstate and flip all
    ! spins
    old_state = this%SubStates(sst)%States(0)
    spin_down = ibits(old_state, 0, this%NBands)
    spin_up = ibits(old_state, this%NBands, this%NBands)
    new_state = ior(spin_up, ishft(spin_down, this%NBands))

    ! the proposed outer superstate is the one containing that state
    sst_flipped = this%StatesSStates(new_state, 2)
end function

!===============================================================================
subroutine dest_States(this)
!===============================================================================
    type(TStates)                       :: this
    integer                             :: iSS

    if(allocated(this%SubStates))then
        do iSS=0,this%NSStates-1
            if(allocated(this%SubStates(iSS)%States))&
                deallocate(this%SubStates(iSS)%States)
        enddo
    endif
    if(allocated(this%StatesSStates))&
        deallocate(this%StatesSStates)
    if(allocated(this%SubStates))&
        deallocate(this%SubStates)
end subroutine dest_States

!===============================================================================
end module MStates
!===============================================================================
