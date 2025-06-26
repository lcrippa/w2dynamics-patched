!===============================================================================
module MParameters
!===============================================================================
    use iso_c_binding, only: c_int64_t, c_double, c_long_double
    private

    public :: c_double, c_int64_t

    integer, parameter :: idlength = 100, linelength = 5000
    integer, parameter, public :: vallength = 500

    type :: TParameter
        character(idlength)              :: ID=' '
        character(vallength)             :: Value=' '
        type(TParameter),pointer         :: next
    end type TParameter

    type(TParameter), pointer :: first=>null()
    type(TParameter), pointer :: last=>null()

    public :: get_Integer_Parameter, get_LongInteger_Parameter
    public :: get_Real_Parameter, get_String_Parameter
    public :: get_Integer_List
    public :: read_ParameterString, dest_Parameters

contains

!===============================================================================
integer function get_Integer_Parameter(ID)
!===============================================================================
!input
    character(*),intent(in)             :: ID
!local
    type(TParameter),pointer            :: tmp

    real :: rval
    tmp=>first
    do while(index(tmp%ID,ID).ne.1)
        if(associated(tmp%next).eqv..false.)then
            write(0,*)ID," not found in parameter list!"
            error stop 'Parameter not found'
        endif
        tmp=>tmp%next
    enddo
    read(tmp%value,*,err=98)get_Integer_Parameter
    return

    ! Attempt with real to support scientific notation
98 read(tmp%value,*,err=99) rval
    if( (rval - int(rval)) > 2*epsilon(1.) ) goto 99
    get_Integer_Parameter = int(rval)
    return

99 write(0,*) ID," is not an integer"
    error stop 'Parameter is not an integer'
end function get_Integer_Parameter

!===============================================================================
integer(c_int64_t) function get_LongInteger_Parameter(ID)
!===============================================================================
!input
    character(*),intent(in)             :: ID
!local
    type(TParameter),pointer            :: tmp
    real(kind=c_long_double)            :: rval

    tmp=>first
    do while(index(tmp%ID,ID).ne.1)
        if(associated(tmp%next).eqv..false.)then
            write(0,*)ID," not found in parameter list!"
            error stop 'Parameter not found'
        endif
        tmp=>tmp%next
    enddo
    read(tmp%value,*,err=99)get_LongInteger_Parameter
    return

    ! Attempt with real to support scientific notation
99 read(tmp%value,*) rval
    if( (rval - int(rval, KIND=c_int64_t)) > 2*epsilon(1.) ) then
        write(0,*) ID," is not an integer"
        error stop 'Parameter not found'
    endif
    get_LongInteger_Parameter = int(rval,KIND=c_int64_t)

end function get_LongInteger_Parameter

!===============================================================================
real(c_double) function get_Real_Parameter(ID)
!===============================================================================
!input
    character(*),intent(in)             :: ID
!local
    type(TParameter),pointer            :: tmp

    tmp=>first
    do while(index(tmp%ID,ID).ne.1)
        if(associated(tmp%next).eqv..false.)then
            write(0,*)ID," not found in parameter list!"
            error stop 'Parameter not found'
        endif
        tmp=>tmp%next
    enddo
    read(tmp%value, *, err=99) get_Real_Parameter
    return

99 write (*,*) 'Invalid parameter:', ID, '=', tmp%value
    error stop 'Parameter is of wrong type'
end function get_Real_Parameter

!===============================================================================
subroutine get_Integer_List(ID,Val)
!===============================================================================
!input
    character(*),intent(in)          :: ID
!output
    integer,intent(out)              :: Val(:)
!local
    type(TParameter),pointer         :: tmp

    tmp=>first
    do while(index(tmp%ID,ID).ne.1)
        if(associated(tmp%next).eqv..false.)then
            write(0,*)ID," not found in parameter list!"
            error stop 'Parameter not found'
        endif
        tmp=>tmp%next
    enddo
    read(tmp%value,*)Val
end subroutine get_Integer_List

!===============================================================================
character(vallength) function get_String_Parameter(ID)
!===============================================================================
!input
    character(*),intent(in)     :: ID
!local
    type(TParameter),pointer            :: tmp

    tmp=>first
    do while(index(tmp%ID,ID).ne.1)
        if(associated(tmp%next).eqv..false.)then
            write(0,*)ID," not found in parameter list!"
            error stop 'Parameter not found'
        endif
        tmp=>tmp%next
    enddo
    get_String_Parameter=tmp%value
end function get_String_Parameter

!===============================================================================
subroutine push_Parameter(ID,Value)
!===============================================================================
!input
    character(*),intent(in)     :: ID,Value

    if(associated(first).eqv..false.)then
        allocate(first)
        last=>first
        first%next=>null()
    else
        allocate(last%next)
        last=>last%next
        last%next=>null()
    endif
    last%ID(:)=' '
    last%ID=ID
    last%Value(:)=' '
    last%Value=Value
end subroutine push_Parameter

!===============================================================================
subroutine read_ParameterString(ParameterString)
!===============================================================================
!input
    character(*),intent(in)             :: ParameterString
!local
    integer                             :: PosOld,PosNew
    character(idlength)                 :: id
    character(vallength)                :: val
    character(linelength)               :: line

    first=>null()
    last=>null()
    line(:)=' '
    id(:)=' '
    val(:)=' '
    PosOld=1
    PosNew=index(ParameterString(PosOld:),char(10))
    line=ParameterString(PosOld:PosNew)
    do while(PosNew.ne.0)
        line=ParameterString(PosOld:PosOld+PosNew-1)
        if(index(line,"#").ne.1)then
            do while(index(line,char(09)).ne.0)
                line(index(line,char(09)):index(line,char(09)))=" "
            enddo
            if(index(line,"#").gt.0)line(index(line,"#"):)=" "
            if(index(line,"=").gt.0)then
                id=trim(adjustl(line(:index(line,"=")-1)))
                val=trim(adjustl(line(index(line,"=")+1:)))
                do while(index(val,",").ne.0)
                    val(index(val,","):index(val,","))=" "
                enddo
                call push_Parameter(id,val)
            endif
        endif
        PosOld=PosOld+PosNew
        PosNew=index(ParameterString(PosOld:),char(10))
    end do
end subroutine read_ParameterString

!===============================================================================
subroutine dest_Parameters()
!===============================================================================
!local
    type(TParameter),pointer            :: tmp1,tmp2

    tmp1=>first
    do while(associated(tmp1))
        tmp2=>tmp1
        tmp1=>tmp1%next
        deallocate(tmp2)
    enddo
    first=>null()
end subroutine dest_Parameters
!===============================================================================
end module MParameters
!===============================================================================
