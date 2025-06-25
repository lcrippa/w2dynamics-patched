!> Module providing C wrappers for the accumulators
module caccumulators
    use MAccumulator
    use iso_c_binding
    implicit none

contains
    function caccumulator_create(iscomplex, typ, ncomp, nlevels) result(this) &
                                        bind(C, name='accumulator_create')
        logical(c_bool), value :: iscomplex
        integer(c_int64_t), value :: typ, ncomp, nlevels
        type(c_ptr) :: this

        type(TAccumulator), pointer :: thisf

        this = c_null_ptr
        if (ncomp < 0 .or. nlevels < 0) &
            return
        if (typ /= ACC_NULL .and. typ /= ACC_MEAN .and. typ /= ACC_BLOCKS) &
            return

        allocate(thisf)
        call accumulator_init(thisf, logical(iscomplex), int(typ), &
                              ncomp, nlevels)
        this = c_loc(thisf)
    end function

    function caccumulator_n(this) result(n) bind(C, name='accumulator_n')
        type(c_ptr), value :: this
        integer(c_int64_t) :: n

        type(TAccumulator), pointer :: thisf

        call c_f_pointer(this, thisf)
        n = accumulator_count(thisf)
    end function

    subroutine caccumulator_add(this) bind(C, name='accumulator_add')
        type(c_ptr), value :: this

        type(TAccumulator), pointer :: thisf

        call c_f_pointer(this, thisf)
        call accumulator_add(thisf)
    end subroutine

    function caccumulator_buffer(this) result(buf) bind(C, name='accumulator_buffer')
        type(c_ptr), value :: this
        type(c_ptr) :: buf

        type(TAccumulator), pointer :: thisf
        real(c_double), pointer :: rbuf(:)
        complex(c_double_complex), pointer :: cbuf(:)

        call c_f_pointer(this, thisf)
        if (accumulator_is_complex(thisf)) then
            call accumulator_buffer(thisf, cbuf)
            buf = c_loc(cbuf)
        else
            call accumulator_buffer(thisf, rbuf)
            buf = c_loc(rbuf)
        endif
    end function

    subroutine caccumulator_reset(this) bind(C, name='accumulator_reset')
        type(c_ptr), value :: this

        type(TAccumulator), pointer :: thisf

        call c_f_pointer(this, thisf)
        call accumulator_reset(thisf)
    end subroutine

    subroutine caccumulator_delete(this) bind(C, name='accumulator_delete')
        type(c_ptr), value :: this

        type(TAccumulator), pointer :: thisf

        call c_f_pointer(this, thisf)
        call accumulator_delete(thisf)
        deallocate(thisf)
    end subroutine

    function caccumulator_mean(this, out_) result(success) &
                                            bind(C, name='accumulator_mean')
        logical(c_bool) :: success
        type(c_ptr), value :: this, out_

        type(TAccumulator), pointer :: thisf
        real(c_double), pointer :: dout(:)
        complex(c_double_complex), pointer :: zout(:)
        integer(c_int64_t) :: ncomp

        call c_f_pointer(this, thisf)
        success = accumulator_has_mean(thisf)
        if (success) then
            ncomp = accumulator_num_comp(thisf)
            if (accumulator_is_complex(thisf)) then
                call c_f_pointer(out_, zout, (/ ncomp /))
                call accumulator_mean(thisf, zout)
            else
                call c_f_pointer(out_, dout, (/ ncomp /))
                call accumulator_mean(thisf, dout)
            endif
        endif
    end function

end module caccumulators
