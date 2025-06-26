program test
   use MLibC
   use MWormState
   use MSignal
   use testing
   implicit none

   type(TWormState) :: worm
   integer :: wormpos_loref(10)
   integer, allocatable :: wormpos_input(:)

   call enter_z_sector(worm)

   ! TEST 1 corresponding to adding operators with hybridization
   ! config: from W3 H1 H2 W1 H3 H4    W5 W4 H5 H6    W2
   !         to   W3 H1 H2 W1 H3 H4 H5 W5 W4 H6 H7 H8 W2
   allocate(wormpos_input, source=[4, 11, 1, 8, 7])
   call assign_pos_from_locpos_typeorder(worm, wormpos_input)

   wormpos_loref(1:5) = [1, 4, 7, 8, 11]
   call check_array_equality( &
            get_worm_oper_positions_time_ordered(worm), wormpos_loref, __LINE__)

   ! insert hybops at [7, 12]; FIXME: test also generation of the diff array
   call wormstate_update_pos(worm, [0, 0, 1, 1, 2])

   wormpos_loref(1:5) = [1, 4, 8, 9, 13]
   call check_array_equality( &
            get_worm_oper_positions_time_ordered(worm), wormpos_loref, __LINE__)

   ! TEST 2 corresponding to removing operators with hybridization
   ! config: from W3 H1 H2 W1 H3 H4 H5 W5 W4 H6 H7 H8 W2
   !         to   W3 H1 H2 W1 H3 H4    W5 W4 H5 H6    W2
   wormpos_input = [4, 13, 1, 9, 8]
   call assign_pos_from_locpos_typeorder(worm, wormpos_input)

   call check_array_equality( &
            get_worm_oper_positions_time_ordered(worm), wormpos_loref, __LINE__)

   call wormstate_remove_hybops(worm, [7, 12])

   wormpos_loref(1:5) = [1, 4, 7, 8, 11]
   call check_array_equality( &
            get_worm_oper_positions_time_ordered(worm), wormpos_loref, __LINE__)

   ! TEST 3 corresponding to a cyclical shift of the configuration
   ! config: from W3 H1 H2 W1 H3 H4 H5 W5 W4 H6 H7 H8 W2
   !         to   W4 H6 H7 H8 W2 W3 H1 H2 W1 H3 H4 H5 W5
   wormpos_input = [4, 13, 1, 9, 8]
   call assign_pos_from_locpos_typeorder(worm, wormpos_input)

   wormpos_loref(1:5) = [1, 4, 8, 9, 13]
   call check_array_equality( &
            get_worm_oper_positions_time_ordered(worm), wormpos_loref, __LINE__)


   call wormstate_taushift(worm, 13, 8)

   wormpos_loref(1:5) = [1, 5, 6, 9, 13]
   call check_array_equality( &
            get_worm_oper_positions_time_ordered(worm), wormpos_loref, __LINE__)

   ! Second set of tests with different initial configs
   ! TEST 4 corresponding to adding operators with hybridization
   ! config: from    H1 W8 H2 W2 H3    H4 W5    W3 H5 W1 H6 H7 W7 W4 H8 W6 H9
   !         to   H1 H2 W8 H3 W2 H4 H5 H6 W5 H7 W3 H8 W1 H9 Ha W7 W4 Hb W6 Hc Hd
   wormpos_input = [10, 4, 8, 14, 7, 16, 13, 2]
   call assign_pos_from_locpos_typeorder(worm, wormpos_input)

   wormpos_loref(1:8) = [2, 4, 7, 8, 10, 13, 14, 16]
   call check_array_equality( &
        get_worm_oper_positions_time_ordered(worm), wormpos_loref, __LINE__)

   ! insert hybops at [1, 7, 10, 21]; FIXME: test also generation of the diff array
   call wormstate_update_pos(worm, [1, 1, 2, 3, 3, 3, 3, 3])

   wormpos_loref(1:8) = [3, 5, 9, 11, 13, 16, 17, 19]
   call check_array_equality( &
            get_worm_oper_positions_time_ordered(worm), wormpos_loref, __LINE__)

   ! TEST 5 corresponding to removing operators with hybridization
   ! config: from H1 H2 W8 H3 W2 H4 H5 H6 W5 H7 W3 H8 W1 H9 Ha W7 W4 Hb W6 Hc Hd
   !         to         W8 H3 W2 H4       W5 H7 W3 H8 W1       W7 W4 Hb W6    Hd
   call wormstate_remove_hybops(worm, [1, 2, 7, 8, 14, 15, 20])

   wormpos_loref(1:8) = [1, 3, 5, 7, 9, 10, 11, 13]
   call check_array_equality( &
            get_worm_oper_positions_time_ordered(worm), wormpos_loref, __LINE__)

   ! TEST 6 corresponding to a cyclical shift of the configuration
   ! config: from W8 H3 W2 H4 W5 H7 W3 H8 W1 W7 W4 Hb W6 Hd
   !         to   Hd W8 H3 W2 H4 W5 H7 W3 H8 W1 W7 W4 Hb W6
   call wormstate_taushift(worm, 14, 13)

   wormpos_loref(1:8) = [2, 4, 6, 8, 10, 11, 12, 14]
   call check_array_equality( &
            get_worm_oper_positions_time_ordered(worm), wormpos_loref, __LINE__)

contains

   subroutine check_array_equality(arr1, arr2, line)
      integer, intent(in) :: arr1(:), arr2(:), line
      integer :: i
      logical :: failed

      failed = .false.
      do i = 1, size(arr1)
         if (arr1(i) /= arr2(i)) failed = .true.
      end do

      if (failed) then
         write(*, *) "First array: ", arr1
         write(*, *) "Second array: ", arr2
         call test_error(line)
      end if
   end subroutine check_array_equality

   subroutine test_error(line)
      use iso_c_binding
      integer :: line

      write(*, "('Test error in line ', I4)") line
      call libc_abort()
   end subroutine test_error

end program test
