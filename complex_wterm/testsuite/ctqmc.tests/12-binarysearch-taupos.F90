#define TEST_ERROR() call terminate(line=__LINE__)

program test
   use MCommon
   use MDiagram
   use testing
   implicit none

   type(TBaseOper) :: a(9)
   type(TLocalOper), allocatable :: b(:)

   a(1) = TBaseOper(1.00d0, 1, 1)
   a(2) = TBaseOper(1.15d0, 1, 2)
   a(3) = TBaseOper(1.15d0, 1, 1)
   a(4) = TBaseOper(2.00d0, 2, 1)
   a(5) = TBaseOper(3.00d0, 1, 2)
   a(6) = TBaseOper(3.00d0, 2, 1)
   a(7) = TBaseOper(4.00d0, 1, 1)
   a(8) = TBaseOper(4.50d0, 2, 2)
   a(9) = TBaseOper(5.00d0, 1, 1)

   if (binarysearch_taupos(a, 0.5d0) /= 1) then
      TEST_ERROR()
   end if

   if (binarysearch_taupos(a, 1.0d0) /= 1) then
      TEST_ERROR()
   end if

   if (binarysearch_taupos(a, 1.5d0) /= 4) then
      TEST_ERROR()
   end if

   if (binarysearch_taupos(a, 3.0d0) /= 5) then
      TEST_ERROR()
   end if

   if (binarysearch_taupos(a, 5.0d0) /= 9) then
      TEST_ERROR()
   end if

   if (binarysearch_taupos(a, 6.0d0) /= 10) then
      TEST_ERROR()
   end if

   if (binarysearch_taupos(a(3:9), 0.5d0) + 2 /= 3) then
      TEST_ERROR()
   end if

   if (binarysearch_taupos(a(3:9), 1.15d0) + 2 /= 3) then
      TEST_ERROR()
   end if

   if (binarysearch_taupos(a(3:9), 3.5d0) + 2 /= 7) then
      TEST_ERROR()
   end if

   if (binarysearch_taupos(a(1:5), 3.0d0) /= 5) then
      TEST_ERROR()
   end if

   if (binarysearch_taupos(a(1:5), 3.1d0) /= 6) then
      TEST_ERROR()
   end if

   allocate(b(9))
   b(1) = TLocalOper(1.00d0, 1, 1, 1)
   b(2) = TLocalOper(1.15d0, 1, 2, 1)
   b(3) = TLocalOper(1.15d0, 1, 1, 1)
   b(4) = TLocalOper(2.00d0, 2, 1, 1)
   b(5) = TLocalOper(3.00d0, 1, 2, 1)
   b(6) = TLocalOper(3.00d0, 2, 1, 1)
   b(7) = TLocalOper(4.00d0, 1, 1, 1)
   b(8) = TLocalOper(4.50d0, 2, 2, 1)
   b(9) = TLocalOper(5.00d0, 1, 1, 1)

   if (binarysearch_taupos(b, 0.5d0) /= 1) then
      TEST_ERROR()
   end if

   if (binarysearch_taupos(b, 1.0d0) /= 1) then
      TEST_ERROR()
   end if

   if (binarysearch_taupos(b, 1.5d0) /= 4) then
      TEST_ERROR()
   end if

   if (binarysearch_taupos(b, 3.0d0) /= 5) then
      TEST_ERROR()
   end if

   if (binarysearch_taupos(b, 5.0d0) /= 9) then
      TEST_ERROR()
   end if

   if (binarysearch_taupos(b, 6.0d0) /= 10) then
      TEST_ERROR()
   end if

   if (binarysearch_taupos(b(3:9), 0.5d0) + 2 /= 3) then
      TEST_ERROR()
   end if

   if (binarysearch_taupos(b(3:9), 1.15d0) + 2 /= 3) then
      TEST_ERROR()
   end if

   if (binarysearch_taupos(b(3:9), 3.5d0) + 2 /= 7) then
      TEST_ERROR()
   end if

   if (binarysearch_taupos(b(1:5), 3.0d0) /= 5) then
      TEST_ERROR()
   end if

   if (binarysearch_taupos(b(1:5), 3.1d0) /= 6) then
      TEST_ERROR()
   end if

   deallocate(b)
end program test
