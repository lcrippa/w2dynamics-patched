#define TEST_ERROR() call terminate(line=__LINE__)

program test
   use MOperator
   use DMOperator
   use ZMOperator
   use MExtRange
   use Testing
   implicit none

   interface set_from_dconst
      procedure :: d_set_from_dconst, z_set_from_dconst
   end interface set_from_dconst

   real(c_double), parameter :: vk3(10) = [ &
          -1044.8542536976215d0, 389.8365509568198d0, 1581.432914594263d0, &
          907.2842517507811d0, -1683.022394551576d0, -1948.0010054324275d0, &
          -392.0350889647234d0, -1473.499789358022d0, -893.5702865971407d0, &
          534.8438106460144d0 ]
   real(c_double), parameter :: vb3(10) = [ &
          -372.0230376288182d0, 1599.8650058490066d0, -221.04067000063446d0, &
          1541.0848067556458d0, 1326.8205389095883d0, -1403.9400550562245d0, &
          -97.33173629675343d0, 1805.1284360605155d0, 1196.6490754790052d0, &
          242.52274025745828d0 ]
   real(c_double), parameter :: vk2(5) = [ &
          35102.33103685331d0, -21183.920194565537d0, 10160.87685952634d0, &
          -59045.55815458095d0, -44650.12774911754d0 ]
   real(c_double), parameter :: vb2(5) = [ &
          74437.73534281671d0, 24396.686280786816d0, -11331.661256534382d0, &
          54471.4331031708d0, -69316.98032080755d0]
   real(c_double), parameter :: vk1(3) = [0.0d0, 0.3565636028997164d0, 1.0d0]
   real(c_double), parameter :: vb1(3) = [1.0d0, 0.0d0, 0.607133050176991d0]

   call d_test_braket()
   call z_test_braket()
   call z_test_braket_complex()

contains

   subroutine d_test_braket()
     type(DTKet) :: k1
     type(DTKet) :: k2
     type(DTKet) :: k3
     type(DTBra) :: b1
     type(DTBra) :: b2
     type(DTBra) :: b3
     type(DExtRange) :: res, cmp

     call set_state(k3, 0, vk3, 7)
     call set_state(b3, 0, vb3, 1)

     call set_from_dconst(cmp, -255598728.4197058d0)
     res = braket_dotprod(b3, k3)
     call assert_close(res, cmp, 1d-14)

     call set_state(k3, 0, vk3(1:7), 7)
     call set_state(b3, 0, vb3(1:7), 1)

     call set_from_dconst(cmp, 665856539.6926042d0)
     res = braket_dotprod(b3, k3)
     call assert_close(res, cmp, 1d-14)

     call set_state(k3, 0, vk3(1:5), 7)
     call set_state(b3, 0, vb3(1:5), 1)

     call set_from_dconst(cmp, -44040188.55432612d0)
     res = braket_dotprod(b3, k3)
     call assert_close(res, cmp, 1d-14)

     call set_state(k2, 0, vk2, 250)
     call set_state(b2, 0, vb2, -17)

     call set_from_dconst(cmp, 2.5670311374140373d+79)
     res = braket_dotprod(b2, k2)
     call assert_close(res, cmp, 1d-13)

     call set_from_dconst(cmp, -7.216494984277578d+83)
     res = braket_dotprod(b3, k2)
     call assert_close(res, cmp, 1d-12)

     call set_from_dconst(cmp, 78024.39176290495d0)
     res = braket_dotprod(b2, k3)
     call assert_close(res, cmp, 1d-13)

     call set_state(k1, 0, vk1, 45)
     call set_state(b1, 0, vb1, 46)

     call set_from_dconst(cmp, 1.503188623975114d+27)
     res = braket_dotprod(b1, k1)
     call assert_close(res, cmp, 1d-14)

     call set_state(k2, 0, vk2(1:3), 250)
     call set_state(b2, 0, vb2(1:3), -17)
     call set_state(k3, 0, vk3(1:3), 7)
     call set_state(b3, 0, vb3(1:3), 1)

     call set_from_dconst(cmp, 5.254449661734455d+93)
     res = braket_dotprod(b1, k2)
     call assert_close(res, cmp, 1d-13)

     call set_from_dconst(cmp, -7.6303645961873d+17)
     res = braket_dotprod(b1, k3)
     call assert_close(res, cmp, 1d-14)

     call set_from_dconst(cmp, -706707581562.1484d0)
     res = braket_dotprod(b2, k1)
     call assert_close(res, cmp, 1d-14)

     call set_from_dconst(cmp, 2.458775123950005d+16)
     res = braket_dotprod(b3, k1)
     call assert_close(res, cmp, 1d-14)
   end subroutine

   subroutine z_test_braket()
      type(ZTKet) :: k1
      type(ZTKet) :: k2
      type(ZTKet) :: k3
      type(ZTBra) :: b1
      type(ZTBra) :: b2
      type(ZTBra) :: b3
      type(ZExtRange) :: res, cmp

      call set_state(k3, 0, vk3(1:7), 7)
      call set_state(b3, 0, vb3(1:7), 1)

      call set_from_dconst(cmp, 665856539.6926042d0)
      res = braket_dotprod(b3, k3)
      call assert_close(res, cmp, 1d-14)

      call set_state(k3, 0, vk3, 7)
      call set_state(b3, 0, vb3, 1)

      call set_from_dconst(cmp, -255598728.4197058d0)
      res = braket_dotprod(b3, k3)
      call assert_close(res, cmp, 1d-14)

      call set_state(k3, 0, vk3(1:5), 7)
      call set_state(b3, 0, vb3(1:5), 1)

      call set_from_dconst(cmp, -44040188.55432612d0)
      res = braket_dotprod(b3, k3)
      call assert_close(res, cmp, 1d-14)

     call set_state(k2, 0, vk2, 250)
     call set_state(b2, 0, vb2, -17)

      call set_from_dconst(cmp, 2.5670311374140373d+79)
      res = braket_dotprod(b2, k2)
      call assert_close(res, cmp, 1d-13)

      call set_from_dconst(cmp, -7.216494984277578d+83)
      res = braket_dotprod(b3, k2)
      call assert_close(res, cmp, 1d-13)

      call set_from_dconst(cmp, 78024.39176290495d0)
      res = braket_dotprod(b2, k3)
      call assert_close(res, cmp, 1d-13)

      call set_state(k1, 0, vk1, 45)
      call set_state(b1, 0, vb1, 46)

      call set_from_dconst(cmp, 1.503188623975114d+27)
      res = braket_dotprod(b1, k1)
      call assert_close(res, cmp, 1d-14)

      call set_state(k2, 0, vk2(1:3), 250)
      call set_state(b2, 0, vb2(1:3), -17)
      call set_state(k3, 0, vk3(1:3), 7)
      call set_state(b3, 0, vb3(1:3), 1)

      call set_from_dconst(cmp, 5.254449661734455d+93)
      res = braket_dotprod(b1, k2)
      call assert_close(res, cmp, 1d-13)

      call set_from_dconst(cmp, -7.6303645961873d+17)
      res = braket_dotprod(b1, k3)
      call assert_close(res, cmp, 1d-14)

      call set_from_dconst(cmp, -706707581562.1484d0)
      res = braket_dotprod(b2, k1)
      call assert_close(res, cmp, 1d-14)

      call set_from_dconst(cmp, 2.458775123950005d+16)
      res = braket_dotprod(b3, k1)
      call assert_close(res, cmp, 1d-14)
   end subroutine

   pure subroutine d_set_from_dconst(e, v)
      type(DExtRange), intent(out) :: e
      real(c_double), intent(in)   :: v

      e = extrange(v)
   end subroutine d_set_from_dconst

   pure subroutine z_set_from_dconst(e, v)
      type(ZExtRange), intent(out) :: e
      real(c_double), intent(in)   :: v

       e = extrange(cmplx(v, kind=c_double_complex))
   end subroutine z_set_from_dconst

   subroutine z_test_braket_complex()
     type(ZTKet) :: k3
     type(ZTBra) :: b3
     type(ZExtRange) :: res, cmp

     call set_state(k3, 0, [&
            cmplx(-1044.8542536976215d0, -386.69842135d0, c_double_complex),&
            cmplx(389.8365509568198d0, -933.33064604d0, c_double_complex),&
            cmplx(1581.432914594263d0, 968.16396594d0, c_double_complex),&
            cmplx(907.2842517507811d0, -737.26991453d0, c_double_complex),&
            cmplx(-1683.022394551576d0, -500.51052071d0, c_double_complex)], &
            27)
     call set_state(b3, 0, [&
            cmplx(-372.0230376288182d0, -532.99219674d0, c_double_complex),&
            cmplx(1599.8650058490066d0, -955.64465096d0, c_double_complex),&
            cmplx(-221.04067000063446d0, 35.28056198d0, c_double_complex),&
            cmplx(1541.0848067556458d0, -993.39865066d0, c_double_complex),&
            cmplx(1326.8205389095883d0, 166.46480745d0, c_double_complex)], &
            -19)

     cmp = extrange(cmplx(411968157.0848851d0, -620113375.5784745d0, c_double_complex))
     res = braket_dotprod(b3, k3)
     call assert_close(res, cmp, rtol=1d-14)
   end subroutine z_test_braket_complex

end program test
