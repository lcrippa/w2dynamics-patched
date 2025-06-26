module MWormMeas
   use iso_c_binding, only: c_double, c_double_complex
   use MDiagram
   use MDiagramMeasurements
   use MImpurity
   use MWormState
   implicit none

   type :: TWormSample
      real(c_double), allocatable :: taus(:)
   end type TWormSample

   ! A type that handles extracting contribution to measured
   ! quantities from the sign and worm part of Monte Carlo
   ! configurations and uses samplers for processing.
   type, abstract, extends(TSamplingMeasurement) :: TSamplingWormMeasurement
   contains
      procedure(worm_to_sampbuf_interface), deferred, nopass :: worm_to_sampbuf
      procedure, pass :: diagram_to_sampbuf => worm_meas_diagram_to_sampbuf
   end type TSamplingWormMeasurement

   ! A type that handles extracting contribution to measured
   ! quantities from the sign and worm part of Monte Carlo
   ! configurations and uses samplers for processing.
   type, abstract, extends(TSamplingMeasurement2D) :: TSamplingWormMeasurement2D
   contains
      procedure(worm_to_sampbuf2d_interface), deferred, nopass :: worm_to_sampbuf2d
      procedure, pass :: diagram_to_sampbuf2d => worm_meas_diagram_to_sampbuf2d
   end type TSamplingWormMeasurement2D

   ! A type that handles extracting contribution to measured
   ! quantities from the sign and worm part of Monte Carlo
   ! configurations and uses samplers for processing.
   type, abstract, extends(TSamplingMeasurement3D) :: TSamplingWormMeasurement3D
   contains
      procedure(worm_to_sampbuf3d_interface), deferred, nopass :: worm_to_sampbuf3d
      procedure, pass :: diagram_to_sampbuf3d => worm_meas_diagram_to_sampbuf3d
   end type TSamplingWormMeasurement3D

   abstract interface
      ! Interface for the concrete measurement implementation that
      ! generates the samples from the worm part w and MC sign sign of
      ! the current configuration and fills them into the passed
      ! buffer sampbuf.
      subroutine worm_to_sampbuf_interface(w, sign, buf)
         import TWormSample, c_double, c_double_complex, ZSampleBuffer
         type(TWormSample), intent(in) :: w
         !real(c_double), intent(in) :: sign
         complex(c_double_complex), intent(in) :: sign
         type(ZSampleBuffer), intent(inout) :: buf
      end subroutine worm_to_sampbuf_interface

      subroutine worm_to_sampbuf2d_interface(w, sign, buf)
         import TWormSample, c_double, c_double_complex, ZSampleBuffer2D
         type(TWormSample), intent(in) :: w
         !real(c_double), intent(in) :: sign
         complex(c_double_complex), intent(in) :: sign
         type(ZSampleBuffer2D), intent(inout) :: buf
      end subroutine worm_to_sampbuf2d_interface

      subroutine worm_to_sampbuf3d_interface(w, sign, buf)
         import TWormSample, c_double, c_double_complex, ZSampleBuffer3D
         type(TWormSample), intent(in) :: w
         !real(c_double), intent(in) :: sign
         complex(c_double_complex), intent(in) :: sign
         type(ZSampleBuffer3D), intent(inout) :: buf
      end subroutine worm_to_sampbuf3d_interface
   end interface

   ! Type for (two-point) susceptibility function <n(tau)n(tau')>
   ! measurement from worm operators.
   type, extends(TSamplingWormMeasurement) :: TWormP2phM
   contains
      procedure, nopass :: worm_to_sampbuf => measure_P2ph_worm
   end type TWormP2phM

   ! Type for (two-point) susceptibility function < (c^\dag
   ! c^\dag)(tau) (cc)(tau') > measurement from worm operators.
   type, extends(TSamplingWormMeasurement) :: TWormP2ppM
   contains
      procedure, nopass :: worm_to_sampbuf => measure_P2pp_worm
   end type TWormP2ppM

   ! Type for symmetric improved estimator measurement from worm
   ! operators.
   type, extends(TSamplingWormMeasurement) :: TWormQQM
   contains
      procedure, nopass :: worm_to_sampbuf => measure_QQ_worm
   end type TWormQQM

   ! Type for Raman diagram <FIXME> measurement from worm operators.
   type, extends(TSamplingWormMeasurement2D) :: TWormRamanM
   contains
      procedure, nopass :: worm_to_sampbuf2d => measure_Raman_worm
   end type TWormRamanM

   ! Type for custom estimator measurement from worm operators.
   type, extends(TSamplingWormMeasurement) :: TWormCustom1DM
   contains
      procedure, nopass :: worm_to_sampbuf => measure_Custom1D_worm
   end type TWormCustom1DM

   ! Type for custom estimator measurement from worm operators.
   type, extends(TSamplingWormMeasurement2D) :: TWormCustom2DM
   contains
      procedure, nopass :: worm_to_sampbuf2d => measure_Custom2D_worm
   end type TWormCustom2DM

   ! Type for custom estimator measurement from worm operators.
   type, extends(TSamplingWormMeasurement3D) :: TWormCustom3DM
   contains
      procedure, nopass :: worm_to_sampbuf3d => measure_Custom3D_worm
   end type TWormCustom3DM
contains

   type(TWormSample) function extract_worm_sample(diagram) result(w)
      type(TDiagram), intent(in) :: diagram
      integer                    :: i
      type(TLocalOper), allocatable :: ops(:)

      allocate(w%taus(get_ntaus(diagram%worm)))
      allocate(ops, source=get_operators(diagram%local))
      do i = 1, size(w%taus)
         w%taus(i) = ops(get_tau_oper_position(diagram%worm, i))%tau
      end do
   end function extract_worm_sample

   !!! FIXME: replace all sectors by python-configured "custom" and
   !!! replace all the rest with one function (possibly merge into the
   !!! one above)

   ! Implementation of diagram_to_sampbuf for subtypes that need only
   ! the worm times and sign of the MC configuration.
   subroutine worm_meas_diagram_to_sampbuf(meas, diag, sampbuf)
      class(TSamplingWormMeasurement), intent(in) :: meas
      type(TDiagram), intent(in) :: diag
      type(ZSampleBuffer), intent(inout) :: sampbuf

      call meas%worm_to_sampbuf( &
                extract_worm_sample(diag), get_diagram_sign(diag), sampbuf)
   end subroutine worm_meas_diagram_to_sampbuf

   subroutine worm_meas_diagram_to_sampbuf2d(meas, diag, sampbuf)
      class(TSamplingWormMeasurement2D), intent(in) :: meas
      type(TDiagram), intent(in) :: diag
      type(ZSampleBuffer2D), intent(inout) :: sampbuf

      call meas%worm_to_sampbuf2d( &
                extract_worm_sample(diag), get_diagram_sign(diag), sampbuf)
   end subroutine worm_meas_diagram_to_sampbuf2d

   subroutine worm_meas_diagram_to_sampbuf3d(meas, diag, sampbuf)
      class(TSamplingWormMeasurement3D), intent(in) :: meas
      type(TDiagram), intent(in) :: diag
      type(ZSampleBuffer3D), intent(inout) :: sampbuf

      call meas%worm_to_sampbuf3d( &
                extract_worm_sample(diag), get_diagram_sign(diag), sampbuf)
   end subroutine worm_meas_diagram_to_sampbuf3d

   ! Get samples of the susceptibility from current times of worm
   ! operators.
   subroutine measure_P2ph_worm(w, sign, buf)
      type(TWormSample), intent(in)      :: w
      complex(c_double_complex), intent(in)         :: sign
      type(ZSampleBuffer), intent(inout) :: buf
      type(ZSamplePoint)                 :: p

      call prepare_buffer(buf, 1)
      p%i = 1
      p%x = w%taus(1) - w%taus(2)
      p%fx = sign
      buf%buffer(1) = p
   end subroutine measure_P2ph_worm

   ! Get samples of the susceptibility from current times of worm
   ! operators.
   subroutine measure_P2pp_worm(w, sign, buf)
      type(TWormSample), intent(in)      :: w
      complex(c_double_complex), intent(in)         :: sign
      type(ZSampleBuffer), intent(inout) :: buf
      type(ZSamplePoint)                 :: p

      call prepare_buffer(buf, 1)
      p%i = 1
      p%x = w%taus(1) - w%taus(2)
      p%fx = sign
      buf%buffer(1) = p
   end subroutine measure_P2pp_worm

   ! Get QQ sample from current times of worm operators.
   subroutine measure_QQ_worm(w, sign, buf)
      type(TWormSample), intent(in)      :: w
      complex(c_double_complex), intent(in)         :: sign
      type(ZSampleBuffer), intent(inout) :: buf
      type(ZSamplePoint)                 :: p

      call prepare_buffer(buf, 1)
      p%i = 1
      p%x = w%taus(2) - w%taus(1)
      p%fx = sign
      buf%buffer(1) = p
   end subroutine measure_QQ_worm

   ! Get samples of the Raman diagram from current times of worm
   ! operators.
   subroutine measure_Raman_worm(w, sign, buf)
      type(TWormSample), intent(in)        :: w
      complex(c_double_complex), intent(in)           :: sign
      type(ZSampleBuffer2D), intent(inout) :: buf
      type(ZSamplePoint2D)                 :: p

      call prepare_buffer(buf, 1)
      p%i = 1
      ! FIXME: was already "check this" in old version
      p%x(1) = w%taus(3) - w%taus(2)
      p%x(2) = w%taus(3) - w%taus(1)
      p%fx = sign
      buf%buffer(1) = p
   end subroutine measure_Raman_worm

   ! Get custom estimator sample from current times of worm operators.
   subroutine measure_Custom1D_worm(w, sign, buf)
      type(TWormSample), intent(in)      :: w
      complex(c_double_complex), intent(in)         :: sign
      type(ZSampleBuffer), intent(inout) :: buf
      type(ZSamplePoint)                 :: p

      call prepare_buffer(buf, 1)
      p%i = 1
      p%x = w%taus(1) - w%taus(2)
      p%fx = sign
      buf%buffer(1) = p
   end subroutine measure_Custom1D_worm

   ! Get custom estimator sample from current times of worm operators.
   subroutine measure_Custom2D_worm(w, sign, buf)
      type(TWormSample), intent(in)      :: w
      complex(c_double_complex), intent(in)         :: sign
      type(ZSampleBuffer2D), intent(inout) :: buf
      type(ZSamplePoint2D)                 :: p

      call prepare_buffer(buf, 1)
      p%i = 1
      p%x(:) = w%taus(1:2) - w%taus(3)
      p%fx = sign
      buf%buffer(1) = p
   end subroutine measure_Custom2D_worm

   ! Get custom estimator sample from current times of worm operators.
   subroutine measure_Custom3D_worm(w, sign, buf)
      type(TWormSample), intent(in)      :: w
      complex(c_double_complex), intent(in)         :: sign
      type(ZSampleBuffer3D), intent(inout) :: buf
      type(ZSamplePoint3D)                 :: p

      call prepare_buffer(buf, 1)
      p%i = 1
      p%x(:) = w%taus(1:3) - w%taus(4)
      p%fx = sign
      buf%buffer(1) = p
   end subroutine measure_Custom3D_worm
end module MWormMeas
