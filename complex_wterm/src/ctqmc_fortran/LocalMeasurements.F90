module MLocalMeasurements
    use iso_c_binding, only: c_double_complex
    use MImpurity
    use MOperator
    use MDiagramMeasurements
    implicit none

    ! A type that handles extracting contribution to measured
    ! quantities from the sign and local part of Monte Carlo
    ! configurations and uses an accumulator for storage and
    ! processing.
    type, abstract, extends(TAccumulatingMeasurement) :: TAccumulatingLocalMeasurement
    contains
        procedure(local_to_array_interface), deferred, nopass :: local_to_array
        procedure, pass :: diagram_to_array => local_meas_diagram_to_array
    end type TAccumulatingLocalMeasurement

    ! A type that handles extracting contribution to measured
    ! quantities from the sign and local part of Monte Carlo
    ! configurations and uses samplers for processing.
    type, abstract, extends(TSamplingMeasurement) :: TSamplingLocalMeasurement
    contains
        procedure(local_to_sampbuf_interface), deferred, nopass :: local_to_sampbuf
        procedure, pass :: diagram_to_sampbuf => local_meas_diagram_to_sampbuf
    end type TSamplingLocalMeasurement

    abstract interface
        ! Interface for the concrete measurement implementation that
        ! generates the correctly shaped contribution to the accumulator
        ! buffer out from the local part l and MC sign sign of the
        ! current configuration.
        subroutine local_to_array_interface(l, sign, out)
            import :: TImpurity, c_double, c_double_complex
            type(TImpurity), intent(in) :: l
            complex(c_double_complex), intent(in) :: sign
            complex(c_double_complex), intent(inout) :: out(:)
        end subroutine local_to_array_interface

        ! Interface for the concrete measurement implementation that
        ! generates the samples from the local part l and MC sign sign of
        ! the current configuration and fills them into the passed
        ! buffer sampbuf.
        subroutine local_to_sampbuf_interface(l, sign, buf)
            import TImpurity, c_double, c_double_complex, ZSampleBuffer
            type(TImpurity), intent(in) :: l
            complex(c_double_complex), intent(in) :: sign
            type(ZSampleBuffer), intent(inout) :: buf
        end subroutine local_to_sampbuf_interface
    end interface

    ! Type for density matrix measurement from the local part of the
    ! MC configuration.
    type, extends(TAccumulatingLocalMeasurement) :: TDensityMatrixM
    contains
        procedure, nopass :: local_to_array => measure_dm
    end type TDensityMatrixM

    ! Type for susceptibility measurement from the local part of the
    ! MC configuration.
    !!! TODOCOMPL: activate again
    !type, extends(TAccumulatingLocalMeasurement) :: TNtauN0M
    !contains
        !procedure, nopass :: local_to_array => measure_ntaun0
    !end type TNtauN0M
contains

    ! Implementation of diagram_to_array for subtypes that need only
    ! the local part and sign of the MC configuration.
    subroutine local_meas_diagram_to_array(meas, diag, out)
        class(TAccumulatingLocalMeasurement), intent(in) :: meas
        type(TDiagram), intent(in) :: diag
        complex(c_double_complex), intent(inout) :: out(:)

        call meas%local_to_array(diag%local, get_diagram_sign(diag), out)
    end subroutine local_meas_diagram_to_array

    ! Implementation of diagram_to_sampbuf for subtypes that need only
    ! the local part and sign of the MC configuration.
    subroutine local_meas_diagram_to_sampbuf(meas, diag, sampbuf)
        class(TSamplingLocalMeasurement), intent(in) :: meas
        type(TDiagram), intent(in) :: diag
        type(ZSampleBuffer), intent(inout) :: sampbuf

        call meas%local_to_sampbuf(diag%local, get_diagram_sign(diag), sampbuf)
    end subroutine local_meas_diagram_to_sampbuf

    ! Calculate density matrix contribution of the current configuration.
    subroutine measure_dm(l, sign, out)
        type(TImpurity), intent(in)          :: l
        complex(c_double_complex), intent(in)       :: sign
        complex(c_double_complex), intent(inout)          :: out(:)

        complex(c_double_complex), allocatable :: out_zmat(:,:)
        real(c_double), allocatable :: out_dmat(:,:)
        type(TOperator) :: out_op

        call compute_density_matrix(l, sign, out_op)

        ! XXX we do not exploit the block structure of the density matrix
        if (is_complex_op(out_op)) then
            allocate(out_zmat, source=getmatrix(out_op%zop))
            out(:) = out(:) + reshape(out_zmat, (/ size(out_zmat) /))
        else
            allocate(out_dmat, source=getmatrix(out_op%dop))
            out(:) = out(:) + reshape(out_dmat, (/ size(out_dmat) /))
        endif
    end subroutine measure_dm

end module MLocalMeasurements
