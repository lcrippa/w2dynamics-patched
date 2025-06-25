module MDiagramMeasurements
    use iso_c_binding, only: c_double_complex
    use MDiagram
    use MAccumulator
    use caccumulators
    use MSamplerD
    use MSamplerZ
    use csamplers
    use MSamplerRegister
    implicit none

    ! A type that handles extracting contribution to measured
    ! quantities from Monte Carlo configurations.
    type, abstract :: TMeasurement
    contains
        ! Perform one measurement of the quantity the concrete
        ! implementation handles on the passed diagram and add it to the
        ! stored data.
        procedure(measure_interface), deferred, pass :: measure
        ! Clear all references to non-owned memory contained by the
        ! concrete implementation if applicable.
        procedure(decouple_interface), deferred, pass :: decouple
    end type TMeasurement

    ! A type that handles extracting contribution to measured
    ! quantities from Monte Carlo configurations and uses an
    ! accumulator for storage and processing.
    type, abstract, extends(TMeasurement) :: TAccumulatingMeasurement
        type(TAccumulator), pointer :: accumulator => null()
    contains
        ! Set the (single) internal accumulator pointer to the passed
        ! accumulator so it will be used to handle new data points from
        ! the next measurement onward.
        procedure, pass :: set_accumulator => accu_meas_set_accumulator
        ! Set internal accumulator pointer as above, but takes a C void*
        ! pointing to an accumulator reference wrapper type.
        procedure, pass :: c_set_accumulator => accu_meas_c_set_accumulator
        ! Concrete measurement implementation that generates the
        ! correctly shaped contribution to the accumulator buffer from
        ! the current diagram.
        procedure(diagram_to_array_interface), deferred, pass :: diagram_to_array
        ! Checks whether we have an accumulator and if so, whether that
        ! has mean or blocks available
        procedure, pass :: has_mean => accu_has_mean
        procedure, pass :: has_blocks => accu_has_blocks
        procedure, pass :: measure => accu_measure
        procedure, pass :: decouple => accu_meas_clear_reference
    end type TAccumulatingMeasurement

    ! A type that handles extracting contribution to measured
    ! quantities from Monte Carlo configurations and uses samplers for
    ! processing.
    type, abstract, extends(TMeasurement) :: TSamplingMeasurement
        type(sampler_register) :: sampler_register
    contains
        ! Add references to a sampler and its backing accumulator to the
        ! list of samplers that will be passed samples new measurements.
        procedure, pass :: add_sampler => samp_meas_add_sampler
        ! Add references to a sampler and its backing accumulator as
        ! above, but takes C void* to reference wrapper types.
        procedure, pass :: c_add_sampler => samp_meas_c_add_sampler
        ! Concrete measurement implementation that generates the samples
        ! from the current diagram and fills them into the passed
        ! buffer.
        procedure(diagram_to_sampbuf_interface), deferred, pass :: diagram_to_sampbuf
        procedure, pass :: measure => samp_measure
        procedure, pass :: decouple => samp_meas_clear_register
    end type TSamplingMeasurement

    ! A type that handles extracting contribution to measured
    ! quantities from Monte Carlo configurations and uses samplers for
    ! processing.
    type, abstract, extends(TMeasurement) :: TSamplingMeasurement2D
        type(sampler2d_register) :: sampler_register
    contains
        ! Add references to a sampler and its backing accumulator to the
        ! list of samplers that will be passed samples new measurements.
        procedure, pass :: add_sampler => samp_meas_add_sampler2d
        ! Add references to a sampler and its backing accumulator as
        ! above, but takes C void* to reference wrapper types.
        procedure, pass :: c_add_sampler => samp_meas_c_add_sampler2d
        ! Concrete measurement implementation that generates the samples
        ! from the current diagram and fills them into the passed
        ! buffer.
        procedure(diagram_to_sampbuf2d_interface), deferred, pass :: diagram_to_sampbuf2d
        procedure, pass :: measure => samp2d_measure
        procedure, pass :: decouple => samp2d_meas_clear_register
    end type TSamplingMeasurement2D

    ! A type that handles extracting contribution to measured
    ! quantities from Monte Carlo configurations and uses samplers for
    ! processing.
    type, abstract, extends(TMeasurement) :: TSamplingMeasurement3D
        type(sampler3d_register) :: sampler_register
    contains
        ! Add references to a sampler and its backing accumulator to the
        ! list of samplers that will be passed samples new measurements.
        procedure, pass :: add_sampler => samp_meas_add_sampler3d
        ! Add references to a sampler and its backing accumulator as
        ! above, but takes C void* to reference wrapper types.
        procedure, pass :: c_add_sampler => samp_meas_c_add_sampler3d
        ! Concrete measurement implementation that generates the samples
        ! from the current diagram and fills them into the passed
        ! buffer.
        procedure(diagram_to_sampbuf3d_interface), deferred, pass :: diagram_to_sampbuf3d
        procedure, pass :: measure => samp3d_measure
        procedure, pass :: decouple => samp3d_meas_clear_register
    end type TSamplingMeasurement3D

    abstract interface
        ! Perform one measurement of the quantity the concrete
        ! implementation meas handles on the passed diagram diag and add
        ! it to the measurement objects backing storage.
        subroutine measure_interface(meas, diag)
            import :: TMeasurement, TDiagram
            class(TMeasurement), intent(inout) :: meas
            type(TDiagram), intent(in) :: diag
        end subroutine measure_interface

        ! Clear all references to non-owned memory contained by the
        ! measurement object meas if applicable.
        subroutine decouple_interface(meas)
            import :: TMeasurement
            class(TMeasurement), intent(inout) :: meas
        end subroutine decouple_interface

        ! Interface for the concrete measurement implementation that
        ! generates the correctly shaped contribution to the accumulator
        ! buffer out from the current diagram diag. The argument meas
        ! may be used in implementations of this subroutine in abstract
        ! subtypes to call methods that process restricted subsets of
        ! the full configuration information (enforced by the interface)
        ! after preprocessing.
        subroutine diagram_to_array_interface(meas, diag, out)
            import :: TAccumulatingMeasurement, TDiagram, c_double, c_double_complex
            class(TAccumulatingMeasurement), intent(in) :: meas
            type(TDiagram), intent(in) :: diag
            complex(c_double_complex), intent(inout) :: out(:)
        end subroutine diagram_to_array_interface

        ! Interface for the concrete measurement implementation that
        ! generates the samples from the current diagram diag and fills
        ! them into the passed buffer sampbuf. The argument meas may be
        ! used in implementations of this subroutine in abstract
        ! subtypes to call methods that process restricted subsets of
        ! the full configuration information (enforced by the interface)
        ! after preprocessing.
        subroutine diagram_to_sampbuf_interface(meas, diag, sampbuf)
            import :: TSamplingMeasurement, TDiagram, ZSampleBuffer
            class(TSamplingMeasurement), intent(in) :: meas
            type(TDiagram), intent(in) :: diag
            type(ZSampleBuffer), intent(inout) :: sampbuf
        end subroutine diagram_to_sampbuf_interface

        subroutine diagram_to_sampbuf2d_interface(meas, diag, sampbuf)
            import :: TSamplingMeasurement2D, TDiagram, ZSampleBuffer2D
            class(TSamplingMeasurement2D), intent(in) :: meas
            type(TDiagram), intent(in) :: diag
            type(ZSampleBuffer2D), intent(inout) :: sampbuf
        end subroutine diagram_to_sampbuf2d_interface

        subroutine diagram_to_sampbuf3d_interface(meas, diag, sampbuf)
            import :: TSamplingMeasurement3D, TDiagram, ZSampleBuffer3D
            class(TSamplingMeasurement3D), intent(in) :: meas
            type(TDiagram), intent(in) :: diag
            type(ZSampleBuffer3D), intent(inout) :: sampbuf
        end subroutine diagram_to_sampbuf3d_interface
    end interface

    ! MC congiguration sign / weight measurement type.
    type, extends(TAccumulatingMeasurement) :: TSignM
    contains
        procedure, pass :: diagram_to_array => measure_sign
    end type TSignM
contains

    ! Implementation of the measurement of the contribution of a single
    ! diagram for a measurement type using an accumulator as backing
    ! storage / for statistical processing. The measurement is
    ! performed by calling the concrete implementation of
    ! diagram_to_array to write to the accumulator's buffer and
    ! signalling the addition of a new contribution by calling the
    ! accumulator's add method.
    subroutine accu_measure(meas, diag)
        class(TAccumulatingMeasurement), intent(inout) :: meas
        type(TDiagram), intent(in) :: diag
        complex(c_double_complex), pointer :: buf(:)

        if (associated(meas%accumulator)) then
            call accumulator_buffer(meas%accumulator, buf)
            call meas%diagram_to_array(diag, buf)
            call accumulator_add(meas%accumulator)
        end if
    end subroutine accu_measure

    ! Set the (single) internal accumulator pointer of meas to the
    ! passed accumulator acc so it will be used to handle new data
    ! points from the next measurement using meas onward.
    subroutine accu_meas_set_accumulator(meas, acc)
        class(TAccumulatingMeasurement), intent(inout) :: meas
        type(TAccumulator), intent(inout), target :: acc

        meas%accumulator => acc
    end subroutine accu_meas_set_accumulator

    ! Set internal accumulator pointer of meas as above, but takes a C
    ! void* cacc pointing to an accumulator reference wrapper type.
    subroutine accu_meas_c_set_accumulator(meas, cacc)
        class(TAccumulatingMeasurement), intent(inout) :: meas
        type(c_ptr), intent(inout) :: cacc
        type(TAccumulator), pointer :: acc

        call c_f_pointer(cacc, acc)
        call meas%set_accumulator(acc)
    end subroutine accu_meas_c_set_accumulator

    ! Disassociates the internal accumulator pointer of meas.
    subroutine accu_meas_clear_reference(meas)
        class(TAccumulatingMeasurement), intent(inout) :: meas

        nullify(meas%accumulator)
    end subroutine accu_meas_clear_reference

    ! Do we have an accumulator with mean?
    logical function accu_has_mean(meas)
        class(TAccumulatingMeasurement), intent(inout) :: meas

        if (associated(meas%accumulator)) then
            accu_has_mean = accumulator_has_mean(meas%accumulator)
        else
            accu_has_mean = .false.
        end if
    end function accu_has_mean

    ! Do we have an accumulator with blocks?
    logical function accu_has_blocks(meas)
        class(TAccumulatingMeasurement), intent(inout) :: meas

        if (associated(meas%accumulator)) then
            accu_has_blocks = accumulator_has_blocks(meas%accumulator)
        else
            accu_has_blocks = .false.
        end if
    end function accu_has_blocks

    ! Implementation of the measurement of the contribution of a single
    ! diagram for a measurement type using a sampler for statistical
    ! processing. The measurement is performed by calling the concrete
    ! implementation of diagram_to_sampbuf to write to a sample buffer
    ! and passing the buffer on to the sampler register to add the
    ! points to all registered samplers.
    subroutine samp_measure(meas, diag)
        class(TSamplingMeasurement), intent(inout) :: meas
        type(TDiagram), intent(in) :: diag
        type(ZSampleBuffer) :: sampbuf  ! FIXME: decide whether to keep buffer single-use and directly here

        if (.not. meas%sampler_register%empty) then  ! to not waste time calling diagram_to_sampbuf
            call meas%diagram_to_sampbuf(diag, sampbuf)
            call meas%sampler_register%add(sampbuf%buffer(1:sampbuf%n))
        end if
    end subroutine samp_measure

    ! Add references to a sampler samp and its backing accumulator acc
    ! to the list of samplers that will be passed samples new
    ! measurements using meas.
    subroutine samp_meas_add_sampler(meas, samp, acc)
        class(TSamplingMeasurement), intent(inout) :: meas
        class(ZSampler), intent(inout), target :: samp
        type(TAccumulator), intent(inout), target :: acc

        call meas%sampler_register%register(samp, acc)
    end subroutine samp_meas_add_sampler

    ! Add references to a sampler and its backing accumulator as above,
    ! but takes C void* csamp and cacc to reference wrapper types.
    subroutine samp_meas_c_add_sampler(meas, csamp, cacc)
        class(TSamplingMeasurement), intent(inout) :: meas
        type(c_ptr), intent(inout) :: csamp, cacc
        type(csampler), pointer :: samp
        type(TAccumulator), pointer :: acc

        call c_f_pointer(csamp, samp)
        call c_f_pointer(cacc, acc)
        call meas%add_sampler(samp%inner, acc)
    end subroutine samp_meas_c_add_sampler

    ! Clear the sampler register of meas, removing all references.
    subroutine samp_meas_clear_register(meas)
        class(TSamplingMeasurement), intent(inout) :: meas

        call meas%sampler_register%clear()
    end subroutine samp_meas_clear_register

    ! Implementation of the measurement of the contribution of a single
    ! diagram for a measurement type using a sampler for statistical
    ! processing. The measurement is performed by calling the concrete
    ! implementation of diagram_to_sampbuf to write to a sample buffer
    ! and passing the buffer on to the sampler register to add the
    ! points to all registered samplers.
    subroutine samp2d_measure(meas, diag)
        class(TSamplingMeasurement2D), intent(inout) :: meas
        type(TDiagram), intent(in) :: diag
        type(ZSampleBuffer2D) :: sampbuf  ! FIXME: decide whether to keep buffer single-use and directly here

        if (.not. meas%sampler_register%empty) then  ! to not waste time calling diagram_to_sampbuf
            call meas%diagram_to_sampbuf2d(diag, sampbuf)
            call meas%sampler_register%add(sampbuf%buffer(1:sampbuf%n))
        end if
    end subroutine samp2d_measure

    ! Add references to a sampler samp and its backing accumulator acc
    ! to the list of samplers that will be passed samples new
    ! measurements using meas.
    subroutine samp_meas_add_sampler2d(meas, samp, acc)
        class(TSamplingMeasurement2D), intent(inout) :: meas
        class(ZSampler2D), intent(inout), target :: samp
        type(TAccumulator), intent(inout), target :: acc

        call meas%sampler_register%register(samp, acc)
    end subroutine samp_meas_add_sampler2d

    ! Add references to a sampler and its backing accumulator as above,
    ! but takes C void* csamp and cacc to reference wrapper types.
    subroutine samp_meas_c_add_sampler2d(meas, csamp, cacc)
        class(TSamplingMeasurement2D), intent(inout) :: meas
        type(c_ptr), intent(inout) :: csamp, cacc
        type(csampler2d), pointer :: samp
        type(TAccumulator), pointer :: acc

        call c_f_pointer(csamp, samp)
        call c_f_pointer(cacc, acc)
        call meas%add_sampler(samp%inner, acc)
    end subroutine samp_meas_c_add_sampler2d

    ! Clear the sampler register of meas, removing all references.
    subroutine samp2d_meas_clear_register(meas)
        class(TSamplingMeasurement2D), intent(inout) :: meas

        call meas%sampler_register%clear()
    end subroutine samp2d_meas_clear_register

    ! Implementation of the measurement of the contribution of a single
    ! diagram for a measurement type using a sampler for statistical
    ! processing. The measurement is performed by calling the concrete
    ! implementation of diagram_to_sampbuf to write to a sample buffer
    ! and passing the buffer on to the sampler register to add the
    ! points to all registered samplers.
    subroutine samp3d_measure(meas, diag)
        class(TSamplingMeasurement3D), intent(inout) :: meas
        type(TDiagram), intent(in) :: diag
        type(ZSampleBuffer3D) :: sampbuf  ! FIXME: decide whether to keep buffer single-use and directly here

        if (.not. meas%sampler_register%empty) then  ! to not waste time calling diagram_to_sampbuf
            call meas%diagram_to_sampbuf3d(diag, sampbuf)
            call meas%sampler_register%add(sampbuf%buffer(1:sampbuf%n))
        end if
    end subroutine samp3d_measure

    ! Add references to a sampler samp and its backing accumulator acc
    ! to the list of samplers that will be passed samples new
    ! measurements using meas.
    subroutine samp_meas_add_sampler3d(meas, samp, acc)
        class(TSamplingMeasurement3D), intent(inout) :: meas
        class(ZSampler3D), intent(inout), target :: samp
        type(TAccumulator), intent(inout), target :: acc

        call meas%sampler_register%register(samp, acc)
    end subroutine samp_meas_add_sampler3d

    ! Add references to a sampler and its backing accumulator as above,
    ! but takes C void* csamp and cacc to reference wrapper types.
    subroutine samp_meas_c_add_sampler3d(meas, csamp, cacc)
        class(TSamplingMeasurement3D), intent(inout) :: meas
        type(c_ptr), intent(inout) :: csamp, cacc
        type(csampler3d), pointer :: samp
        type(TAccumulator), pointer :: acc

        call c_f_pointer(csamp, samp)
        call c_f_pointer(cacc, acc)
        call meas%add_sampler(samp%inner, acc)
    end subroutine samp_meas_c_add_sampler3d

    ! Clear the sampler register of meas, removing all references.
    subroutine samp3d_meas_clear_register(meas)
        class(TSamplingMeasurement3D), intent(inout) :: meas

        call meas%sampler_register%clear()
    end subroutine samp3d_meas_clear_register

    ! Measures the MC configuration sign from the current configuration
    ! diag and adds it to buffer out.
    subroutine measure_sign(meas, diag, out)
        class(TSignM), intent(in) :: meas
        type(TDiagram), intent(in) :: diag
        complex(c_double_complex), intent(inout) :: out(:)

        if (size(out) /= 1) &
          error stop 'Wrong size'

        out(:) = out(:) + get_diagram_sign(diag)
    end subroutine measure_sign
end module MDiagramMeasurements
