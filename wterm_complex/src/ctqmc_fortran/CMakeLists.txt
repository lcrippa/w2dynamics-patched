cmake_minimum_required(VERSION 3.9...3.31)

# Templated files
foreach (value_type double zcomplex)
    set("ctqmc_src_${value_type}"
        ${value_type}/Accumulator.F90
        ${value_type}/Bath.F90
        ${value_type}/BlockMatrix.F90
        ${value_type}/Buffer.F90
        ${value_type}/ExtRange.F90
        ${value_type}/Hybridization.F90
        ${value_type}/Impurity.F90
        ${value_type}/Operator.F90
        ${value_type}/Psi.F90
        ${value_type}/QRDecomposition.F90
        ${value_type}/Sampler.F90
        ${value_type}/UpdatableMatrix.F90
        )
endforeach()

add_library(CTQMC SHARED
  CTQMC.F90
  libraries/BLAS.F90
  libraries/FFTW.F90
  libraries/LAPACK.F90
  libraries/LibC.F90
  Parameters.F90
  AccumulatorBase.F90
  RandomNumbers.F90
  Accumulator.F90
  c-interface/Accumulator.F90
  SamplerRegister.F90
  c-interface/Sampler.F90
  Common.F90
  LogDiscretization.F90
  CompoundIndex.F90
  WormState.F90
  Diagram.F90
  States.F90
  BathBase.F90
  Bath.F90
  BathMeasurements.F90
  LocalBase.F90
  LocalMeasurements.F90
  Impurity.F90
  ComplexNumbers.F90
  MoveStatistics.F90
  DiagramMoves.F90
  DiagramMeasurements.F90
  ExtRange.F90
  MeasurementInterface.F90
  Operator.F90
  Psi.F90
  Permutation.F90
  Printing.F90
  Progress.F90
  Signals.F90
  integer/Sorting.F90
  double/Sorting.F90
  Sorting.F90
  Utilities.F90
  SymmetryMoves.F90
  Testing.F90
  Window.F90
  WormMeas.F90
  integer/Buffer.F90
  ${ctqmc_src_double}
  ${ctqmc_src_zcomplex}
  )
set_target_properties(CTQMC PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/w2dyn/auxiliaries")
target_link_libraries(CTQMC PRIVATE
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    FFTW::FFTW
    )
