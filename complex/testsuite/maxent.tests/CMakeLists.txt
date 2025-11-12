cmake_minimum_required(VERSION 3.9...3.31)

function(declare_test binary_file src_file)
    # Add the executable for the test case
    add_executable(${binary_file} ${src_file})

    # Link the proper libraries
    set_target_properties(${binary_file} PROPERTIES LINKER_LANGUAGE Fortran)
    target_link_libraries(${binary_file} MAXENT)
    target_link_libraries(${binary_file} ${Python3_LIBRARIES})
    target_link_libraries(${binary_file} ${LAPACK_LIBRARIES})

    # Add the test case
    add_test(${binary_file} ${binary_file})
endfunction()

function(add_test_aux_file filename)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/${filename}
        ${CMAKE_CURRENT_BINARY_DIR}/${filename}
        COPYONLY)
endfunction()

enable_testing()
declare_test(1-getgridunit 1-getgridunit.F90)
