# Find FFTW3 library
#
# Understands the following environment variables
#
#  - FFTWDIR           : path to FFTW
#  - FFTW_ROOT         : path to FFTW
#
# Allows the following cache variables (detected by default):
#
#  - FFTW_ROOT_DIR     : where to find
#  - FFTW_INCLUDE_DIR  : directory to include file
#  - FFTW_LIBRARY      : path to library
#  - FFTWF_LIBRARY     : path to float library
#  - FFTWL_LIBRARY     : path to long double library
#
# Sets the following variables:
#
#  - FFTW_FOUND        : if FFTW has been found
#  - FFTW_INCLUDE_DIRS : header files
#  - FFTW_LIBRARIES    : libraries
#  - FFTW_VERSION      : version info
#
# Imported targets (if library has been found):
#
#  - FFTW::FFTW        : FFTW for double precision
#  - FFTW::FFTWf       : FFTW for singe precision
#  - FFTW::FFTWl       : FFTW for long double precision
#
# Copyright (C) 2024 Markus Wallerberger and others
# SPDX-License-Identifier: MIT
#
cmake_minimum_required(VERSION 3.9...3.31)
include(FindPackageHandleStandardArgs)

find_package(PkgConfig)
if (PKG_CONFIG_FOUND)
    pkg_check_modules(PKG_FFTW QUIET "fftw3")
endif()

set(_fftw_standard_paths ${FFTW_ROOT_DIR} ENV FFTWDIR ENV FFTW_ROOT)

# --- find include directory

if (NOT FFTW_INCLUDE_DIR)
    find_path(FFTW_INCLUDE_DIR
        NAMES "fftw3.h"
        HINTS ${PKG_FFTW_INCLUDE_DIRS} ${_fftw_standard_paths}
        PATH_SUFFIXES include
        DOC "FFTW3 - include directory"
        )
    mark_as_advanced(FFTW_INCLUDE_DIR)
endif()

# --- find libraries

function(_find_fftw_library varname libname)
    if (NOT "${varname}")
        find_library("${varname}"
            NAMES "${libname}"
            HINTS ${PKG_FFTW_LIBRARY_DIRS} ${_fftw_standard_paths}
            PATH_SUFFIXES lib lib64
            DOC "FFTW3 - library"
            )
        mark_as_advanced("${varname}")
    endif()
endfunction()

_find_fftw_library(FFTW_LIBRARY fftw3)
_find_fftw_library(FFTWF_LIBRARY fftw3f)
_find_fftw_library(FFTWL_LIBRARY fftw3l)

# --- find version

function(_extract_fftw_version versionvar)
    set(testsrc "${CMAKE_BINARY_DIR}/fftw_check.c")
    string(JOIN "\n" filetext
        "/* Test program to extract fftw version from library itself */"
        "#include <stdio.h>"
        "#include <fftw3.h>"
        ""
        "int main()"
        "{"
        "    printf(\"%s\", fftw_version);"
        "}\n"
        )
    file(WRITE "${testsrc}" "${filetext}")
    try_run(run_result compile_result "${CMAKE_BINARY_DIR}" "${testsrc}"
        CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${FFTW_INCLUDE_DIR}"
        LINK_LIBRARIES ${FFTW_LIBRARY}
        COMPILE_OUTPUT_VARIABLE compile_output
        RUN_OUTPUT_VARIABLE run_output
        )
    if (NOT compile_result)
        message(WARNING "FFTW test - error compiling: ${compile_output}")
        return()
    elseif (run_result)
        message(WARNING "FFTW test - error running: ${run_output}")
        return()
    endif()

    message(VERBOSE "FFTW test - outputs: ${run_output}")
    if (run_output MATCHES "fftw-([0-9][0-9.]+[0-9])")
        set(version "${CMAKE_MATCH_1}")
        message(VERBOSE "FFTW extracted version: ${version}")

        # export back to parent scope ("pass by reference")
        set("${versionvar}" "${version}" PARENT_SCOPE)
    else()
        message(WARNING "FFTW test - unable to extract version")
    endif()
endfunction()

if (FFTW_INCLUDE_DIR)
    if (NOT FFTW_VERSION)
        _extract_fftw_version(FFTW_VERSION)
    endif()
    if (NOT FFTW_VERSION AND PKG_FFTW_VERSION)
        set(FFTW_VERSION "${PKG_FFTW_VERSION}")
        message(VERBOSE "FFTW packagekit-detected version: ${FFTW_VERSION}")
    endif()
endif()

# --- set outputs

find_package_handle_standard_args(FFTW
    REQUIRED_VARS FFTW_LIBRARY FFTW_INCLUDE_DIR
    VERSION_VAR FFTW_VERSION
    )

function(_add_fftw_target targetname libname)
    if (DEFINED "${libname}" AND NOT TARGET "${targetname}")
        add_library("${targetname}" UNKNOWN IMPORTED)
        set_target_properties("${targetname}" PROPERTIES
            IMPORTED_LOCATION "${${libname}}"
            INTERFACE_COMPILE_OPTIONS "${FFTW_DEFINITIONS}"
            INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIR}"
            )
        message(VERBOSE "FFTW target ${targetname} linking to ${libname}")
    endif()
endfunction()

if (FFTW_FOUND)
    set(FFTW_INCLUDE_DIRS "${FFTW_INCLUDE_DIR}")
    set(FFTW_DEFINITIONS ${PKG_FFTW_CFLAGS_OTHER})
    set(FFTW_LIBRARIES "${FFTW_LIBRARY}" "${FFTWF_LIBRARY}" "${FFTWL_LIBRARY}")

    _add_fftw_target(FFTW::FFTW FFTW_LIBRARY)
    _add_fftw_target(FFTW::FFTWf FFTWF_LIBRARY)
    _add_fftw_target(FFTW::FFTWl FFTWL_LIBRARY)
endif()
