# CMake project file for ctqmc_fortran

##################################################
# Define the project and the dependencies that it has
##################################################

cmake_minimum_required(VERSION 3.9...3.31)
PROJECT(ctqmc Fortran C)

# Set the module path, so that CMake also considers our modules
SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules/")

# Uncomment if it is required that Fortran 90 is supported
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Set some options the user may choose
option(OPT_HOST_ARCH "Target host machine instruction set" ON)
option(TRY_PIP_INSTALL "Try to install Python package dependencies using pip" ON)

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, and TESTING.  You should  review this file and make sure the flags
# are to your liking.
include(SetFortranFlags)

find_package(BLAS)
find_package(LAPACK)

# There is an error in CMAKE with this flag for pgf90.  Unset it
GET_FILENAME_COMPONENT(FCNAME ${CMAKE_Fortran_COMPILER} NAME)
IF(FCNAME STREQUAL "pgf90")
    UNSET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
ENDIF(FCNAME STREQUAL "pgf90")

############################################################
# Define the actual files and folders that make up the build
############################################################

include(PythonPolyfill)
find_python_interpreter()
find_package(PackageHandleStandardArgs)

if (TRY_PIP_INSTALL)
    include(installviapip)
else()
    function(install_via_pip dummy result)
        set(${result} 0 PARENT_SCOPE)
    endfunction()
endif()

find_package(NUMPY 1.10)
IF(NOT NUMPY_FOUND)
    install_via_pip("numpy>=1.10" ERRORCODE)
    if (${ERRORCODE})
    MESSAGE(WARNING "Failed to find or install test and runtime dependency numpy")
    endif()
ENDIF()


find_package(H5PY)
if(NOT H5PY_FOUND)
    install_via_pip("h5py" ERRORCODE)
    if (${ERRORCODE})
    MESSAGE(WARNING "Failed to find or install h5py required for some tests and at runtime")
    endif()
endif()

find_package(FFTW)
if(FFTW_FOUND)
  SET( CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -DHAVE_FFTW" )
ENDIF()

find_package(MPI4PY)
if(NOT MPI4PY_FOUND)
    install_via_pip("mpi4py" ERRORCODE)
    if (${ERRORCODE})
    MESSAGE(WARNING "Failed to find or install mpi4py required for some tests and at runtime")
    endif()
endif()

find_package(SCIPY 0.10)
IF(NOT SCIPY_FOUND)
    install_via_pip("scipy" ERRORCODE)
    if (${ERRORCODE})
    MESSAGE(WARNING "Failed to find or install scipy required for some tests and at runtime")
    endif()
ENDIF()

find_package(CONFIGOBJ)
if(NOT CONFIGOBJ_FOUND)
    install_via_pip("configobj" ERRORCODE)
    if (${ERRORCODE})
    MESSAGE(WARNING "Failed to find or install configobj required for some tests and at runtime")
    endif()
endif()

# Define the executable name
SET(CTQMCEXE ctqmc)

# Define some directories
SET(SRC ${PROJECT_SOURCE_DIR}/src)
SET(LIB ${CMAKE_CURRENT_BINARY_DIR}/lib)
#SET(BIN ${PROJECT_SOURCE_DIR}/bin)

SET(SRCCTQMC ${SRC}/ctqmc_fortran)
SET(SRCMAXENT ${SRC}/maxent)

# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${LIB})

# The source for the ctqmc binary and have it placed in the bin folder
ADD_SUBDIRECTORY(${SRCCTQMC})
ADD_SUBDIRECTORY(${SRCMAXENT})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean
    COMMAND ${CMAKE_COMMAND} -P ${PROJECT_SOURCE_DIR}/distclean.cmake
)

ADD_SUBDIRECTORY(testsuite/maxent.tests)
ADD_SUBDIRECTORY(testsuite/ctqmc.tests)
ADD_SUBDIRECTORY(testsuite/w2dyn.tests)
enable_testing()


# INSTALLATION

message(STATUS "Configuration for installation (optional):")

# default / fallback paths relative to CMAKE_INSTALL_PREFIX
set(INSTALL_PATH_SITE_PACKAGES "lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages")
set(INSTALL_PATH_SCRIPTS "bin")

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_path('platlib'));"
    RESULT_VARIABLE PY_SYSCONFIG_PLATLIB_FAILURE
    OUTPUT_VARIABLE PY_SYSCONFIG_PLATLIB_PATH
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if (NOT PY_SYSCONFIG_PLATLIB_FAILURE)
    message(STATUS "Using path obtained from Python sysconfig module for package installation")
    message(STATUS "Package installation path: ${PY_SYSCONFIG_PLATLIB_PATH}")
    set(INSTALL_PATH_SITE_PACKAGES "${PY_SYSCONFIG_PLATLIB_PATH}")
  endif()

  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_path('scripts'));"
    RESULT_VARIABLE PY_SYSCONFIG_SCRIPTS_FAILURE
    OUTPUT_VARIABLE PY_SYSCONFIG_SCRIPTS_PATH
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if (NOT PY_SYSCONFIG_SCRIPTS_FAILURE)
    message(STATUS "Using path obtained from Python sysconfig module for script installation")
    message(STATUS "Script installation path: ${PY_SYSCONFIG_SCRIPTS_PATH}")
    set(INSTALL_PATH_SCRIPTS "${PY_SYSCONFIG_SCRIPTS_PATH}")
  endif()

  if ((NOT PY_SYSCONFIG_PLATLIB_FAILURE) OR (NOT PY_SYSCONFIG_SCRIPTS_FAILURE))
    message(STATUS "Set CMAKE_INSTALL_PREFIX to override default installation paths")
  endif()

else()

  message(STATUS "Installation prefix directory (variable CMAKE_INSTALL_PREFIX): ${CMAKE_INSTALL_PREFIX}")

endif()

option(INSTALL_EXECUTABLES "Install executable scripts DMFT.py, cthyb, Maxent.py and hgrep to the executable directory (depends on prefix, e.g. /usr/bin)" ON)

message(STATUS "Installing executable scripts DMFT.py, cthyb, Maxent.py, hgrep (option INSTALL_EXECUTABLES): ${INSTALL_EXECUTABLES}")

install(DIRECTORY w2dyn
  DESTINATION "${INSTALL_PATH_SITE_PACKAGES}"
  COMPONENT Python_Package)
if(INSTALL_EXECUTABLES)
  set(SCRIPT_FILES DMFT.py cthyb hgrep Maxent.py)
  foreach(SCRIPT_FILE ${SCRIPT_FILES})
    install(
      CODE  "file(READ \"${CMAKE_CURRENT_SOURCE_DIR}/${SCRIPT_FILE}\" SCRIPT_CONTENT)
             file(WRITE \"${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT_FILE}.cmake.out\" \"#!${Python3_EXECUTABLE}\n# Preceding path inserted by cmake during installation\n\")
             file(APPEND \"${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT_FILE}.cmake.out\" \"\${SCRIPT_CONTENT}\")"
      COMPONENT Executable_Scripts
    )
    install(
      PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT_FILE}.cmake.out"
      DESTINATION "${INSTALL_PATH_SCRIPTS}"
      COMPONENT Executable_Scripts
      RENAME "${SCRIPT_FILE}"
    )
    install(CODE  "file(REMOVE \"${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT_FILE}.cmake.out\")"
            COMPONENT Executable_Scripts)
  endforeach()
endif()
