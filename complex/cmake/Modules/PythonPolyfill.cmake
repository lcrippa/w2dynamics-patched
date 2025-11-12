# Polyfills for old python versions to work around CMP0148

macro(find_python_dev)
    if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.12)
        find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
    else()
        find_python_interpreter()
        find_package(PythonLibs REQUIRED)

        set(Python3_Development_FOUND "${PYTHONLIBS_FOUND}")
        set(Python3_LIBRARIES         "${PYTHON_LIBRARIES}")
        set(Python3_INCLUDE_DIRS      "${PYTHON_INCLUDE_DIRS}")
    endif()
endmacro()

macro(find_python_interpreter)
    if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.12)
        find_package(Python3 COMPONENTS Interpreter REQUIRED)
    else()
        if (NOT Python3_Interpreter_FOUND)
            message(WARNING
                "Outdated cmake version ${CMAKE_VERSION}, using hacks to make"
                " Python work. An update to cmake >= 3.12 is recommended.")
        endif()

        find_package(PythonInterp REQUIRED)
        if (PYTHON_VERSION_MAJOR LESS 3)
            message(ERROR "Python 3 is required")
        endif()

        set(Python3_Interpreter_FOUND "${PYTHONINTERP_FOUND}")
        set(Python3_EXECUTABLE        "${PYTHON_EXECUTABLE}")
        set(Python3_VERSION           "${PYTHON_VERSION_STRING}")
        set(Python3_VERSION_MAJOR     "${PYTHON_VERSION_MAJOR}")
        set(Python3_VERSION_MINOR     "${PYTHON_VERSION_MINOR}")
        set(Python3_VERSION_PATCH     "${PYTHON_VERSION_PATCH}")
    endif()
endmacro()
