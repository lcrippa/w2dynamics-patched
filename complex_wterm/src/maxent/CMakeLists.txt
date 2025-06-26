cmake_minimum_required(VERSION 3.9...3.31)

add_library(MAXENT SHARED
            MaximumEntropy.F90
	    MersenneTwister.F90)
set_target_properties(MAXENT PROPERTIES
                      LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/w2dyn/maxent")
target_link_libraries(MAXENT PRIVATE
                      ${BLAS_LIBRARIES}
                      ${LAPACK_LIBRARIES}
                      ${CMAKE_THREAD_LIBS_INIT}
)
