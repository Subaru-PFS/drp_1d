# - Try to find FFTW.
# Usage: find_package(FFTW [COMPONENTS [single double (long-double|long) threads (omp|openmp)]])
#
# Variables used by this module:
#  FFTW_ROOT_DIR                   - FFTW root directory
# Variables defined by this module:
#  FFTW_FOUND                      - system has FFTW
#  FFTW_INCLUDE_DIR                - the FFTW include directory (cached)
#  FFTW_INCLUDE_DIRS (not cached)  - the FFTW include directories
#                                    (identical to FFTW_INCLUDE_DIR)
#  FFTW[FL]?_LIBRARY               - the FFTW library - double, single(F),
#                                    long-double(L) precision (cached)
#  FFTW[FL]?_THREADS_LIBRARY       - the threaded FFTW library - double, single(F),
#                                    long-double(L) precision (cached)
#  FFTW[FL]?_OMP_LIBRARY           - the omped FFTW library - double, single(F),
#                                    long-double(L) precision (cached)
#  FFTW_LIBRARIES                  - list of all FFTW libraries found


if(NOT FFTW_FOUND)

  if(FFTW_FIND_COMPONENTS MATCHES "^$")
    set(_components double)
  else()
    set(_components ${FFTW_FIND_COMPONENTS})
  endif()

  # Loop over each component.
  set(_libraries)
  foreach(_comp ${_components})
    if(_comp STREQUAL "single")
      list(APPEND _libraries fftw3f)
    elseif(_comp STREQUAL "double")
      list(APPEND _libraries fftw3)
    elseif(_comp STREQUAL "long-double" OR _comp STREQUAL "long")
      list(APPEND _libraries fftw3l)
    elseif(_comp STREQUAL "quad")
      list(APPEND _libraries fftw3q)
    elseif(_comp STREQUAL "threads")
      set(_use_threads ON)
    elseif(_comp STREQUAL "omp" OR _comp STREQUAL "openmp")
      set(_use_omp ON)
    else()
      message(FATAL_ERROR "FindFFTW: unknown component `${_comp}' specified. "
        "Valid components are `single', `double', `long-double|long', and `threads|omp|openmp'.")
    endif()
  endforeach(_comp ${_components})

# If using threads, we need to link against threaded libraries as well.
  if(_use_threads)
    set(_thread_libs)
    foreach(_lib ${_libraries})
      list(APPEND _thread_libs ${_lib}_threads)
    endforeach()
    set(_libraries ${_thread_libs} ${_libraries})
  endif()


# If using openmp, we need to link against omped libraries as well.
  if(_use_omp)
    set(_omp_libs)
    foreach(_lib ${_libraries})
      list(APPEND _omp_libs ${_lib}_omp)
    endforeach()
    set(_libraries ${_omp_libs} ${_libraries})
  endif()


# Keep a list of variable names that we need to pass on to
# find_package_handle_standard_args().
  set(_check_list)

# Search for all requested libraries.
  foreach(_lib ${_libraries})
    string(TOUPPER ${_lib} _LIB)
    find_library(${_LIB}_LIBRARY ${_lib}
                 HINTS ENV FFTW_ROOT_DIR
                 PATH_SUFFIXES lib)
    mark_as_advanced(${_LIB}_LIBRARY)
    list(APPEND FFTW_LIBRARIES ${${_LIB}_LIBRARY})
    list(APPEND _check_list ${_LIB}_LIBRARY)
  endforeach()

# Search for the header file.
  find_path(FFTW_INCLUDE_DIR fftw3.h
            HINTS ENV FFTW_ROOT_DIR
            PATH_SUFFIXES include)
  mark_as_advanced(FFTW_INCLUDE_DIR)
  set(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
  list(APPEND _check_list FFTW_INCLUDE_DIR)


# Handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(FFTW DEFAULT_MSG ${_check_list})


  mark_as_advanced(FFTW_FOUND FFTW_INCLUDE_DIRS FFTW_LIBRARIES)

  list(REMOVE_DUPLICATES FFTW_LIBRARIES)
  list(REMOVE_DUPLICATES FFTW_INCLUDE_DIRS)

endif()
