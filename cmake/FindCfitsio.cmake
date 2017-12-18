# -*- cmake -*-
#
# FindCfitsio.cmake
#
# Created on: Aug 19, 2013
#     Author: Pavel Binko
#
# -- Check for CFITSIO library --
#
# Defines:
#
#  CFITSIO_FOUND - system has the CFITSIO library
#  CFITSIO_INCLUDE_DIR - the CFITSIO include directory
#  CFITSIO_LIBRARY - The library needed to use CFITSIO
#
# Hints for non-standard locations:
#
#   CFITSIO_ROOT_DIR (used for finding cfitsio)
#
# Source of information :
#
#   https://github.com/wangd/SciDB-HDF5/blob/master/cmake/Modules/FindCfitsio.cmake
#

if(NOT CFITSIO_FOUND)

  FIND_PATH(CFITSIO_INCLUDE_DIR fitsio.h
            HINTS ENV CFITSIO_INSTALL_DIR
            PATH_SUFFIXES include cfitsio include/cfitsio
          )

  FIND_LIBRARY(CFITSIO_LIBRARY cfitsio
               HINTS ENV CFITSIO_INSTALL_DIR
               PATH_SUFFIXES lib
            )


  set(CFITSIO_LIBRARIES ${CFITSIO_LIBRARY})
  set(CFITSIO_INCLUDE_DIRS ${CFITSIO_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set CFITSIO_FOUND to TRUE if
# all listed variables are TRUE
  INCLUDE(FindPackageHandleStandardArgs)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(Cfitsio DEFAULT_MSG CFITSIO_INCLUDE_DIRS CFITSIO_LIBRARIES)


  if(NOT CMAKE_CROSSCOMPILING)
    INCLUDE(CheckCXXSourceRuns)
    set(CMAKE_REQUIRED_INCLUDES ${CFITSIO_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_LIBRARIES ${CFITSIO_LIBRARIES})
    CHECK_CXX_SOURCE_RUNS("
#include <fitsio.h>

int main(int, char**) {

  int answer=fits_is_reentrant();
  if (answer) {
    return 0;
  } else {
    return 1;
  }


}"
    CFITSIO_IS_REENTRANT)
    mark_as_advanced(CFITSIO_IS_REENTRANT)
  endif()

  mark_as_advanced(CFITSIO_FOUND CFITSIO_INCLUDE_DIRS CFITSIO_LIBRARIES)

endif()