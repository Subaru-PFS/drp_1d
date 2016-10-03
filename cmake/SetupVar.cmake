#-*-cmake-*-

# CACHE variables macro
macro( typed_cache_set type doc var )
    set ( ${var} ${ARGN} CACHE ${type} ${doc} FORCE )
endmacro()
#--------------------------------------------------------
# Compiler options
#--------------------------------------------------------
# Set a default build type if none is given
if ( NOT CMAKE_BUILD_TYPE ) # Debug default
  typed_cache_set ( STRING "Build type: Release or Debug" CMAKE_BUILD_TYPE Release   )
endif()

SET( CMAKE_CXX_FLAGS_RELEASE "-O3")

SET( CMAKE_CXX_FLAGS_DEBUG " -g -DDEBUG_BUILD")

add_definitions("-std=c++11 -Wno-deprecated -Werror=return-type")

MESSAGE (STATUS "CMAKE_BUILD_TYPE =" ${CMAKE_BUILD_TYPE} )
SET(EXT "-debug")
IF (${CMAKE_BUILD_TYPE}  STREQUAL "Release")
  MESSAGE( STATUS "CMAKE_CXX_FLAGS  = " ${CMAKE_CXX_FLAGS_RELEASE} )
  SET(EXT "")
ELSE()
  MESSAGE( STATUS "CMAKE_CXX_FLAGS  = " ${CMAKE_CXX_FLAGS_DEBUG} )
ENDIF()
