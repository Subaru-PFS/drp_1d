#-*-cmake-*-


# Force build type to release if not defined
if(NOT CMAKE_BUILD_TYPE AND NOT CODE_COVERAGE)
  set(CMAKE_BUILD_TYPE Release)
elseif(CODE_COVERAGE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

message( STATUS "CMAKE_BUILD_TYPE = " ${CMAKE_BUILD_TYPE} )

#--------------------------------------------------------
# Compiler options
#--------------------------------------------------------

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG_BUILD")

add_definitions("-std=c++11 -Wno-deprecated -Werror=return-type")

if(${CMAKE_BUILD_TYPE}  STREQUAL "Release" AND NOT CODE_COVERAGE)
  message( STATUS "CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS_RELEASE} )
  set(EXT "")
else()
  message( STATUS "CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS_DEBUG} )
  set(EXT "-debug")
endif()
