#-*-cmake-*-


# Force build type to release if not defined
if(NOT CMAKE_BUILD_TYPE AND NOT CODE_COVERAGE)
  set(CMAKE_BUILD_TYPE Release)
elseif(CODE_COVERAGE)
  set(CMAKE_BUILD_TYPE Coverage)
endif()

message( STATUS "CMAKE_BUILD_TYPE = " ${CMAKE_BUILD_TYPE} )

#--------------------------------------------------------
# Compiler options
#--------------------------------------------------------

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG_BUILD")

set(CMAKE_CXX_FLAGS_COVERAGE "-O0 -g -fprofile-arcs -ftest-coverage --coverage")

mark_as_advanced(CMAKE_CXX_FLAGS_COVERAGE)

add_definitions("-std=c++11 -Wno-deprecated -Werror=return-type")

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  message( STATUS "CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS_RELEASE} )
  set(EXT "")
elseif(CMAKE_BUILD_TYPE STREQUAL "Coverage")
  message( STATUS "CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS_COVERAGE} )
  set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE ON)
  set(EXT "-coverage")
else()
  message( STATUS "CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS_DEBUG} )
  set(EXT "-debug")
endif()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-sign-compare -Wno-parentheses")
