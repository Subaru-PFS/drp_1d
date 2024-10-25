#-*-cmake-*-


# Force build type to release if not defined
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message( STATUS "CMAKE_BUILD_TYPE = " ${CMAKE_BUILD_TYPE} )

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)
message(STATUS "BUILD_SHARED_LIBS = " ${BUILD_SHARED_LIBS})

set(BUILD_STATIC_LIBS OFF)
message(STATUS "BUILD_STATIC_LIBS = " ${BUILD_STATIC_LIBS})

option(BUILD_TESTING "Build testings" ON)
message(STATUS "BUILD_TESTING = " ${BUILD_TESTING})

option(COVERAGE "Test code coverage" OFF)
message(STATUS "COVERAGE = " ${COVERAGE})

#--------------------------------------------------------
# Compiler options
#--------------------------------------------------------

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -D_FORTIFY_SOURCE=2")

set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG_BUILD -Wl,--no-as-needed")

# set(CMAKE_CXX_FLAGS_COVERAGE "-O0 -g -fprofile-arcs -ftest-coverage --coverage")
set(CMAKE_CXX_FLAGS_COVERAGE "-g -O0 --coverage")

mark_as_advanced(CMAKE_CXX_FLAGS_COVERAGE)

add_definitions("-std=c++17 -Wno-deprecated -Werror=return-type")

string(TOLOWER ${CMAKE_BUILD_TYPE} BUILD_TYPE)

if(BUILD_TYPE STREQUAL "debug")
  message( STATUS "CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS_DEBUG} )
  set(EXT "-debug")
elseif(BUILD_TYPE STREQUAL "coverage")
  message( STATUS "CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS_COVERAGE} )
  set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE ON)
  set(EXT "-coverage")
else()
  message( STATUS "CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS_RELEASE} )
  set(EXT "")
endif()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-sign-compare -Wno-parentheses")
