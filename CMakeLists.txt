CMAKE_MINIMUM_REQUIRED (VERSION 2.8)
PROJECT( cpf-redshift )

# extra path for modules
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
include(SetupVar)
include(gitrevision)

# Force build type to release if not defined
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

message("Using git revision : " ${GIT_REVISION})

IF(BUILD_SHARED_LIBS)
  SET(Boost_USE_STATIC_LIBS OFF)
ELSE()
  SET(Boost_USE_STATIC_LIBS ON)
ENDIF()

FIND_PACKAGE(Boost COMPONENTS filesystem system thread timer chrono program_options unit_test_framework)
FIND_PACKAGE(Cfitsio)
FIND_PACKAGE(GSL 2.1)
FIND_PACKAGE(FFTW)

IF(NOT GSL_FOUND)
  # older versions of cmake don't have FindGSL
  FIND_PACKAGE(PkgConfig)
  pkg_check_modules(GSL gsl>=2.1)
ENDIF()

#--------------------------------------------------------
# Setup shortcut access directory
#--------------------------------------------------------
SET( ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}" )

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/usr" CACHE PATH "" FORCE)
endif()

SET( cpf-redshift_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include )
SET( cpf-redshift_LINK_DIR ${CMAKE_INSTALL_PREFIX}/lib )

# Set up THIRDPARTY_LIBS
IF (Boost_FOUND)
  IF (CMAKE_BUILD_TYPE STREQUAL "Release")
    SET(BOOST_TARGET "RELEASE")
  ELSE()
    SET(BOOST_TARGET "DEBUG")
  ENDIF()

  SET(boost_libs ${Boost_FILESYSTEM_LIBRARY_${BOOST_TARGET}}
    ${Boost_SYSTEM_LIBRARY_${BOOST_TARGET}}
    ${Boost_THREAD_LIBRARY_${BOOST_TARGET}}
    ${Boost_TIMER_LIBRARY_${BOOST_TARGET}}
    ${Boost_CHRONO_LIBRARY_${BOOST_TARGET}}
    ${Boost_PROGRAM_OPTIONS_LIBRARY_${BOOST_TARGET}}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_${BOOST_TARGET}})
  LIST(APPEND cpf-redshift_THIRDPARTY_LIBS ${boost_libs})
  LIST(APPEND cpf-redshift_INCLUDE_DIR ${Boost_INCLUDE_DIRS})
  LIST(APPEND cpf-redshift_LINK_DIR ${Boost_LIBRARY_DIRS})
ELSE()
  LIST(APPEND MISSING_THIRDPARTIES "boost")
  LIST(APPEND cpf-redshift_THIRDPARTY_LIBS pthread boost_filesystem
    boost_system boost_thread boost_timer boost_chrono boost_program_options)
ENDIF()

IF (CFITSIO_FOUND)
  LIST(APPEND cpf-redshift_THIRDPARTY_LIBS ${CFITSIO_LIBRARIES})
  LIST(APPEND cpf-redshift_INCLUDE_DIR ${CFITSIO_INCLUDE_DIRS})
ELSE()
  LIST(APPEND MISSING_THIRDPARTIES "cfitsio")
  LIST(APPEND cpf-redshift_THIRDPARTY_LIBS cfitsio)
ENDIF()

IF (FFTW_FOUND)
  LIST(APPEND cpf-redshift_THIRDPARTY_LIBS ${FFTW_LIBRARIES})
  LIST(APPEND cpf-redshift_INCLUDE_DIR ${FFTW_INCLUDE_DIR})
ELSE()
  LIST(APPEND MISSING_THIRDPARTIES "fftw")
  LIST(APPEND cpf-redshift_THIRDPARTY_LIBS fftw3)
ENDIF()

IF (GSL_FOUND)
  LIST(APPEND cpf-redshift_THIRDPARTY_LIBS ${GSL_LIBRARIES})
  LIST(APPEND cpf-redshift_INCLUDE_DIR ${GSL_INCLUDE_DIR})
ELSE()
  LIST(APPEND MISSING_THIRDPARTIES "gsl")
  LIST(APPEND cpf-redshift_THIRDPARTY_LIBS gsl gslcblas)
ENDIF()

# use ccache if any
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
  set(ENV{CC} "ccache ${CMAKE_C_COMPILER}")
  set(ENV{CXX} "ccache ${CMAKE_CXX_COMPILER}")
endif(CCACHE_FOUND)

#--------------------------------------------------------
# Parse version, this will be used later to configure
# cpf-core-redshift-version.cmake with the appropriate version number
#--------------------------------------------------------

FILE (STRINGS ${ROOT_DIR}/VERSION VERSION)
MESSAGE( STATUS "Version used: ${VERSION}" )

#SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wno-sign-compare -Wno-parentheses")
IF(BUILD_SHARED_LIBS)
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto")
ENDIF()
SET(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -g")

IF(TEST_COVERAGE)
  SET(GCC_COVERAGE_COMPILE_FLAGS "--coverage")
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
  SET(cpf-redshift_COVERAGE "True")
ELSE()
  SET(cpf-redshift_COVERAGE "False")
ENDIF()

IF(GPROF)
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Og -pg --coverage")
ENDIF()

#--------------------------------------------------------
# Preprocessor and linker path
#--------------------------------------------------------

INCLUDE_DIRECTORIES( ${ROOT_DIR}/RedshiftLibrary ${cpf-redshift_INCLUDE_DIR} )
LINK_DIRECTORIES( ${cpf-redshift_LINK_DIR} )

#--------------------------------------------------------
# Exported configuration variables
#--------------------------------------------------------

IF( NOT (CFITSIO_FOUND AND FFTW_FOUND AND GSL_FOUND AND Boost_FOUND) )
  MESSAGE( STATUS "Some library not found on system. Using thirdparty/" )
  set(PARALLEL_BUILD_FLAG "-j4")
  ADD_SUBDIRECTORY( tools )
ELSE()
  MESSAGE( STATUS "Using system third-party libraries." )
ENDIF()

MESSAGE("Using third party libraries: " "${cpf-redshift_THIRDPARTY_LIBS}")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#-------------------------  -------------------------------
# Extract version number
#--------------------------------------------------------
FILE (STRINGS ${ROOT_DIR}/VERSION VERSION_NUMBER)

#--------------------------------------------------------
# Compile library
#--------------------------------------------------------
SET( LIB_NAME cpf-redshift )

FILE( GLOB_RECURSE LIB_SOURCE_FILES ${ROOT_DIR}/RedshiftLibrary/src/lib/*.cpp )

# install with rpath to find third-parties
SET(CMAKE_MACOSX_RPATH 1)
SET(CMAKE_SKIP_BUILD_RPATH FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY( ${LIB_NAME} SHARED ${LIB_SOURCE_FILES} )
  TARGET_LINK_LIBRARIES(${LIB_NAME} ${cpf-redshift_THIRDPARTY_LIBS})
  INSTALL( TARGETS ${LIB_NAME}
    LIBRARY DESTINATION lib)
ELSE()
  ADD_LIBRARY( ${LIB_NAME}-static STATIC ${LIB_SOURCE_FILES} )
  SET_TARGET_PROPERTIES(${LIB_NAME}-static PROPERTIES OUTPUT_NAME ${LIB_NAME})
  TARGET_LINK_LIBRARIES(${LIB_NAME}-static ${cpf-redshift_THIRDPARTY_LIBS})
  INSTALL( TARGETS ${LIB_NAME}-static
    ARCHIVE DESTINATION lib)
  SET(LIB_NAME ${LIB_NAME}-static)
ENDIF()

#--------------------------------------------------------
# Compile test
#--------------------------------------------------------
IF(BUILD_TESTS AND BUILD_SHARED_LIBS)
ENABLE_TESTING()

FILE( GLOB_RECURSE TEST_SOURCE_FILES ${ROOT_DIR}/RedshiftLibrary/tests/src/*_test.cpp )
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK)

CONFIGURE_FILE(${ROOT_DIR}/RedshiftLibrary/tests/src/test-config.h.in
  ${ROOT_DIR}/RedshiftLibrary/tests/src/test-config.h)

FOREACH(TEST_SOURCE ${TEST_SOURCE_FILES})
  GET_FILENAME_COMPONENT(TEST_NAME ${TEST_SOURCE} NAME_WE)
  ADD_EXECUTABLE( ${TEST_NAME} ${TEST_SOURCE} )
  TARGET_LINK_LIBRARIES( ${TEST_NAME} ${LIB_NAME} ${cpf-redshift_THIRDPARTY_LIBS} ${tests_THIRDPARTY_LIBS})
  ADD_TEST(NAME ${TEST_NAME} COMMAND ${TEST_NAME} )
  SET_PROPERTY(SOURCE ${TEST_SOURCE} APPEND PROPERTY COMPILE_DEFINITIONS BOOST_TEST_MODULE=${TEST_NAME} )
ENDFOREACH(TEST_SOURCE)
ENDIF()

#--------------------------------------------------------
# Install directive
#--------------------------------------------------------

INSTALL( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/RedshiftLibrary/RedshiftLibrary DESTINATION include )
INSTALL( FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/cpf-redshift-config.cmake DESTINATION share/cmake )
INSTALL( FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/cpf-redshift-config-version.cmake DESTINATION share/cmake )

#--------------------------------------------------------
# Build python wrapper
#--------------------------------------------------------

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/swig_ext.py.in
  ${CMAKE_CURRENT_SOURCE_DIR}/swig_ext.py
  )
SET (GENERATED_FILES_LIST "swig_ext.py" "pyamazed/redshift_wrap.cpp"
  "pyamazed/redshift_wrap.h" "pyamazed/redshift.py" "pyamazed.egg-info")

SET_DIRECTORY_PROPERTIES(PROPERTIES
  ADDITIONAL_MAKE_CLEAN_FILES "${GENERATED_FILES_LIST}")

#--------------------------------------------------------
# Package deployement and config
#--------------------------------------------------------

INCLUDE ( CMakePackageConfigHelpers )

# Configure config file (this is used by find_package to find where the library is)
CONFIGURE_PACKAGE_CONFIG_FILE(
       ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cpf-redshift-config.cmake.in
       ${CMAKE_CURRENT_BINARY_DIR}/cmake/cpf-redshift-config.cmake
       INSTALL_DESTINATION share/cmake
       PATH_VARS cpf-redshift_INCLUDE_DIR cpf-redshift_LINK_DIR
)

# Configure version file (this is used by find_package to determine a matching version)
CONFIGURE_FILE(
       ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cpf-redshift-config-version.cmake.in
       ${CMAKE_CURRENT_BINARY_DIR}/cmake/cpf-redshift-config-version.cmake
       @ONLY
)

CONFIGURE_FILE(
	${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.h.in
	${CMAKE_CURRENT_SOURCE_DIR}/RedshiftLibrary/RedshiftLibrary/version.h
)

# Generate TAGS file if etags is installed
find_program(ETAGS NAMES etags ctags)
if(ETAGS)
  execute_process(COMMAND "find ${ROOT_DIR}/RedshiftLibrary/ ( -name *.cpp -o -name *.h -o -name *.hpp -o -name *.py -o -name *.c ) -exec ${ETAGS} -a ${ROOT_DIR}/TAGS {} +"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif(ETAGS)
