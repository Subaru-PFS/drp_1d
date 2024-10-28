CMAKE_MINIMUM_REQUIRED (VERSION 3.15)
PROJECT( cpf-redshift )

# Dependencies requirements
set(BOOST_MINIMUM_REQUIRED 1.74)
set(CFISTIO_MINIMUM_REQUIRED 3.36)
set(FFTW_MINIMUM_REQUIRED 3.3.8)
set(OPENBLASS_MINIMUM_REQUIRED 0.3)
set(GSL_MINIMUM_REQUIRED 2.5)

# Environment settings
set( ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}" )
message(STATUS "ROOT_DIR = " ${ROOT_DIR})
message(STATUS "CMAKE_CURRENT_SOURCE_DIR = " ${CMAKE_CURRENT_SOURCE_DIR})
message(STATUS "CMAKE_BINARY_DIR = " ${CMAKE_BINARY_DIR})
MESSAGE( STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )

# Regular cmake tools
FIND_PACKAGE(Python3 COMPONENTS Interpreter REQUIRED)
set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
FIND_PACKAGE(PkgConfig REQUIRED)

# package required to build python wrapper
FIND_PACKAGE(SWIG REQUIRED)
FIND_PACKAGE(Python REQUIRED Interpreter Development.Module)

# Install directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local/pylibamazed" CACHE PATH "" FORCE)
endif()
message(STATUS "CMAKE_INSTALL_PREFIX=" ${CMAKE_INSTALL_PREFIX})

# Extra path for modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# Define build flags
include(SetupVar)

# Get library revision
include(gitrevision)
message(STATUS "Using git revision : " ${GIT_REVISION})
FILE (STRINGS ${ROOT_DIR}/VERSION PROJECT_VERSION)
MESSAGE( STATUS "pylibamazed version used: ${PROJECT_VERSION}" )

# Building environment
LIST(APPEND cpf-redshift_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/RedshiftLibrary)
LIST(APPEND cpf-redshift_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include)
LIST(APPEND cpf-redshift_LINK_DIR ${CMAKE_INSTALL_PREFIX}/lib)
if(CMAKE_PREFIX_PATH)
  LIST(APPEND cpf-redshift_INCLUDE_DIR ${CMAKE_PREFIX_PATH}/include)
  LIST(APPEND cpf-redshift_LINK_DIR ${CMAKE_PREFIX_PATH}/lib)
endif()
MESSAGE(STATUS "cpf-redshift_INCLUDE_DIR=" "${cpf-redshift_INCLUDE_DIR}")
MESSAGE(STATUS "cpf-redshift_LINK_DIR=" "${cpf-redshift_LINK_DIR}")


#--------------------------------------------------------
# Setup thirdparties
#--------------------------------------------------------

# install with rpath to find third-parties
SET(CMAKE_MACOSX_RPATH 1)
IF(CMAKE_PREFIX_PATH)
  LIST(APPEND CMAKE_INSTALL_RPATH ${CMAKE_PREFIX_PATH}/lib)
ENDIF()

# Add specific PREFIX PATH for users on Appel osx managed by brew
list(APPEND CMAKE_PREFIX_PATH /usr/local/opt/openblas/lib)
LIST(APPEND CMAKE_INSTALL_RPATH /usr/local/opt/openblas/lib)

message(STATUS "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}")
message(STATUS "CMAKE_INSTALL_RPATH =  ${CMAKE_INSTALL_RPATH}")

# Setup thirdparty directory : Boost
find_package(Boost ${BOOST_MINIMUM_REQUIRED} COMPONENTS filesystem system thread timer chrono program_options unit_test_framework   REQUIRED)
set(Boost_USE_STATIC_LIBS OFF)
SET(BOOST_TARGET "RELEASE")
SET(boost_libs ${Boost_FILESYSTEM_LIBRARY_${BOOST_TARGET}}
  ${Boost_SYSTEM_LIBRARY_${BOOST_TARGET}}
  ${Boost_THREAD_LIBRARY_${BOOST_TARGET}}
  ${Boost_TIMER_LIBRARY_${BOOST_TARGET}}
  ${Boost_CHRONO_LIBRARY_${BOOST_TARGET}}
  ${Boost_PROGRAM_OPTIONS_LIBRARY_${BOOST_TARGET}}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_${BOOST_TARGET}})
message(STATUS "Boost_FILESYSTEM_LIBRARY_${BOOST_TARGET} = ${Boost_FILESYSTEM_LIBRARY_${BOOST_TARGET}}")
message(STATUS "Boost_SYSTEM_LIBRARY_${BOOST_TARGET} = ${Boost_SYSTEM_LIBRARY_${BOOST_TARGET}}")
message(STATUS "Boost_THREAD_LIBRARY_${BOOST_TARGET} = ${Boost_THREAD_LIBRARY_${BOOST_TARGET}}")
message(STATUS "Boost_TIMER_LIBRARY_${BOOST_TARGET} = ${Boost_TIMER_LIBRARY_${BOOST_TARGET}}")
message(STATUS "Boost_CHRONO_LIBRARY_${BOOST_TARGET} = ${Boost_CHRONO_LIBRARY_${BOOST_TARGET}}")
message(STATUS "Boost_PROGRAM_OPTIONS_LIBRARY_${BOOST_TARGET} = ${Boost_PROGRAM_OPTIONS_LIBRARY_${BOOST_TARGET}}")
message(STATUS "Boost_UNIT_TEST_FRAMEWORK_LIBRARY_${BOOST_TARGET} = ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_${BOOST_TARGET}}")
message(STATUS "Boost_LIBRARY_DIRS = ${Boost_LIBRARY_DIRS}")
message(STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")
LIST(APPEND cpf-redshift_THIRDPARTY_LIBS ${boost_libs})
LIST(APPEND cpf-redshift_LINK_DIR ${Boost_LIBRARY_DIRS})
LIST(APPEND cpf-redshift_INCLUDE_DIR ${Boost_INCLUDE_DIRS})

# Setup thirdparty directory : CFITSIO
pkg_check_modules(CFITSIO cfitsio>=${CFISTIO_MINIMUM_REQUIRED} REQUIRED)
message(STATUS "Found cfitsio: ${CFITSIO_LINK_LIBRARIES}")
message(STATUS "CFITSIO_LIBRARIES = ${CFITSIO_LIBRARIES}")
message(STATUS "CFITSIO_LIBRARY_DIRS = ${CFITSIO_LIBRARY_DIRS}")
message(STATUS "CFITSIO_INCLUDE_DIRS = ${CFITSIO_INCLUDE_DIRS}")
LIST(APPEND cpf-redshift_THIRDPARTY_LIBS ${CFITSIO_LIBRARIES})
LIST(APPEND cpf-redshift_LINK_DIR ${CFITSIO_LIBRARY_DIRS})
LIST(APPEND cpf-redshift_INCLUDE_DIR ${CFITSIO_INCLUDE_DIRS})

# Setup thirdparty directory : FFTW
pkg_check_modules(FFTW fftw3>=${FFTW_MINIMUM_REQUIRED} REQUIRED)
message(STATUS "Found fftw3: ${FFTW_LINK_LIBRARIES}")
message(STATUS "FFTW_LIBRARIES = ${FFTW_LIBRARIES}")
message(STATUS "FFTW_LIBRARY_DIRS = ${FFTW_LIBRARY_DIRS}")
message(STATUS "FFTW_INCLUDE_DIRS = ${FFTW_INCLUDE_DIRS}")
LIST(APPEND cpf-redshift_THIRDPARTY_LIBS ${FFTW_LIBRARIES})
LIST(APPEND cpf-redshift_LINK_DIR ${FFTW_LIBRARY_DIRS})
LIST(APPEND cpf-redshift_INCLUDE_DIR ${FFTW_INCLUDE_DIRS})

# Setup thirdparty directory : OpenBLAS
set(BLA_VENDOR OpenBLAS)
set(BLA_STATIC OFF)
set(BLA_F95 OFF)
find_package(BLAS REQUIRED)
message(STATUS "BLAS_LIBRARIES = ${BLAS_LIBRARIES}")
LIST(APPEND cpf-redshift_THIRDPARTY_LIBS ${BLAS_LIBRARIES})

# Setup thirdparty directory : GSL
pkg_check_modules(GSL gsl>=${GSL_MINIMUM_REQUIRED} REQUIRED)
message(STATUS "Found gsl: ${GSL_LINK_LIBRARIES}")
message(STATUS "GSL_LIBRARY_DIRS = ${GSL_LIBRARY_DIRS}")
message(STATUS "GSL_INCLUDE_DIRS = ${GSL_INCLUDE_DIRS}")
LIST(APPEND cpf-redshift_THIRDPARTY_LIBS gsl)
LIST(APPEND cpf-redshift_LINK_DIR ${GSL_LIBRARY_DIRS})
LIST(APPEND cpf-redshift_INCLUDE_DIR ${GSL_INCLUDE_DIRS})

# setup thirdparty directory : eigen3
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
message(STATUS "Found eigen")
message(STATUS "EIGEN3_INCLUDE_DIRS = ${EIGEN3_INCLUDE_DIRS}")
LIST(APPEND cpf-redshift_THIRDPARTY_HEADERONLY Eigen3::Eigen)

# setup thirdparty directory : lbfgspp
find_package (lbfgspp 0.3.0 REQUIRED NO_MODULE)
get_target_property(LBFGSPP_INCLUDE_DIRS lbfgspp INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "Found lbfgspp:")
message(STATUS "LBFGSPP_INCLUDE_DIRS = ${LBFGSPP_INCLUDE_DIRS}")
LIST(APPEND cpf-redshift_THIRDPARTY_HEADERONLY lbfgspp)

# Thirdparties setup summary
MESSAGE("Using third party libraries: " "${cpf-redshift_THIRDPARTY_LIBS}")
MESSAGE("Using third party header only libraries: " "${cpf-redshift_THIRDPARTY_HEADERONLY}")

#--------------------------------------------------------
# Setting ccache
#--------------------------------------------------------

find_program(CCACHE_PROGRAM ccache)
if (CCACHE_PROGRAM)
    message(NOTICE "Selected CCache: ${CCACHE_PROGRAM}")
    MESSAGE (STATUS "CCACHE_BASEDIR=${PROJECT_SOURCE_DIR}")
    MESSAGE (STATUS "CCACHE_DIR=$ENV{CCACHE_DIR}")
    set(ccacheEnv
            CCACHE_BASEDIR=${PROJECT_SOURCE_DIR}
            CCACHE_DIR=$ENV{CCACHE_DIR}
            CCACHE_SLOPPINESS="pch_defines,file_macro,include_file_ctime"
            CCACHE_NOHASHDIR=true
    )
    set(CMAKE_C_COMPILER_LAUNCHER ${CMAKE_COMMAND} -E env ${ccacheEnv} ${CCACHE_PROGRAM})
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CMAKE_COMMAND} -E env ${ccacheEnv} ${CCACHE_PROGRAM})
endif ()

#--------------------------------------------------------
# Preprocessor and linker path
#--------------------------------------------------------

list(REMOVE_DUPLICATES cpf-redshift_INCLUDE_DIR)
INCLUDE_DIRECTORIES( ${cpf-redshift_INCLUDE_DIR} )
list(REMOVE_DUPLICATES cpf-redshift_LINK_DIR)
LINK_DIRECTORIES( ${cpf-redshift_LINK_DIR} )

#--------------------------------------------------------
# Compile library
#--------------------------------------------------------
SET( LIB_NAME cpf-redshift )

FILE( GLOB_RECURSE LIB_SOURCE_FILES ${ROOT_DIR}/RedshiftLibrary/src/lib/*.cpp )

ADD_DEFINITIONS(-DLBFGSBFITTER)

ADD_LIBRARY( ${LIB_NAME} SHARED ${LIB_SOURCE_FILES} )
TARGET_LINK_LIBRARIES(${LIB_NAME} ${cpf-redshift_THIRDPARTY_LIBS} ${cpf-redshift_THIRDPARTY_HEADERONLY})
INSTALL( TARGETS ${LIB_NAME} LIBRARY DESTINATION lib)

#--------------------------------------------------------
# Compile test
#--------------------------------------------------------
if(BUILD_TESTING)
  ENABLE_TESTING()

  FILE( GLOB_RECURSE TEST_SOURCE_FILES ${ROOT_DIR}/RedshiftLibrary/tests/src/*_test.cpp )
  ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK)

  CONFIGURE_FILE(${ROOT_DIR}/RedshiftLibrary/tests/src/test-config.h.in
    ${ROOT_DIR}/RedshiftLibrary/tests/src/test-config.h)

  FOREACH(TEST_SOURCE ${TEST_SOURCE_FILES})
    GET_FILENAME_COMPONENT(TEST_NAME ${TEST_SOURCE} NAME_WE)
    ADD_EXECUTABLE( ${TEST_NAME} ${TEST_SOURCE} )
    TARGET_LINK_LIBRARIES( ${TEST_NAME} ${LIB_NAME} ${cpf-redshift_THIRDPARTY_LIBS} ${CMAKE_DL_LIBS})
    ADD_TEST(NAME ${TEST_NAME} COMMAND ${TEST_NAME} )
    SET_PROPERTY(SOURCE ${TEST_SOURCE} APPEND PROPERTY COMPILE_DEFINITIONS BOOST_TEST_MODULE=${TEST_NAME} )
  ENDFOREACH(TEST_SOURCE)
ENDIF()

#--------------------------------------------------------
# Generate coverage data : HTML and XML output
#--------------------------------------------------------
if(BUILD_TYPE STREQUAL "coverage")

  find_program(GCOVR_EXECUTABLE gcovr REQUIRED)

  # Create a coverage report in xml format
  set(GCOVR_XML_COMMAND ${GCOVR_EXECUTABLE} --xml
      --filter ../RedshiftLibrary/src/lib/ --filter ../RedshiftLibrary/RedshiftLibrary/
      --exclude-unreachable-branches --print-summary -o ${CMAKE_CURRENT_SOURCE_DIR}/coverage/coverage.xml
      --root ${CMAKE_BINARY_DIR}
  )

  # Create a coverage report in html format
  set(GCOVR_HTML_COMMAND ${GCOVR_EXECUTABLE} --html --html-details
      --filter ../RedshiftLibrary/src/lib/ --filter ../RedshiftLibrary/RedshiftLibrary/
      --exclude-unreachable-branches --print-summary -o ${CMAKE_CURRENT_SOURCE_DIR}/coverage/coverage.html
      --root ${CMAKE_BINARY_DIR}
  )

  # Target to run coverage tests with 'make coverage_xml' command
  add_custom_target(coverage_xml
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/coverage
      COMMAND ${CMAKE_MAKE_PROGRAM} test
      COMMAND ${GCOVR_XML_COMMAND}
      COMMENT "Coverage files have been output to ${CMAKE_CURRENT_SOURCE_DIR}/coverage"
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      VERBATIM
  )

  # Target to run coverage tests with 'make coverage_html' command
  add_custom_target(coverage_html
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/coverage
      COMMAND ${CMAKE_MAKE_PROGRAM} test
      COMMAND ${GCOVR_HTML_COMMAND}
      COMMENT "Coverage files have been output to ${CMAKE_CURRENT_SOURCE_DIR}/coverage"
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      VERBATIM
  )

  set(COVERAGE_REPORT_FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/coverage/coverage.html"
  "${CMAKE_CURRENT_SOURCE_DIR}/coverage/coverage.xml"
  )

  set_property(DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${COVERAGE_REPORT_FILES})

  # Target to clean the coverage repository as well as generated .gcda and .gcno files
  add_custom_target(clean_coverage
      COMMAND rm -rf ${COVERAGE_REPORT_FILES}
      COMMAND find . -name "*.gcda" -exec rm -f {} \;
      COMMAND find . -name "*.gcno" -exec rm -f {} \;
      COMMENT "Removing generated coverage files"
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      VERBATIM
  )
endif()

#--------------------------------------------------------
# Install directive
#--------------------------------------------------------

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/pylibamazed/auxdir/pylibamazed/ DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/pylibamazed/python/pylibamazed/resources)


#--------------------------------------------------------
# Build python wrapper
#--------------------------------------------------------
IF (SKBUILD_PROJECT_NAME)
  MESSAGE(STATUS "==== python wrapper ====")
  MESSAGE(STATUS "SWIG_EXECUTABLE =" ${SWIG_EXECUTABLE})
  MESSAGE(STATUS "Python_EXECUTABLE =" ${Python_EXECUTABLE})
  MESSAGE(STATUS "SKBUILD_PROJECT_NAME =" ${SKBUILD_PROJECT_NAME})

  MESSAGE(STATUS "CMAKE_INSTALL_RPATH=" ${CMAKE_INSTALL_RPATH})
  # Grab the variables from a local Python installation F2PY headers
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -c
            "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  MESSAGE(STATUS "NUMPY_INCLUDE_DIR =" ${NUMPY_INCLUDE_DIR})

  # Generate the wrapper files
  add_custom_command(
    OUTPUT redshift.py redshift_wrap.cpp
    COMMAND ${SWIG_EXECUTABLE} -c++ -python -I${CMAKE_CURRENT_SOURCE_DIR}/RedshiftLibrary
            -outdir ${CMAKE_CURRENT_BINARY_DIR}
            -o ${CMAKE_CURRENT_BINARY_DIR}/redshift_wrap.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/RedshiftLibrary/RedshiftLibrary/redshift.i
    VERBATIM
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/RedshiftLibrary/RedshiftLibrary/redshift.i
  )

  python_add_library(
    _redshift
    MODULE
    WITH_SOABI
    ${CMAKE_CURRENT_BINARY_DIR}/redshift_wrap.cpp
  )

  target_include_directories(_redshift PUBLIC  ${NUMPY_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../src)
  
  MESSAGE (STATUS "cpf-redshift_LINK_DIR :" ${cpf-redshift_LINK_DIR})
  if (APPLE)
      set_target_properties(_redshift PROPERTIES INSTALL_RPATH "@loader_path")
  else()
      set_target_properties(_redshift PROPERTIES INSTALL_RPATH "$ORIGIN;${cpf-redshift_LINK_DIR}")
  endif()
  target_link_libraries(_redshift PUBLIC cpf-redshift ${cpf-redshift_THIRDPARTY_LIBS} ${cpf-redshift_THIRDPARTY_HEADERONLY} )

  install(TARGETS _redshift DESTINATION pylibamazed)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/redshift.py DESTINATION pylibamazed)
  install(FILES ${cpf-redshift_BINARY_DIR}/libcpf-redshift.so DESTINATION pylibamazed)
  #install(FILES ${cpf-redshift_BINARY_DIR}/libcpf-redshift.${LIBEXT} DESTINATION pylibamazed)
ENDIF (SKBUILD_PROJECT_NAME)
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


IF (NOT SKBUILD_PROJECT_NAME)
# Configure version file (this is used by find_package to determine a matching version)
WRITE_BASIC_PACKAGE_VERSION_FILE(
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/cpf-redshift-config-version.cmake
  VERSION ${VERSION}
  COMPATIBILITY SameMajorVersion
)

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/RedshiftLibrary/RedshiftLibrary/version.h.in
  ${CMAKE_CURRENT_SOURCE_DIR}/RedshiftLibrary/RedshiftLibrary/version.h
)
ENDIF (NOT SKBUILD_PROJECT_NAME)
