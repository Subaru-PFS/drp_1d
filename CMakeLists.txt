CMAKE_MINIMUM_REQUIRED (VERSION 2.8)
PROJECT( cpf-redshift )

# extra path for modules
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
include(SetupVar)
include(gitrevision)

message("Using git revision : " ${GIT_REVISION})

FIND_PACKAGE(Boost COMPONENTS filesystem system thread timer chrono program_options unit_test_framework)
FIND_PACKAGE(Cfitsio)
FIND_PACKAGE(GSL)
FIND_PACKAGE(FFTW)

# use ccache if any
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
  set(ENV{CC} "ccache cc")
  set(ENV{CXX} "ccache c++")
endif(CCACHE_FOUND)

#--------------------------------------------------------
# Setup shortcut acces directory
#--------------------------------------------------------
SET( ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}" )

SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${ROOT_DIR}/lib )
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${ROOT_DIR}/lib )
SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${ROOT_DIR}/bin )

#--------------------------------------------------------
# Parse version, this willbe used later to configure
# cpf-core-redshift-version.cmake with the appropriate version number
#--------------------------------------------------------
FILE (STRINGS ${ROOT_DIR}/VERSION VERSION)
MESSAGE( STATUS "Version used: ${VERSION}" )

#SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wno-sign-compare -Wno-parentheses")

IF(TEST_COVERAGE)
  SET(GCC_COVERAGE_COMPILE_FLAGS "--coverage")
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
  SET(cpf-redshift_COVERAGE "True")
ELSE()
  SET(cpf-redshift_COVERAGE "False")
ENDIF()

#--------------------------------------------------------
# Exported configuration variables
#--------------------------------------------------------
SET( cpf-redshift_INCLUDE_DIR ${ROOT_DIR}/RedshiftLibrary/ )
SET( cpf-redshift_LINK_DIR ${ROOT_DIR}/lib/ )

IF( NOT (CFITSIO_FOUND AND FFTW_FOUND AND GSL_FOUND AND BOOST_FOUND) )
  MESSAGE( STATUS "Some library not found on system. Using thirdparty/" )
  SET( cpf-redshift_THIRDPARTY_INCLUDE_DIR ${ROOT_DIR}/thirdparty/include/ )
  SET( cpf-redshift_THIRDPARTY_LINK_DIR ${ROOT_DIR}/thirdparty/lib/)
  SET( cpf-redshift_THIRDPARTY_LIBS
    cfitsio fftw3 gsl gslcblas pthread boost_filesystem boost_system boost_thread boost_timer
    boost_chrono boost_program_options boost_regex)
  SET( cpf-redshift_THIRDPARTY_STATIC_LIBS
    libcfitsio.a libfftw3.a libgsl.a libgslcblas.a libboost_filesystem.a
    libboost_system.a libboost_thread.a libboost_timer.a libboost_chrono.a libboost_program_options.a
    libboost_regex.a libboost_unit_test_framework.a)
  IF(BUILD_SHARED_LIBS)
    SET( tests_THIRDPARTY_LIBS boost_unit_test_framework)
  ELSE()
    SET( tests_THIRDPARTY_LIBS libboost_unit_test_framework.a)
  ENDIF()
  INCLUDE_DIRECTORIES( ${cpf-redshift_THIRDPARTY_INCLUDE_DIR} )
  set(PARALLEL_BUILD_FLAG "-j4")
  ADD_SUBDIRECTORY( tools )
ELSE()
  MESSAGE( STATUS "Using system third-party libraries." )
  SET( cpf-redshift_THIRDPARTY_INCLUDE_DIR "" )
  SET( cpf-redshift_THIRDPARTY_LINK_DIR "" )
  LIST(APPEND cpf-redshift_THIRDPARTY_LIBS ${CFITSIO_LIBRARIES} ${FFTW_LIBRARIES} ${GSL_LIBRARIES}
    ${Boost_FILESYSTEM_LIBRARY_RELEASE}
    ${Boost_SYSTEM_LIBRARY_RELEASE}
    ${Boost_THREAD_LIBRARY_RELEASE}
    ${Boost_TIMER_LIBRARY_RELEASE}
    ${Boost_CHRONO_LIBRARY_RELEASE}
    ${Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE})
  SET(cpf-redshift_THIRDPARTY_STATIC_LIBS ${cpf-redshift_THIRDPARTY_LIBS})
  SET(tests_THIRDPARTY_LIBS ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE})

  INCLUDE_DIRECTORIES( ${CFITSIO_INCLUDE_DIR} ${FFTW_INCLUDE_DIR} ${GSL_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS})
ENDIF()

MESSAGE("Using third party libraries: " "${cpf-redshift_THIRDPARTY_LIBS}")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#--------------------------------------------------------
# Preprocessor and linker path
#--------------------------------------------------------
INCLUDE_DIRECTORIES( ${cpf-redshift_INCLUDE_DIR} )

LINK_DIRECTORIES( ${cpf-redshift_LINK_DIR} ${cpf-redshift_THIRDPARTY_LINK_DIR}  )

#-------------------------  -------------------------------
# Extract version number
#--------------------------------------------------------
FILE (STRINGS ${ROOT_DIR}/VERSION VERSION_NUMBER)

#--------------------------------------------------------
# Compile library
#--------------------------------------------------------
SET( LIB_NAME cpf-redshift )

FILE( GLOB_RECURSE LIB_SOURCE_FILES ${ROOT_DIR}/RedshiftLibrary/src/lib/*.cpp )

IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY( ${LIB_NAME} SHARED ${LIB_SOURCE_FILES} )
  TARGET_LINK_LIBRARIES(${LIB_NAME} ${cpf-redshift_THIRDPARTY_LIBS})
  INSTALL( TARGETS ${LIB_NAME} ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
  EXPORT( PACKAGE ${LIB_NAME} )
  EXPORT(TARGETS ${LIB_NAME} FILE FindCPFRedshift.cmake)
ELSE()
  ADD_LIBRARY( ${LIB_NAME}-static STATIC ${LIB_SOURCE_FILES} )
  SET_TARGET_PROPERTIES(${LIB_NAME}-static PROPERTIES OUTPUT_NAME ${LIB_NAME})
  TARGET_LINK_LIBRARIES(${LIB_NAME}-static ${cpf-redshift_THIRDPARTY_LIBS})
  INSTALL( TARGETS ${LIB_NAME}-static ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
  EXPORT( PACKAGE ${LIB_NAME}-static )
  EXPORT(TARGETS ${LIB_NAME}-static FILE FindCPFRedshift.cmake)
  SET(LIB_NAME ${LIB_NAME}-static)
ENDIF()


#--------------------------------------------------------
# Compile test
#--------------------------------------------------------
IF(BUILD_SHARED_LIBS)
ENABLE_TESTING()

FILE( GLOB_RECURSE TEST_SOURCE_FILES ${ROOT_DIR}/RedshiftLibrary/tests/src/*_test.cpp )
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK)

CONFIGURE_FILE(${ROOT_DIR}/RedshiftLibrary/tests/src/test-config.h.in
  ${ROOT_DIR}/RedshiftLibrary/tests/src/test-config.h)

FOREACH(TEST_SOURCE ${TEST_SOURCE_FILES})
  GET_FILENAME_COMPONENT(TEST_NAME ${TEST_SOURCE} NAME_WE)
  ADD_EXECUTABLE( ${TEST_NAME} ${TEST_SOURCE} )
  TARGET_LINK_LIBRARIES( ${TEST_NAME} ${LIB_NAME} ${cpf-redshift_THIRDPARTY_LIBS} ${tests_THIRDPARTY_LIBS})
  ADD_TEST(NAME ${TEST_NAME}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/bin/${TEST_NAME} )
  SET_PROPERTY(SOURCE ${TEST_SOURCE} APPEND PROPERTY COMPILE_DEFINITIONS BOOST_TEST_MODULE=${TEST_NAME} )
ENDFOREACH(TEST_SOURCE)
ENDIF()
#--------------------------------------------------------
# Install directive
#--------------------------------------------------------
INSTALL( DIRECTORY ${ROOT_DIR}/RedshiftLibrary/ DESTINATION RedshiftLibrary )
INSTALL( FILES  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cpf-redshift-config.cmake DESTINATION cmake )

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
       INSTALL_DESTINATION ${LIB_INSTALL_DIR}/cmake
       PATH_VARS cpf-redshift_INCLUDE_DIR cpf-redshift_THIRDPARTY_INCLUDE_DIR cpf-redshift_LINK_DIR cpf-redshift_THIRDPARTY_LINK_DIR
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
  execute_process(COMMAND "find RedshiftLibrary/ ( -name *.cpp -o -name *.h -o -name *.hpp -o -name *.py -o -name *.c ) -exec ${ETAGS} -a TAGS {} +"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif(ETAGS)

