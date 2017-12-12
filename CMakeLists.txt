CMAKE_MINIMUM_REQUIRED (VERSION 2.8)
PROJECT( cpf-redshift )

# extra path for modules
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
include(SetupVar)

FIND_PACKAGE(Boost REQUIRED COMPONENTS filesystem system thread timer chrono program_options unit_test_framework)
FIND_PACKAGE(Cfitsio)
FIND_PACKAGE(GSL)
FIND_PACKAGE(FFTW)

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

#--------------------------------------------------------
# Exported configuration variables
#--------------------------------------------------------
SET( cpf-redshift_INCLUDE_DIR ${ROOT_DIR}/RedshiftLibrary/ )
SET( cpf-redshift_LINK_DIR ${ROOT_DIR}/lib/ )

SET( cpf-redshift_THIRDPARTY_LIBS "")

IF( NOT DEFINED CFITSIO_INCLUDE_DIR )
  LIST(APPEND cpf-redshift_THIRDPARTY_LIBS cfitsio )
ENDIF()

IF( NOT DEFINED FFTW_INCLUDE_DIR )
  LIST(APPEND cpf-redshift_THIRDPARTY_LIBS fftw3 )
ENDIF()

IF( NOT DEFINED GSL_INCLUDE_DIRS )
  LIST(APPEND cpf-redshift_THIRDPARTY_LIBS gsl gslcblas)
ENDIF()

IF( NOT DEFINED Boost_INCLUDE_DIRS )
  LIST(APPEND cpf-redshift_THIRDPARTY_LIBS pthread boost_filesystem boost_system boost_thread boost_timer boost_chrono boost_program_options boost_regex)
  SET(tests_THIRDPARTY_LIBS boost_unit_test_framework)
ENDIF()

IF( cpf-redshift_THIRDPARTY_LIBS )
  # Set THIRDPARTY_{INCLUDE,LNIK}_DIR if using any third-party
  SET( cpf-redshift_THIRDPARTY_INCLUDE_DIR ${ROOT_DIR}/thirdparty/include/ )
  SET( cpf-redshift_THIRDPARTY_LINK_DIR ${ROOT_DIR}/thirdparty/lib/)
ELSE()
  SET( cpf-redshift_THIRDPARTY_INCLUDE_DIR ${cpf-redshift_INCLUDE_DIR} )
  SET( cpf-redshift_THIRDPARTY_LINK_DIR ${cpf-redshift_LINK_DIR} )
ENDIF()

IF( DEFINED BUILD_THIRDPARTY )
  STRING( TOLOWER ${BUILD_THIRDPARTY} __BUILD3RD )
  IF( __BUILD3RD STREQUAL "true")
    ADD_SUBDIRECTORY( tools )
  ENDIF ()
ENDIF ()

MESSAGE("Using third party libraries: " ${cpf-redshift_THIRDPARTY_LIBS})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#--------------------------------------------------------
# Preprocessor and linker path
#--------------------------------------------------------
INCLUDE_DIRECTORIES( ${cpf-redshift_INCLUDE_DIR} ${cpf-redshift_THIRDPARTY_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} ${FFTW_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR})

LINK_DIRECTORIES( ${cpf-redshift_LINK_DIR} ${cpf-redshift_THIRDPARTY_LINK_DIR}  )

#-------------------------  -------------------------------
# Extract version number
#--------------------------------------------------------
FILE (STRINGS ${ROOT_DIR}/VERSION VERSION_NUMBER)

#--------------------------------------------------------
# Compile library
#--------------------------------------------------------
SET( LIB_NAME cpf-redshift${EXT} )

FILE( GLOB_RECURSE LIB_SOURCE_FILES ${ROOT_DIR}/RedshiftLibrary/src/lib/*.cpp )

ADD_LIBRARY( ${LIB_NAME} SHARED ${LIB_SOURCE_FILES} )
TARGET_LINK_LIBRARIES(${LIB_NAME}
  ${cpf-redshift_THIRDPARTY_LIBS}
  ${Boost_FILESYSTEM_LIBRARY_RELEASE}
  ${Boost_SYSTEM_LIBRARY_RELEASE}
  ${Boost_THREAD_LIBRARY_RELEASE}
  ${Boost_TIMER_LIBRARY_RELEASE}
  ${Boost_CHRONO_LIBRARY_RELEASE}
  ${Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE}
  ${CFITSIO_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES})

#--------------------------------------------------------
# Compile test
#--------------------------------------------------------
ENABLE_TESTING()

SET( TEST_NAME cpf-redshift-test${EXT} )

FILE( GLOB_RECURSE TEST_SOURCE_FILES ${ROOT_DIR}/RedshiftLibrary/tests/src/*.cpp )

FOREACH(TEST_SOURCE ${TEST_SOURCE_FILES})
  GET_FILENAME_COMPONENT(TEST_NAME ${TEST_SOURCE} NAME_WE)
  ADD_EXECUTABLE( ${TEST_NAME} ${TEST_SOURCE} )
  TARGET_LINK_LIBRARIES( ${TEST_NAME} ${LIB_NAME} ${cpf-redshift_THIRDPARTY_LIBS} ${tests_THIRDPARTY_LIBS}
  ${Boost_LIBRARIES} ${CFITSIO_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES})
  ADD_TEST(NAME ${TEST_NAME}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/bin/${TEST_NAME} )
ENDFOREACH(TEST_SOURCE)

#--------------------------------------------------------
# Install directive
#--------------------------------------------------------
#INSTALL( TARGETS ${LIB_NAME} ARCHIVE DESTINATION lib )
INSTALL( TARGETS ${LIB_NAME} LIBRARY DESTINATION lib )
INSTALL( DIRECTORY ${ROOT_DIR}/RedshiftLibrary/ DESTINATION RedshiftLibrary )
INSTALL( FILES  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cpf-redshift-config.cmake DESTINATION cmake )


#--------------------------------------------------------
# Build python wrapper
#--------------------------------------------------------
FIND_PACKAGE(SWIG)
INCLUDE(${SWIG_USE_FILE})
FIND_PACKAGE(PythonLibs)
INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_DIRS} )

SET_PROPERTY(SOURCE RedshiftLibrary/RedshiftLibrary/redshift.i PROPERTY CPLUSPLUS ON)
SWIG_ADD_LIBRARY(redshift LANGUAGE python SOURCES RedshiftLibrary/RedshiftLibrary/redshift.i)
SWIG_LINK_LIBRARIES( redshift ${PYTHON_LIBRARIES}
  ${cpf-redshift_THIRDPARTY_LIBS}
  ${Boost_LIBRARIES}
  ${CFITSIO_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LIB_NAME} )

#--------------------------------------------------------
# Package deployement and config
#--------------------------------------------------------
INCLUDE ( CMakePackageConfigHelpers )

# Configure config file (this is used by find_package to find where the library is)
CONFIGURE_PACKAGE_CONFIG_FILE(
       ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cpf-redshift-config.cmake.in
       ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cpf-redshift-config.cmake
       INSTALL_DESTINATION cmake
        )

# Configure version file (this is used by find_package to determine a matching version)
CONFIGURE_FILE(
       ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cpf-redshift-config-version.cmake.in
       ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cpf-redshift-config-version.cmake
       @ONLY
        )

EXPORT( PACKAGE ${LIB_NAME} )
