CMAKE_MINIMUM_REQUIRED (VERSION 2.8)
PROJECT( cpf-redshift )

# extra path for modules
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
include(SetupVar)

#--------------------------------------------------------
# Setup shortcut acces directory
#--------------------------------------------------------
SET( ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}" )

SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${ROOT_DIR}/lib )
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${ROOT_DIR}/lib )
SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${ROOT_DIR}/bin )

IF( DEFINED BUILD_THIRDPARTY )
  STRING( TOLOWER ${BUILD_THIRDPARTY} __BUILD3RD )
  IF( __BUILD3RD STREQUAL "true")
    ADD_SUBDIRECTORY( tools )
  ENDIF ()
ENDIF ()


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
SET( cpf-redshift_THIRDPARTY_INCLUDE_DIR ${ROOT_DIR}/thirdparty/include/ )

SET( cpf-redshift_LINK_DIR ${ROOT_DIR}/lib/ )
SET( cpf-redshift_THIRDPARTY_LINK_DIR ${ROOT_DIR}/thirdparty/lib/ )

SET( cpf-redshift_THIRDPARTY_LIBS cfitsio pthread boost_filesystem boost_system boost_thread boost_timer boost_chrono  boost_program_options boost_regex  libgsl.a libgslcblas.a fftw3.a)

#--------------------------------------------------------
# Preprocessor and linker path
#--------------------------------------------------------
INCLUDE_DIRECTORIES(    ${cpf-redshift_INCLUDE_DIR}  ${cpf-redshift_THIRDPARTY_INCLUDE_DIR} )

LINK_DIRECTORIES(       ${cpf-redshift_LINK_DIR} ${cpf-redshift_THIRDPARTY_LINK_DIR}  )

#-------------------------  -------------------------------
# Extract version number
#--------------------------------------------------------
FILE (STRINGS ${ROOT_DIR}/VERSION VERSION_NUMBER)

#--------------------------------------------------------
# Compile library
#--------------------------------------------------------
SET( LIB_NAME cpf-redshift${EXT} )

FILE( GLOB_RECURSE LIB_SOURCE_FILES ${ROOT_DIR}/RedshiftLibrary/src/lib/*.cpp )

ADD_LIBRARY( ${LIB_NAME} STATIC ${LIB_SOURCE_FILES} )

#--------------------------------------------------------
# Compile test
#--------------------------------------------------------
SET( TEST_NAME cpf-redshift-test${EXT} )

FILE( GLOB_RECURSE TEST_SOURCE_FILES ${ROOT_DIR}/RedshiftLibrary/tests/src/*.cpp )

ADD_EXECUTABLE( ${TEST_NAME} ${TEST_SOURCE_FILES} )

# Supply link librairies
TARGET_LINK_LIBRARIES( ${TEST_NAME} ${LIB_NAME}  ${cpf-redshift_THIRDPARTY_LIBS} boost_unit_test_framework )

#--------------------------------------------------------
# Install directive
#--------------------------------------------------------
INSTALL( TARGETS ${LIB_NAME} ARCHIVE DESTINATION lib )
INSTALL( DIRECTORY ${ROOT_DIR}/RedshiftLibrary/ DESTINATION RedshiftLibrary )
INSTALL( FILES  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cpf-redshift-config.cmake DESTINATION cmake )

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
