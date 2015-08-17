CMAKE_MINIMUM_REQUIRED (VERSION 2.8) 
PROJECT( cpf-redshift ) 

#--------------------------------------------------------
# Setup shortcut acces directory
#--------------------------------------------------------
SET( ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.." )

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
# Find Package
#--------------------------------------------------------
FIND_PACKAGE( cpf-core ${VERSION} REQUIRED )

#--------------------------------------------------------
# Compiler options
#--------------------------------------------------------
SET( CMAKE_CXX_FLAGS "-O3 ")

SET( CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG_BUILD")

#--------------------------------------------------------
# Exported configuration variables
#--------------------------------------------------------
SET( cpf-redshift_INCLUDE_DIR ${ROOT_DIR}/include/ )
SET( cpf-redshift_THIRDPARTY_INCLUDE_DIR ${ROOT_DIR}/thirdparty/include/ )

SET( cpf-redshift_LINK_DIR ${ROOT_DIR}/lib/ )
SET( cpf-redshift_THIRDPARTY_LINK_DIR ${ROOT_DIR}/thirdparty/lib/ )

SET( cpf-redshift_THIRDPARTY_LIBS cfitsio)

#--------------------------------------------------------
# Preprocessor and linker path
#--------------------------------------------------------
INCLUDE_DIRECTORIES(    ${cpf-redshift_INCLUDE_DIR}  ${cpf-redshift_THIRDPARTY_INCLUDE_DIR} 
                        ${cpf-core_INCLUDE_DIR} ${cpf-core_THIRDPARTY_INCLUDE_DIR} )

LINK_DIRECTORIES(       ${cpf-redshift_LINK_DIR} ${cpf-redshift_THIRDPARTY_LINK_DIR} 
                        ${cpf-core_LINK_DIR} ${cpf-core_THIRDPARTY_LINK_DIR})

#-------------------------  -------------------------------
# Extract version number
#--------------------------------------------------------
FILE (STRINGS ${ROOT_DIR}/VERSION VERSION_NUMBER)

#--------------------------------------------------------
# Compile library
#--------------------------------------------------------
SET( LIB_NAME cpf-redshift )

FILE( GLOB_RECURSE LIB_SOURCE_FILES ${ROOT_DIR}/src/*.cpp )

ADD_LIBRARY( ${LIB_NAME} STATIC ${LIB_SOURCE_FILES} )

#--------------------------------------------------------
# Compile test
#--------------------------------------------------------
SET( TEST_NAME cpf-redshift-test )
  
FILE( GLOB_RECURSE TEST_SOURCE_FILES ${ROOT_DIR}/test/*.cpp )

ADD_EXECUTABLE( ${TEST_NAME} ${TEST_SOURCE_FILES} )

# Supply link librairies
TARGET_LINK_LIBRARIES( ${TEST_NAME} ${LIB_NAME}   ${cpf-core_LIB} ${cpf-core_THIRDPARTY_LIBS} ${cpf-redshift_THIRDPARTY_LIBS} cppunit.a)

#--------------------------------------------------------
# Install directive
#--------------------------------------------------------
INSTALL( TARGETS ${LIB_NAME} ARCHIVE DESTINATION lib )
INSTALL( DIRECTORY ${ROOT_DIR}/include/ DESTINATION include )
INSTALL( FILES  ${CMAKE_CURRENT_SOURCE_DIR}/cpf-redshift-config.cmake DESTINATION cmake )

#--------------------------------------------------------
# Package deployement and config
#--------------------------------------------------------
INCLUDE ( CMakePackageConfigHelpers )

# Configure config file (this is used by find_package to find where the library is)
CONFIGURE_PACKAGE_CONFIG_FILE( 
       ${CMAKE_CURRENT_SOURCE_DIR}/cpf-redshift-config.cmake.in
       ${CMAKE_CURRENT_SOURCE_DIR}/cpf-redshift-config.cmake
       INSTALL_DESTINATION cmake
        )

# Configure version file (this is used by find_package to determine a matching version)
CONFIGURE_FILE( 
       ${CMAKE_CURRENT_SOURCE_DIR}/cpf-redshift-config-version.cmake.in
       ${CMAKE_CURRENT_SOURCE_DIR}/cpf-redshift-config-version.cmake
       @ONLY
        )

EXPORT( PACKAGE ${LIB_NAME} )
