CMAKE_MINIMUM_REQUIRED (VERSION 2.6) 
PROJECT( epic_redshift ) 

#--------------------------------------------------------
# Find Package
#--------------------------------------------------------
FIND_PACKAGE( epic_core REQUIRED )

#--------------------------------------------------------
# Setup shortcut acces directory
#--------------------------------------------------------
SET( EPIC_REDSHIFT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.." )

MESSAGE( STATUS "------------------------------------------------------" )
MESSAGE( STATUS "EPIC_REDSHIFT:" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = " ${CMAKE_BUILD_TYPE} )
MESSAGE( STATUS "CMAKE_SYSTEM_NAME = " ${CMAKE_SYSTEM_NAME} )
MESSAGE( STATUS "CMAKE_SYSTEM_PROCESSOR = " ${CMAKE_SYSTEM_PROCESSOR} )
MESSAGE( STATUS "CMAKE_CURRENT_SOURCE_DIR = " ${CMAKE_CURRENT_SOURCE_DIR} )
MESSAGE( STATUS "CMAKE_MODULE_PATH = " ${CMAKE_MODULE_PATH} )
MESSAGE( STATUS "CMAKE_SOURCE_DIR = " ${CMAKE_SOURCE_DIR} )

SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${EPIC_REDSHIFT_DIR}/lib )
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${EPIC_REDSHIFT_DIR}/lib )
SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EPIC_REDSHIFT_DIR}/bin )

#--------------------------------------------------------
# Compiler options
#--------------------------------------------------------
SET( CMAKE_CXX_FLAGS "-O3 ")

IF( DEFINED ENABLE_PROFILING )
  SET( CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG_BUILD -pg")
ELSE()
  SET( CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG_BUILD")
ENDIF()

#--------------------------------------------------------
# Exported configuration variables
#--------------------------------------------------------
SET( epic_redshift_INCLUDE_DIR ${EPIC_REDSHIFT_DIR}/include/ )
SET( epic_redshift_thirdparty_INCLUDE_DIR ${EPIC_REDSHIFT_DIR}/thirdparty/include/ )

SET( epic_redshift_LINK_DIR ${EPIC_REDSHIFT_DIR}/lib/ )
SET( epic_redshift_thirdparty_LINK_DIR ${EPIC_REDSHIFT_DIR}/thirdparty/lib/ )

SET( epic_redshift_LINK_LIBS cfitsio)

#--------------------------------------------------------
# Preprocessor and linker path
#--------------------------------------------------------
INCLUDE_DIRECTORIES(    ${epic_redshift_INCLUDE_DIR}  ${epic_redshift_thirdparty_INCLUDE_DIR} 
                        ${epic_core_INCLUDE_DIR} ${epic_core_thirdparty_INCLUDE_DIR} )

LINK_DIRECTORIES(       ${epic_redshift_LINK_DIR} ${epic_redshift_thirdparty_LINK_DIR} 
                        ${epic_core_LINK_DIR} ${epic_core_thirdparty_LINK_DIR})

#-------------------------  -------------------------------
# Extract version number
#--------------------------------------------------------
FILE (STRINGS ${EPIC_REDSHIFT_DIR}/VERSION VERSION_NUMBER)
MESSAGE( STATUS "VERSION_NUMBER: " ${VERSION_NUMBER} )

#--------------------------------------------------------
# Compile library
#--------------------------------------------------------
SET( EPIC_REDSHIFT epic_redshift )

FILE( GLOB_RECURSE LIB_REDSHIFT_SOURCE_FILE_LIST ${EPIC_REDSHIFT_DIR}/src/*.cpp )
ADD_LIBRARY( ${EPIC_REDSHIFT} STATIC ${LIB_REDSHIFT_SOURCE_FILE_LIST} )

#--------------------------------------------------------
# Compile test
#--------------------------------------------------------
SET( EPIC_REDSHIFT_TEST epic_redshift_test )
  
FILE( GLOB_RECURSE TEST_SOURCE_FILE_LIST ${EPIC_REDSHIFT_DIR}/test/*.cpp )
ADD_EXECUTABLE( ${EPIC_REDSHIFT_TEST} ${TEST_SOURCE_FILE_LIST} )

# Supply link librairies
TARGET_LINK_LIBRARIES( ${EPIC_REDSHIFT_TEST} ${EPIC_REDSHIFT} cppunit.a ${epic_core_LINK_LIBS} ${epic_redshift_LINK_LIBS})

ADD_CUSTOM_TARGET( redshift_runtest COMMAND cd ${EPIC_REDSHIFT_DIR}/bin/ && ./${EPIC_REDSHIFT_TEST} DEPENDS ${EPIC_REDSHIFT_TEST} )
  
#--------------------------------------------------------
# Install directive
#--------------------------------------------------------
INSTALL( TARGETS ${EPIC_REDSHIFT} ARCHIVE DESTINATION lib )
INSTALL( DIRECTORY  ${EPIC_REDSHIFT_DIR}/include/ DESTINATION include )
INSTALL( FILES  ${CMAKE_CURRENT_SOURCE_DIR}/epic_redshift-config.cmake DESTINATION cmake )

#--------------------------------------------------------
# Package deployement and config
#--------------------------------------------------------
INCLUDE ( CMakePackageConfigHelpers )

CONFIGURE_PACKAGE_CONFIG_FILE( 
       ${CMAKE_CURRENT_SOURCE_DIR}/epic_redshift-config.cmake.in
       ${CMAKE_CURRENT_SOURCE_DIR}/epic_redshift-config.cmake
       INSTALL_DESTINATION cmake
       
        )
       
EXPORT( PACKAGE ${EPIC_REDSHIFT} )

#--------------------------------------------------------
# Compile docs
#--------------------------------------------------------
ADD_CUSTOM_TARGET( redshift_makedoc COMMAND cd ${EPIC_REDSHIFT_DIR}/doc/ && ../thirdparty/bin/doxygen ./doxygen.conf )
