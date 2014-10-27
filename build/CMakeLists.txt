CMAKE_MINIMUM_REQUIRED (VERSION 2.6) 
PROJECT( epic_redshift ) 

SET( COMPILE_TEST "" )
#SET( ENABLE_PROFILING "" )

#--------------------------------------------------------
# Setup shortcut acces directory
#--------------------------------------------------------
SET( EPIC_REDSHIFT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../" )

MESSAGE( STATUS "CMAKE_BUILD_TYPE = " ${CMAKE_BUILD_TYPE} )
MESSAGE( STATUS "CMAKE_SYSTEM_NAME = " ${CMAKE_SYSTEM_NAME} )
MESSAGE( STATUS "CMAKE_SYSTEM_PROCESSOR = " ${CMAKE_SYSTEM_PROCESSOR} )

SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${EPIC_REDSHIFT_DIR}/lib )
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${EPIC_REDSHIFT_DIR}/lib )
SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EPIC_REDSHIFT_DIR}/bin )

#--------------------------------------------------------
# Determine where Epic_Core is
#--------------------------------------------------------
IF( NOT DEFINED EPIC_CORE_DIR )
    SET( EPIC_CORE_DIR "${EPIC_REDSHIFT_DIR}../core/" )
ENDIF()


MESSAGE( STATUS "EPIC_CORE_DIR = " ${EPIC_CORE_DIR} )

INCLUDE( ${EPIC_CORE_DIR}/build/CMakeCommon.in )

#--------------------------------------------------------
# Setup bin sufix
# Save architecture and system info used by the package_release.py script
#--------------------------------------------------------
SET( BIN_SUFFIX "" )
UPDATE_SYSTEM_NAME_AND_PROCESSOR( ${EPIC_REDSHIFT_DIR} BIN_SUFFIX )

#--------------------------------------------------------
# Compiler options
#--------------------------------------------------------
SET( CMAKE_CXX_FLAGS "-O3 ")

IF( DEFINED ENABLE_PROFILING )
  SET( CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG_BUILD -pg")
ELSE()
  SET( CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG_BUILD")
ENDIF()

INCLUDE_DIRECTORIES( "${EPIC_REDSHIFT_DIR}/include/" )
INCLUDE_DIRECTORIES( "${EPIC_REDSHIFT_DIR}/thirdparty/include/" )
INCLUDE_DIRECTORIES( "${EPIC_CORE_DIR}/include/" )
INCLUDE_DIRECTORIES( "${EPIC_CORE_DIR}/thirdparty/include/" )

LINK_DIRECTORIES( "${EPIC_REDSHIFT_DIR}/lib" )
LINK_DIRECTORIES( "${EPIC_REDSHIFT_DIR}/thirdparty/lib/" )
LINK_DIRECTORIES( "${EPIC_CORE_DIR}/thirdparty/lib/" )
LINK_DIRECTORIES( "${EPIC_CORE_DIR}/lib/" )

#-------------------------  -------------------------------
# Extract version number
#--------------------------------------------------------
FILE (STRINGS "${EPIC_REDSHIFT_DIR}/VERSION" VERSION_NUMBER)


MESSAGE( STATUS "VERSION_NUMBER: " ${VERSION_NUMBER} )

#--------------------------------------------------------
# Compile library
#--------------------------------------------------------

FILE( GLOB_RECURSE LIB_REDSHIFT_SOURCE_FILE_LIST ${EPIC_REDSHIFT_DIR}/src/*.cpp )

SET( EPIC_REDSHIFT_STATIC epic_redshift_static )

ADD_LIBRARY( ${EPIC_REDSHIFT_STATIC} STATIC ${LIB_REDSHIFT_SOURCE_FILE_LIST} )


#--------------------------------------------------------
# Compile test
#--------------------------------------------------------
IF( DEFINED COMPILE_TEST )

  SET( EPIC_REDSHIFT_TEST epic_redshift_test )
  
  FILE( GLOB_RECURSE TEST_SOURCE_FILE_LIST ${EPIC_REDSHIFT_DIR}/test/*.cpp )
  
  ADD_EXECUTABLE( ${EPIC_REDSHIFT_TEST} ${TEST_SOURCE_FILE_LIST} )

  # Supply link librairies
  TARGET_LINK_LIBRARIES( ${EPIC_REDSHIFT_TEST} ${EPIC_REDSHIFT_STATIC} )

  TARGET_LINK_LIBRARIES( ${EPIC_REDSHIFT_TEST} pthread )
  
  IF( ${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Darwin" )
    
  ELSEIF( ${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Linux" )
    TARGET_LINK_LIBRARIES( ${EPIC_REDSHIFT_TEST} rt )
  ENDIF()
  
  TARGET_LINK_LIBRARIES( ${EPIC_REDSHIFT_TEST} epic_core_static.a )
  TARGET_LINK_LIBRARIES( ${EPIC_REDSHIFT_TEST} cppunit.a )
  TARGET_LINK_LIBRARIES( ${EPIC_REDSHIFT_TEST} cfitsio.a )
  TARGET_LINK_LIBRARIES( ${EPIC_REDSHIFT_TEST} boost_filesystem.a )
  TARGET_LINK_LIBRARIES( ${EPIC_REDSHIFT_TEST} boost_system.a )
  TARGET_LINK_LIBRARIES( ${EPIC_REDSHIFT_TEST} boost_thread.a )
  TARGET_LINK_LIBRARIES( ${EPIC_REDSHIFT_TEST} zmq.a )
  
  
  ADD_CUSTOM_TARGET( runtest COMMAND cd ${EPIC_REDSHIFT_DIR}/bin/ && ./${EPIC_REDSHIFT_TEST} DEPENDS ${EPIC_REDSHIFT_TEST} )
  
ENDIF()




#--------------------------------------------------------
# Compile docs
#--------------------------------------------------------
ADD_CUSTOM_TARGET( makedoc COMMAND cd ${EPIC_REDSHIFT_DIR}/doc/ && ../thirdparty/bin/doxygen ./doxygen.conf )
