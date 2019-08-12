# Build all third party dependencies
IF (BUILD_SHARED_LIBS)
  SET(SHARED_FLAG "--shared")
ELSE()
  SET(SHARED_FLAG "--static")
ENDIF()

FIND_PACKAGE(Python3 COMPONENTS Interpreter)
IF (NOT Python3_Interpreter_FOUND)
  MESSAGE(FATAL_ERROR "Can't find Python3 interpreter")
ENDIF()

EXECUTE_PROCESS( COMMAND ${Python3_EXECUTABLE} buildandinstallthirdparty.py ${PARALLEL_BUILD_FLAG} ${SHARED_FLAG} --prefix ${CMAKE_INSTALL_PREFIX} ${MISSING_THIRDPARTIES}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  RESULT_VARIABLE result )
