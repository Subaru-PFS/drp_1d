# Build all third party dependencies
IF (BUILD_SHARED_LIBS)
  SET(SHARED_FLAG "--shared")
ELSE()
  SET(SHARED_FLAG "--static")
ENDIF()

EXECUTE_PROCESS( COMMAND ${PYTHON_EXECUTABLE} buildandinstallthirdparty.py ${PARALLEL_BUILD_FLAG} ${SHARED_FLAG} --prefix ${CMAKE_INSTALL_PREFIX} ${MISSING_THIRDPARTIES}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  RESULT_VARIABLE result )
