find_program(GIT git)
if(GIT)
  set(GIT_DIR "${PROJECT_SOURCE_DIR}/.git")
  if(NOT IS_DIRECTORY ${GIT_DIR})
    file(READ ${PROJECT_SOURCE_DIR}/VERSION version)
    string(STRIP ${version} GIT_REVISION)
    message(STATUS ".git repository not found. Default version used: ${GIT_REVISION}")
  else()
    execute_process(COMMAND ${GIT} rev-parse --short HEAD OUTPUT_VARIABLE GIT_REVISION_UNSTRIPPED)
    string(STRIP ${GIT_REVISION_UNSTRIPPED} GIT_REVISION)
    execute_process(COMMAND ${GIT} status --untracked-files=no --short OUTPUT_VARIABLE GIT_UNCOMMITED)
    if(GIT_UNCOMMITED)
      message(STATUS "You have uncommited changes. Using -untracked flag.")
      set(GIT_REVISION "${GIT_REVISION}-untracked")
    endif()
  endif()
else()
  file(READ ${PROJECT_SOURCE_DIR}/VERSION version)
  string(STRIP ${version} GIT_REVISION)
  message(STATUS "Git not found. Default version used: ${GIT_REVISION}")
endif()
