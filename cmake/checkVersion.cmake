#-*-cmake-*-

# The purpose of the following script is to
# reduce compilation time when in developpement mode
# we lauch many times the command :
# time pip install -v -e . -C build-dir=build-pip  -C cmake.define.CMAKE_PREFIX_PATH=$SHARED --no-build-isolation -C cmake.define.BUILD_TESTING=ON
# it prevents to rebuild un-necessary files
execute_process(COMMAND "python" "${PROJECT_SOURCE_DIR}/cmake/check_version_h.py")
