# I. Prerequisite :

See Prerequisite in Epic Core package README:

http://gitlab.oamp.dev/epic/epic_core/tree/develop/README


# II. Build and install (Unix and unix like):

## Clone Git repository:

  `git clone git@gitlab.oamp.dev:epic/epic_redshift.git $ROOT_DIR`

## Build and install third party library by running the buildandinstall.sh script:

`
 cd $ROOT_DIR/tools/
 python ./buildandinstallthirdparty.py`

Please note that each tool has its own dependency list. If a tool fails to build, remove the relevant files from $ROOT_DIR/thirdparty, install the missing dependency and try again.

The "doxygen" tool requires at least flex and bison to work.

## Create makefiles:

Build process uses cmake tool

you can build either in *Debug* or *Release* mode

### example : build in Debug mode
  `
  cd ROOT_AMAZED
  mkdir build-debug
  cd build-debug
  cmake .. -DCMAKE_BUILD_TYPE=Debug
  make
  make install
  make package`

### example : build in Release  mode
  `
  cd ROOT_AMAZED
  mkdir build
  cd build
  cmake .. -DCMAKE_BUILD_TYPE=Release
  make
  make install
  make package`

Note :
if you don't specify any *-DCMAKE_BUILD_TYPE=xxxxx* , it will build by default in Release mode


# III. Additional documentation

Detailed documentation about this software can be found by building the provided documentation:

* Build documentation:

  `
  cd $ROOT_DIR/tools/
  python ./builddoc.py`

* Then open in your web browser:

  `$ROOT_DIR/docs/html/index.html`
