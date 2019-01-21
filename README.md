# cpf-redshift


## Per-distribution installation

### Install dependencies on CentOS7

As root:

    yum install -y epel-release
    yum -y install https://centos7.iuscommunity.org/ius-release.rpm
    rpm --import /etc/pki/rpm-gpg/IUS-COMMUNITY-GPG-KEY
    yum install -y git gcc-c++ make cmake swig boost-devel cfitsio-devel fftw-devel \
                   python36-numpy python36u-pip python-virtualenv python36-devel

### Install dependencies on Debian/Ubuntu

As root:

    apt-get install -y git cmake ccache build-essential swig python3-pip \
                       libboost-filesystem-dev libboost-system-dev libboost-thread-dev \
                       libboost-timer-dev libboost-chrono-dev libboost-program-options-dev \
                       libboost-regex-dev libboost-test-dev \
                       libcfitsio-dev libgsl-dev libfftw3-dev \
                       pkg-config \
                       python3-numpy python3-astropy

### Install depencies on MacOS

Use `brew` as packet manager on MacOS:

    brew install boost cfitsio gsl fftw

### Download, build and install

As a user, in `$HOME`:

    git clone -b develop git@gitlab.lam.fr:CPF/cpf-redshift.git
    mkdir cpf-redshift/build
    cd cpf-redshift/build
    cmake .. -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release
    make -j4
    make install
    virtualenv -p python3.6 --system-site-packages $HOME/venv
    source $HOME/venv/bin/activate
    pip3.6 install astropy # on CentOS7 and MacOSX

#### For linux users

    pip3.6 install -e $HOME/cpf-redshift

#### For MacOSX users

    MACOSX_DEPLOYMENT_TARGET=10.13 CC=clang CXX=clang++ pip3.6 install -e $HOME/cpf-redshift

Set `MACOSX_DEPLOYMENT_TARGET` variable to properly MacOSX version (`sw_vers` command on terminal).

### Build options

Build process uses cmake tool

#### -DCMAKE_BUILD_TYPE

You can build either in `Debug` or `Release` mode (defaults to `Release`).

##### example : build in Debug mode

    cd $ROOT_DIR
    mkdir build-debug
    cd build-debug
    cmake .. -DCMAKE_BUILD_TYPE=Debug
    make
    make install

##### example : build in Release  mode

    cd $ROOT_DIR
    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make
    make install


#### -DBUILD_SHARED_LIBS

You can build either static or shared library (defaults to `ON`).

    cmake .. -DBUILD_SHARED_LIBS=ON
or

    cmake .. -DBUILD_SHARED_LIBS=OFF


#### -DCMAKE_INSTALL_PREFIX

You can specify install directory (defaults to `$HOME/usr`).

    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local


#### -DBUILD_TESTS

In order to build tests, shared libs version must be enabled :

    cmake .. -DBUILD_SHARED_LIBS=ON -DBUILD_TESTS=ON

Run tests with :

    make test

#### -DTEST_COVERAGE

Tests can output line coverage statistics if you define `TEST_COVERAGE` :

    cmake .. -DBUILD_SHARED_LIBS=ON -DBUILD_TESTS=ON -DTEST_COVERAGE=ON

Run tests with :

    make test

Create coverage reports with :

     GCOV_PREFIX=$HOME/src/cpf-redshift/RedshiftLibrary/ GCOV_PREFIX_STRIP=4 lcov -q -c -t "result" -o tests.cov --no-external -b $HOME/src/cpf-redshift/RedshiftLibrary/ -d CMakeFiles
     lcov -q -r tests.cov '*/tests/src/*' -o coverage.info
     genhtml -o coverage coverage.info

### Usage in client code

#### CMakeLists.txt

You client `CMakeLists.txt` must include :

    FIND_PACKAGE( cpf-redshift )
    INCLUDE_DIRECTORIES( ${cpf-redshift_INCLUDE_DIR} )
    LINK_DIRECTORIES( ${cpf-redshift_LINK_DIR} )

    TARGET_LINK_LIBRARIES( YOUR_TARGET ${cpf-redshift_LIB})

##### Building

Say you have built cpf-redshift in `$HOME/usr` directory.

You can build the client with :

     CMAKE_MODULE_PATH=$HOME/usr/share/cmake cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON


### Build and install python interface and client with setuptools

#### With python setup.py install

Create a virtualenv and install amazed :

     virtualenv venv
     source venv/bin/activate

Install amazed with setup.py:

     python setup.py install

Or, installwith pip :

     pip install -e .

Run it with :

     amazed --help

## Build wheel package

Once cpf-redshift package has been compiled :

     pip install numpy auditwheel

     python setup.py bdist
     pip wheel /build/cpf-redshift/ -w wheelhouse
     auditwheel repair \
          wheelhouse/pyamazed-0.0.1-cp${MAJ}${MIN}-cp${MAJ}${MIN}m-linux_x86_64.whl \
          -w wheel
     pip install twine
     twine upload --repository-url https://test.pypi.org/legacy/ \
          wheel/pyamazed-0.0.1-cp36-cp36m-manylinux1_x86_64.whl

## Additional documentation

Detailed documentation about this software can be found by building the provided documentation:

Build documentation:

    cd $ROOT_DIR/tools/
    python ./builddoc.py`

Then open in your web browser:

    $ROOT_DIR/docs/html/index.html

## Libraries

+ GSL : GNU Scientific Library

[More informations here](https://www.gnu.org/software/gsl/)

## Contacts

Please send your bug reports or questions to amazed-support@lam.fr
