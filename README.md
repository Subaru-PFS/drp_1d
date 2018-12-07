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
	virtualenv -p python3.6 --system-site-packages $HOME/venv
	source $HOME/venv/bin/activate
	pip3.6 install astropy # on CentOS7 only
	pip3.6 install -e $HOME/cpf-redshift/

## Build and install (Unix and unix like):

#### 1. Clone Git repository:

Create a project folder :

    mkdir Projects
    cd Projects

Clone the cpf-redshift repository :

	git clone git@gitlab.lam.fr:CPF/cpf-redshift.git

#### 2. Build and install third party library by running the buildandinstall.sh script:

Execute the script in cpf-redshift/tools :

	cd cpf-redshift/tools/
	./buildandinstallthirdparty.py

Please note that each tool has its own dependency list. If a tool fails to build, remove the relevant files from $ROOT_DIR/thirdparty, install the missing dependency and try again.

The "doxygen" tool requires at least flex and bison to work.

This step can be automatically done with cmake, by adding -DBUILD_THIRDPARTY=true parameter. e.g. :

     cmake . -DBUILD_THIRDPARTY=true

#### 3. Build libcpf-redshift

Build process uses cmake tool

You can build either in **Debug** or **Release** mode

##### example : build in Debug mode

	cd $ROOT_DIR
	mkdir build-debug
	cd build-debug
	cmake .. -DCMAKE_BUILD_TYPE=Debug
	make
	make install
	make package

##### example : build in Release  mode

	cd $ROOT_DIR
	mkdir build
	cd build
	cmake .. -DCMAKE_BUILD_TYPE=Release
	make
	make install
	make package

Note :
If you don't specify any **-DCMAKE_BUILD_TYPE=xxxxx** , it will build by default in Release mode


##### shared and static

You can build either static or shared library :

	cmake .. -DBUILD_SHARED_LIBS=ON
or

	cmake .. -DBUILD_SHARED_LIBS=OFF

##### building and running tests

In order to build tests, shared libs version must be enabled :

    cmake .. -DBUILD_SHARED_LIBS=ON

Tests can output line coverage statistics if you define `TEST_COVERAGE` :

    cmake .. -DBUILD_SHARED_LIBS=ON -DTEST_COVERAGE=ON

Run tests with :

    make test

Create coverage reports with :

     GCOV_PREFIX=$HOME/src/cpf-redshift/RedshiftLibrary/ GCOV_PREFIX_STRIP=4 lcov -q -c -t "result" -o tests.cov --no-external -b $HOME/src/cpf-redshift/RedshiftLibrary/ -d CMakeFiles
     lcov -q -r tests.cov '*/tests/src/*' -o coverage.info
	 genhtml -o coverage coverage.info

#### 4. Usage in client code

##### CMakeLists.txt

You client `CMakeLists.txt` must include :

    FIND_PACKAGE( cpf-redshift )
    INCLUDE_DIRECTORIES( ${cpf-redshift_INCLUDE_DIR} ${cpf-redshift_THIRDPARTY_INCLUDE_DIR}  )
	LINK_DIRECTORIES( ${cpf-redshift_LINK_DIR} ${cpf-redshift_THIRDPARTY_LINK_DIR} )

    TARGET_LINK_LIBRARIES( YOUR_TARGET ${cpf-redshift_THIRDPARTY_LIBS} ${cpf-redshift_LIB})
    # or alternatively, if you want to link to boost-tests :
    TARGET_LINK_LIBRARIES( YOUR_TARGET ${cpf-redshift_THIRDPARTY_LIBS} ${tests_THIRDPARTY_LIBS} ${cpf-redshift_LIB})

###### Building

Say you have built a Release and a Debug version of cpf-redshift in `~/src/cpf-redshift/build` and
`~/src/cpf-redshift/build-debug` directories.

You can build each version with :

     CMAKE_MODULE_PATH=$HOME/src/cpf-redshift/build/ cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON

or :

     CMAKE_MODULE_PATH=$HOME/src/cpf-redshift/build-debug/ cmake .. -DCMAKE_BUILD_TYPE=Debug -DBUILD_SHARED_LIBS=ON

#### 5. Build and install python interface and client with setuptools

##### With python setup.py install

Create a virtualenv and install amazed :

     virtualenv venv
	 source venv/bin/activate

Install amazed with setup.py:

     python setup.py install

Or, installwith pip :

     pip install -e .

Run it with :

     amazed --help


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

Please send your bug reports or questions to alain DOT schmitt AT lam DOT fr
