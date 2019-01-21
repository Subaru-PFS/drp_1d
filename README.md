# drp_1d

## Requirements

drp_1d has the following strict requirements:
* [gcc](https://gcc.gnu.org/)
* [python](https://www.python.org/) 3.5
* [cmake](https://cmake.org/)
* [swig](http://www.swig.org/)
* [boost](https://www.boost.org/)
* [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/)
* [gsl](https://www.gnu.org/software/gsl/) 2.1
* [fftw](http://www.fftw.org/)

drp_1d also depends on other python packages
* [numpy](http://www.numpy.org/)
* [astropy](http://www.astropy.org/)

### Installing dependencies on CentOS7

As root:

    yum install -y epel-release
    yum -y install https://centos7.iuscommunity.org/ius-release.rpm
	rpm --import /etc/pki/rpm-gpg/IUS-COMMUNITY-GPG-KEY
    yum install -y git gcc-c++ make cmake swig boost-devel cfitsio-devel fftw-devel \
	               python36-numpy python36u-pip python-virtualenv python36-devel

### Installing dependencies on Debian/Ubuntu

As root:

    apt-get install -y git cmake ccache build-essential swig python3-pip \
                       libboost-filesystem-dev libboost-system-dev libboost-thread-dev \
					   libboost-timer-dev libboost-chrono-dev libboost-program-options-dev \
					   libboost-regex-dev libboost-test-dev \
					   libcfitsio-dev libgsl-dev libfftw3-dev \
					   pkg-config \
					   python3-numpy python3-astropy

### Installing depencies on MacOS

Use `brew` as packet manager on MacOS:

    brew install gcc cmake swig boost cfitsio gsl fftw

Use [Anaconda](https://www.anaconda.com/) as python3 provider and then install python dependencies with `pip`:

    pip3 install astropy


## Installing drp_1d

### Building C++ code from source

As a user:

    git clone git@github.com:Subaru-PFS/drp_1d.git
    mkdir drp_1d/build
    cd drp_1d/build
    cmake .. -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=ON
    make -j4
	make install

You can specify install directory with `CMAKE_INSTALL_PREFIX` (defaults to `$HOME/usr`).

    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local

### Installing drp_1d python module from pip

From `drp_1d` root directory:

#### For linux users

    pip install .

#### For MacOSX users

    MACOSX_DEPLOYMENT_TARGET=10.13 CC=clang CXX=clang++ pip3.6 install .

Set `MACOSX_DEPLOYMENT_TARGET` variable to properly MacOSX version (`sw_vers` command on terminal).

### Testing an installed drp_1d

From `drp_1d` root directory:

    cd build/
    make test

## Contacts

Please send your bug reports or questions to amazed-support@lam.fr
