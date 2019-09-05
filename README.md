# drp_1d
The **pylibamazed** library for Subaru-PFS project.

## Requirements

drp_1d has the following strict requirements:
* [gcc](https://gcc.gnu.org/)
* [python](https://www.python.org/) >=3.6
* [cmake](https://cmake.org/) >=3.6
* [swig](http://www.swig.org/) >=2.0
* [boost](https://www.boost.org/) >=1.53
* [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/) >=3.36
* [gsl](https://www.gnu.org/software/gsl/) >=2.5
* [fftw](http://www.fftw.org/) >=3.3.8

drp_1d also depends on other python packages
* [numpy](http://www.numpy.org/) >=1.16.0
* [astropy](http://www.astropy.org/) >=3.1.1

### Installing dependencies on CentOS7

As root:

    yum -y install https://centos7.iuscommunity.org/ius-release.rpm
    yum install -y git gcc-c++ make cmake3 swig boost-devel cfitsio-devel \
      patchelf python36u python36u-libs python36u-devel python36u-pip

### Installing dependencies on Debian - buster

As root:

    apt-get update
    apt-get install -y git g++ cmake swig pkg-config libboost-all-dev libgsl-dev \
      libcfitsio-dev libfftw3-dev python3 python3-pip


### Installing depencies on MacOS

Use `brew` as packet manager on MacOS:

    brew install gcc cmake swig boost cfitsio gsl fftw

Use [Anaconda](https://www.anaconda.com/) as python3 provider and then install python dependencies with `pip`:

    pip3 install numpy
    pip3 install astropy


## Installing pylibamazed

### Building C++ code from source

As a user:

    git clone git@github.com:Subaru-PFS/drp_1d.git
    mkdir drp_1d/build
    cd drp_1d/build
    cmake ..
    make -j 4
    make install

You can specify install directory with `CMAKE_INSTALL_PREFIX` (defaults to `$HOME/local`).

    cmake .. -DCMAKE_INSTALL_PREFIX=/my/path/local

### Installing pylibamazed python module from pip

From `drp_1d` root directory:

#### For linux users

    pip install .

#### For MacOSX users

    MACOSX_DEPLOYMENT_TARGET=10.13 CC=clang CXX=clang++ pip3 install .

Set `MACOSX_DEPLOYMENT_TARGET` variable to properly MacOSX version (`sw_vers` command on terminal).

### Testing an installed pylibamazed

From `drp_1d` root directory:

    cd build/
    make test

## Contacts

Please send your bug reports or questions to amazed-support@lam.fr
