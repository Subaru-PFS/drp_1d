# drp_1d
The **pylibamazed** library for Subaru-PFS project.

## Requirements

`drp_1d` has the following strict requirements:
* [gcc](https://gcc.gnu.org/)
* [python](https://www.python.org/) >=3.6
* [cmake](https://cmake.org/) >=3.12
* [swig](http://www.swig.org/) >=2.0

## Dependencies

`pylibamazed` depends on following third parties:
* [boost](https://www.boost.org/) >=1.57
* [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/) >=3.36
* [gsl](https://www.gnu.org/software/gsl/) >=2.5
* [fftw](http://www.fftw.org/) >=3.3.8

`drp_1d` also depends on following python packages
* [numpy](https://www.numpy.org/) >=1.16.0
* [astropy](https://www.astropy.org/) >=3.1.1
* [cython](https://cython.org/) >=0.17.0
* [pandas](https://pandas.pydata.org/) >=1.0.0
* [h5py](https://www.h5py.org/) >=2.9

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

### Installing python dependencies

Activate your virtual environment as needed then install python dependencies with `pip`:

    pip install numpy
    pip install astropy


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

    pip install .

### Testing an installed pylibamazed

From `drp_1d` root directory:

    cd build/
    make test

## Contacts

Please send your bug reports or questions to amazed-support@lam.fr
