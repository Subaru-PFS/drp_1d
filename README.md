# drp_1d
The **pylibamazed** library for Subaru-PFS project.

## Requirements

`drp_1d` has the following strict requirements:
* [gcc](https://gcc.gnu.org/)
* [python](https://www.python.org/) >=3.6
* [cmake](https://cmake.org/) >=3.12
* [swig](http://www.swig.org/) >=4.0

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

## Copyright & License

Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM

https://www.lam.fr/

This software is a computer program whose purpose is to estimate the
spectrocopic redshift of astronomical sources (galaxy/quasar/star)
from there 1D spectrum.

This software is governed by the CeCILL-C license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL-C
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-C license and that you accept its terms.
