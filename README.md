# pylibamazed

## About ?

pylibamazed is a Python package wrapping numerical algorithms for the analysis of 1D spectroscopic data of astrophysical sources.

## Main features

* Estimate source category (galaxy, star or QSO)
* Classify source type (sub-classification into the source category, e.g. the spectral type of a star for the star category)
* Estimate redshift
* Provide redshift reliability
* Estimate radial velocity
* Measure fluxes of emission lines

## Playing pylibamazed
pylibamazed should work on Linux and OS X.

### Dependencies

The pylibamazed algorithms are mainly developped in C++ and wrapped in Python. These algorithms are dependant on the thirdparties listed below. It is recommended to install third parties on your system using your own package manager. However, pylibamazed provides a python script to install these thirdparties. To install third parties using pylibamazed internal script, refer to the related [third parties](#Third-parties-install-guide) section.

Required third parties:
* [boost](https://www.boost.org/) >=1.57
* [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/) >=3.36
* [gsl](https://www.gnu.org/software/gsl/) >=2.5
* [fftw](http://www.fftw.org/) >=3.3.8
* [openblas](https://www.openblas.net/) >= 0.3.19
* [eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) >= 3.4.0
* [lbfgspp](https://lbfgspp.statr.me) == 0.2.0-13-g563106b

Required python packages:
* [numpy](https://www.numpy.org/) >=1.16.0
* [astropy](https://www.astropy.org/) >=3.1.1
* [cython](https://cython.org/) >=0.17.0
* [pandas](https://pandas.pydata.org/) >=1.0.0
* [h5py](https://www.h5py.org/) >=2.9

### Build and install

Prerequisites to install pylibamazed from source code:

* [gcc](https://gcc.gnu.org/)
* [python](https://www.python.org/) >=3.6
* [cmake](https://cmake.org/) >=3.12
* [swig](http://www.swig.org/) >=4.0

To build the C++ part, in `pylibamazed` root directory, run:

    mkdir build
    cd build
    cmake ..
    make
    make install

To run tests in `build` directory, run:

    make test

To build and install `pylibamazed` python module, in `pylibamazed` root directory run:

    pip install .

### Build options

The Build process uses the cmake tool

#### -DCMAKE_BUILD_TYPE

You can build the library either in `Release`, `Debug` or `Coverage` mode (default to `Release`).

For instance to build pylibamazed in Debug mode, run: 

    cmake .. -DCMAKE_BUILD_TYPE=Debug

#### -DCMAKE_INSTALL_PREFIX

You can specify the absolute path to your installation directory (default to `$HOME/.local`). If you are working with several versions of pylibamazed, it is hightly recommended to specify a proper installation directory for each version. 

    cmake .. -DCMAKE_INSTALL_PREFIX=/my/own/directory

#### -DCMAKE_PREFIX_PATH

If the thirdparties are not installed in a regular directory, you can specify the path to find thirdparties. If some thirdparties have been installed with the internal pylibamazed script, you must specify the  corresponding directory as follows

    cmake .. -DCMAKE_PREFIX_PATH=/my/thirdparty/directory


## Additional documentation

The python API documentation could be generated as follows:

Build the documentation:

    cd $ROOT_DIR/pylibamazed/doc
    make html

Then open in your web browser:

    $ROOT_DIR/pylibamazed/build/html/index.html

## Third parties installation guide

As stated earlier pylibamazed depends on several third parties (refer to [this section](#dependencies) for the complete list). It is recommended to install third parties on your system using your own package manager. However, pylibamazed provides a python script to install theses thirdparties.


### Installing with package manager

On ubuntu :
```sh
sudo apt install -y \
libboost-dev libboost-filesystem-dev libboost-thread-dev libboost-timer-dev libboost-program-options-dev libboost-test-dev \
libcfitsio-dev \
libgsl-dev \
libfftw3-dev \
libopenblas-dev \
libeigen3-dev \
liblbfgspp-dev
```

### Installing from source

    buildandinstallthirdparty.py [-h] [--workdir WORKDIR] [--prefix PREFIX] [-j PARALLEL] [--extra_flags EXTRA_FLAGS] [--force] [name1 ...]

Name argument corresponds to the third party name and could take the following values:  [`boost` | `cfitsio` | `gsl` | `fftw` | `openblas` | `eigen` | `lbfgspp`].

For instance, to install the fftw and cfitsio third parties into the `thirdparty` directory, execute:

    python tools/buildandinstallthirdparty.py fftw cfitsio

Other command line options:

`--workdir`: specifies the working directory for the third party building (absolute path)

    python tools/buildandinstallthirdparty.py fftw cfitsio --workdir=/tmp

`--prefix`: specifies the installation directory for third parties (absolute path)

    python tools/buildandinstallthirdparty.py fftw cfitsio --prefix=/usr/local

`-j`: specifies the number of make jobs to run simultaneously 

    python tools/buildandinstallthirdparty.py fftw cfitsio -j 4

`--extra_flags`: specifies extra_flag to give to the build stage of third party

    python tools/buildandinstallthirdparty.py fftw cfitsio --extra_flags=

`--force`: forces the library building and overwrites existing built library 

    python tools/buildandinstallthirdparty.py fftw cfitsio --force

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
