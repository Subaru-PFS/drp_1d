# cpf-redshift


## Build and install (Unix and unix like):

1. Clone Git repository:

Create a project folder :

    mkdir Projects
    cd Projects

Clone the cpf-redshift repository :

	git@gitlab.lam.fr:CPF/cpf-redshift.git

2. Build and install third party library by running the buildandinstall.sh script:

Execute the script in cpf-redshift/tools :

	cd cpf-redshift/tools/
	./buildandinstallthirdparty.py

Please note that each tool has its own dependency list. If a tool fails to build, remove the relevant files from $ROOT_DIR/thirdparty, install the missing dependency and try again.

The "doxygen" tool requires at least flex and bison to work.

3. Create makefiles:

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
If you don't specify any **-DCMAKE_BUILD_TYPE=xxxxx** , it will build  by default in Release mode

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
