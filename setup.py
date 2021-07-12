import os
from setuptools import setup
from distutils.command.build import build
from swig_ext import swig_ext

root_dir = os.path.dirname(__file__)
version_file = os.path.join(root_dir, 'VERSION')
version = open(version_file).read().strip()

gitrevision_file = os.path.join(root_dir, 'RedshiftLibrary/RedshiftLibrary/version.h')
gitrevision = open(gitrevision_file).readline().split()[2].replace('"','').replace("-",".").strip()
if gitrevision: 
    gitrevision_python="+git"+gitrevision 
else:
    gitrevision_python=''

__version__ = version


class CustomBuild(build):
    def run(self):
        self.run_command('build_ext')
        build.run(self)


setup(
    name="pylibamazed",
    version=__version__, 
    author="LAM - Laboratoire d'Astrophysique de Marseille",
    author_email="amazed-support@lam.fr",
    description=("AMAZED python library."),
    license="GPLv3+",
    url="http://www.lam.fr",
    packages=['pylibamazed'],
    package_dir = {'':'pylibamazed/python'},
    include_package_data=True,
    long_description=open(os.path.join(os.path.dirname(__file__),
                                       'README.md')).read(),
    setup_requires=['setuptools-git-versioning>=1.2.0', 'numpy>=1.16.0', 'astropy>=3.1.1', 'cython>=0.17.0', 'pandas>=1.0.0', 'h5py>=2.9'],
    install_requires=['setuptools-git-versioning>=1.2.0', 'numpy>=1.16.0', 'astropy>=3.1.1', 'cython>=0.17.0', 'pandas>=1.0.0', 'h5py>=2.9'],
    tests_require=['pytest-runner', 'pytest', ],
    py_modules=['pylibamazed/redshift'],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Topic :: Scientific/Engineering",
        "Programming Language :: Python",
        "License :: OSI Approved :: GNU General Public License v3 or "
        "later (GPLv3+)",
    ],
    cmdclass={'build': CustomBuild},
    ext_modules=[swig_ext],
)
