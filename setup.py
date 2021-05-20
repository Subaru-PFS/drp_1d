import os
from setuptools import setup
from distutils.command.build import build
from swig_ext import swig_ext

root_dir = os.path.dirname(__file__)
version_file = os.path.join(root_dir, 'VERSION')
version = open(version_file).read().strip()

class CustomBuild(build):
    def run(self):
        self.run_command('build_ext')
        build.run(self)

setup(
    name="pylibamazed",
    version_config={
        "count_commits_from_version_file": True,
        "template": version,
        "dev_template": version+".dev{sha}",
        "dirty_template": version+".dev{sha}.dirty",
        "version_file": version_file
    },
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
