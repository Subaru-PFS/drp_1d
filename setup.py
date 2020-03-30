import os
from setuptools import setup
from distutils.command.build import build
from swig_ext import swig_ext

class CustomBuild(build):
    def run(self):
        self.run_command('build_ext')
        build.run(self)

setup(
    name="pylibamazed",
    version="0.9.0.dev5701",
    author="LAM - Laboratoire d'Astrophysique de Marseille",
    author_email="amazed-support@lam.fr",
    description=("AMAZED python library."),
    license="GPLv3+",
    url="http://www.lam.fr",
    packages=['pylibamazed'],
    long_description=open(os.path.join(os.path.dirname(__file__),
                                       'README.md')).read(),
    setup_requires=['numpy>=1.16.0', 'astropy>=3.1.1'],
    install_requires=['numpy>=1.16.0', 'astropy>=3.1.1'],
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
