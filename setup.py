import os
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
from setuptools import setup, Extension
from pathlib import Path
import numpy

try:
    from swig_ext import swig_ext
except ImportError:
    print('swig_ext.py not found. Please run cmake before setup !')
    exit()

setup(
    name="pyamazed",
    version="0.0.1",
    author="LAM - Laboratoire d'Astrophysique de Marseille",
    author_email="amazed-support@lam.fr",
    description=("CPF-redshift python client."),
    license="GPLv3+",
    url="http://www.lam.fr",
    packages=['pyamazed'],
    long_description=open(os.path.join(os.path.dirname(__file__),
                                       'README.md')).read(),

    setup_requires=['pytest-runner', 'numpy', 'astropy'],
    tests_require=['pytest', ],
    install_requires=['numpy', 'astropy'],
    py_modules=['pyamazed/redshift'],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Topic :: Scientific/Engineering",
        "Programming Language :: Python",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
    ],
    ext_modules=[swig_ext],
    entry_points={
        'console_scripts': [
            'amazed=pyamazed.amazed:amazed'
        ],
    },

)
