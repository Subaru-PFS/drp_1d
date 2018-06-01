import os
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
from setuptools import setup, Extension
from pathlib import Path

try:
    from swig_ext import swig_ext
except ImportError:
    print('swig_ext.py not found. Please run cmake before setup !')
    exit()

class build_cmake_ext(build_ext):

    def run(self):
        for ext in self.extensions:
            self.build_cmake(ext)
        super().run()

    def build_cmake(self, ext):
        cwd = Path().absolute()
        print('Building in {}'.format(cwd))

        # these dirs will be created in build_py, so if you don't have
        # any python sources to bundle, the dirs will be missing
        build_temp = Path(self.build_temp)
        if not build_temp.exists():
            build_temp.mkdir(parents=True)
        extdir = Path(self.get_ext_fullpath(ext.name))
        if not extdir.exists():
            extdir.mkdir(parents=True)

        # example of cmake args
        config = 'Debug' if self.debug else 'Release'
        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + str(extdir.parent.absolute()),
            '-DCMAKE_BUILD_TYPE=' + config,
            '-DBUILD_SHARED_LIBS=ON'
        ]

        # example of build args
        build_args = [
            #'--config', config,
            '--', '-j4'
        ]

        os.chdir(str(build_temp))
        self.spawn(['cmake', str(cwd)] + cmake_args)
        if not self.dry_run:
            self.spawn(['cmake', '--build', '.'] + build_args)
        os.chdir(str(cwd))

setup(
    name = "pyamazed",
    version = "0.0.1",
    author = "CeSAM",
    author_email = "alain DOT schmitt AT lam DOT fr",
    description = ("CPF-redshift python client."),
    license = "GPLv3+",
    url = "http://www.lam.fr",

    packages = ['pyamazed'],

    long_description=open(os.path.join(os.path.dirname(__file__), 'README.md')).read(),

    setup_requires=['pytest-runner', ],
    tests_require=['pytest', ],

    classifiers=[
        "Development Status :: 1 - Planning",
        "Topic :: Scientific/Engineering",
        "Programming Language :: Python",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
    ],

    ext_modules=[ swig_ext ],
    #cmdclass={
    #    'build_ext': build_cmake_ext,
    #},

    entry_points={
        'console_scripts': [
            'amazed=pyamazed.amazed:amazed'
        ],
    },

)
