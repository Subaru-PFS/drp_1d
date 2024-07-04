import os
from setuptools import setup
from setuptools.command.build import build
from swig_ext import swig_ext

root_dir = os.path.dirname(__file__)
version_file = os.path.join(root_dir, 'VERSION')
version = open(version_file).read().strip()

gitrevision_file = os.path.join(root_dir, 'RedshiftLibrary/RedshiftLibrary/version.h')
gitrevision = open(gitrevision_file).readline().split()[2].replace('"', '').replace("-", ".").strip()
if gitrevision:
    gitrevision_python = "+git" + gitrevision
else:
    gitrevision_python = ''

__version__ = gitrevision


class CustomBuild(build):
    def run(self):
        self.run_command('build_ext')
        build.run(self)


setup(
    version=__version__,
    include_package_data=True,
    py_modules=['pylibamazed/redshift'],
    cmdclass={'build': CustomBuild},
    ext_modules=[swig_ext]
)
