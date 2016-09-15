from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need
# fine tuning.
buildOptions = dict(packages = ["PyQt5"], excludes = [])

import sys
base = 'Win32GUI' if sys.platform=='win32' else None

executables = [
    Executable('aprocessgui.py', base=base)
]

setup(name='amazed-gui',
      version = 'O.1',
      description = 'gui for amazed',
      options = dict(build_exe = buildOptions),
      executables = executables)
