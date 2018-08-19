
from distutils.core import setup
import os
import sys

# First the ParmedTools packages:
packages = ['mdoutanalyzer']

modules = []

# Scripts
scripts = ['cpinutil.py', 'ceinutil.py', 'charmmlipid2amber.py', 'mdout_analyzer.py',
           'finddgref.py', 'fitpkaeo.py', 'fixremdcouts.py', 'genremdinputs.py', 'softcore_setup.py']

if __name__ == '__main__':
    try:
        from distutils.command.build_py import build_py_2to3 as build_py
        from distutils.command.build_scripts import build_scripts_2to3 as build_scripts
        PY3 = True
    except ImportError:
        from distutils.command.build_py import build_py
        from distutils.command.build_scripts import build_scripts
        PY3 = False

    setup(name='AmberTools',
          version='17.0',
          description='Various modules needed for AmberTools Python programs',
          author='Jason M. Swails, Ben Madej, Thomas T. Joseph, and Vinicius Wilian D. Cruzeiro',
          author_email='jason.swails -at- gmail.com',
          url='http://ambermd.org',
          license='GPL v2 or later',
          packages=packages,
          py_modules=modules,
          cmdclass={'build_py': build_py, 'build_scripts': build_scripts},
          scripts=scripts)
