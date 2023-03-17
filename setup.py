#
# Building the extension.
#

from setuptools import setup
from mebuex import MesonExtension, build_ext

ext = MesonExtension('pybalonor.bayes', builddir='builddir')

setup(ext_modules=[ext], cmdclass={'build_ext' : build_ext})
