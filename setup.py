import setuptools
from setuptools import setup
from distutils.core import setup, Extension

# Available at setup time due to pyproject.toml
#from pybind11.setup_helpers import Pybind11Extension, build_ext
#from pybind11 import get_cmake_dir

import sys

__version__ = "0.1.0"

# The main interface is through Pybind11Extension.
# * You can add cxx_std=11/14/17, and then build_ext can be removed.
# * You can set include_pybind11=false to add the include directory yourself,
#   say from a submodule.
#
# Note:
#   Sort input source files if you glob sources to ensure bit-for-bit
#   reproducible builds (https://github.com/pybind/python_example/pull/53)

ext_modules = [
#  Pybind11Extension("_jnpy", ["ext/main.cpp"], define_macros = [('VERSION_INFO', __version__)])
  Extension("_jnpy", ["ext/main.cpp"], extra_compile_args=['-std=c++11'])
]

setup(
  name="jnpy",
  version=__version__,
  author="Jian Wang",
  author_email="jianopt@gmail.com",
  url="https://github.com/hust220/jnpy",
  description="Extract binding pocket from protein-ligand complexes",
  long_description="",
  packages=setuptools.find_packages(),
  ext_modules=ext_modules,
  extras_require={"test": "pytest"},
  # Currently, build_ext only provides an optional "highest supported C++
  # level" feature, but in the future it may provide more features.
#  cmdclass={"build_ext": build_ext},
  zip_safe=False,
)

