import setuptools
from setuptools import setup
from distutils.core import setup, Extension
import glob

# Available at setup time due to pyproject.toml
#from pybind11.setup_helpers import Pybind11Extension, build_ext
#from pybind11 import get_cmake_dir

import sys

__version__ = "0.1.12"

cppfiles = glob.glob("ext/*.cpp")
print(cppfiles)

ext_modules = [
    Extension("_cobio", cppfiles, extra_compile_args=['-std=c++20'])
]

setup(
    name="cobio",
    version=__version__,
    author="Jian Wang",
    author_email="jianopt@gmail.com",
    url="https://github.com/hust220/cobio",
    description="Jian's Python Library",
    long_description="",
    include_package_data=True,
    # package_data={
    #     '': ['medusa_parameter/*']
    # },
    packages=setuptools.find_packages(),
    ext_modules=ext_modules,
    extras_require={"test": "pytest"},
    # Currently, build_ext only provides an optional "highest supported C++
    # level" feature, but in the future it may provide more features.
    #  cmdclass={"build_ext": build_ext},
    zip_safe=False,
)
