# setup.py
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy, pathlib

ext_modules = [
    Extension("fp_tools.utils.sequences", [str(pathlib.Path("src/fp_tools/utils/sequences.pyx"))],
              include_dirs=[numpy.get_include()]),
    Extension("fp_tools.utils.ngs",       [str(pathlib.Path("src/fp_tools/utils/ngs.pyx"))],
              include_dirs=[numpy.get_include()]),
    Extension("fp_tools.utils.signals",   [str(pathlib.Path("src/fp_tools/utils/signals.pyx"))],
              include_dirs=[numpy.get_include()]),
]

setup(
    ext_modules=cythonize(
        ext_modules,
        language_level="3",
        compiler_directives={"boundscheck": False, "wraparound": False},
    ),
)
