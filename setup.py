from setuptools import setup
# from setuptools.extension import Extension
from Cython.Build import cythonize

extensions = cythonize(
    "src/orthofinder/tools/ete4/*.pyx",  
    compiler_directives={"language_level": "3"},
)

setup(
    ext_modules=extensions
)

