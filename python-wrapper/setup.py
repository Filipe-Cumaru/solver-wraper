from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

moab_root = "/home/facsa/MOAB"
moab_include = moab_root + "/include"
moab_lib = moab_root + "/lib"
trilinos_root = "/home/facsa/Trilinos-Serial"
trilinos_include = trilinos_root + "/include"
trilinos_lib = trilinos_root + "/lib"

include_paths = [moab_include, trilinos_include, np.get_include()]
lib_paths = [moab_lib, trilinos_lib]

setup(
    name="TPFASolver",
    ext_modules=cythonize("pytpfa/tpfa.pyx", language='c++'),
    version="0.0.1",
    author="Filipe Cumaru",
    author_email="facsa@cin.ufpe.br"
)
