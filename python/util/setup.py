# from setuptools import setup
# from Cython.Build import cythonize
# import os

# # print(os.path.join(os.path.dirname(os.path.realpath(__file__)),"exp.pyx"))
# setup(
#     name='Experiment Utils',
#     ext_modules=cythonize(
#         os.path.join(os.path.dirname(os.path.realpath(__file__)),"exp.pyx"),
#         compiler_directives={'language_level' : "3"}
#         ),
#     zip_safe=False,    
# )

from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

# setup(
#     ext_modules=[
#         Extension("exp", ["exp.c"],
#                   include_dirs=[numpy.get_include()]),
#     ],
# )

# # Or, if you use cythonize() to make the ext_modules list,
# # include_dirs can be passed to setup()

setup(
    ext_modules=cythonize("exp.pyx", compiler_directives={'language_level' : "3"}),
    include_dirs=[numpy.get_include()],
    package_dir={'util': ''},
)    
