import os, sys
import subprocess
import json

from setuptools import setup, Extension


# Call cmake to generate json file with include and link information about dolfin and pybind11
if not os.path.isfile(os.path.join("build", "config.json")) :
    if not os.path.exists("build") :
        os.mkdir("build")
    subprocess.check_call(["cmake", os.getcwd()], cwd=os.path.abspath("build"))

with open(os.path.join("build", "config.json"), 'r') as infile :
    config = json.load(infile)

include_dirs = config["pybind11"]["include_dirs"].split(";") + \
               config["dolfin"]["include_dirs"].split(";")   + \
               config["mshr"]["include_dirs"].split(";")


mshr_ext = Extension('mshr.cpp',
                     ['src/mshr.cpp'],
                     include_dirs=include_dirs,
                     library_dirs=config['mshr']['lib_dirs'].split(";"),
                     libraries=config['mshr']['libs'].split(";"),
                     extra_compile_args=['-std=c++11'],
                     language='c++11')


setup(name             = 'mshr',
      version          = '2019.1.0',
      author           = 'FEniCS Project',
      description      = 'mshr python interface (via pybind11)',
      long_description = '',
      packages         = ["mshr",],
      ext_modules      = [mshr_ext],
      install_requires = ["numpy", "fenics-dolfin"]
      #zip_safe         = False)
      )
