# Schroedinger Equation
***
This is a library to solve the stationary schrodinger equation. As a eigenvalue solver the jacobi method is used.

## Pre-Requisits
* `eigen-3.4.0` (not implemented)
* `CMAKE` 
* `C++17 Compiler`
* `python3`
* `python3-pip`
* `pybind11`

## Supported Platforms
Currently, only Linux is supported and tested for. If you find problems for Windows/Mac please let us know how you solved them!
***

## Building the python-module 

* `mkdir build` -create a build directory
* `cd build` -change directory to build folder
* `cmake .. ` -create makefiles
* `make` -build python library
* `cd ../lib` - change to lib directory
* `pip install -e .` - install python library

***

## Python example with viennaio library

Interactive plot of the the wavefunction in the infinite potential

* `cd test` - change to example dir
* `python3 plot.py` - run testfile
***

