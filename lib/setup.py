#!/usr/bin/python3

from setuptools import setup,find_packages

setup(
    name='Schroedinger',
    version='0.1.0',
    description="ViennaIO",
    long_description='Semiconductor Monte Carlo simulation package',
    packages=['viennaio'], #get_active_folders(), #find_packages(),
    #packages=find_packages(),
    package_dir={'Schroedinger': 'Schroedinger'},
    #install_requires=['numpy', 'scipy', 'matplotlib'], 
    url='https://tcl-git.iue.tuwien.ac.at/TCAD/viennamc',
)