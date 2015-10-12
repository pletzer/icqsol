#!/usr/bin/env python

"""
Set up script for icqsol
pletzer@psu.edu
"""

from distutils.core import setup

setup(name='icqsol',
      version='0.1',
      description='Solving engineering problems on the web',
      author='Alex Pletzer',
      author_email='alexander@gokliya.net',
      url='https://github.com/pletzer/icqsol/wiki',
      package_dir = {'icqsol': ''}, # the present working directory maps to icqsol below
      packages=['icqsol',
                'icqsol.color', 
                'icqsol.discretization',
                'icqsol.shapes',
                'icqsol.util'],
      requires = ['numpy', 'vtk', 'csg', 'triangle',],
     )
