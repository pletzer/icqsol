#!/usr/bin/env python

"""
Set up script for icqsol
pletzer@psu.edu
"""

from setuptools import setup

setup(name='icqsol',
      version='0.3.1',
      description='Solving engineering problems on the web',
      author='Alex Pletzer',
      author_email='alexander@gokliya.net',
      url='https://github.com/pletzer/icqsol/wiki',
      install_requires=['pytriangle>=1.0.6', 'pycsg>=0.2.1'],
      dependency_links=['http://github.com/pletzer/pytriangle/tarball/master#egg=pytriangle-1.0.3',
                        'http://github.com/pletzer/pycsg/tarball/master#egg=pycsg-0.2.1',
      ],
      package_dir = {'icqsol': ''}, # the present working directory maps to icqsol below
      data_files = [('icqsol/textures', ['textures/Swietenia_macrophylla_wood.jpg',
                                  'textures/220px-COnglomerate-sandstone_layers_Nerriga.jpg',
                                  'textures/checkerboard.png',])],
      packages=['icqsol',
                'icqsol.color', 
                'icqsol.bem',
                'icqsol.discretization',
                'icqsol.shapes',
                'icqsol.util'],
      requires = ['numpy', 'vtk',],
     )
