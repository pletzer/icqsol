#!/usr/bin/env python

"""
Set up script for ICQSol
pletzer@psu.edu
"""

from distutils.core import setup

setup(name='ICQSol',
      version='0.1',
      description='Solving engineering problems on the web',
      author='Alex Pletzer',
      author_email='pletzer@psu.edu',
      url='https://github.com/pletzer/icqsol/wiki',
      package_dir = {'icqsol': ''}, # the present working directory maps to icqsol below
      packages=['icqsol', 'icqsol.tools', 
                'icqsol.tools.geometry', 'icqsol.tools.common'],
      requires = ['numpy', 'vtk'],
     )
