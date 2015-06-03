#!/usr/bin/env python

from distutils.core import setup

setup(name='ICQSol',
      version='0.1',
      description='Solving engineering problems on the web',
      author='Alex Pletzer',
      author_email='pletzer@psu.edu',
      url='https://github.com/pletzer/icqsol/wiki',
      package_dir = {'icqsol': ''}, # the present working directory maps to icqsol below
      packages=['icqsol.tools', 'icqsol.tools.geometry'],
     )
