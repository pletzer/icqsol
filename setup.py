#!/usr/bin/env python

"""
Set up script for icqsol
pletzer@psu.edu
"""

import os
from setuptools import setup, Extension


def getVtkLibraryDirs():
    """
    Guess the path of the VTK libraries 
    @return list of paths
    """
    
    # Sniff LD_LIBRARY_PATH (or equivalent)
    ld_library_path_name = 'LD_LIBRARY_PATH'
    separator = ':'
    shared_library_suffix = 'so'
    uname = os.uname()[0]
    
    if uname == 'Darwin':
        ld_library_path_name = 'DYLD_LIBRARY_PATH'
        shared_library_suffix = 'dylib'
    elif uname == 'Windows':
        # untested!
        ld_library_path_name = 'PATH'
        separator = ';'
        shared_library_suffix = 'dll'

    library_dirs = []
    for lib_dir in os.environ.get(ld_library_path_name, '').split(separator):
        if os.path.exists(lib_dir + '/libvtkCommon.' + shared_library_suffix):
            library_dirs.append(lib_dir)

    return library_dirs

def getVtkIncludeDirs():
    """
    Guess the path of the VTK include files
    @return list of paths
    """

    library_dirs = getVtkLibraryDirs()
    include_dirs = []
    for path in library_dirs:

        if len(path) > 0 and path[-1] == '/': 
            path = path[:-1]

        bn = os.path.basename(path)
        if os.path.exists(path + '/../include'):
            include_dirs.append(path + '/../include')
        elif os.path.exists(path + '/../../include/' + bn):
            include_dirs.append(path + '/../../include/' + bn)
    return include_dirs

setup(name='icqsol',
      version='0.3.8',
      description='Solving engineering problems on the web',
      author='Alex Pletzer and Greg Von Kuster',
      author_email='alexander@gokliya.net',
      url='https://github.com/pletzer/icqsol/wiki',
      install_requires=['pytriangle>=1.0.6', 'pycsg>=0.2.3'], # should really be 'pycsg>=0.3.3'],
      dependency_links=['http://github.com/pletzer/pytriangle/tarball/master#egg=pytriangle-1.0.6',
                        'http://github.com/pletzer/pycsg/tarball/master#egg=pycsg-0.3.3',
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
      ext_modules = [Extension('icqsol.icqQuadratureCpp', # name of the shared lib
                               ['bem/icqFunctor.cpp',
                                'bem/icqLaplaceFunctor.cpp',
                                'bem/icqQuadrature.cpp'],
                               include_dirs=['bem'],
                               ),
                    Extension('icqsol.icqLaplaceMatricesCpp', 
                              ['bem/icqLaplaceMatrices.cpp'],
                              include_dirs=['bem'] + getVtkIncludeDirs(),
                              library_dirs=getVtkLibraryDirs(),
                              libraries=['vtkCommon', 'vtkFiltering'],
                              )
      ],
      requires = ['numpy', 'vtk',],
     )
