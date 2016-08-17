#!/usr/bin/env python

"""
Set up script for icqsol
alexander@gokliya.net
"""

import os
import glob
from setuptools import setup, Extension
import re
import vtk # for version number
import __init__ # for version number

def anacondaLibraryNameFix(vtkLibraries, vtkLibraryDir):
    """
    On Anaconda, some of the libraries only exist with the vtk version 
    appended.
    """
    vtkVersion = '-{0}.{1}'.format(vtk.VTK_MAJOR_VERSION, vtk.VTK_MINOR_VERSION)
    suffix = '.so'
    if os.uname()[0] == 'Darwin':
        suffix = '.dylib'
    elif os.uname()[0] == 'Windows':
        # I don't have a windows machine to test...
        # do we have to worry about capitalization on Windows?
        suffix = '.dll'
    count = 0
    for lib in vtkLibraries:
        if not os.path.exists(vtkLibraryDir + '/lib' + lib + suffix):
            newName = lib + vtkVersion
            li = vtkLibraryDir + '/lib' + newName + suffix
            if os.path.exists(vtkLibraryDir + '/lib' + newName + suffix) or \
               os.path.exists(vtkLibraryDir + '/lib' + newName + '.a') or \
               os.path.exists(vtkLibraryDir + '/' + newName + suffix):
                # set the loibrary name to name-MAJOR_VERSION.MINOR_VERSION
                vtkLibraries[count] = newName
        count += 1

# Because we're linking C++ code to the VTK library, we 
# need to know where VTK was installed. 
VTK_INCLUDE_DIRS = os.environ.get('VTK_INCLUDE_DIRS', '').split(';')
VTK_LIBRARIES = os.environ.get('VTK_LIBRARIES', '').split(';')
VTK_RUNTIME_LIBRARY_DIRS = os.environ.get('VTK_RUNTIME_LIBRARY_DIRS', '').split(';')
# Override the libraries if one provides VTK_LIBRARY_DIR
vtk_lib_dir = os.environ.get('VTK_LIBRARY_DIR', '')
if vtk_lib_dir:
    libfiles = glob.glob(vtk_lib_dir + '/libvtk*')
    VTK_LIBRARIES = []
    VTK_RUNTIME_LIBRARY_DIRS = [vtk_lib_dir]
    for libfile in libfiles:
        # Skip libraries ending with version number
        if re.search(r'\.\d+$', libfile): continue
        lib = os.path.basename(libfile)
        lib = re.sub(r'^lib', '', lib)
        lib = re.sub(r'\.\w+$', '', lib) 
        # exclude python libraries
        if re.search('Python', lib): continue
        VTK_LIBRARIES.append(lib) 

# Anaconda renamed the libraries
anacondaLibraryNameFix(VTK_LIBRARIES, VTK_RUNTIME_LIBRARY_DIRS[0])

# remove all the Python libraries since they are causing a link issue
# on CentOS 7 and we don't really need them
indicesToDelete = []
for i in range(len(VTK_LIBRARIES)):
  lib = VTK_LIBRARIES[i]
  if lib.find('Python') >= 0 or lib.find('verdict') >=0: 
    indicesToDelete.append(i)
indicesToDelete.reverse()

for i in indicesToDelete:
  del VTK_LIBRARIES[i]

print('VTK_INCLUDE_DIRS         = {0}\n'.format(VTK_INCLUDE_DIRS))
print('VTK_LIBRARIES            = {0}\n'.format(VTK_LIBRARIES))
print('VTK_RUNTIME_LIBRARY_DIRS = {0}\n'.format(VTK_RUNTIME_LIBRARY_DIRS))

setup(name='icqsol',
      version=__init__.__version__, 
      description='Solving engineering problems on the web',
      author='Alex Pletzer and Greg Von Kuster',
      author_email='alexander@gokliya.net',
      url='https://github.com/pletzer/icqsol/wiki',
      install_requires=['pytriangle>=1.0.9', 'pycsg>=0.3.12'],
      dependency_links=['http://github.com/pletzer/pytriangle/tarball/master#egg=pytriangle-1.0.9',
                        'http://github.com/pletzer/pycsg/tarball/master#egg=pycsg-0.3.12',
      ],
      package_dir = {'icqsol': ''}, # the present working directory maps to icqsol below
      data_files = [('icqsol/textures', ['textures/Swietenia_macrophylla_wood.jpg',
                                  'textures/220px-COnglomerate-sandstone_layers_Nerriga.jpg',
                                  'textures/checkerboard.png',])],
      packages=['icqsol',
                'icqsol.color', 
                'icqsol.bem',
                'icqsol.solvers',
                'icqsol.discretization',
                'icqsol.shapes',
                'icqsol.util'],
      ext_modules = [Extension('icqsol.icqLaplaceMatricesCpp', # name of the shared lib
                               ['bem/icqFunctor.cpp',
                                'bem/icqLaplaceFunctor.cpp',
                                'bem/icqQuadrature.cpp',
                                'bem/icqLaplaceMatrices.cpp'] + glob.glob('bem/*.h'),
                               include_dirs=['bem'] + VTK_INCLUDE_DIRS,
                               library_dirs=VTK_RUNTIME_LIBRARY_DIRS,
                               libraries=VTK_LIBRARIES,
                               ),
                      Extension('icqsol.icqInsideLocatorCpp', 
                                ['csg/icqInsideLocator.cpp'] + glob.glob('csg/*.h'),
                                include_dirs=['csg'] + VTK_INCLUDE_DIRS,
                                library_dirs=VTK_RUNTIME_LIBRARY_DIRS,
                                libraries=VTK_LIBRARIES,
                                ),
      ],
      requires = ['numpy', 'vtk',],
     )
