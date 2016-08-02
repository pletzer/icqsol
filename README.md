icqsol
======

icqsol is a collection of tools for constructing complex geometries from 
primitve shapes on the web with the goal of solving various engineering
problems including thermal and electrostatic. 

How to build icqsol 
-------------------

### Prerequisites


* python. Version 2.7 or 3.5 will work
* numpy. Version 1.9.2 will work
* VTK with python wrapping enabled. Version 6.3.0 or 7.0 will work
* CMake

### Installing icqsol

Icqsol uses a CMake build system, which will call "python setyp.py"
to install the package. You have three choices:

a) Let python choose the installation directory, typically 
   under /usr/lib/python<version>/site-packages. May require
   sudo privilege

```bash
cmake .
```

(The period at the end is critical.)

b) Install under ~/.local/lib/python<version>/site-packages/

```bash
cmake -DINSTALL_USER=ON .
```

c) Install in a user specified directory. Note that when using this 
option, the user will need to set the PYTHONPATH environment variable 
to point to <myDirectory>/lib/site-packages.

```bash
cmake -DINSTALL_PREFIX=ON -DCMAKE_INSTALL_PREFIX=<myDirectory> .
```

You may need to help icqsol find the VTK package by setting VTK_DIR

```bash
cmake -DVTK_DIR=<vtk_install_dir>/lib/cmake/vtk [other options] .
```

to point to the directory where VTKCOnfig.cmake is located.

How to test icqsol
------------------

You should be able to import icqsol in the python interpreter

```bash
python
Python 2.7.5 (default, Feb 11 2014, 07:46:25) 
[GCC 4.8.2 20140120 (Red Hat 4.8.2-13)] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> import icqsol
>>> 
```

More extensive tests can be run by typing

```bash
ctest
```


