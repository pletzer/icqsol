cmake_minimum_required (VERSION 2.8)
project(icqsol)

# need python
find_package(PythonInterp REQUIRED)

# make sure the numpy python module can be imported
execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c "import numpy" ERROR_VARIABLE ERROR_NUMPY)
if (NOT ERROR_NUMPY STREQUAL "")
  message(FATAL_ERROR "You must have python-numpy installed!")
endif()

# make sure the vtk python module can be imported
execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c "import vtk" ERROR_VARIABLE ERROR_VTK)
if (NOT ERROR_VTK STREQUAL "")
  message(FATAL_ERROR "You must have python-vtk installed!")
endif()

option(INSTALL_USER "Install in user directory" OFF)
option(INSTALL_PREFIX "Install under CMAKE_INSTALL_PREFIX" OFF)

# run python setup.py install
if (NOT INSTALL_USER AND NOT INSTALL_PREFIX)
  # let python choose where to install, may require sudo access
  execute_process(COMMAND "${PYTHON_EXECUTABLE}" setup.py install)
elseif (INSTALL_USER)
  # install in user directory
  execute_process(COMMAND "${PYTHON_EXECUTABLE}" setup.py install --user)
else ()
  # user specifies where to install, must set PYTHONPATH
  message("\n\n*** Installing under ${CMAKE_INSTALL_PREFIX}. Be sure to set PYTHONPATH, e.g. ***")
  message("export PYTHONPATH=${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages:$PYTHONPATH\n\n")
  execute_process(COMMAND "${PYTHON_EXECUTABLE}" setup.py install --prefix=${CMAKE_INSTALL_PREFIX})
endif()

# enable ctest
enable_testing()

#
# class tests
#

set(TESTS_DIR "${CMAKE_SOURCE_DIR}/tests/tools")
add_test(NAME testRing COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/geometry/testRing.py")
add_test(NAME testScrew COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/geometry/testScrew.py")

#
# tool tests
#

# create a simple spherical object
add_test(NAME testBuildGeometry1 COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tools/geometry/buildGeometry.py"
         --create "s = Sphere(radius=0.5, origin=(0.,0.,0.))"
         --assemble "s"
         --output testBuildGeometry1.ply
         )

# a cube with a notch and a cylindrical hole
add_test(NAME testBuildGeometry2 COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tools/geometry/buildGeometry.py"
         --create "cube = Box(loBound=(0.,0.,0.), hiBound=(1.,1.,1.))"
                  "sphere = Sphere(radius=0.2, origin=(1., 1., 1.))"
                  "cyl = Cylinder(radius=0.1, origin=(0.1, 0.2, 0.3), length=1.1)"
         --assemble "cube - sphere - cyl"
         --output testBuildGeometry2.ply
         )




