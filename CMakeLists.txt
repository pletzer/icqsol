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

# run python setup.py install
execute_process(COMMAND "${PYTHON_EXECUTABLE}" ${CMAKE_SOURCE_DIR}/setup.py install --user)

# enable ctest
enable_testing()

#
# tests
#
set(TESTS_DIR "${CMAKE_SOURCE_DIR}/tests/tools")
add_test(NAME testRing COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/geometry/testRing.py")
add_test(NAME testScrew COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/geometry/testScrew.py")






