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

set(TESTS_DIR "${CMAKE_SOURCE_DIR}/tests")
add_test(NAME testCreateCylinder 
         COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/geometry/testCreateCylinder.py" --output testCreateCylinder.vtk)
add_test(NAME testCreateBox 
         COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/geometry/testCreateBox.py" --output testCreateBox.vtk)
add_test(NAME testCreateSphere 
         COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/geometry/testCreateSphere.py" --output testCreateSphere.vtk)
add_test(NAME testCreateCone
         COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/geometry/testCreateCone.py" --output testCreateCone.vtk)
add_test(NAME testSaveLoad COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/geometry/testSaveLoad.py")
add_test(NAME testToPolygons COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/geometry/testToPolygons.py")
add_test(NAME testFromPolygons COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/geometry/testFromPolygons.py")
add_test(NAME testUnion COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/geometry/testUnion.py")
add_test(NAME testRing COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/geometry/testRing.py")
add_test(NAME testBolt COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/geometry/testBolt.py")
add_test(NAME testSurfaceInsideOtherShape COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/geometry/testSurfaceInsideOtherShape.py")
add_test(NAME testPointInsideShape COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/geometry/testPointInsideShape.py" --n_theta 6 --n_phi 4)
add_test(NAME testPointInsideShape2 COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/geometry/testPointInsideShape.py" --n_theta 32 --n_phi 16)

add_test(NAME testColorMap
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/color/testColorMap.py")
add_test(NAME testPrimitiveSphere 
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/createPrimitiveShape.py" 
         --output testPrimitiveSphere.vtk 
         --type "sphere" --ascii 
         --options "radius=1" "origin=(0, 1, 2)")
add_test(NAME testPrimitiveSphere2x4 
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/createPrimitiveShape.py" 
         --output testPrimitiveSphere4x4.vtk 
         --type "sphere" --ascii 
         --options "radius=1.0" "origin=(0., 0., 0.)" "n_theta=4" "n_phi=4")
add_test(NAME testCylOuter
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/createPrimitiveShape.py" 
         --output cylOuter.ply 
         --type "cylinder"  
         --options "radius=1.0" 
                   "origin=(0., 0., 0.)" 
                   "lengths=(1., 0., 0.)"
                   "n_theta=8")
add_test(NAME testCylInner
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/createPrimitiveShape.py" 
         --output cylInner.ply 
         --type "cylinder"  
         --options "radius=0.5" 
                   "origin=(-0.1, 0., 0.)" 
                   "lengths=(2., 0., 0.)"
                   "n_theta=3")
add_test(NAME testPrimitiveCylinder
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/createPrimitiveShape.py" 
         --output testPrimitiveCylinder.ply
         --type "cylinder"  
         --options "radius=0.5" 
                   "origin=(0.5, 0.2, 0.3)" 
                   "lengths=(0.35, 0., 0.)"
                   "n_theta=16")
add_test(NAME testPrimitiveCone
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/createPrimitiveShape.py" 
         --output testPrimitiveCone.vtk 
         --type "cone" --ascii 
         --options "radius=0.5" 
                   "origin=(0.5, 0.2, 0.3)" 
                   "lengths=(0.35, 0.0, 0.0)")
add_test(NAME testPrimitiveBox
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/createPrimitiveShape.py" 
         --output testPrimitiveBox.vtk 
         --type "box" --ascii 
         --options "origin=(0.5, 0.2, 0.3)" 
                   "lengths=(0.7, 0.8, 0.9)")
add_test(NAME testCompositeUnion
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/composeShapes.py" 
         --shapeTuples "s1,testPrimitiveSphere.vtk" "s2,testPrimitiveCylinder.ply"
         --compose "s1 + s2"
         --output testCompositeUnion.vtk)
add_test(NAME testCompositeIntersection
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/composeShapes.py"
         --shapeTuples "s1,testPrimitiveSphere.vtk" "s2,testPrimitiveCylinder.ply"
         --compose "s1 * s2"
         --output testCompositeIntersection.vtk)
add_test(NAME testRingCoarse
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/composeShapes.py" 
         --shapeTuples "s1, cylOuter.ply" "s2, cylInner.ply"
         --compose "s1 - s2"
         --ascii
         --output ringCoarse.vtk)
add_test(NAME testCompositeRemoval
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/composeShapes.py" 
         --shapeTuples "a, testPrimitiveSphere.vtk" "b, testPrimitiveCylinder.ply"
         --compose "a - b"
         --output testCompositeRemoval.vtk)
add_test(NAME testTranslate
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/translateShape.py" 
         --input testPrimitiveBox.vtk
         --translate "1, 2, 3"
         --output testTranslatedBox.vtk)
add_test(NAME testRotate
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/rotateShape.py" 
         --input testPrimitiveBox.vtk
         --angle 45.0 --axis "1,0,0"
         --output testRotatedBox.vtk)
add_test(NAME testRotate2
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/rotateShape.py" 
         --input testPrimitiveCylinder.ply
         --angle 45.0 --axis "1,0,0"
         --output testRotatedCylinder.ply)
add_test(NAME testRefine
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/refineShape.py"
         --input testPrimitiveSphere4x4.vtk
         --refine 2
         --output testPrimitiveSphere8x8.vtk
         --ascii)
add_test(NAME testAddSurfaceFieldFromExpression
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/addSurfaceFieldFromExpression.py"
         --expression "1 + x**2*sin(-log(0.1+abs(x*y)*pi))"
         --name "my_field"
         --input testPrimitiveCylinder.ply
         --output testPrimitiveCylinderWithField.vtk
         --ascii)
add_test(NAME testAddSurfaceFieldFromExpressionTimeDependent
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/addSurfaceFieldFromExpression.py"
         --expression "1 + x**2*sin(-log(0.1+abs(x*y)*pi))*sin(pi*t + 0.3)"
         --name "my_field"
         --times "linspace(0, 1., 4)"
         --input testPrimitiveCylinder.ply
         --output testPrimitiveCylinderWithTimeDependentField.vtk
         --ascii)
add_test(NAME testColorSurfaceFieldGnu
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/colorSurfaceField.py"
         --colormap "gnu"
         --name "my_field"
         --input testPrimitiveCylinderWithField.vtk
         --output testPrimitiveCylinderWithFieldColored.vtk
         --ascii)




