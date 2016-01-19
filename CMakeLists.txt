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
option(VTK_INCLUDE_DIR "Directory where VTK include files are located" "")
option(VTK_LIBRARY_DIR "Directory where VTK libraries are located" "")

# run python setup.py install
if (NOT INSTALL_USER AND NOT INSTALL_PREFIX)
  # let python choose where to install, may require sudo access
  execute_process(COMMAND VTK_INCLUDE_DIR="${VTK_INCLUDE_DIR}" VTK_LIBRARY_DIR="${VTK_LIBRARY_DIR}"
                          "${PYTHON_EXECUTABLE}" setup.py install)
elseif (INSTALL_USER)
  # install in user directory
  execute_process(COMMAND VTK_INCLUDE_DIR="${VTK_INCLUDE_DIR}" VTK_LIBRARY_DIR="${VTK_LIBRARY_DIR}" 
                          "${PYTHON_EXECUTABLE}" setup.py install --user)
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
add_test(NAME testInsideLocatorCpp
         COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/testInsideLocatorCpp.py")
add_test(NAME testQuadrature 
         COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/testQuadrature.py")
add_test(NAME testQuadratureCpp
         COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/testQuadratureCpp.py")
add_test(NAME testLaplace 
         COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/testLaplace.py")
add_test(NAME testCreateCylinder 
         COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/testCreateCylinder.py" --output testCreateCylinder.vtk)
add_test(NAME testCreateBox 
         COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/testCreateBox.py" --output testCreateBox.vtk)
add_test(NAME testCreateSphere 
         COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/testCreateSphere.py" --output testCreateSphere.vtk)
add_test(NAME testCreateCone
         COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/testCreateCone.py" --output testCreateCone.vtk)
add_test(NAME testSaveLoad COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/testSaveLoad.py")
add_test(NAME testToPolygons COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/testToPolygons.py")
add_test(NAME testFromPolygons COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/testFromPolygons.py")
add_test(NAME testUnion COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/testUnion.py")
add_test(NAME testRing COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/testRing.py")
add_test(NAME testBolt COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/testBolt.py")
add_test(NAME testSurfaceInsideOtherShape COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/testSurfaceInsideOtherShape.py")
add_test(NAME testPointInsideShape COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/testPointInsideShape.py" --n_theta 6 --n_phi 4)
add_test(NAME testPointInsideShape2 COMMAND "${PYTHON_EXECUTABLE}" "${TESTS_DIR}/testPointInsideShape.py" --n_theta 32 --n_phi 16)
add_test(NAME testStructuredGridRead
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testStructuredGridRead.py")
add_test(NAME testUnstructuredGridRead
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testUnstructuredGridRead.py")
add_test(NAME testColormaps
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testColormaps.py")

add_test(NAME testColorMap
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testColorMap.py")
add_test(NAME testPrimitiveSphereOptionList 
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/createPrimitiveShape.py" 
         --type "sphere" --list)
add_test(NAME testPrimitiveConeOptionList 
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/createPrimitiveShape.py" 
         --type "cone" --list)
add_test(NAME testPrimitiveBoxOptionList 
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/createPrimitiveShape.py" 
         --type "box" --list)
add_test(NAME testPrimitiveCylinderOptionList 
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/createPrimitiveShape.py" 
         --type "cylinder" --list)
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
add_test(NAME textAddTextureCheckerboard
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/addTexture.py" 
         --texture "${CMAKE_SOURCE_DIR}/textures/checkerboard.png"
         --maxedge 0.1
         --input testPrimitiveBox.vtk 
         --output testPrimitiveBoxCheckerboard.vtk
         --ascii)
add_test(NAME textAddTextureWood
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/addTexture.py" 
         --texture "${CMAKE_SOURCE_DIR}/textures/Swietenia_macrophylla_wood.jpg"
         --maxedge 0.1
         --input testPrimitiveBox.vtk 
         --output testPrimitiveBoxWood.vtk
         --ascii)
add_test(NAME textAddTextureWoodSpecifyFormat
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/addTexture.py" 
         --texture "${CMAKE_SOURCE_DIR}/textures/Swietenia_macrophylla_wood.jpg"
         --texture_file_format jpeg
         --maxedge 0.1
         --input testPrimitiveBox.vtk 
         --output testPrimitiveBoxWood2.vtk
         --ascii)
add_test(NAME textAddTextureStone
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/addTexture.py" 
         --texture "${CMAKE_SOURCE_DIR}/textures/220px-COnglomerate-sandstone_layers_Nerriga.jpg"
         --maxedge 0.01
         --input testPrimitiveBox.vtk 
         --output testPrimitiveBoxStone.vtk
         --ascii)
# should fail because we don’t support gif format
add_test(NAME textAddTextureGif
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/addTexture.py" 
         --texture "${CMAKE_SOURCE_DIR}/textures/sample.gif"
         --maxedge 0.01
         --input testPrimitiveBox.vtk 
         --output testPrimitiveBoxGif.vtk
         --ascii)
set_tests_properties(textAddTextureGif PROPERTIES WILL_FAIL TRUE)
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
add_test(NAME testRotate3
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/rotateShape.py" 
         --input testRemove2.vtk
         --angle 45.0 --axis "0, 0, 1"
         --output testRotatedRemoved.vtk)
add_test(NAME testScale
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/scaleShape.py" 
         --input testRemove2.vtk
         --scale "1.3, 2.1, 4.5"
         --output testScaledRemoved.vtk)
add_test(NAME testRefine
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/refineShape.py"
         --input testPrimitiveSphere4x4.vtk
         --refine 2
         --output testPrimitiveSphere8x8.vtk
         --ascii)
add_test(NAME testRefineSurface
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/refineSurface.py"
         --input testBolt.vtk
         --maxedge 0.05
         --output testBoltFine.vtk
         --ascii)
add_test(NAME testRefineSurface2
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/refineSurface.py"
         --input testBoltFine.vtk
         --maxedge 0.05
         --output testBoltFine2.vtk
         --ascii)
add_test(NAME testRefineSurfaceInfiniteEdgeCondition
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/refineSurface.py"
         --input testBoltFine.vtk
         --output testBoltFine2.vtk
         --ascii)
add_test(NAME testAddSurfaceFieldFromExpressionBolt
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/addSurfaceFieldFromExpression.py"
         --expression "exp(-0.1*t)*cos(6*pi*x)*cos(10*pi*sqrt(y**2 + z**2))/sqrt(y**2 + z**2 + 0.1)"
         --times "linspace(0, 1.0, 10)"
         --name "wave"
         --input testBoltFine.vtk
         --output testBoltFineWithField.vtk)
add_test(NAME testAddSurfaceFieldFromExpressionBoltCell
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/addSurfaceFieldFromExpression.py"
         --expression "exp(-0.1*t)*cos(6*pi*x)*cos(10*pi*sqrt(y**2 + z**2))/sqrt(y**2 + z**2 + 0.1)"
         --times "linspace(0, 1.0, 10)"
         --name "wave"
         --location "cell"
         --input testBoltFine.vtk
         --output testBoltFineWithFieldCell.vtk)
add_test(NAME testColorBoltCell
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/colorSurfaceField.py"
         --input testBoltFineWithFieldCell.vtk 
         --name "wave"
         --output testColorBoltCell.vtk --ascii)
add_test(NAME testIntegrateCellField
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/integrateSurfaceField.py"
         --input testBoltFineWithFieldCell.vtk 
         --name "wave")
add_test(NAME testIntegratePointField
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/integrateSurfaceField.py"
         --input testBoltFineWithField.vtk 
         --name "wave")
add_test(NAME testAddSurfaceFieldFromExpressionBolt2
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/addSurfaceFieldFromExpression.py"
         --expression "1.0"
         --name "ones"
         --input testBoltFineWithField.vtk
         --output testBoltFineWithField2.vtk
         --ascii)
add_test(NAME testAddSurfaceFieldFromExpressionBolt3
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/addSurfaceFieldFromExpression.py"
         --expression "2.0"
         --name "twos"
         --input testBoltFineWithField2.vtk
         --output testBoltFineWithField3.vtk
         --ascii)
add_test(NAME testAddSurfaceFieldFromExpressionSmallPolygon
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/addSurfaceFieldFromExpression.py"
         --expression "exp(-0.1*t)*cos(6*pi*x)*cos(10*pi*sqrt(y**2 + z**2))/sqrt(y**2 + z**2 + 0.1)"
         --times "linspace(0, 1.0, 10)"
         --name "wave"
         --input "${CMAKE_SOURCE_DIR}/data/surface_field.vtk"
         --output testAddSurfaceFieldFromExpressionSmallPolygon.vtk)
add_test(NAME testColorBoltComponent0
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/colorSurfaceField.py"
         --name "wave"
         --component 0
         --colormap gnu
         --input testBoltFineWithField.vtk
         --ascii 
         --output testBoltFineColored_0.vtk)
add_test(NAME testColorBoltComponent1
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/colorSurfaceField.py"
         --name "wave"
         --component 1
         --colormap gnu
         --input testBoltFineWithField.vtk
         --ascii
         --output testBoltFineColored_1.vtk)
add_test(NAME testAddSurfaceFieldFromExpression
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/addSurfaceFieldFromExpression.py"
         --expression "1 + x**2*sin(-log(0.1+abs(x*y)*pi))"
         --name "my_field"
         --input testPrimitiveCylinder.ply
         --output testPrimitiveCylinderWithField.vtk
         --ascii)
add_test(NAME testAddSurfaceFieldFromBooleanExpression
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/addSurfaceFieldFromExpression.py"
         --expression "(z/4.)*(sqrt(x**2+y**2) __gt__ 1 or (z __le__ 2))"
         --name "my_field"
         --input testPrimitiveCylinder.ply
         --output testPrimitiveCylinderWithBooleanField.vtk
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
add_test(NAME testAddSurfaceFieldFromExpressionPlate
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/addSurfaceFieldFromExpression.py"
         --expression "sin(pi*x/0.1)*cos(pi*y/0.15)*exp(-(z-0.1)**2)"
         --name "my_field"
         --refine 0.02
         --input "${CMAKE_SOURCE_DIR}/data/plate.vtk"
         --output testAddSurfaceFieldFromExpressionPlate.vtk
         --ascii)
add_test(NAME testAddSurfaceFieldFromExpressionPlate2
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/addSurfaceFieldFromExpression.py"
         --expression "sin(pi*x/0.1)*cos(pi*y/0.15)*exp(-(z-0.1)**2)"
         --name "my_field"
         --refine 0.01
         --input "${CMAKE_SOURCE_DIR}/data/plate.vtk"
         --output testAddSurfaceFieldFromExpressionPlate.vtk
         --ascii)
add_test(NAME testColorSurfaceFieldGnu
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/examples/colorSurfaceField.py"
         --colormap "gnu"
         --name "my_field"
         --input testPrimitiveCylinderWithField.vtk
         --output testPrimitiveCylinderWithFieldColored.vtk
         --ascii)
add_test(NAME testRemove2
         COMMAND "${PYTHON_EXECUTABLE}"
         "${CMAKE_SOURCE_DIR}/examples/composeShapes.py"
         --shapeTuples "a1,${CMAKE_SOURCE_DIR}/data/a1.vtk" "b1,${CMAKE_SOURCE_DIR}/data/b1.vtk" 
         --compose "a1 - b1" --output testRemove2.vtk --ascii)
add_test(NAME testRefineWithField
         COMMAND "${PYTHON_EXECUTABLE}"
         "${CMAKE_SOURCE_DIR}/examples/refineSurface.py"
         --input testPrimitiveCylinderWithTimeDependentField.vtk
         --output testRefineWithField.vtk
         --maxedge 0.1
         --ascii)
add_test(NAME testRefineWithCellField
         COMMAND "${PYTHON_EXECUTABLE}"
         "${CMAKE_SOURCE_DIR}/examples/refineSurface.py"
         --input testBoltFineWithFieldCell.vtk
         --output testRefineWithCellField.vtk
         --maxedge 0.1
         --ascii)
add_test(NAME testPrimitiveBoxSolveLaplaceConst
         COMMAND time "${PYTHON_EXECUTABLE}"
         "${CMAKE_SOURCE_DIR}/examples/solveLaplace.py"
         --input testPrimitiveBox.vtk 
         --dirichlet "1.0" 
         --refine 0.5 
         --input_name "v" --output_name "surf_charge" --verbose
         --output testPrimitiveBoxSolveLaplaceConst.vtk --ascii)
add_test(NAME testPrimitiveSphereSolveLaplaceConst
         COMMAND time "${PYTHON_EXECUTABLE}"
         "${CMAKE_SOURCE_DIR}/examples/solveLaplace.py"
         --input testPrimitiveSphere.vtk 
         --dirichlet "1.0" 
         --refine 0.5
         --input_name "v" --output_name "surf_charge" --verbose
         --output testPrimitiveSphereSolveLaplaceConst.vtk --ascii)
add_test(NAME testPrimitiveSphereSolveLaplace1
         COMMAND time "${PYTHON_EXECUTABLE}"
         "${CMAKE_SOURCE_DIR}/examples/solveLaplace.py"
         --input testPrimitiveSphere.vtk 
         --dirichlet "1./sqrt(x**2 + (y-2.)**2 + (z-4.)**2)" 
         --refine 0.5 
         --input_name "v" --output_name "surf_charge" --verbose
         --output testPrimitiveSphereSolveLaplace1.vtk --ascii)
add_test(NAME testPrimitiveSphereSolveLaplaceOffset
         COMMAND time "${PYTHON_EXECUTABLE}"
         "${CMAKE_SOURCE_DIR}/examples/solveLaplace.py"
         --input testPrimitiveSphere.vtk 
         --dirichlet "1./sqrt((x+10.)**2 + (y-2.)**2 + (z-4.)**2)" 
         --refine 0.5 
         --input_name "v" --output_name "surf_charge" --verbose
         --output testPrimitiveSphereSolveLaplaceOffset.vtk --ascii)
add_test(NAME testPrimitiveSphereSolveLaplace2
         COMMAND time "${PYTHON_EXECUTABLE}"
         "${CMAKE_SOURCE_DIR}/examples/solveLaplace.py"
         --input testPrimitiveSphere.vtk 
         --dirichlet "1.0" 
         --refine 0.5 
         --input_name "v" --output_name "surf_charge" --verbose
         --output testPrimitiveSphereSolveLaplace2.vtk --ascii)

# this should be the last test
add_test(NAME testOutputFiles
         COMMAND "sh" "${CMAKE_SOURCE_DIR}/checkOutput.sh")

