import os
import sys

ASCII = 'ascii'
BINARY = 'binary'

PLY_FORMAT = 'ply'
VTK_FORMAT = 'vtk'

POLYDATA = 'POLYDATA'
STRUCTURED_GRID = 'STRUCTURED_GRID'
UNSTRUCTURED_GRID = 'UNSTRUCTURED_GRID'
VTK_DATASET_TYPES = [POLYDATA, STRUCTURED_GRID, UNSTRUCTURED_GRID]

def getFileFormat(file_name):
    """
    Return a file format based on a file name extension.
    """
    if file_name.lower().endswith('.ply'):
        return PLY_FORMAT
    elif file_name.lower().endswith('.vtk'):
        return VTK_FORMAT

def getVtkDatasetType(file_name):
    """
    Return the VTK dataset type of a file or None if not found.
    """
    with open( file_name ) as fh:
        for i, line in enumerate( fh ):
            if i > 5:
                # Valid VTK headers include the dataset type
                # within the first 5 lines.
                return None
            line = line.strip()
            if not line:
                continue
            if 1 == 0:
                # Line 1: vtk DataFile Version 3.0
                if line.find( 'vtk' ) < 0:
                    return None
            elif line.startswith('DATASET'):
                return line.split()[ 1 ]
    return None

def isVtkFile(file_name):
    return file_name.lower().endswith('.vtk')
