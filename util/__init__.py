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
    # Open the file in binary mode (works for both ASCII and binary files)
    with open(file_name, 'rb') as fh:
        line = fh.readline()
        lineNumber = 1
        # First line: # vtk DataFile Version 3.0
        if line.find(b'vtk') < 0:
            # not a vtk file
            return None
        # Read the first 4 lines
        while lineNumber < 4:
            line = fh.readline()
            lineNumber += 1
        if line.startswith(b'DATASET'):
                return line.split()[1].decode('utf-8')
    return None

def isVtkFile(file_name):
    return file_name.lower().endswith('.vtk')
