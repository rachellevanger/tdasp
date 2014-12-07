# TDA Statistical Package : Project
#
# Copyright (C) 2014 Rachel Levanger
# Author: Rachel Levanger <rachel.levanger@gmail.com>
#                         <rachel@math.rutgers.edu>
# URL: <http://www.rachellevanger.com/math/TDASP>

import numpy as np
import collections as c
from lib import persistenceDiagram as pd
from numpy.lib import recfunctions as rf


######################################################################
## Table of Contents
######################################################################
## - Constants
## - Data Classes
##   - Project
## - 

######################################################################
## Constants
######################################################################



######################################################################
## Data Classes
######################################################################

class ProjectError(Exception):
    """An exception class for project-related errors."""

class Project(object):
    """
    The top-level data container for a TDA Statistical Project.

    Create a Project by supplying a name and a description.

    Project attributes:

    - name: The name of this project.
    - description: A detailed description of this project.
    - data: Collections of data associated to the project.

    Project methods:

    A Project only has methods pertaining to loading and retrieving
    data from either the datapoints or the distance matrices. All
    statistical analyses are done at the data level.

    loadDistanceMatrix - Attaches a distance matrix to a project.
    
    """

    ######################################################################
    ## top-level methods
    ######################################################################
    def __init__(self, name="Name", description="Description", dirOut="", dirPD="", dirLD=""):
        self.name = name
        self.description = description
        self._numDataPoints = 0
        self._dataPoints = []
        self.dirPD = dirPD
        self.fileFormatPD = ''
        self.dirLD = dirLD
        self._persistenceDiagrams = [] # Array of persistence diagram objects, in order of point idx field
        self.hom_max = None
        self.dirOut = dirOut

    def __repr__(self):
        tup = type(self).__name__, self.name, self.description
        return "%s('%s : %s')" % tup

    ######################################################################
    ## accessor methods 
    ######################################################################

    def numDataPoints(self):
        return self._numDataPoints

    def dataPoints(self):
        return self._dataPoints

    def distanceMatrix(self):
        return self._distanceMatrix

    def distanceMatrixDesc(self):
        return self._distanceMatrixDesc

    def persistenceDiagrams(self):
        return self._persistenceDiagrams

    ######################################################################
    ## functional methods
    ######################################################################

    def loadDistanceMatrix(self, pathToDistanceMatrix, description="Description", delim=" "):
        # Associate a distance matrix and its description to the project. Initialize the
        # data points record array with the number of dimensions in the matrix.
        tmpDistanceMatrix = np.genfromtxt(pathToDistanceMatrix, delimiter=delim)

        if self._checkSizeOkay(tmpDistanceMatrix.shape[0]):
            self._distanceMatrix = tmpDistanceMatrix
            self._distanceMatrixDesc = description
            if self._numDataPoints == 0:
                self._numDataPoints = self._distanceMatrix.shape[0]
                self._dataPoints = self._initDataPoints(self._numDataPoints)
        else:
            raise ProjectError( ('Current project size is %d points, but distance matrix is '
                           '%d points. Dimensions must match to load new distance matrix.') % \
                               (self._numDataPoints, tmpDistanceMatrix.shape[0]) )

    def loadDataAttributesFromFile(self, pathToDataFile, label, datatype):
        # Currently only handles importing a vector of data
        dataAttributes = np.genfromtxt(pathToDataFile, dtype=datatype)

        if self._checkSizeOkay(dataAttributes.shape[0]):
            if self._numDataPoints == 0:
                self._numDataPoints = dataAttributes.shape[0]
                self._dataPoints = self._initDataPoints(self._numDataPoints)
            self._dataPoints = rf.append_fields(self._dataPoints, label, dataAttributes, dtypes=datatype)
        else:
            raise ProjectError( ('Current project size is %d points, but attribute array is '
                           '%d points. Dimensions must match to load data attributes.') % \
                               (self._numDataPoints, dataAttributes.shape[0]) )

    def loadDataAttributesFromArray(self, dataAttributes, label, datatype):
        # Currently only handles importing a vector of data
        if self._checkSizeOkay(dataAttributes.shape[0]):
            dataAttributes.tolist()
            if self._numDataPoints == 0:
                self._numDataPoints = dataAttributes.shape[0]
                self._dataPoints = self._initDataPoints(self._numDataPoints)
            self._dataPoints = rf.append_fields(self._dataPoints, label, dataAttributes, dtypes=datatype)
        else:
            raise ProjectError( ('Current project size is %d points, but attribute array is '
                           '%d points. Dimensions must match to load data attributes.') % \
                                (self._numDataPoints, dataAttributes.shape[0]) )

    def loadDataAttributesFromRecordArray(self, dataAttributes, label):
        # Currently only handles importing a vector of data
        if self._checkSizeOkay(dataAttributes.shape[0]):
            if self._numDataPoints == 0:
                self._numDataPoints = dataAttributes.shape[0]
                self._dataPoints = self._initDataPoints(self._numDataPoints)
            self._dataPoints = rf.append_fields(self._dataPoints, label, dataAttributes)
        else:
            raise ProjectError( ('Current project size is %d points, but attribute array is '
                           '%d points. Dimensions must match to load data attributes.') % \
                                (self._numDataPoints, dataAttributes.shape[0]) )

    def removeDataAttributes(self, labels):
        # Remove the fields passed as a parameter to the function
        self._dataPoints = rf.drop_fields(self._dataPoints, labels)

    def loadPersistenceDiagrams(self):
        if self.hom_max == None:
            raise ProjectError( 'Either hom_min or hom_max is empty. Both variables must be set to load persistence diagrams.')
        if self.dirPD == "":
            raise ProjectError( 'Path to persistence diagrams not set.' )
        if self.fileFormatPD == "":
            raise ProjectError( 'File format for names of persistence diagrams not set.' )

        # Remove any existing persistence diagrams
        self._persistenceDiagrams = []

        for hom in range(0,self.hom_max+1):
            print('Loading diagrams for dimension ' + str(hom) + '.')
            hom_diagrams = []
            for j in range(0,self._numDataPoints):
                pathToPersistenceDiagram = self.dirPD + (self.fileFormatPD % (self._dataPoints[j]['fileName'], hom))
                hom_diagrams.append(pd.Diagram(pathToPersistenceDiagram))
            self._persistenceDiagrams.append(hom_diagrams)

        print('Diagrams loaded.')

    ######################################################################
    ## utility methods
    ######################################################################

    def _checkSizeOkay(self, size):
        if self._numDataPoints == 0:
            return 1
        elif self._numDataPoints == size:
            return 1
        else:
            return 0

    def _initDataPoints(self, size):
        return np.array(np.arange(0,size,1), dtype=[('idx','int')])



