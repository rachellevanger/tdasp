# TDA Statistical Package : Persistence
#
# Copyright (C) 2014 Rachel Levanger
# Author: Rachel Levanger <rachel.levanger@gmail.com>
#                         <rachel@math.rutgers.edu>
# URL: <http://www.rachellevanger.com/math/TDASP>

import numpy as np
from numpy.lib import recfunctions as rf

######################################################################
## Table of Contents
######################################################################
## - Constants
## - Data Classes
##   - Diagram

######################################################################
## Constants
######################################################################



######################################################################
## Data Classes
######################################################################

class DiagramError(Exception):
    """An exception class for project-related errors."""

class Diagram(object):
    """
    An object representation of a persistence diagram. Currently takes
    in file formats that are output from the Perseus application.

    """


    ######################################################################
    ## top-level methods
    ######################################################################

    def __init__(self, pathToFile=""):
        self.path = pathToFile
        self._points = np.genfromtxt(pathToFile, delimiter=' ', names='birth, death', dtype='f8, f8')

        # Generate commonly-used derived fields: lifespan, avg_coord
        self._points = rf.append_fields(self._points, 'lifespan', self._points['death'] - self._points['birth'], dtypes='f8')
        self._points = rf.append_fields(self._points, 'avg_coord', (self._points['death'] + self._points['birth'])/2, dtypes='f8')


    def __repr__(self):
        tup = type(self).__name__, self.path
        return "%s('path=%s')" % tup

    ######################################################################
    ## accessor methods 
    ######################################################################

    def points(self):
        return self._points

    ######################################################################
    ## functional methods
    ######################################################################



    ######################################################################
    ## utility methods
    ######################################################################


