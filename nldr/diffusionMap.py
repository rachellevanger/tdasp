# TDA Statistical Package : NLDR : diffusionMap
#
# Copyright (C) 2014 Rachel Levanger
# Author: Rachel Levanger <rachel.levanger@gmail.com>
#                         <rachel@math.rutgers.edu>
# URL: <http://www.rachellevanger.com/math/TDASP>

import numpy as np
#from scipy.sparse.linalg.eigen.arpack import eigsh as largest_eigsh
from scipy.sparse.linalg.eigen.arpack import eigen_symmetric as largest_eigsh
import warnings as w

######################################################################
## Table of Contents
######################################################################
## - Constants
## - Data Classes
##   - DiffusionMap
## - 

######################################################################
## Constants
######################################################################



######################################################################
## Data Classes
######################################################################


class DiffusionMapError(Exception):
    """An exception class for diffusion map-related errors."""

class DiffusionMap(object):
    """
    The parameter-level data container for a diffusion map embedding.

    Create a Diffusion Map by supplying an epsilon parameter.
    
    """

    ######################################################################
    ## top-level methods
    ######################################################################
    def __init__(self, epsilon="1"):
        self.epsilon = epsilon

    def __repr__(self):
        tup = type(self).__name__, self.epsilon
        return "%s('parameter=%f')" % tup

    ######################################################################
    ## accessor methods 
    ######################################################################



    ######################################################################
    ## functional methods
    ######################################################################

    def generateEmbeddedCoordinates(self, distanceMatrix, numCoords=0):
        """ Computes the diffusion coordinates from the distance matrix supplied.
        The diffusion kernel is $k(x,y) = exp{ distance(x,y)^2 / epsilon }
        SVD is used to compuet the singular values of a matrix that is conjugated
        to the kernel. Then the corrsponding vectors are rescaled to get the
        eigen vectors. (Original version: Miro Kramar, 2013, in MATLAB)
        
        NOTES: - The first eigenvalue is always +1, so it is not returned.

        Rachel Levanger, 2014."""
        if numCoords==0:
            numCoords=distanceMatrix.shape[0] - 2

        min_coords = min(numCoords, distanceMatrix.shape[1] - 2)
        if min_coords < numCoords:
            w.warn('Number of coordinates truncated to max value of [ndim(distmat)-2].')
            numCoords = min_coords

        numCoords = numCoords + 1 # Make up for zero-based indexing

        matrix_size = distanceMatrix.shape
        K = np.zeros(matrix_size)

        # Construct the kernel
        for i in range(0,matrix_size[0]):
            for j in range(0,matrix_size[1]):
                K[i,j] = np.exp(-((distanceMatrix[i,j])/self.epsilon))
                
        # Construct the symmetric conjugate (by \sqrt(pi) ) of the kernel 
        dx = sum(K)
        dx_square_root = np.sqrt(dx)
        
        P = np.zeros(matrix_size)
        for i in range(0,matrix_size[0]):
            for j in range(0,matrix_size[1]):
                P[i,j] = np.divide(K[i,j], dx_square_root[i] * dx_square_root[j] )
                
        # Get the largest eigenvalues and corresponding eigenvectors
        EigenValues, EigenVectors = largest_eigsh(P, numCoords, which="LM")

        # Transform the eigen vectors to get the (right) eigenvectors of
        # k(x,y)/dx(x) 

        d = sum(dx)
        pi_sqrt = np.sqrt(np.divide(dx, d))

        for i in range(0,numCoords):
            for j in range(0,EigenVectors.shape[0]):
                EigenVectors[j,i] = np.divide(EigenVectors[j,i], pi_sqrt[j])

        # Sort by eigenvalue descending
        I = EigenValues.argsort()[::-1]
        EigenValues = EigenValues[I]
        EigenVectors = EigenVectors[:,I]

        # Account for sign changes in eigenvectors and scale by eigenvalues
        for i in range(0,EigenValues.shape[0]):
            EigenVectors[:,i] = EigenVectors[:,i]*EigenValues[i]
            if EigenVectors[0,i] < 0:
                EigenVectors[:,i] = EigenVectors[:,i]*(-1)

        return EigenVectors[:,1:], EigenValues[1:]

    ######################################################################
    ## utility methods
    ######################################################################


