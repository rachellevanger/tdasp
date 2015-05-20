


import tdasp
from nldr import diffusionMap as dm
from lib import persistenceDiagram as pd
from scipy.stats import kendalltau
import numpy as np

import sys
from optparse import OptionParser


parser = OptionParser('usage: -o //path/options.txt')

parser.add_option("-o", dest="optionsfile",
                  help="path to the file contatining all of the variable settings")

(options, args) = parser.parse_args()


options_file = open(options.optionsfile).read()
options_file_vars = options_file.split("\n")

for var in options_file_vars:
    var_pair = var.split("=")
    var_name = var_pair[0]
    var_value = var_pair[1]

    if var_name == 'project_name':
        opt_project_name = var_value
    elif var_name == 'project_desc':
        opt_project_desc = var_value
    elif var_name == 'project_dirOut':
        opt_project_dirOut = var_value
    elif var_name == 'project_distmat':
        opt_project_distmat = var_value
    elif var_name == 'pd_dir':
        opt_pd_dir = var_value
    elif var_name == 'order_of_files':
        opt_order_of_files = var_value
    elif var_name == 'enstrophy_values':
        opt_enstrophy_values = var_value
    elif var_name == 'eps_low':
        opt_eps_low = float(var_value)
    elif var_name == 'eps_high':
        opt_eps_high = float(var_value)
    elif var_name == 'eps_num_steps':
        opt_eps_num_steps = int(var_value)
    elif var_name == 'num_eigenvecs':
        opt_num_eigenvecs = int(var_value)
    elif var_name == 'eps_kendall_tau':
        opt_eps_kendall_tau = var_value
    elif var_name == 'hom_max':
        opt_hom_max = int(var_value)
    elif var_name == 'param_max':
        opt_param_max = int(var_value)
    elif var_name == 'param_skip':
        opt_param_skip = int(var_value)



# Number of eigenvectors to compare
numCoords = opt_num_eigenvecs

# Initialize the project.
myProject = tdasp.Project(opt_project_name, opt_project_desc, opt_project_dirOut)

# Import the associated distance matrix into the project.
myProject.loadDistanceMatrix(opt_project_distmat, '', ' ')

# Load field-level data into the project.
myProject.loadDataAttributesFromFile(opt_order_of_files,'fileName','a20')
myProject.loadDataAttributesFromFile(opt_enstrophy_values,'enstrophy','f8')



######### RUN SINGLE EPSILON ANALYSIS ###########

if opt_eps_kendall_tau == 'mean':
    eps = myProject.distanceMatrix().mean()
else:
    eps = float(opt_eps_kendall_tau)

# Generate Diffusion Map Coordinates.
myMap = dm.DiffusionMap(eps)
EVecs, EVals = myMap.generateEmbeddedCoordinates(myProject.distanceMatrix(),numCoords)

for j in range(0,EVals.shape[0]):
    myProject.loadDataAttributesFromArray(EVecs[:,j], 'diffMap_mean_' + str(j+1), 'f8')

# Output diffusion map coordinates to file.
f = open(opt_project_dirOut+"diffusion_projections_eps_{0}.txt".format(eps), 'w')
for j in range(0,myProject.numDataPoints()):
    f.write("{0}, {1}, {2}\n".format(myProject.dataPoints()['fileName'][j], myProject.dataPoints()['diffMap_mean_1'][j], myProject.dataPoints()['diffMap_mean_2'][j]))

f.close()


