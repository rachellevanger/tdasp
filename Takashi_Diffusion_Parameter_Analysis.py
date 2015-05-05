


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


# How to sweep through values of epsilon
eps_low = opt_eps_low
eps_high = opt_eps_high
eps_num_steps = opt_eps_num_steps

# Number of eigenvectors to compare
numCoords = opt_num_eigenvecs

# Initialize the project.
myProject = tdasp.Project(opt_project_name, opt_project_desc, opt_project_dirOut)

# Import the associated distance matrix into the project.
myProject.loadDistanceMatrix(opt_project_distmat, '', ' ')


# Load field-level data into the project.
myProject.loadDataAttributesFromFile(opt_order_of_files,'fileName','a20')
myProject.loadDataAttributesFromFile(opt_enstrophy_values,'enstrophy','f8')

# Set single epsilon step
eps_step = (eps_high - eps_low)/float(eps_num_steps)

# Initialize structured array to hold the data, data built as follows:
# eps = (scalar) epsilon parameter for diffusion map
# evals = (float vector length numcoords) vector of eigenvalues for associated epsilon parameter
# ktau = (float vector length numcoords) kendall tau correlation coefficients for each eigenvalue
# pvals = (float vector length numcoords) p-values for correlation coefficients for eahch eigenvalue
diff_descriptor = {'names': ('eps','evals','ktau', 'pvals'), 'formats': ('f4', str(numCoords)+'f8', str(numCoords)+'f8', str(numCoords)+'f8', str(numCoords)+'f8')}
diffusion_array = np.zeros(eps_num_steps, dtype=diff_descriptor)


####### RUN DIFFUSION MAP ANALYSIS ######

print('Running analysis for diffusion map parameter...')

# Fill the structured array
for k in range(0,eps_num_steps):

    # Set epsilon
    eps = eps_low + k*eps_step
    diffusion_array['eps'][k] = eps
    
    # Generate Diffusion Map Coordinates.
    myMap = dm.DiffusionMap(eps)
    EVecs, EVals = myMap.generateEmbeddedCoordinates(myProject.distanceMatrix(),numCoords)
    diffusion_array['evals'][k] = EVals
    for j in range(0,EVals.shape[0]):
        myProject.loadDataAttributesFromArray(EVecs[:,j], 'diffMap_' + str(j+1), 'f8')
        diffusion_array['ktau'][k][j], diffusion_array['pvals'][k][j] = kendalltau(myProject.dataPoints()['enstrophy'], myProject.dataPoints()['diffMap_' + str(j+1)])
        myProject.removeDataAttributes('diffMap_' + str(j+1))
    
diffusion_array['ktau'] = np.absolute(diffusion_array['ktau'])


# Diffusion Parameter Analysis Output
# Eigenvalues
f = open(opt_project_dirOut+'diffusion_parameter_analysis_eigenvalues.txt', 'w')
x = diffusion_array['eps']
y = diffusion_array['evals']
for i in range(0,len(x)):
    f.write("{0}, {1}\n".format(str(x[i])," ".join([str(val) for val in y[i]])))
f.close()

# Kendall-Tau Correlation
f = open(opt_project_dirOut+'diffusion_parameter_analysis_kendall_tau.txt', 'w')
x = diffusion_array['eps']
y = diffusion_array['ktau']
for i in range(0,len(x)):
    f.write("{0}, {1}\n".format(str(x[i])," ".join([str(val) for val in y[i]])))
f.close()

# Kendall-Tau p-values
f = open(opt_project_dirOut+'diffusion_parameter_analysis_pvals.txt', 'w')
x = diffusion_array['eps']
y = diffusion_array['pvals']
for i in range(0,len(x)):
    f.write("{0}, {1}\n".format(str(x[i])," ".join([str(val) for val in y[i]])))
f.close()

print('...analysis done!')






