


import tdasp
from nldr import diffusionMap as dm
from lib import persistenceDiagram as pd
from scipy.stats import kendalltau
import numpy as np

from pylab import *
import matplotlib.pyplot as plt

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

# plt.subplots(1,3,figsize=(20,5))

# subplot(1,3,1)
# x = diffusion_array['eps']
# y = diffusion_array['evals']
# semilogy(x, y)
# xlabel('epsilon')
# ylabel('Eigenvalues')
# title('Top 10 eigenvalues as function of epsilon')

# subplot(1,3,2)
# x = diffusion_array['eps']
# y = diffusion_array['ktau']
# plot(x, y)
# xlabel('epsilon')
# ylabel('Kendall-tau')
# title('Kendall-tau Correlation for top 10 eigenvectors')

# subplot(1,3,3)
# x = diffusion_array['eps']
# y = diffusion_array['pvals']
# semilogy(x, y)
# xlabel('epsilon')
# ylabel('p-value')
# title('p-values for Kendall-tau correlations')

# savefig(myProject.dirOut + 'diffusion_parameter_analysis.png')

print('...analysis done!')


# ######### RUN SINGLE EPSILON ANALYSIS ###########

# if opt_eps_kendall_tau == 'mean':
#     eps = myProject.distanceMatrix().mean()
# else:
#     eps = float(opt_eps_kendall_tau)

# # Generate Diffusion Map Coordinates.
# myMap = dm.DiffusionMap(eps)
# EVecs, EVals = myMap.generateEmbeddedCoordinates(myProject.distanceMatrix(),numCoords)

# for j in range(0,EVals.shape[0]):
#     myProject.loadDataAttributesFromArray(EVecs[:,j], 'diffMap_mean_' + str(j+1), 'f8')

# print('Loading persistence diagrams into project...')

# # Set persistence data directory for project
# myProject.dirPD = opt_pd_dir
# myProject.fileFormatPD = '%s_persistence_%d'


# # Set homology range
# myProject.hom_max = opt_hom_max

# # Load the Persistence Diagrams
# myProject.loadPersistenceDiagrams()

# print('Capturing persistence diagram statistics for correlations...')

# # (Avg) birth, death, lifespan, avg_coord for top n points, sorted by those values, descending
# fields = ['birth', 'death', 'lifespan', 'avg_coord']

# # Set parameters
# stat_param_max = opt_param_max
# stat_param_skip = opt_param_skip
# print("Number of steps: " + str(stat_param_max/stat_param_skip+1))

# all_stats = []
# for hom in range(0,myProject.hom_max+1):
#     print("Hom dim: " + str(hom))
#     hom_stats = []
#     for f in range(0,np.size(fields)):
#         print("Field: " + fields[f])
#         field_stats = []
#         for j in range(0,myProject.numDataPoints()):
#             point_stats = []
#             for n in range(0,stat_param_max/stat_param_skip+1):
#                 data = myProject.persistenceDiagrams()[hom][j].points()[fields[f]]
#                 data.sort()
#                 data = data[::-1]
#                 m = min(n*stat_param_skip, np.size(data))
#                 if m==0:
#                     m=1
#                 data = data[0:m]
#                 stat = [m, j, data.mean()] # n-value, point index, average
#                 point_stats.append(stat)
#             field_stats.append(point_stats)
#         hom_stats.append(field_stats)
#     all_stats.append(hom_stats)

# avg_stats = np.array(all_stats)
    
# print('...statistics captured!')

# # for charts
# x = avg_stats[0,0,0,:,0]
# field_labels = ['Birth', 'Death', 'Lifespan', '(Birth + Death)/2']

# # array to hold correlations
# corr_stats = np.zeros((np.shape(fields)[0], stat_param_max/stat_param_skip+1, myProject.hom_max+1, 2))

# # Loop through each statistic collected
# for field in fields:
#     # gather correlations at different parameter values
#     for hom in range(0, myProject.hom_max+1):
#         for n in range(0,stat_param_max/stat_param_skip+1):
#             corr_stats[fields.index(field),n,hom,0], corr_stats[fields.index(field),n,hom,1] = kendalltau(myProject.dataPoints()['diffMap_mean_1'], avg_stats[hom,fields.index(field), :, n,2])
    
#     # make chart
#     plt.subplots(1,2,figsize=(15,5))
    
#     subplot(1,2,1)
#     y = corr_stats[fields.index(field), :,:,0]
#     plot(x, y)
#     xlabel('n')
#     ylabel('Kendall-tau correlation coefficient')
#     title('Kendall-tau correlations: Average(' + field_labels[fields.index(field)] + ')')
    
#     subplot(1,2,2)
#     y = corr_stats[fields.index(field), :,:,1]
#     semilogy(x, y)
#     xlabel('n')
#     ylabel('Kendall-tau p-value')
#     title('Kendall-tau correlations: Average(' + field_labels[fields.index(field)] + ')')
    
#     savefig(myProject.dirOut + 'KT_corr_avg_' + field + '.png')
    
    



