


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

print('Loading persistence diagrams into project...')

# Set persistence data directory for project
myProject.dirPD = opt_pd_dir
myProject.fileFormatPD = '%s_persistence_%d'

# Set homology range
myProject.hom_max = opt_hom_max

# Load the Persistence Diagrams
myProject.loadPersistenceDiagrams()

print('Capturing persistence diagram statistics for correlations...')

# (Avg) birth, death, lifespan, avg_coord for top n points, sorted by those values, descending
fields = ['birth', 'death', 'lifespan', 'avg_coord']

# Set parameters
stat_param_max = opt_param_max
stat_param_skip = opt_param_skip
print("Number of steps: " + str(stat_param_max/stat_param_skip+1))

all_stats = []
for hom in range(0,myProject.hom_max+1):
    print("Hom dim: " + str(hom))
    hom_stats = []
    for f in range(0,np.size(fields)):
        print("Field: " + fields[f])
        field_stats = []
        for j in range(0,myProject.numDataPoints()):
            point_stats = []
            for n in range(0,stat_param_max/stat_param_skip+1):
                data = myProject.persistenceDiagrams()[hom][j].points()[fields[f]]
                data.sort()
                data = data[::-1]
                m = min(n*stat_param_skip, np.size(data))
                if m==0:
                    m=1
                data = data[0:m]
                stat = [m, j, data.mean()] # n-value, point index, average
                point_stats.append(stat)
            field_stats.append(point_stats)
        hom_stats.append(field_stats)
    all_stats.append(hom_stats)

avg_stats = np.array(all_stats)
    
print('...statistics captured!')

# for data output
x = avg_stats[0,0,0,:,0]
field_labels = ['Birth', 'Death', 'Lifespan', '(Birth + Death)/2']

# array to hold correlations
corr_stats = np.zeros((np.shape(fields)[0], stat_param_max/stat_param_skip+1, myProject.hom_max+1, 2))

# Loop through each statistic collected
for field in fields:
    # gather correlations at different parameter values
    for hom in range(0, myProject.hom_max+1):
        for n in range(0,stat_param_max/stat_param_skip+1):
            corr_stats[fields.index(field),n,hom,0], corr_stats[fields.index(field),n,hom,1] = kendalltau(myProject.dataPoints()['diffMap_mean_1'], avg_stats[hom,fields.index(field), :, n,2])
    
    # Output data
    # Correlations
    f = open(opt_project_dirOut+"KT_corr_avg_{0}.txt".format(field), 'w')
    y = corr_stats[fields.index(field), :,:,0]
    for i in range(0,len(x)):
        f.write("{0}, {1}\n".format(str(x[i])," ".join([str(abs(val)) for val in y[i]])))
    f.close()

    # P-Values
    f = open(opt_project_dirOut+"KT_pval_avg_{0}.txt".format(field), 'w')
    y = corr_stats[fields.index(field), :,:,1]
    for i in range(0,len(x)):
        f.write("{0}, {1}\n".format(str(x[i])," ".join([str(val) for val in y[i]])))
    f.close()




