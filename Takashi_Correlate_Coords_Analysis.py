


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
    elif var_name == 'coords_epsilon':
        opt_coords_epsilon = float(var_value)
    elif var_name == 'coords_pts_min':
        opt_coords_pts_min = int(var_value)
    elif var_name == 'coords_pts_max':
        opt_coords_pts_max = int(var_value)
    elif var_name == 'coords_p_min':
        opt_coords_p_min = int(var_value)
    elif var_name == 'coords_p_max':
        opt_coords_p_max = int(var_value)



# Initialize the project.
myProject = tdasp.Project(opt_project_name, opt_project_desc, opt_project_dirOut)

# Import the associated distance matrix into the project.
myProject.loadDistanceMatrix(opt_project_distmat, '', ' ')

# Load field-level data into the project.
myProject.loadDataAttributesFromFile(opt_order_of_files,'fileName','a20')
myProject.loadDataAttributesFromFile(opt_enstrophy_values,'enstrophy','f8')



print('Loading persistence diagrams into project...')

# Set persistence data directory for project
myProject.dirPD = opt_pd_dir
myProject.fileFormatPD = '%s_persistence_%d'

# Set homology range
myProject.hom_max = opt_hom_max

# Load the Persistence Diagrams
myProject.loadPersistenceDiagrams()


print('Generating coordinates for persistence diagrams...')

# Generate the coordinates
myProject.generateCoordinates(opt_coords_epsilon, opt_coords_pts_min, opt_coords_pts_max, opt_coords_p_min, opt_coords_p_max)

print('Collecting coordinates for persistence diagrams...')

# Set parameters
# (Avg) birth, death, lifespan, avg_coord for top n points, sorted by those values, descending
coord_names = ['lifespan', 'birth']

all_stats = []
for hom in range(0,myProject.hom_max+1):
    print("Hom dim: " + str(hom))
    hom_stats = []
    for c in range(0,np.size(coord_names)):
        print("Coordinate: " + coord_names[c])
        coord_stats = []
        for j in range(0,myProject.numDataPoints()):
            if coord_names[c] == 'lifespan':
                coord_stats.append(myProject.persistenceDiagrams()[hom][j].lifespan_coords)
            elif coord_names[c] == 'birth':
                coord_stats.append(myProject.persistenceDiagrams()[hom][j].birth_coords)
        hom_stats.append(coord_stats)
    all_stats.append(hom_stats)

coord_stats = np.array(all_stats)
    
print('...coordinates captured!')

num_pts = opt_coords_pts_max - opt_coords_pts_min + 1
num_ps = opt_coords_p_max - opt_coords_p_min + 1

# for data output
x = range(opt_coords_pts_min, opt_coords_pts_max + 1) # n's are on x-axis.
chart_labels = ['Lifespan Coordinates', 'Birth Coordinates']

# array to hold correlations
corr_stats = np.zeros((np.shape(coord_names)[0], num_pts, num_ps, myProject.hom_max+1, 2))

print('Correlating coordinates...')

# Loop through each coordinate
for c in coord_names:
    # gather correlations at different parameter values
    for hom in range(0, myProject.hom_max+1):
        for n in range(0, num_pts):
            for p in range(0, num_ps):
                corr_stats[coord_names.index(c),n,p,hom,0], corr_stats[coord_names.index(c),n,p,hom,1] = kendalltau(myProject.dataPoints()['enstrophy'], coord_stats[hom,coord_names.index(c), :, p, n])
                corr_stats[coord_names.index(c),n,p,hom,0] = np.absolute(corr_stats[coord_names.index(c),n,p,hom,0])

        # Output data
        # Correlations
        f = open(opt_project_dirOut+"KT_corr_coords_{0}_hom_{1}.txt".format(c,str(hom)), 'w')
        print(opt_project_dirOut+"KT_corr_coords_{0}_hom_{1}.txt".format(c,str(hom)))
        y = corr_stats[coord_names.index(c), :,:,hom,0]
        for i in range(0,len(x)):
            f.write("{0}, {1}\n".format(str(x[i])," ".join([str(abs(val)) for val in y[i]])))
        f.close()

        # P-Values
        f = open(opt_project_dirOut+"KT_pval_coords_{0}_hom_{1}.txt".format(c,str(hom)), 'w')
        y = corr_stats[coord_names.index(c), :,:,hom,1]
        for i in range(0,len(x)):
            f.write("{0}, {1}\n".format(str(x[i])," ".join([str(val) for val in y[i]])))
        f.close()
    
print('...correlations captured!')



        # # make chart
        # plt.subplots(1,2,figsize=(15,5))
        
        # subplot(1,2,1)
        # y = corr_stats[coord_names.index(c), :,:,hom,0]
        # plot(x, y)
        # xlabel('n')
        # ylabel('Kendall-tau correlation coefficient')
        # title('Kendall-tau correlations: Average(' + chart_labels[coord_names.index(c)] + ') hom=' + str(hom))
        
        # subplot(1,2,2)
        # y = corr_stats[coord_names.index(c), :,:,hom,1]
        # semilogy(x, y)
        # xlabel('n')
        # ylabel('Kendall-tau p-value')
        # title('Kendall-tau correlations: Average(' + chart_labels[coord_names.index(c)] + ') hom=' + str(hom))
        
        # savefig(myProject.dirOut + 'KT_corr_avg_' + c + '_hom_'+str(hom)+'.png')


