


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
# Set homology range
myProject.hom_max = opt_hom_max

# Import the associated distance matrix into the project.
myProject.loadDistanceMatrix(opt_project_distmat, '', ' ')

# Set parameters
# (Avg) birth, death, lifespan, avg_coord for top n points, sorted by those values, descending
coord_names = ['lifespan', 'birth']

# for data output
x = range(opt_coords_pts_min, opt_coords_pts_max + 1) # n's are on x-axis.
chart_labels = ['Lifespan Coordinates', 'Birth Coordinates']

# Loop through each coordinate
for c in coord_names:
    # gather correlations at different parameter values
    for hom in range(0, myProject.hom_max+1):

        # make charts
        plt.subplots(1,2,figsize=(15,5))
        
        subplot(1,2,1)
        f = open(opt_project_dirOut+"KT_corr_coords_{0}_hom_{1}.txt".format(c,str(hom))).read()
        f_lines = f.split('\n')
        x = []
        y = []
        for line in f_lines:
            if line == "":
                break
            line_pair = line.split(", ")
            x.append(float(line_pair[0]))
            y.append(line_pair[1].split(" "))
        plot(x, y)
        xlabel('n')
        ylabel('Kendall-tau correlation coefficient')
        title('Kendall-tau correlations: ' + chart_labels[coord_names.index(c)])
        
        subplot(1,2,2)
        f = open(opt_project_dirOut+"KT_pval_coords_{0}_hom_{1}.txt".format(c,str(hom))).read()
        f_lines = f.split('\n')
        x = []
        y = []
        for line in f_lines:
            if line == "":
                break
            line_pair = line.split(", ")
            x.append(float(line_pair[0]))
            y.append(line_pair[1].split(" "))
        semilogy(x, y)
        xlabel('n')
        ylabel('Kendall-tau p-value')
        title('Kendall-tau p-values: ' + chart_labels[coord_names.index(c)])
        
        savefig(myProject.dirOut + 'KT_corr_' + c + '_hom_'+str(hom)+'.png')
        
        

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


