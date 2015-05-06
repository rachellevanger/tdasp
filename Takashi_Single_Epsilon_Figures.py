


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



# Initialize the project.
myProject = tdasp.Project(opt_project_name, opt_project_desc, opt_project_dirOut)

# (Avg) birth, death, lifespan, avg_coord for top n points, sorted by those values, descending
fields = ['birth', 'death', 'lifespan', 'avg_coord']
field_labels = ['Birth', 'Death', 'Lifespan', '(Birth + Death)/2']


# Loop through each statistic collected
for field in fields:


    # make charts
    plt.subplots(1,2,figsize=(15,5))
    
    subplot(1,2,1)
    f = open(opt_project_dirOut+"KT_corr_avg_{0}.txt".format(field)).read()
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
    title('Kendall-tau correlations: Average(' + field_labels[fields.index(field)] + ')')
    
    subplot(1,2,2)
    f = open(opt_project_dirOut+"KT_pval_avg_{0}.txt".format(field)).read()
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
    title('Kendall-tau correlations: Average(' + field_labels[fields.index(field)] + ')')
    
    savefig(myProject.dirOut + 'KT_corr_avg_' + field + '.png')
    
    



