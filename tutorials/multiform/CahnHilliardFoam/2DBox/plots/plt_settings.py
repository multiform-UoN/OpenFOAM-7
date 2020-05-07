import itertools
import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import sys

#sns.set(style='white')
#sns.set()

rcParams.update({'figure.autolayout': True})
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True



# Dims given in [width, height]
fig_dims = [7.33, 5] # height in inches, single column figure
double_dims = [6.66, 2.058] # double width fig
si_dims = [6.66, 4.116] # double size figure

figsize = {'paper':fig_dims,
           'si':si_dims, 
           'double':double_dims}
prefix = {'paper':'',
          'si':'SI_',
          'double':'d_'}
fig_locs = ['paper', 'si']

cols = sns.color_palette('Paired')

beltv_all = ['1mm','3mm','5mm','10mm','25mm','50mm','100mm']
beltv_lim = ['1mm','3mm','5mm','10mm','25mm','50mm']

#colours = {'1mm':cols[1],
#           '3mm':cols[3],
#           '5mm':cols[5],
#           '10mm':cols[7],
#           '25mm':cols[9],
#           '50mm':cols[11],
#           '100mm':cols[0],   
#       }
         

