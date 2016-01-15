"""
Created on Wed Oct 28 10:05:37 2015

@author: pieter
"""
import time 
import tables
import scipy
import sys
import numpy as np
import matplotlib.pylab as plt
from numpy import *

from scipy.linalg import norm
from scipy import sqrt

import matplotlib.pylab as plt

from pylab import *
import os
# set the interactive mode of pylab to ON
ion()
# opens a new figure to plot into
fig_hndl = figure()
# make an empty list into which we'll append
# the filenames of the PNGs that compose each
# frame.
files=[]    
# filename for the name of the resulting movie
filename = 'animation'

os.system("mencoder 'mf://movieplots/*.jpg' -mf type=jpg:fps=30 -ovc lavc  -lavcopts vcodec=wmv2 -oac copy -o output.avi")
# cleanup

#mencoder 'mf://*.png' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=wmv2 -oac copy  -o output.avi
#'movieplots/plot_rho_t*.png'



#mencoder "mf://movieplots/*.png" -mf type=png:fps=10 -ovc lavc  -lavcopts vcodec=wmv2 -oac copy -o output.avi