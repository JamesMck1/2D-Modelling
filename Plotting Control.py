#Load Packages
import numpy as np
import Scheme
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from celluloid import Camera
from datetime import datetime
import Plotting

######################################################################################################################################################################
# Inputs
######################################################################################################################################################################

rows = 25
cols = 100
init_depth = np.full((rows,cols), 0.1)
#init_depth[0] = 1
#init_depth[:,:1] = 1
init_depth[40:61,40:61] = 1

init_momentum_x = np.zeros((rows,cols))
init_momentum_y = np.zeros((rows,cols))
boundary_conditions = [np.repeat(1,cols)],[np.repeat(0,rows)],[np.repeat(1,cols)],[np.repeat(0,rows)]
t_end = 10
cell_width = 0.01

######################################################################################################################################################################
# Output
######################################################################################################################################################################

plot = Plotting.Animate(rows, cols, init_depth, init_momentum_x, init_momentum_y, boundary_conditions, t_end, cell_width)