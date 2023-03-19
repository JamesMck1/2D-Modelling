#Load Packages
import numpy as np
import Scheme

###################################################################################################################################################################################
# Inputs
##################################################################################################################################################################################

rows = 10
cols = 10
init_depth = np.full((rows,cols), 0.1)
init_depth[0] = 1

init_momentum_x = np.zeros((rows,cols))
init_momentum_y = np.zeros((rows,cols))
boundary_conditions = [np.repeat(1,rows)],[np.repeat(1,cols)],[np.repeat(1,rows)],[np.repeat(1,cols)]
t_end = 100
cell_width = 0.1

#Test = Scheme.scheme(rows, cols, init_depth, init_momentum_x, init_momentum_y, boundary_conditions, t_end, cell_width)

Test = Scheme.domain(rows, cols, init_depth, init_momentum_x, init_momentum_y, boundary_conditions, cell_width)
Test.generate_cells_and_interfaces()