#Load Packages
import numpy as np
import Scheme
import matplotlib.pyplot as plt
from celluloid import Camera
from datetime import datetime

###################################################################################################################################################################################
# Inputs
##################################################################################################################################################################################

rows = 10
cols = 10
init_depth = np.full((rows,cols), 0.25)
#init_depth[0] = 1
#init_depth[:,:1] = 1
init_depth[4:5,4:5] = 1.5

init_momentum_x = np.zeros((rows,cols))
init_momentum_y = np.zeros((rows,cols))
boundary_conditions = [np.repeat(1,rows)],[np.repeat(1,cols)],[np.repeat(1,rows)],[np.repeat(1,cols)]
t_end = 2.5
cell_width = 0.1

#Test = Scheme.scheme(rows, cols, init_depth, init_momentum_x, init_momentum_y, boundary_conditions, t_end, cell_width)

dt_string = datetime.now().strftime("%d%m%Y %H%M%S")
    
Domain = Scheme.domain(rows, cols, init_depth, init_momentum_x, init_momentum_y, boundary_conditions, cell_width)
Domain.generate_cells_and_interfaces() 

fig, axs = plt.subplots(3, 1)
cmaps = ['Blues', 'RdBu']
camera = Camera(fig)
#fig.tight_layout()

#depth
nx, ny = (rows+1, cols+1)
x = np.linspace(0, 1, nx)
y = np.linspace(0, 1, ny)
xv, yv = np.meshgrid(x, y)
zmin, zmax = -1,1
#Depth
ax = axs[0]
pcm = ax.pcolormesh(x, y, Domain.init_depth, cmap=cmaps[0], vmin=zmin, vmax=zmax)
ax.set_title('Depth')
fig.colorbar(pcm, ax=ax)
#x-momentum
ax = axs[1]
pcm = ax.pcolormesh(x,y, init_momentum_x, cmap=cmaps[1], vmin=zmin, vmax=zmax)
ax.set_title('x-Momentum')
fig.colorbar(pcm, ax=ax)
#y-momentum
ax = axs[2]
pcm = ax.pcolormesh(x,y, init_momentum_y, cmap=cmaps[1], vmin=zmin, vmax=zmax)
ax.set_title('y-Momentum')
fig.colorbar(pcm, ax=ax) 

camera.snap()

while Domain.sim_time <= t_end:
    Domain.calculate_fluxes()
    Domain.update_cells()
    Domain.mass_conservation_check()

    #Depth
    ax = axs[0]
    pcm = ax.pcolormesh(x,y, Domain.depths, cmap=cmaps[0], vmin=zmin, vmax=zmax)
    #ax.set_title('Depth')
    fig.colorbar(pcm, ax=ax)
    #x-momentum
    ax = axs[1]
    pcm = ax.pcolormesh(x,y, Domain.x_momentums, cmap=cmaps[1], vmin=zmin, vmax=zmax)
    #ax.set_title('x-Momentum')
    fig.colorbar(pcm, ax=ax)
    #y-momentum
    ax = axs[2]
    pcm = ax.pcolormesh(x,y, Domain.y_momentums, cmap=cmaps[1], vmin=zmin, vmax=zmax)
    #ax.set_title('y-Momentum')
    fig.colorbar(pcm, ax=ax) 
    
    camera.snap()

Animation = camera.animate() #Animate snapshots by combining into a .gif
Animation.save(fr"C:\Users\User\Documents\PhD\GitHub\2D Outputs\Animated Solution - {dt_string}.gif", writer='pillow')














