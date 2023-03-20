#Load Packages
import numpy as np
import Scheme
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from celluloid import Camera
from datetime import datetime

def Animate(rows, cols, init_depth, init_momentum_x, init_momentum_y, boundary_conditions, t_end, cell_width):

    dt_string = datetime.now().strftime("%d%m%Y %H%M%S")
        
    Domain = Scheme.domain(rows, cols, init_depth, init_momentum_x, init_momentum_y, boundary_conditions, cell_width)
    Domain.generate_cells_and_interfaces() 
    
    fig, axs = plt.subplots(3, 1)
    cmaps = ['Blues', 'RdBu']
    camera = Camera(fig)
    #fig.tight_layout()
    
    nx, ny = (rows+1, cols+1)
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    xv, yv = np.meshgrid(x, y)
    dmin, dmax = 0, np.max(init_depth)
    xmin, xmax = np.min(init_momentum_x), np.max(init_momentum_x)
    ymin, ymax = np.min(init_momentum_y), np.max(init_momentum_y)
    
    #Depth
    ax = axs[0]
    pcm = ax.pcolormesh(y, x, Domain.init_depth, cmap=cmaps[0], vmin = dmin, vmax = dmax)
    ax.set_title('Depth')
    cbar_1 = fig.colorbar(pcm, ax=ax)
    #x-momentum
    ax = axs[1]
    pcm = ax.pcolormesh(y, x, init_momentum_x, cmap=cmaps[1], vmin = xmin, vmax = xmax)
    ax.set_title('x-Momentum')
    cbar_2 = fig.colorbar(pcm, ax=ax)
    #y-momentum
    ax = axs[2]
    pcm = ax.pcolormesh(y, x, init_momentum_y, cmap=cmaps[1], vmin = ymin, vmax = ymax)
    ax.set_title('y-Momentum')
    cbar_3 = fig.colorbar(pcm, ax=ax)
    #cbar_2 = fig.colorbar(pcm, ax=axs[1:3]) 
    
    camera.snap()
    
    while Domain.sim_time <= t_end:
        Domain.calculate_fluxes()
        Domain.update_cells()
        #Domain.mass_conservation_check()
        
        #Depth
        cbar_1.remove()
        cbar_2.remove()
        cbar_3.remove()
        ax = axs[0]
        pcm = ax.pcolormesh(y, x, Domain.depths, cmap=cmaps[0], vmin = dmin, vmax = dmax)
        #ax.set_title('Depth')
        cbar_1 = fig.colorbar(pcm, ax=ax)
        #x-momentum
        ax = axs[1]
        pcm = ax.pcolormesh(y, x, Domain.x_momentums, cmap=cmaps[1], vmin = xmin, vmax = xmax)
        #ax.set_title('x-Momentum')
        cbar_2 = fig.colorbar(pcm, ax=ax)
        #y-momentum
        ax = axs[2]
        pcm = ax.pcolormesh(y, x, Domain.y_momentums, cmap=cmaps[1], vmin = ymin, vmax = ymax)
        #ax.set_title('y-Momentum')
        cbar_3 = fig.colorbar(pcm, ax=ax) 
        
        camera.snap()
    
    Animation = camera.animate() #Animate snapshots by combining into a .gif
    Animation.save(fr"C:\Users\User\Documents\PhD\GitHub\2D Outputs\Animated Solution - {dt_string}.gif", writer='pillow')










