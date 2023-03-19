# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 16:42:54 2023

@author: James Mckenna
"""

#Load Packages
import numpy as np
import math
import os.path
from datetime import date
from datetime import datetime
import HLL

#Cell Objects

class cell():
    
    #Class variables
    g = 9.81 #gravitational acceleration (m/s^2)
    
    def __init__(self, depth, x_momentum, y_momentum, cell_width, mannings_n):
        self.depth = depth
        self.x_momentum = x_momentum
        self.y_momentum = y_momentum
        self.x_velocity = np.nan_to_num(x_momentum/depth)
        self.y_velocity = np.nan_to_num(y_momentum/depth)
        self.cell_width = cell_width
        self.mannings_n = mannings_n
    
class interior_cell(cell):
    
    def __init__(self, depth, x_momentum, y_momentum, cell_width):
        super().__init__(depth, x_momentum, y_momentum, cell_width)
        
    def assign_interfaces(self, x_interface_L, x_interface_R, y_interface_L, y_interface_R): #assign interface objects
        self.x_interface_L = x_interface_L #-ve x = left
        self.x_interface_R = x_interface_R #+ve x = right
        self.y_interface_L = y_interface_L #-ve y = left
        self.y_interface_R = y_interface_R #+ve y = right
        
    def update(self, stable_tstep):
        """
        Use a strang splitting scheme to update the conserved variables in two dimensions
        
        0) Using stable timestep delta t = C*(delta x/S max), S max = max(|max(u,v)|+c)
        1) Update in x-direction using initial condition U^n to give U^n+1/2
        2) Update in y-direction using initial condition U^n+1/2 to give U^n+1
        
        """
        
        #x-sweep
        self.depth = self.depth - (stable_tstep/self.cell_width)*(self.x_interface_L.depth_flux-self.x_interface_R.depth_flux) #updated depth = depth - dt/dx(F_L-F_R)
        self.x_momentum = self.x_momentum - (stable_tstep/self.cell_width)*(self.x_interface_L.momentum_flux-self.x_interface_R.momentum_flux) #updated depth = depth - dt/dx(F_L-F_R)
        
        #y-sweep
        self.depth = self.depth - (stable_tstep/self.cell_width)*(self.y_interface_L.depth_flux-self.y_interface_R.depth_flux) #updated depth = depth - dt/dx(F_L-F_R)
        self.y_momentum = self.y_momentum - (stable_tstep/self.cell_width)*(self.y_interface_L.momentum_flux-self.y_interface_R.momentum_flux) #updated depth = depth - dt/dx(F_L-F_R)
        
class exterior_cell(cell): #ghost cell
    
    def __init__(self, depth, x_momentum, y_momentum, cell_width, mannings_n):
        super().__init__(cell_width, mannings_n)
        
class reflective_cell(exterior_cell): #reflective ghost cell
    
    def __init__(self, cell_width, mannings_n, int_cell):
        super().__init__(cell_width, mannings_n)
        self.depth = int_cell.depth
        self.x_momentum = -int_cell.x_momentum
        self.y_momentum = -int_cell.y_momentum
        self.x_velocity = -int_cell.x_velocity
        self.y_velocity = -int_cell.y_velocity
        
class transmissive_cell(exterior_cell): #transmissive ghost cell
    
    def __init__(self, cell_width, mannings_n, int_cell):
        super().__init__(cell_width, mannings_n)
        self.depth = int_cell.depth
        self.x_momentum = int_cell.x_momentum
        self.y_momentum = int_cell.y_momentum
        self.x_velocity = int_cell.x_velocity
        self.y_velocity = int_cell.y_velocity
        
#Interface Objects

class interface():
    
    #Class variables
    g = 9.81 #gravitational acceleration (m/s^2)
    
    def __init__(self, left_cell, right_cell):
        self.left_cell = left_cell #cell object assigned as left state
        self.right_cell = right_cell #cell object assigned as right state
        
class x_interface(interface): #interface in the x-plane
    
    def __init__(self, left_cell, right_cell):
        super().__init__(left_cell, right_cell) #+ve = right, -ve = left
        
    def wavespeeds(self): #calculate wavespeeds
        
        wavespeeds = HLL.Wavespeeds(self.g, self.left_cell.x_velocity, self.right_cell.x_velocity, self.left_cell.depth, self.right_cell.depth)
    
        self.left_wavespeed = wavespeeds[0]
        self.right_wavespeed = wavespeeds[1]
    
    def HLL_Solver(self): #use a HLL approximate Riemann solver to resolve numerical fluxes
        
        fluxes = HLL.HLL_Solver(self.g, self.left_cell.x_velocity, self.right_cell.x_velocity, self.left_cell.depth, self.right_cell.depth,
                                self.left_wavespeed, self.right_wavespeed)
        
        self.depth_flux = fluxes[0]
        self.momentum_flux = fluxes[1]
        
class y_interface(interface): #interface in the y-plane
    
    def __init__(self, left_cell, right_cell):
        super().__init__(left_cell, right_cell) #+ve = right, -ve = left
        
    def wavespeeds(self): #calculate wavespeeds
        
        wavespeeds = HLL.Wavespeeds(self.g, self.left_cell.y_velocity, self.right_cell.y_velocity, self.left_cell.depth, self.right_cell.depth)
    
        self.left_wavespeed = wavespeeds[0]
        self.right_wavespeed = wavespeeds[1]
    
    def HLL_Solver(self): #use a HLL approximate Riemann solver to resolve numerical fluxes
        
        fluxes = HLL.HLL_Solver(self.g, self.left_cell.y_velocity, self.right_cell.y_velocity, self.left_cell.depth, self.right_cell.depth,
                                self.left_wavespeed, self.right_wavespeed)
        
        self.depth_flux = fluxes[0]
        self.momentum_flux = fluxes[1]
        
#Domain Object

class Domain():
    
    def __init__(self, rows, cols, init_depth, init_momentum_x, init_momentum_y):
        self.cells = [] #initialise list to store cell objects
        self.interfaces = [] #initialise list to store interface objects
        self.grid = (rows, cols) #number of cells in each row and column
        
        self.init_depth = init_depth
        self.init_momentum_x = init_momentum_x
        self.init_momentum_y = init_momentum_y
        
    def generate_cells_and_interfaces(self):
        
        #Generate internal cells
        """
        Generate cells and store in a list
        
        """
        for iy in np.arange(0, self.grid[1]):
                row = [] #initialise an array to store cells in a row
                for ix in np.arange(0, self.grid[0]):
                    Cell = cell()
                
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        