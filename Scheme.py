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
    
class interior_cell(cell):
    
    def __init__(self, depth, x_momentum, y_momentum, cell_width):
        self.cell_width = cell_width
        self.depth = depth
        self.x_momentum = x_momentum
        self.y_momentum = y_momentum
        self.x_velocity = np.nan_to_num(x_momentum/depth)
        self.y_velocity = np.nan_to_num(y_momentum/depth)
        
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
        
        #x-sweep (update variables in the x-direction across perpendicular interfaces aligned with the y-axis)
        self.depth = self.depth + (stable_tstep/self.cell_width)*(self.y_interface_L.depth_flux-self.y_interface_R.depth_flux) #updated depth = depth - dt/dx(F_L-F_R)
        self.x_momentum = self.x_momentum + (stable_tstep/self.cell_width)*(self.y_interface_L.momentum_flux-self.y_interface_R.momentum_flux) #updated depth = depth - dt/dx(F_L-F_R)
        self.x_velocity = self.x_momentum/self.depth
        
        #y-sweep (update variables in the y-direction across perpendicular interfaces aligned with the x-axis)
        self.depth = self.depth + (stable_tstep/self.cell_width)*(self.x_interface_L.depth_flux-self.x_interface_R.depth_flux) #updated depth = depth - dt/dx(F_L-F_R)
        self.y_momentum = self.y_momentum + (stable_tstep/self.cell_width)*(self.x_interface_L.momentum_flux-self.x_interface_R.momentum_flux) #updated depth = depth - dt/dx(F_L-F_R)
        self.y_velocity = self.y_momentum/self.depth
        
class reflective_cell(cell): #reflective ghost cell
    
    def __init__(self, Cell):
        self.Cell = Cell
        self.depth = self.Cell.depth
        self.x_momentum = -self.Cell.x_momentum
        self.y_momentum = -self.Cell.y_momentum
        self.x_velocity = -self.Cell.x_velocity
        self.y_velocity = -self.Cell.y_velocity
    
    def update(self):
        self.depth = self.Cell.depth
        self.x_momentum = -self.Cell.x_momentum
        self.y_momentum = -self.Cell.y_momentum
        self.x_velocity = -self.Cell.x_velocity
        self.y_velocity = -self.Cell.y_velocity
        
class transmissive_cell(cell): #transmissive ghost cell
    
    def __init__(self, Cell):
        self.Cell = Cell
        self.depth = self.Cell.depth
        self.x_momentum = self.Cell.x_momentum
        self.y_momentum = self.Cell.y_momentum
        self.x_velocity = self.Cell.x_velocity
        self.y_velocity = self.Cell.y_velocity
    
    def update(self):
        self.depth = self.Cell.depth
        self.x_momentum = self.Cell.x_momentum
        self.y_momentum = self.Cell.y_momentum
        self.x_velocity = self.Cell.x_velocity
        self.y_velocity = self.Cell.y_velocity
        
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
        
        wavespeeds = HLL.Wavespeeds(self.g, self.left_cell.depth, self.right_cell.depth, self.left_cell.y_velocity, self.right_cell.y_velocity)
    
        self.left_wavespeed = wavespeeds[0]
        self.right_wavespeed = wavespeeds[1]
    
    def HLL_Solver(self): #use a HLL approximate Riemann solver to resolve numerical fluxes
        
        fluxes = HLL.HLL_Solver(self.g, self.left_cell.y_velocity, self.right_cell.y_velocity, self.left_cell.depth, self.right_cell.depth,
                                self.left_wavespeed, self.right_wavespeed)
        
        self.depth_flux = fluxes[0]
        self.momentum_flux = fluxes[1]
        
class y_interface(interface): #interface in the y-plane
    
    def __init__(self, left_cell, right_cell):
        super().__init__(left_cell, right_cell) #+ve = right, -ve = left
        
    def wavespeeds(self): #calculate wavespeeds
        
        wavespeeds = HLL.Wavespeeds(self.g, self.left_cell.depth, self.right_cell.depth, self.left_cell.x_velocity, self.right_cell.x_velocity)
    
        self.left_wavespeed = wavespeeds[0]
        self.right_wavespeed = wavespeeds[1]
    
    def HLL_Solver(self): #use a HLL approximate Riemann solver to resolve numerical fluxes
        
        fluxes = HLL.HLL_Solver(self.g, self.left_cell.x_velocity, self.right_cell.x_velocity, self.left_cell.depth, self.right_cell.depth,
                                self.left_wavespeed, self.right_wavespeed)
        
        self.depth_flux = fluxes[0]
        self.momentum_flux = fluxes[1]
        
#Domain Object

class domain():
    
    def __init__(self, rows, cols, init_depth, init_momentum_x, init_momentum_y, boundary_conditions, cell_width):
        self.cells = [] #initialise list to store cell objects
        self.x_interfaces = [] #initialise list to store interface objects aligned with the x-axis
        self.y_interfaces = [] #initialise list to store interface objects aligned with the y-axis
        self.grid = (rows, cols) #number of cells in each row and column
        self.cell_width = cell_width
        
        #initial condition list = list of array where each array represents the initial value for each row
        self.init_depth = init_depth
        self.init_momentum_x = init_momentum_x
        self.init_momentum_y = init_momentum_y
        
        #boundary conditions list = list of four arrays. Index i of the array represents the boundary condition for cell i. Boundary condition represented by a numerical value
        self.BCs_N = boundary_conditions[0] #northern boundary conditions
        self.BCs_E = boundary_conditions[1] #eastern boundary conditions
        self.BCs_S = boundary_conditions[2] #southern boundary conditions
        self.BCs_W = boundary_conditions[3] #western boundary conditions
        
        self.tstep = 0 #initialise stable timestep
        self.sim_time = 0 #initialise counter for the simulation duration
        
        #auditing
        self.x_depth_fluxes = np.zeros((rows+1, cols))
        self.y_depth_fluxes = np.zeros((rows, cols+1))
        self.x_momentum_fluxes = np.zeros((rows+1, cols))
        self.y_momentum_fluxes = np.zeros((rows, cols+1))
        
        self.depths = np.zeros((rows, cols)) 
        self.x_momentums = np.zeros((rows, cols))
        self.y_momentums = np.zeros((rows, cols))
        
        self.masses = [np.sum(self.init_depth)]
        
    def generate_cells_and_interfaces(self):
        
        #Generate internal cells
        """
        Generate cells and store in a list, cell_i,j = cell list index [i][j]
        
        """
        for ix in np.arange(0, self.grid[0]):
                row = [] #Initialise an array to store cells in a row
                for iy in np.arange(0, self.grid[1]):
                    Cell = interior_cell(self.init_depth[ix][iy], self.init_momentum_x[ix][iy], self.init_momentum_y[ix][iy], self.cell_width) #Generate cell object
                    row.append(Cell)
                    
                self.cells.append(row)
        
        #Generate external ghost cells
        """
        Generate external ghost cells and link to their internal neighbours
        
        External ghost cells generated in clockwise order (north boundary -> east boundary -> south boundary -> west boundary)
        
        """ 
        north = [] #northern ghost cells
        for iy in np.arange(0, self.grid[1]):
            if self.BCs_N[0][iy] == 0: #transmissive
                Ghost_cell = transmissive_cell(self.cells[0][iy]) #northern
                north.append(Ghost_cell)
                
            elif self.BCs_N[0][iy] == 1: #reflective
                Ghost_cell = reflective_cell(self.cells[0][iy]) #northern
                north.append(Ghost_cell)
            
            else: #invalid boundary condition
                print('Invalid Boundary Condition, Boundary condition should be formatted as an integer')
        
        south = [] #southern ghost cells
        for iy in np.arange(0, self.grid[1]):
            if self.BCs_S[0][iy] == 0: #transmissive
                Ghost_cell = transmissive_cell(self.cells[-1][iy]) #northern
                south.append(Ghost_cell)
                
            elif self.BCs_S[0][iy] == 1: #reflective
                Ghost_cell = reflective_cell(self.cells[-1][iy]) #northern
                south.append(Ghost_cell)
            
            else: #invalid boundary condition
                print('Invalid Boundary Condition, Boundary condition should be formatted as an integer')
        
        east = [] #eastern ghost cells
        for ix in np.arange(0, self.grid[0]):
            if self.BCs_E[0][ix] == 0: #transmissive
                Ghost_cell = transmissive_cell(self.cells[ix][-1]) #northern
                east.append(Ghost_cell)
                
            elif self.BCs_E[0][ix] == 1: #reflective
                Ghost_cell = reflective_cell(self.cells[ix][-1]) #northern
                east.append(Ghost_cell)
            
            else: #invalid boundary condition
                print('Invalid Boundary Condition, Boundary condition should be formatted as an integer') 
        
        west = [] #western ghost cells
        for ix in np.arange(0, self.grid[0]):
            if self.BCs_W[0][ix] == 0: #transmissive
                Ghost_cell = transmissive_cell(self.cells[ix][0]) #northern
                west.append(Ghost_cell)
                
            elif self.BCs_W[0][ix] == 1: #reflective
                Ghost_cell = reflective_cell(self.cells[ix][0]) #northern
                west.append(Ghost_cell)
            
            else: #invalid boundary condition
                print('Invalid Boundary Condition, Boundary condition should be formatted as an integer') 
        
        self.BCs = [north, east, south, west]
        
        #Generate internal interfaces
        """
        Generate interfaces by first performing a x-sweep (generate all interfaces aligned with the x-axis), then repeat for the y-axis
        
        """
        #x-sweep
        #first row
        row = []
        for iy in np.arange(0, self.grid[1]): #for each cell in the first row of interfaces (northern boundary)
            Interface = x_interface(self.BCs[0][iy], self.cells[0][iy]) #left cell = northern ghost cell, right cell = northern internal cell
            row.append(Interface)
        
        self.x_interfaces.append(row)
        
        #middle rows
        for ix in np.arange(1, self.grid[0]): #for each row of interfaces
            row = []
            for iy in np.arange(0, self.grid[1]): #for each cell in the row of interfaces
                Interface = x_interface(self.cells[ix-1][iy], self.cells[ix][iy]) 
                row.append(Interface)
                
            self.x_interfaces.append(row)
            
        #last row
        row = []
        for iy in np.arange(0, self.grid[1]): #for each cell in the last row of interfaces (southern boundary)
            Interface = x_interface(self.cells[-1][iy], self.BCs[2][iy]) #left cell = southern internal cell, right cell = southern ghost cell
            row.append(Interface) 
            
        self.x_interfaces.append(row)    
            
        #y-sweep    
        for ix in np.arange(0, self.grid[0]):
            row = []
            Interface = y_interface(self.BCs[3][ix], self.cells[ix][0]) #left cell = western ghost cell, right cell = first internal cell
            row.append(Interface)
            
            for iy in np.arange(1, self.grid[1]):
                Interface = y_interface(self.cells[ix][iy-1], self.cells[ix][iy])
                row.append(Interface)
            
            Interface = y_interface(self.cells[ix][-1], self.BCs[1][ix]) #left cell = last internal cell, right cell = eastern ghost cell
            row.append(Interface)
            
            self.y_interfaces.append(row)
            
        #Connect cells to interfaces
        for ix in np.arange(0, self.grid[0]):
            for iy in np.arange(0, self.grid[1]):
                self.cells[ix][iy].assign_interfaces(self.x_interfaces[ix][iy], self.x_interfaces[ix+1][iy], self.y_interfaces[ix][iy], self.y_interfaces[ix][iy+1])
        
    def calculate_fluxes(self):     
        
        self.Sy_max = 0
        
        #across interfaces aligned with the x-axis
        for ix in np.arange(0, len(self.x_interfaces)):
            row = self.x_interfaces[ix]
            for iy in np.arange(0, len(row)):
                row[iy].wavespeeds()
                self.Sy_max = max(abs(row[iy].left_wavespeed), abs(row[iy].right_wavespeed), self.Sy_max)
                row[iy].HLL_Solver()
                
                self.x_depth_fluxes[ix][iy] = row[iy].depth_flux
                self.x_momentum_fluxes[ix][iy] = row[iy].momentum_flux
        
        self.Sx_max = 0
        
        #across interfaces aligned with the y-axis                
        for ix in np.arange(0, len(self.y_interfaces)):
            row = self.y_interfaces[ix]
            for iy in np.arange(0, len(row)):
                row[iy].wavespeeds()
                self.Sx_max = max(abs(row[iy].left_wavespeed), abs(row[iy].right_wavespeed), self.Sx_max)
                row[iy].HLL_Solver()
                
                self.y_depth_fluxes[ix][iy] = row[iy].depth_flux
                self.y_momentum_fluxes[ix][iy] = row[iy].momentum_flux
    
    def update_cells(self):       
        
        self.tstep = (0.95*self.cell_width)/(self.Sx_max+self.Sy_max)
        self.sim_time += self.tstep
        
        for ix in np.arange(0, len(self.cells)):
            row = self.cells[ix]
            for iy in np.arange(0, len(row)):
                row[iy].update(self.tstep)
                
                self.depths[ix][iy] = row[iy].depth
                self.x_momentums[ix][iy] = row[iy].x_momentum
                self.y_momentums[ix][iy] = row[iy].y_momentum
                
        for ix in np.arange(0, len(self.BCs)):
            row = self.BCs[ix]
            for iy in np.arange(0, len(row)):
                row[iy].update()
        
        print(f't = {self.sim_time}')
        
    def mass_conservation_check(self):
        
        self.masses.append(np.sum(self.depths)) 
        
            
###################################################################################################################################################################################
# Numerical Scheme
##################################################################################################################################################################################

def scheme(rows, cols, init_depth, init_momentum_x, init_momentum_y, boundary_conditions, t_end, cell_width):
    
    Domain = domain(rows, cols, init_depth, init_momentum_x, init_momentum_y, boundary_conditions, cell_width)
    Domain.generate_cells_and_interfaces()
    
    while Domain.sim_time <= t_end:
        Domain.calculate_fluxes()
        Domain.update_cells()
            
            
            
            
            
            
            
            
            
            
            
            
        
        