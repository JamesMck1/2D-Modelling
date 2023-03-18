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

class exterior_cell(cell):
    
    def __init__(self, depth, x_momentum, y_momentum, cell_width, mannings_n):
        super().__init__(cell_width, mannings_n)
        
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
    
    def __init__(self):
        pass