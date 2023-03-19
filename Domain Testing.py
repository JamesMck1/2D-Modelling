import numpy as np


class obj():
    
    def __init__(self):
        self.h = 1
        self.u = 0
        
    def update(self):
        self.h += 1
        self.u += 5

cells = []


for ix in np.arange(0,5):
    row = []
    for i in np.arange(0,5): #cells in a row
        cell = obj()
        row.append(cell)
        
    cells.append(row)
    
