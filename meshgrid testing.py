import numpy as np

test = np.mgrid[0:10,0:10]

dx = 0.1
rows, cols = 10, 10
nx, ny = (rows+1, cols+1)
x = np.linspace(0, 1, nx)
y = np.linspace(0, 1, ny)
xv, yv = np.meshgrid(x, y)