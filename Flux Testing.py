import numpy as np
import HLL

g = 9.81
h_L = 0.25
h_R = 1.5
u_L = 0
u_R = 0
S = HLL.Wavespeeds(g, h_L, h_R, u_L, u_R)
x_flux_L = HLL.HLL_Solver(g, u_L, u_R, h_L, h_R, S[0], S[1])

g = 9.81
h_L = 1.5
h_R = 0.25
u_L = 0
u_R = 0
S = HLL.Wavespeeds(g, h_L, h_R, u_L, u_R)
x_flux_R = HLL.HLL_Solver(g, u_L, u_R, h_L, h_R, S[0], S[1])

g = 9.81
h_L = 0.25
h_R = 1.5
u_L = 0
u_R = 0
S = HLL.Wavespeeds(g, h_L, h_R, u_L, u_R)
y_flux_L = HLL.HLL_Solver(g, u_L, u_R, h_L, h_R, S[0], S[1])

g = 9.81
h_L = 1.5
h_R = 0.25
u_L = 0
u_R = 0
S = HLL.Wavespeeds(g, h_L, h_R, u_L, u_R)
y_flux_R = HLL.HLL_Solver(g, u_L, u_R, h_L, h_R, S[0], S[1])