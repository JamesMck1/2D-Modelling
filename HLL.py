# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 19:00:46 2023

@author: James Mckenna
"""

import math

def q_k(depth_k, h_star): #determine whether the wave is a shock or rarefaction
    if h_star > depth_k: #shock wave
        q_k = math.sqrt(0.5*(((h_star + depth_k)*h_star)/(depth_k**2)))
    else: #rarefaction wave
        q_k = 1
        
    return(q_k)

def Wavespeeds(g, h_L, h_R, u_L, u_R):
    c_L = math.sqrt(g*h_L) #left wave celerity
    c_R = math.sqrt(g*h_R) #right wave celerity
    
    h_0 = (1/g)*((0.5*(c_L+c_R)) + 0.25*(u_L-u_R))**2 #Initial estimate using two-rarefaction approximation
    
    if h_0 <= min(h_L, h_R):
        h_star = h_0 #Use two-rarefaction approximation
    elif h_0 > min(h_L, h_R):
        p_L = math.sqrt(((g)*(h_0+h_L))/(2*h_0*h_L))
        p_R = math.sqrt(((g)*(h_0+h_R))/(2*h_0*h_R))
        h_star = (p_L*h_L + p_R*h_R + u_L - u_R)/(p_L + p_R)
    else:
        print('Error Calculating h_star depth')
    
    S_L = u_L - c_L*q_k(h_L, h_star) #Left extreme wavespeed
    S_R = u_R + c_R*q_k(h_R, h_star) #Right extreme wavespeed
    
    return(S_L, S_R)

def SWE_Flux(g, h, u): #Shallow water fluxes
    f1 = h*u #depth flux
    f2 = h*(u**2) + 0.5*g*(h**2) #momentum flux
    
    return(f1, f2)

def HLL_Flux(g, u_L, u_R, h_L, h_R, S_L, S_R): #HLL flux for shallow water equations
    F_L = SWE_Flux(g, h_L, u_L) #Left flux
    F_R = SWE_Flux(g, h_R, u_R) #Right flux
    
    f1 = (S_R*F_L[0] - S_L*F_R[0] + S_R*S_L*(h_R - h_L))/(S_R - S_L) #depth flux
    f2 = (S_R*F_L[1] - S_L*F_R[1] + S_R*S_L*(h_R*u_R - h_L*u_L))/(S_R - S_L) #momentum flux
    
    return(f1, f2)
    
def HLL_Solver(g, u_L, u_R, h_L, h_R, S_L, S_R): #Determine depth and momentum fluxes for two adjacent cells
    if S_L <= 0 <= S_R: #interface lies within the star region
        Fluxes = HLL_Flux(g, u_L, u_R, h_L, h_R, S_L, S_R) #HLL Flux
        
    elif S_R <= 0: #interface lies within the right state
        Fluxes = SWE_Flux(g, h_R, u_R) 
    
    else: # S_L >= 0, interface lies within the left state
        Fluxes = SWE_Flux(g, h_L, u_L)

    f1 = Fluxes[0] #depth flux
    f2 = Fluxes[1] #momentum flux
    
    return(f1, f2)