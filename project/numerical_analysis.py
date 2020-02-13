# Main Structural Analysis model we need to develop

# Imports
import numpy as np


# CODE...


# Spanwise station
Nx = 41


la = 2.691

def th_xi(i, Nx=Nx):
    return (i-1)/Nx * np.pi

def xi(i):
    return .5* (la/2*(1-np.cos(th_xi(i))) + la/2*(1-np.cos(th_xi(i+1)) ) )

print(xi(41))