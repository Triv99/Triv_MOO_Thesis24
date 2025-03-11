'******************************************************************************************************************************************************************************'
import numpy as np
from Model_Parameters import *

#Cluster Size Functions

def mass2diameter(mass):
    "Mass in grams, diameter in cm"
    return ((6 / (np.pi * rho_p)) * mass) ** (1/3)

def diameter2mass(diameter):
    "Diameters in cm, mass in grams"
    return  (1/6) * np.pi * rho_p * (diameter**3)

'******************************************************************************************************************************************************************************'
# Generating Functions

#Defining the diameter generating function
def generate_diameters(n, d_min, d_max, units):
    "Generate n diameters between d_min and d_max"
    if units == 'um':
        return np.linspace(d_min, d_max, n) * 10**-4 #Converting micrometers to centimeters
    elif units == 'cm':
        return np.linspace(d_min, d_max, n)
    
#Creating The mass domain function
def generate_mass_domain(diameters):
    "Generate n masses between d_min and d_max"
    "Units of diameters in cm, mass in grams"
    return diameter2mass(diameters)

'******************************************************************************************************************************************************************************'
# Kernel Functions

# Define K1, K2, and P
K1 = (3 / 4) / (np.pi * rho_p)
K2 = (rho_p - rho_f) / mu

def P(mi, mj):
    return min(mi**(1/3), mj**(1/3)) / max(mi**(1/3), mj**(1/3))

def Settling_Velocity(m):
    return (20/9) * a_g * K1**(2/3) * K2 * m**(2/3)

# Coagulation Kernels
def b_br(mi, mj):
    """Brownian Motion Kernel"""
    return ((2 * 10**6)/ 3) * (K_B * T_th / mu) * ((mi**(1/3) + mj**(1/3))**2) / ((mi * mj)**(1/3))

def ls_r(mi, mj):
    """Laminar Shear Rectilinear Kernel"""

    return (4 / 3) * gamma * K1 * ((mi**(1/3) + mj**(1/3))**3)

def ts_r(mi, mj):
    """Turbulent Shear Rectilinear Kernel"""
    return 1.3 * ((epsilon / nu)**0.5) * K1 * ((mi**(1/3) + mj**(1/3))**3)

def ts_c(mi, mj):
    """Turbulent Shear Curvilinear Kernel"""
    # The Square is now applicable to the whole (1 + 2p) term instead of just p. 
    p = P(mi, mj)
    return 9.8 * ((p**2) / ((1 + 2 * p)**2)) * ((epsilon / nu)**0.5) * K1 * ((mi**(1/3) + mj**(1/3))**3)

def ds_r(mi, mj):
    """Differential Settling Rectilinear Kernel"""
    return (20/ 9) * np.pi * a_g * K2 * K1**(4/3) * ((mi**(1/3) + mj**(1/3))**2) * abs(mj**(2/3) - mi**(2/3))

def ds_c(mi, mj):
    """Differential Settling Curvilinear Kernel"""
    # The curvilinear correction factor has been added.
    p = P(mi, mj)
    return (10 / 9) * ((p**2) / ((1 + p)**2)) * np.pi * a_g * K2 * K1**(4/3) * mi**(2/3) * abs(mi**(2/3) - mj**(2/3))

# Main Kernel Function
def beta(mi, mj, Kernel, Type, Flow):
    """Main Kernel Function to Determine Coagulation Kernel Based on Type and Kernel"""
    if Kernel == "all":
        if Type == "rec":
            if Flow == "Laminar":
                return b_br(mi, mj) + ls_r(mi, mj) + ds_r(mi, mj)
            elif Flow == "Turbulent":
                return b_br(mi, mj) + ts_r(mi, mj) + ds_r(mi, mj)            
        elif Type == "cur":
            if Flow == "Laminar":
                return b_br(mi, mj) + ls_r(mi, mj) + ds_c(mi, mj)
            elif Flow == "Turbulent":
                return b_br(mi, mj) + ts_c(mi, mj) + ds_c(mi, mj)
    elif Kernel == "br":
        return b_br(mi, mj)
    elif Kernel == "ls":
        return ls_r(mi, mj)
    elif Kernel == "ts":
        if Type == "rec":
            return ts_r(mi, mj)
        elif Type == "cur":
            return ts_c(mi, mj)
    elif Kernel == "ds":
        if Type == "rec":
            return ds_r(mi, mj)
        elif Type == "cur":
            return ds_c(mi, mj)
        
# Naming the Kernels        
b_br.__name__ = "Brownian Motion Kernel"
ls_r.__name__ = "Laminar Shear Rectilinear Kernel"
ts_r.__name__ = "Turbulent Shear Rectilinear Kernel"
ts_c.__name__ = "Turbulent Shear Curvilinear Kernel"
ds_r.__name__ = "Differential Settling Rectilinear Kernel"
ds_c.__name__ = "Differential Settling Curvilinear Kernel"

'******************************************************************************************************************************************************************************'
