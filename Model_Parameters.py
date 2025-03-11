
'*****************************************************************************************************************************************************************'
import pandas as pd

'*****************************************************************************************************************************************************************'
#Modelling Parameters Smoluchowski Coagulation Model.

'*****************************************************************************************************************************************************************'
#Control Volume / Sample Parameters / Flow Regime

Sample_Volume = 25 #Size of the Sample Container in the LISST 200X in mL which is equivalent cubic centimeters.
mld = 100  # Mixed Layer Depth/Fluid Thickness, Units: cm
Flow ='Turbulent' #Flow Regime: 'Laminar' or 'Turbulent

'*****************************************************************************************************************************************************************'
#Kernel Function Parameters.

K_B = 1.381e-23  # Boltzmann Constant, Units: J/K
T_Celcius = 20  # Input Temperature in Celcius
T_th = T_Celcius + 273.15   # Thermodynamic Temperature, Units: K
rho_f = 1.025  # Fluid Density, Units: g/cm^3, Note: Dependent on Thermodynamic Temperature
rho_p = 1.200  # Diatom/Cocolithospheres  # Particle Density, Units: g/cm^3 #balance between biological/chemical factors.
mu = 10**-3 # To be set  # Dynamic Viscosity, Units: kg/(m*s), Note: Dependent on Thermodynamic Temperature
nu = mu/(rho_f *1000) # To be set  # Kinematic Viscosity, Units: m^2/s, Note: Dependent on Thermodynamic Temperature
gamma = 1  # Shear Gradient, Units: 1/s
epsilon = 10**-3 # Turbulent Energy Dissipation Rate, Units: m^2/s^3
a_g = 9.81  # Acceleration Due to Gravity, Units: m/s^2

'*****************************************************************************************************************************************************************'
#Coagulation Solver Parameters

alpha = 1 # Probability of Coagulation, Units: Dimensionless
Kernel = 'all' # Coagulation Kernel Type: 'all' for all kernels, 'br' for Brownian, 'ls' for laminar shear, 'ts' for turbulent shear,'fs' for total Fluid Shear, 'ds' for Total differential settling.
Type = 'rec' #Coagulation Kenrel type: 'rec' for Rectilinear, 'cur' for Curvilinear.
t_span = (0, 3600) # Time span for the simulation (start and end time) in seconds.
set_time_step = False # True for manual time steping, False for automatic time steping.
dt = 0.1  # Time step for the simulation in seconds for manual time stepping.
include_settling =  False # True to include settling in the model, False to exclude settling in the model.



'*****************************************************************************************************************************************************************'

'*****************************************************************************************************************************************************************'
# Update Descriptions Here

descriptions = {
    "Sample_Volume": "Size of the Sample Container in the LISST 200X in mL which is equivalent cubic centimeters.",
    "mld": "Mixed Layer Depth/Fluid Thickness, Units: cm",
    "Flow": "Flow Regime: 'Laminar' or 'Turbulent",
    "K_B": "Boltzmann Constant, Units: J/K",
    "T_Celcius": "Input Temperature in Celcius",
    "T_th": "Thermodynamic Temperature, Units: K",
    "rho_f": "Fluid Density, Units: g/cm^3, Note: Dependent on Thermodynamic Temperature",
    "rho_p": "Diatom/Cocolithospheres Particle Density, Units: g/cm^3, balance between biological/chemical factors.",
    "mu": "Dynamic Viscosity, Units: kg/(m*s), Note: Dependent on Thermodynamic Temperature",
    "nu": "Kinematic Viscosity, Units: m^2/s, Note: Dependent on Thermodynamic Temperature",
    "gamma": "Shear Gradient, Units: 1/s",
    "epsilon": "Turbulent Energy Dissipation Rate, Units: m^2/s^3",
    "a_g": "Acceleration Due to Gravity, Units: m/s^2",
    "alpha": "Stickiness / Probability of Coagulation, Units: Dimensionless",
    "Kernel": "Coagulation Kernel",
    "Type": "Coagulation Kernel type",
    "t_span": "Time span for the simulation (start and end time) in seconds",
    "set_time_step": "True for manual time stepping, False for automatic time stepping",
    "dt": "Time step for the simulation in seconds for manual time stepping",
    "include_settling": "True to include settling in the model, False to exclude settling in the model"
}

# Functions to update parameter table
def capture_parameters():
    # Use the global scope to capture variables
    global_vars = globals()

    # Collect the values and descriptions
    parameters = []
    for name, description in descriptions.items():
        if name in global_vars:
            value = global_vars[name]
            parameters.append((name, value, description))

    # Create a DataFrame
    df = pd.DataFrame(parameters, columns=["Parameter", "Value", "Description"])
    return df

# Capture and print the parameters in a table format
Params_df = capture_parameters()