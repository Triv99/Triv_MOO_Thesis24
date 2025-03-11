import pandas as pd
import numpy as np
import itertools
from scipy.interpolate import interp1d
from scipy.interpolate import PchipInterpolator
from scipy import integrate
from scipy.integrate import solve_ivp
from sklearn.metrics import mean_squared_error
from scipy.stats import norm
from Model_Parameters import *
from Model_Functions import *

def SCE_Algorithm(m_values, Initial_Concentrations, Kernel, Type, t_span, include_settling):
    """Function to Solve the System of ODEs for Coagulation using the SCE Algorithm"""
    """m_values: The mass domain of the system of ODEs"""
    """Initial_Concentrations: The initial concentrations of the system of ODEs"""
    """Kernel: The type of coagulation kernel to use"""
    """Type: The type of coagulation kernel to use"""
    """t_span: The time span for the simulation (start and end time) in seconds"""
    """include_settling: True to include settling in the model, False to exclude settling in the model"""

    # Define the domain of the system of ODEs
    mass_domain = m_values  # Use the masses as the domain

    # Define the initial conditions for the system of ODEs
    initial_conditions = Initial_Concentrations  # Get the initial concentrations

    def coagulation_system(t, C):
        dCdt = np.zeros_like(C)  # Initialize the derivative array

        # Add a check for the Type and Kernel before calling beta
        valid_kernels = {'all', 'br', 'ls', 'ts', 'ds'}
        valid_types = {'rec', 'cur'}
        if Kernel not in valid_kernels or Type not in valid_types:
            raise ValueError(f"Invalid Type or Kernel: Type={Type}, Kernel={Kernel}")
        
        # Create an interpolator dynamically based on current C and mass_domain
        f = interp1d(mass_domain, C, kind='cubic', fill_value='extrapolate')

        # Loop over the size bins
        for i, m in enumerate(mass_domain):
            # Compute the gain term for smaller particles combining
            gain_integral = integrate.quad(
                lambda mi: 0.5 * beta(m - mi, mi, Kernel, Type, Flow) * f(m - mi) * f(mi),
                0, m
            )[0] if m > 0 else 0

            # Compute the loss term for all particles of size m
            loss_integral = integrate.quad(
                lambda mi: beta(m, mi, Kernel, Type, Flow) * f(mi),
                m, mass_domain[-1]
            )[0]

            # Calculate the rate of change of concentration for this bin
            dCdt[i] = gain_integral - (f(m) * loss_integral)
            if include_settling:
                dCdt[i] -= f(m) * (Settling_Velocity(m) / mld)

        return dCdt

    # Define the time span and solve the system of ODEs
    t_eval = np.arange(t_span[0], t_span[1] + dt, dt) if set_time_step else None
    solution = solve_ivp(coagulation_system, t_span, initial_conditions, method='BDF', t_eval=t_eval)

    # Set negative concentrations to zero
    solution.y[solution.y < 0] = 0

    # Prepare the solution dataframe
    sol_df = pd.DataFrame(solution.y.T, index=solution.t, columns=mass_domain)
    sol_df.columns = [f'Number Concentration (1/cm3) at Mass {m:.3e} g' for m in mass_domain]
    sol_df.index.name = 'Time (s)'

    return sol_df
