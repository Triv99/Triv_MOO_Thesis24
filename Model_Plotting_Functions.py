'*********************************************************************************************************************************************************************************'
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import numpy as np

from Model_Parameters import *
from Model_Functions import *

# Plotting Code Functions For the Coagulation Kernels

'*********************************************************************************************************************************************************************************'

def KernelHeatMaps(P1Diameters, P2Diameters):
    "Function to Generate Heatmaps for all Coagulation Kernels"
    "P1Diameters and P2Diameters are the diameters of the particles in cm"
    "The Output is a 2x3 grid of heatmaps for each kernel type, with the color being the rate of coagulation in cm³/s"

    # Calculate values for m1 and m2
    m1_values = diameter2mass(P1Diameters)
    m2_values = diameter2mass(P2Diameters)

    # Define beta functions, assuming they are already defined
    beta_functions = [b_br, ls_r, ts_r, ts_c, ds_r, ds_c]
    labels = ['Brownian Motion', 'Laminar Shear Rectilinear', 'Turbulent Shear Rectilinear', 
              'Turbulent Shear Curvilinear', 'Differential Settling Rectilinear', 'Differential Settling Curvilinear']
    results = []

    for beta in beta_functions:
        grid = np.array([[beta(m1, m2) for m1 in m1_values] for m2 in m2_values])
        results.append(grid)

    # Convert diameters from cm to µm
    P1Diameters_um = [diameter * 1e4 for diameter in P1Diameters]
    P2Diameters_um = [diameter * 1e4 for diameter in P2Diameters]

    # Plotting
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    tick_spacing = 10  # This determines the spacing between ticks, change it as needed
    for i, ax in enumerate(axes.flat):
        heatmap = sns.heatmap(results[i], ax=ax, cmap='viridis', cbar_kws={'label': 'Rate (cm³/s)'})
        ax.set_title(labels[i])
        ax.set_xlabel('Particle 1 Diameter (µm)')
        ax.set_ylabel('Particle 2 Diameter (µm)')
        
        # Set tick locations and labels with reduced frequency
        ax.set_xticks(np.arange(0, len(P1Diameters_um), tick_spacing))
        ax.set_yticks(np.arange(0, len(P2Diameters_um), tick_spacing))
        ax.set_xticklabels(['{:.1f}'.format(P1Diameters_um[idx]) for idx in np.arange(0, len(P1Diameters_um), tick_spacing)], rotation=45)
        ax.set_yticklabels(['{:.1f}'.format(P2Diameters_um[idx]) for idx in np.arange(0, len(P2Diameters_um), tick_spacing)], rotation=0)

    plt.tight_layout()
    return plt.show()

'**********************************************************************************************'

def SpecificKernelPlot(P1Diameters, P2Diameters, kernel_function):
    "Function to Generate a Heatmap and a rotating 3D Plot for a Specific Coagulation Kernel"
    "P1Diameters and P2Diameters are the diameters of the particles in cm"

    # Calculate values for m1 and m2
    m1_values = diameter2mass(P1Diameters)
    m2_values = diameter2mass(P2Diameters)

    # Calculate grid for the kernel function
    grid = np.array([[kernel_function(m1, m2) for m1 in m1_values] for m2 in m2_values])

    # Get the kernel name from the function's __name__ attribute
    kernel_name = kernel_function.__name__.replace('_', ' ').title()

    # Plotting
    fig = plt.figure(figsize=(14, 7))

    # Heatmap
    ax1 = fig.add_subplot(1, 2, 1)
    sns.heatmap(grid, ax=ax1, cmap='viridis', cbar_kws={'label': 'Rate (cm³/s)'})
    ax1.set_title(f'Heatmap of {kernel_name}')
    ax1.set_xlabel('P1 Diameter (cm)')
    ax1.set_ylabel('P2 Diameter (cm)')
    ax1.set_xticks(np.linspace(0, len(P1Diameters) - 1, num=10))
    ax1.set_yticks(np.linspace(0, len(P2Diameters) - 1, num=10))
    ax1.set_xticklabels(['{:.4f}'.format(d) for d in P1Diameters[::len(P1Diameters)//10]], rotation=45)
    ax1.set_yticklabels(['{:.4f}'.format(d) for d in P2Diameters[::len(P2Diameters)//10]], rotation=0)

    # 3D spinning plot
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    ax2.set_title(f'3D Plot of {kernel_name}')
    X, Y = np.meshgrid(m1_values, m2_values)
    surf = ax2.plot_surface(X, Y, grid, cmap='viridis', edgecolor='none')

    # Animation function
    def update(frame):
        ax2.view_init(elev=30., azim=frame)
        return fig,

    # Animate
    anim = FuncAnimation(fig, update, frames=range(0, 360, 2), interval=150, blit=False)

    plt.tight_layout()
    plt.close()

    # Return animation object for IPython display
    return HTML(anim.to_html5_video())

'*********************************************************************************************************************************************************************************'
# Plotting Code Functions For the Model Solutions

'*********************************************************************************************************************************************************************************'
def InitialConcentrationPlot(diameters, Initial_Concentrations, logscalex = True):
    "Function to plot the initial number concentration of particles against particle diameter"
    "Diameters are in cm, Initial Concentrations are in 1/cm³"
    "logscalex is a boolean to determine if the x-axis should be in log scale"
    
    # Plotting
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(diameters, Initial_Concentrations, marker='o')
    ax.set_title('Initial Number Concentration of Particles vs. Particle Diameter')
    ax.set_xlabel('Particle Diameter (cm)')
    ax.set_ylabel('Initial Number Concentration (1/cm³)')
    if logscalex:
        ax.set_xscale('log')
    plt.grid()
    
    return plt.show()

'**********************************************************************************************'

def AnimatedPSD(sol_df, diameters, logscalex= False):

    "Function to generate an animated plot of the number concentration of particles over time"
    "sol_df is the solution dataframe, diameters are in cm"
    "logscalex is a boolean to determine if the x-axis should be in log scale"
    
    # Create a new figure and axis
    fig, ax = plt.subplots(figsize=(10, 6))
    line, = ax.plot([], [], lw=2)  # Initialize a line object, note the comma
    ax.set_xlabel('Particle Diameter (cm)')
    if logscalex:
        ax.set_xscale('log')
    ax.set_xlim(min(diameters), max(diameters))
    ax.set_ylim(min(sol_df.min()), max(sol_df.max()))
    ax.set_ylabel('Number Concentration (1/cm³)')
    ax.set_title('Number Concentration of Particles Over Time')

    # Initialize function for the animation
    def init():
        line.set_data([], [])
        return line,

    # Update function for animation
    def update(frame):
        x_data =  diameters 
        y_data = sol_df.iloc[frame]  # Get the row corresponding to the current frame by index
        line.set_data(x_data, y_data)
        ax.set_title(f'Number Concentration of Particles at Time {sol_df.index[frame]:.2f}s')
        return line,

    # Create animation
    ani = FuncAnimation(fig, update, frames=len(sol_df.index), init_func=init, blit=True, interval=500)

    plt.close()

    # Display the animation in the notebook
    return HTML(ani.to_html5_video())

'**********************************************************************************************'

def plot_quartile_data(sol_df, diameters, logscalex=False):

    "Function to plot the quartile data of the number concentration of particles over time"
    "sol_df is the solution dataframe, diameters are in cm"
    "logscalex is a boolean to determine if the x-axis should be in log scale"

    # Determine the indices for 0%, 25%, 50%, and 75%
    quartile_indices = [int(len(sol_df) * q) for q in [0, 0.25, 0.5, 0.75]]

    # Extract the quartile data
    quartile_data = sol_df.iloc[quartile_indices]
    quartile_times = sol_df.index[quartile_indices]

    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 8))
    for i, (index, time) in enumerate(zip(quartile_indices, quartile_times)):
        ax.plot(diameters, quartile_data.iloc[i], label=f'Time: {time:.2f} s')
        
    ax.set_xlabel('Particle Diameter (cm)')
    if isinstance(logscalex, bool) and logscalex:
        ax.set_xscale('log')
    ax.set_xlim(min(diameters), max(diameters))
    ax.set_ylim(sol_df.min().min(), sol_df.max().max())
    
    ax.set_ylabel('Number Concentration (1/cm3)')
    ax.set_title('Number Concentration at Different Time Points')
    ax.legend()
    ax.grid(True)

    return plt.show()

'**********************************************************************************************'

def Timeseriesplots(sol_df):

    "Function to plot the analytical solution, sum of number concentration, and total mass of the system over time"
    "sol_df is the solution dataframe"

    # Create 1x3 subplots
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 6))

    # Plot the analytical solution
    def concentration_at_time(t):
        pi = np.pi
        C0 = sol_df.iloc[0].sum()
        return C0 * np.exp(-4 * alpha * gamma * t / pi)

    t = np.linspace(0, 1000, 1000)
    conc_t = concentration_at_time(t)
    conc_t
    ax1.plot(t, conc_t)
    ax1.set_xlabel('Time (s)')
    ax1.set_xscale('log')
    ax1.set_ylim(0, 50000)
    ax1.set_ylabel('Number Concentration (1/cm³)')
    ax1.set_title('Analytical Solution for Number Concentration Over Time')

    # Plot the row sums against the time index column
    # Sum each row of the sol_df DataFrame
    row_sums = sol_df.sum(axis=1)
    ax2.plot(sol_df.index, row_sums)
    ax2.set_xlabel('Time (s)')
    ax2.set_xscale('log')
    ax2.set_ylabel('Sum of Number Concentration (1/cm³)')
    ax2.set_title('Sum of Number Concentration Over Time')

    #Plot the Total Mass of the system over time
    masses = sol_df.columns.str.extract(r'(\d+\.\d+e-\d+)').astype(float)[0]
    masses = masses.values
    mass_sums = ((sol_df * masses).sum(axis=1)) * Sample_Volume
    ax3.plot(sol_df.index, mass_sums)
    ax3.set_xlabel('Time (s)')
    ax3.set_xscale('log')
    ax3.set_ylabel('Total Mass (g)')
    ax3.set_title('Total Mass of the System Over Time')

    # Tight layout to prevent overlap
    plt.tight_layout()

    return plt.show()

'*********************************************************************************************************************************************************************************'
