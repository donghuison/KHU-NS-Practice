import numpy as np
from scipy.io import FortranFile
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


def read_fortran_data(filename):
    """Read data from Fortran unformatted file."""
    try:
        with FortranFile(filename, 'r') as f:
            # Read jmax 
            jmax = f.read_ints()[0]
            
            # Read x coordinates 
            x = f.read_reals(dtype=np.float32)
            
            # Read time steps and temperature data
            times = []
            temperatures = []
            
            while True:
                try:
                    record = f.read_reals(dtype=np.float32)
                    times.append(record[1])        # time
                    temperatures.append(record[2:]) # temperature array
                except Exception:
                    break
                    
        return x, np.array(times), np.array(temperatures, order='F')
    
    except FileNotFoundError:
        raise FileNotFoundError(f"Could not find file: {filename}")
    except Exception as e:
        raise Exception(f"Error reading data: {str(e)}")

def plot_temperature_profiles(x, times, temperatures, output_path):
    """Create and save temperature profile plot."""
    # Set up the plot style
    sns.set_theme(context='notebook', style='darkgrid')
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    
    # Color palette for better visualization
    colors = plt.cm.viridis(np.linspace(0, 1, len(times)//5 + 1))
    
    # Plot temperature profiles
    for i, color in zip(range(0, len(times), max(1, len(times)//5)), colors):
        ax.plot(x, temperatures[i], 
                label=f'$t = {times[i]:.0f}$',
                color=color,
                linewidth=2)
    
    # Customize plot
    ax.set_xlabel(r'$x$', fontsize=18, labelpad=10)
    ax.set_ylabel(r'$T$', fontsize=18, labelpad=10)
    ax.legend(loc='best', fontsize=14)
    
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.tick_params(axis='both', which='major', labelsize=14)
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(output_path, format='pdf', bbox_inches='tight', dpi=500)
    plt.close()

def main():

    DATA_DIR = Path('data')
    OUTPUT_DIR = Path('figs')
    
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    input_file = DATA_DIR / 'diff.da'
    output_file = OUTPUT_DIR / 'temperature_profiles.pdf'
    
    try:
        x, times, temperatures = read_fortran_data(input_file)
        plot_temperature_profiles(x, times, temperatures, output_file)
        print(f"Successfully created plot: {output_file}")
        
    except Exception as e:
        print(f"Error: {str(e)}")
        
if __name__ == "__main__":
    main()