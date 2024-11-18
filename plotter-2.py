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
            numerical_temps = []  
            analytical_temps = [] 
            
            while True:
                try:
                    # Read numerical solution
                    record1 = f.read_reals(dtype=np.float32)
                    times.append(record1[1])        
                    numerical_temps.append(record1[2:]) 
                    
                    # Read analytical solution
                    record2 = f.read_reals(dtype=np.float32)
                    analytical_temps.append(record2[2:]) 
                    
                except Exception:
                    break
                    
        return x, np.array(times), np.array(numerical_temps, order='F'), np.array(analytical_temps, order='F')
    
    except FileNotFoundError:
        raise FileNotFoundError(f"Could not find file: {filename}")
    except Exception as e:
        raise Exception(f"Error reading data: {str(e)}")

def plot_temperature_profiles(x, times, numerical_temps, analytical_temps, output_path):
    """Create and save temperature profile plot."""

    sns.set_theme(context='notebook', style='darkgrid')

    fig, ax = plt.subplots(figsize=(10, 6))
    
    colors = plt.cm.Set2(np.linspace(0, 1, len(times)//5 + 1))
     
    num_line = ax.plot([], [], '-', color='gray', linewidth=2.5, label='Numerical')
    ana_line = ax.plot([], [], '--', color='gray', linewidth=2, label='Analytical')
    
    # Plot temperature profiles
    for i, color in zip(range(0, len(times), max(1, len(times)//5)), colors):
        ax.plot(x, numerical_temps[i], '-', color=color, linewidth=2.5, alpha=0.9)
        ax.plot(x, analytical_temps[i], '--', color=color, linewidth=2, alpha=0.9)
        
    # Add time annotations
    time_labels = [f'$t = {t:.0f}$' for t in times[::max(1, len(times)//5)]]
    legend1 = ax.legend(num_line + ana_line, ['Numerical', 'Analytical'], 
                       loc='upper left', fontsize=14)
    legend2 = ax.legend(handles=[plt.Line2D([], [], color=c, label=l) 
                                for c, l in zip(colors, time_labels)],
                       loc='upper right', fontsize=14)
    
    ax.add_artist(legend1)
    ax.add_artist(legend2)

    ax.set_xlabel(r'$x$', fontsize=18, labelpad=10)
    ax.set_ylabel(r'$T$', fontsize=18, labelpad=10)
    
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(True, linestyle='--', alpha=0.7)


    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.5)

    plt.tight_layout()
    
    plt.savefig(output_path, 
                bbox_inches='tight', 
                dpi=500,
                format='pdf',
                edgecolor='none')
    
    plt.close()

def main():

    DATA_DIR = Path('data')
    OUTPUT_DIR = Path('figs')
     
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    input_file = DATA_DIR / 'diff2.da'
    output_file = OUTPUT_DIR / 'temperature_comparison.pdf'

    try:
        x, times, numerical_temps, analytical_temps = read_fortran_data(input_file)
        plot_temperature_profiles(x, times, numerical_temps, analytical_temps, output_file)
        print(f"Successfully created plot: {output_file}")
        
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    main()