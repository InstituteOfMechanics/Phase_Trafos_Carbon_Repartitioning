"""
Create plots for the homogeneous cooling problem.

"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Enable LaTeX rendering
plt.rc('text', usetex=True)

# Set font to DejaVuSans
plt.rc('font', **{'family': 'serif', 
               'size': 9
              })

plt.rcParams['text.latex.preamble'] = "\n".join([r"\usepackage[utf8]{inputenc}",
                                                  r"\usepackage[T1]{fontenc}",
                                                 r'\usepackage{siunitx}'])

RESULT_DIR = "simulations_homogeneous_cooling"
PLOT_DIR = "plots"

if __name__ == "__main__":
    result_path = Path(RESULT_DIR)
    
    plot_path = Path(PLOT_DIR)
    plot_path.mkdir(exist_ok=True)
    
    files = ["homogeneous_100Cr6_T0_500_htc_1-0000.npz",
             "homogeneous_100CrMnSi64_T0_500_htc_1-0000.npz",
             ]
    
    labels = ["100Cr6", "100CrMnSi6-4"]
    
    styles = ["-", "--",]
    
    plt.close("all")
    fig1, (ax1, ax2, ax3) = plt.subplots(3,1, sharex=True, figsize=[10/2.45,10/2.45])
    fig2, ax4 = plt.subplots(1,1, figsize=[10/2.45, 4/2.45])
    
    for i, file_path in enumerate(files):
        label=labels[i]
        style=styles[i]
    
        npz_path = result_path / file_path
    
        with np.load(npz_path) as file:
            time = file["times"]
            sdv = file["sdvs"].T
            temp = file["temps"]
            strain = file["strain"]

        # conversion to s
        time *= 1e-3
        
        
        ax1.plot(time, temp, style, label=label)
        
        # subplot: microstructure
        ax2.plot(time, sdv[:,0], style, label=r"$\beta_\mathrm{M}$ " + label)
        ax2.plot(time, sdv[:,1], style, label=r"$\beta_\mathrm{B}$ " + label)
                
       # subplot: carbon fraction
        ax3.plot(time, sdv[:,3], style, label=label)
        
        ax4.plot(temp, strain[:,0], style, label=label)
 
        
    ax1.set_ylabel(r"Temperature $T$ / \si{\celsius}")
    ax1.legend()
    ax1.grid()
    
    ax2.set_ylabel(r"Volume fraction $\beta$")
    ax2.legend()
    ax2.grid()
    
    ax3.set_ylabel(r"Carbon content $x_{C,\mathrm{A}}$")
    ax3.legend()
    ax3.grid()

    ax3.set_xlabel(r"Time $t$  / \si{\second}")
    
    ax4.invert_xaxis()
    ax4.set_xlabel(r"Temperature $T$  / \si{\celsius}")
    ax4.set_ylabel(r"Strain $\varepsilon_{11}$")
    ax4.grid()
    ax4.legend()
    
    fig1.tight_layout()
    fig1.align_labels()

    png_path1 = plot_path / "results_homogeneous_cooling.png"
    plt.figure(1)
    plt.savefig(png_path1, dpi=600)
    
    fig2.tight_layout()  
    
    png_path2 = plot_path / "results_homogeneous_cooling_strain.png"
    plt.figure(2)
    plt.savefig(png_path2, dpi=600)
    
