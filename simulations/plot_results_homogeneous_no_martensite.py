"""
Create plots for the homogeneous problem without martensite evolution.

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


RESULT_DIR = "simulations_homogeneous_no_martensite"
PLOT_DIR = "plots"

if __name__ == "__main__":
    result_path = Path(RESULT_DIR)
    
    plot_path = Path(PLOT_DIR)
    plot_path.mkdir(exist_ok=True)
    
    files = ["homogeneous_100Cr6_T0_300_htc_10-0000_isothermal_no_martensite.npz",
             "homogeneous_100CrMnSi64_T0_300_htc_10-0000_isothermal_no_martensite.npz",
             "homogeneous_100Cr6_T0_300_htc_10-0000_isothermal_fixed_temp_no_martensite.npz",
             "homogeneous_100CrMnSi64_T0_300_htc_10-0000_isothermal_fixed_temp_no_martensite.npz",
             ]
    
    labels = ["100Cr6 Convection", "100CrMnSi6-4 Convection", "100Cr6 Dirichlet", "100CrMnSi6-4 Dirichlet"]
    
    styles = ["-", "-", "--", "--"]
    
    plt.close("all")
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1, sharex=True, figsize=[10/2.45, 15/2.45])
    
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
        ax2.plot(time, sdv[:,1], style)

       # subplot: carbon fraction
        ax3.plot(time, sdv[:,3], style)

        ax4.plot(time[1:], sdv[1:,9], style)
 
        
    
    ax1.set_ylabel(r"Temperature $T$ / \si{\celsius}")
    ax1.grid()

    ax2.set_ylabel(r"Bainite fraction $\beta_\mathrm{B}$")
    ax2.grid()
    
    ax3.set_ylabel(r"Carbon content $x_{C,\mathrm{A}}$")
    ax3.grid()

    ax4.set_ylabel(r"Trafo factor $f$")
    ax4.grid()
    ax4.set_xlabel(r"Time $t$  / \si{\second}")    
    
    ax1.legend(bbox_to_anchor=(0.4, 1.50), loc='upper center', borderaxespad=0., ncol=2)
    
    plt.tight_layout()    
    fig.align_ylabels()
    
    png_path = plot_path / "results_homogeneous_no_martensite.png"

    plt.savefig(png_path, dpi=600)
