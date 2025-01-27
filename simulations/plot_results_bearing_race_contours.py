"""
Plot the deformation of the inner contour for the bearing race example.

"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

RESULT_DIR = "simulations_bearing_race"
PLOT_DIR = "plots"

# Enable LaTeX rendering
plt.rc('text', usetex=True)

# Set font to DejaVuSans
plt.rc('font', **{'family': 'serif', 
               'size': 9
              })

plt.rcParams['text.latex.preamble'] = "\n".join([r"\usepackage[utf8]{inputenc}",
                                                  r"\usepackage[T1]{fontenc}",
                                                 r'\usepackage{siunitx}'])

if __name__ == "__main__":

    result_path = Path(RESULT_DIR)
    
    plot_path = Path(PLOT_DIR)
    plot_path.mkdir(exist_ok=True)
    
    files = ["bearing_race_quenching_100Cr6_Tinf_20_displacements.npz",
             "bearing_race_quenching_100CrMnSi64_Tinf_20_displacements.npz",
             "bearing_race_quenching_100Cr6_Tinf_60_displacements.npz",
             "bearing_race_quenching_100CrMnSi64_Tinf_60_displacements.npz",
             ]
    
    labels = [r"100Cr6, $T_\infty=20$",
              r"100CrMnSi6-4, $T_\infty=20$",
              r"100Cr6, $T_\infty=60$",
              r"100CrMnSi6-4, $T_\infty=60$",
              ]
    
    styles = ["-", "--", "-", "--"]
    
    plt.close("all")
    fig, ax1 = plt.subplots(1,1, sharex=True, figsize=[10/2.45,6/2.45])

    for i, file_path in enumerate(files):
        label=labels[i]
        style=styles[i]
        
        npz_path = result_path / file_path
        
        with np.load(npz_path) as file:
            x = file["y"]
            x = x - x.min()
            uy = file["ux"]
    
        ax1.plot(x, uy, label=label)
        
    
    
    ax1.set_xlabel(r"Axial coordinate $z$ / mm")
    ax1.set_ylabel(r"Radial displacement $u_r$ / mm")
    
    ax1.grid()
    ax1.legend(loc="upper center", bbox_to_anchor=(0.5,0.85))
    
    plt.tight_layout()    
    
    
    png_path = plot_path / "results_bearing_race_contours.png"
    plt.savefig(png_path, dpi=600)
