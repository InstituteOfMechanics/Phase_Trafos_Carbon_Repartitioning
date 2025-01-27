#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot selected results for the quenching of the bearing race.

Note: this script creates plots for 3 points, but only point 1 and 3 are
used in the paper, where point 3 is called point 2.

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
    
    files = ["bearing_race_quenching_100Cr6_Tinf_20.npz",
             "bearing_race_quenching_100Cr6_Tinf_60.npz",
             "bearing_race_quenching_100CrMnSi64_Tinf_20.npz",
             "bearing_race_quenching_100CrMnSi64_Tinf_60.npz",
             ]
    
    labels = [r"100Cr6, $T_\infty=20$",
              r"100Cr6, $T_\infty=60$",
              r"100CrMnSi6-4, $T_\infty=20$",
              r"100CrMnSi6-4, $T_\infty=60$",
              ]
    
    styles = ["-", "-", "--", "--"]
    
    plt.close("all")
    
    # change order (model definition 1,2,3, paper uses only 1,3 as 1,2)
    point_numbers = [1, 3]
    point_labels = ["P1", "P2"]
    
    for point_num, point_label in zip(point_numbers, point_labels):
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1, sharex=True, figsize=[10/2.45,15/2.45])
    
        for i, file_path in enumerate(files):
            # label=f"P{point_label}, {labels[i]}"
            label=labels[i]
            style=styles[i]
            
            npz_path = result_path / file_path
            
            with np.load(npz_path) as file:
                time = file["times"]
                
                sdv_n = file[f"sdvs_n{point_num}"].T
                temp_n = file[f"temps_n{point_num}"]
                
                sdv_e = file[f"sdvs_e{point_num}"].T
                temp_e = file[f"temps_e{point_num}"]
                
            # conversion to s
            time *= 1e-3    
            
            selector = time <= 180
            
        
            ax1.plot(time[selector], temp_e[selector], style, label=label)
        
            ax2.plot(time[selector], sdv_e[selector,0], style, label=label)
        
            ax3.plot(time[selector], sdv_e[selector,1], style, label=label)
        
            ax4.plot(time[selector], sdv_e[selector,3], style, label=label)
        
        
        
        ax1.set_ylabel(r"Temperature $T$ / \si{\celsius}")
        
        ax1.grid()
        # ax1.legend()
        ax1.legend(bbox_to_anchor=(0.4, 1.60), loc='upper center', borderaxespad=0., ncol=2)
    
        ax2.set_ylabel(r"Martensite fraction $\beta_\mathrm{M}$")
        ax2.grid()
        # ax2.legend()
    
        ax3.set_ylabel(r"Bainite fraction $\beta_\mathrm{B}$")
        ax3.grid()
        # ax3.legend()
        
        ax4.set_xlabel(r"Time $t$ / \si{\second}")
        ax4.set_ylabel(r"Carbon content $x_{C,\mathrm{A}}$")
        ax4.grid()
        # ax4.legend()
        
        plt.tight_layout()    
        
        fig.align_labels()
        
        png_path = plot_path / f"results_bearing_race_{point_label}.png"

        plt.savefig(png_path, dpi=600)
        
        
