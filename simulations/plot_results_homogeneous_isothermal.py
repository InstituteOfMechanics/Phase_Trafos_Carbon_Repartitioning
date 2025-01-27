#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot selected results for the one element test.

"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


RESULT_DIR = "simulations_homogeneous_isothermal"
PLOT_DIR = "plots"

REFERENCE_TIMES = {(320, 0.01): 10**1.6799144609163803,
                   (320, 0.99): 10**2.556993526859674,
                   (350, 0.99): 10**2.1265851347960307,
                   (350, 0.01): 10**1.126996879172153,
                   (380, 0.99): 10**1.888071668161266,
                   (380, 0.01): 10**0.8852309561741873,
                   }


if __name__ == "__main__":
    result_path = Path(RESULT_DIR)
    plot_path = Path(PLOT_DIR)
    
    plt.close("all")
    fig, ax1 = plt.subplots(1,1, sharex=True, figsize=[6,6])
    
    for T in (320, 350, 380):
        file_name = f"homogeneous_100Cr6_T0_{T}_htc_1-0000_isothermal_fixed_temp.npz"
        npz_path = result_path / Path(file_name)
    
        with np.load(npz_path) as file:
            time = file["times"]
            sdv = file["sdvs"].T
            temp = file["temps"]
            strain = file["strain"]
            betaB = sdv[:,1]
            
        # conversion to s
        time *= 1e-3
        
        t01, t99 = np.interp((0.01, 0.99), betaB, time)
        t01ref = REFERENCE_TIMES[(T, 0.01)]
        t99ref = REFERENCE_TIMES[(T, 0.99)]
        
        print(f"time for 1% trafo: {t01} s, reference: {t01ref}")
        print(f"time for 99% trafo: {t99} s, reference: {t99ref}")
        
        # subplot: temperature over time
        ax1.plot(time, sdv[:,1], label=f"T={T}")
        ax1.set_xlabel(r"time $t$")
        ax1.set_ylabel(r"bainite fraction $\beta_B$")
        ax1.grid()
    
    plt.legend()
    png_path = plot_path / "results_homogeneous_isothermal.png"
    plt.savefig(png_path, dpi=600)
