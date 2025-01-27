#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Export selected results from the bearing race example.

"""

import argparse
import numpy as np

from abaqus import session
from abaqusConstants import NODAL

INNER_RADIUS = 150

def export_displacements(job_name):

    odb_path = job_name + ".odb"
    npz_path = job_name + "_displacements" + ".npz"
    
    odb = session.openOdb(odb_path)
    
    step = odb.steps.values()[-1]
    frame = step.frames[-1]    
    
    coords = read_node_field(odb, frame, "COORD")
    x, y = coords.T
    
    displacements = read_node_field(odb, frame, "U")
    ux, uy = displacements.T
    
    selector = np.isclose(x, min(x))
    
    y_inner = y[selector]
    ux_inner = ux[selector]
    
    # sort from top to bottom
    sorter = np.flip(np.argsort(y_inner))
    
    y_sorted = y_inner[sorter]
    ux_sorted = ux_inner[sorter]
    
    np.savez(npz_path,
             y=y_sorted,
             ux=ux_sorted,
             )
    
        
def parse_arguments():
    """Get a dict of parsed argument values."""
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--job_name",
                        type=str,
                        help="Name of the job without file extensions",
                        required=True,
                        )
    
    args, _ = parser.parse_known_args()
    
    return vars(args)


def read_node_field(odb, frame, field_name):
    
    fieldOut = frame.fieldOutputs[field_name].getSubset(position=NODAL)
    vals = fieldOut.values

    data = []    
    for val in vals:
        try:
            data.append(val.dataDouble)
        except:
            data.append(val.data)

    data = np.array(data)
    
    return data


if __name__ == "__main__":
    args = parse_arguments()
        
    export_displacements(**args)
        
        