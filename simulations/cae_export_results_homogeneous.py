#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Export selected results from the homogeneous state example.

"""

import argparse
import numpy as np

from abaqus import session
from abaqusConstants import CENTROID

FIELD_NAMES = ["TEMP", "SDV1", "SDV2", "SDV3", "SDV4", "SDV5",
               "SDV6", "SDV7", "SDV8", "SDV9", "SDV10", "E"]

def export_results(job_name):

    odb_path = job_name + ".odb"
    npz_path = job_name + ".npz"
    
    odb = session.openOdb(odb_path)
    
    step = odb.steps.values()[-1]
    
    field_data = {}
    
    
    for field_name in FIELD_NAMES:
        times = []
        
        for frame in step.frames:
            times.append(frame.frameValue)
            
        times = np.array(times)
    
    
    for field_name in FIELD_NAMES:
        data = []
        
        for frame in step.frames:
            fieldOut = frame.fieldOutputs[field_name].getSubset(position=CENTROID)
            vals = fieldOut.values
            assert len(vals) == 1
            val = vals[0]
            
            try:
                data.append(val.dataDouble)
            except:
                data.append(val.data)
            
        field_data[field_name] = np.array(data)
    
    sdv_field_names = [fname for fname in FIELD_NAMES if fname.startswith("SDV")]
    sdv_vec = np.vstack((field_data[fname] for fname in sdv_field_names))    
    
    np.savez(npz_path,
             times=times,
             sdvs=sdv_vec,
             temps=field_data["TEMP"],
             strain=field_data["E"])
    
        
def parse_arguments():
    """Get a dict of parsed argument values."""
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--job_name",
                        type=str,
                        help="Name of the job without file extensions",
                        )
    
    args, _ = parser.parse_known_args()
    
    return vars(args)

if __name__ == "__main__":
    args = parse_arguments()
        
    export_results(**args)

    