#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Export selected results from the bearing race example.

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
    
    steps = odb.steps.values()    
    
    times = []
    for step in steps:
        for frame in step.frames:
            times.append(step.totalTime + frame.frameValue)
            
        
    times = np.array(times)
    
    field_data_e1 = read_element_centroid_data(steps, FIELD_NAMES, (155.0, -12.5, 0.0))
    field_data_e2 = read_element_centroid_data(steps, FIELD_NAMES, (163.75, -75.0, 0.0))
    field_data_e3 = read_element_centroid_data(steps, FIELD_NAMES, (177.5, -137.5, 0.0))
    
    field_data_n1 = read_node_point_data(odb, steps, FIELD_NAMES, "P1")
    field_data_n2 = read_node_point_data(odb, steps, FIELD_NAMES, "P2")
    field_data_n3 = read_node_point_data(odb, steps, FIELD_NAMES, "P3")
    
    
    sdv_field_names = [fname for fname in FIELD_NAMES if fname.startswith("SDV")]
    
    sdvvec_n1 = np.vstack((field_data_n1[fname] for fname in sdv_field_names))
    sdvvec_n2 = np.vstack((field_data_n2[fname] for fname in sdv_field_names))
    sdvvec_n3 = np.vstack((field_data_n3[fname] for fname in sdv_field_names))
    
    sdvvec_e1 = np.vstack((field_data_e1[fname] for fname in sdv_field_names))
    sdvvec_e2 = np.vstack((field_data_e2[fname] for fname in sdv_field_names))
    sdvvec_e3 = np.vstack((field_data_e3[fname] for fname in sdv_field_names))
    
    np.savez(npz_path,
             times=times,
             
             sdvs_n1=sdvvec_n1,
             temps_n1=field_data_n1["TEMP"],
             
             sdvs_n2=sdvvec_n2,
             temps_n2=field_data_n2["TEMP"],
             
             sdvs_n3=sdvvec_n3,
             temps_n3=field_data_n3["TEMP"],
             
             sdvs_e1=sdvvec_e1,
             temps_e1=field_data_e1["TEMP"],
             
             sdvs_e2=sdvvec_e2,
             temps_e2=field_data_e2["TEMP"],
             
             sdvs_e3=sdvvec_e3,
             temps_e3=field_data_e3["TEMP"],
             )
    
        
def parse_arguments():
    """Get a dict of parsed argument values."""
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--job_name",
                        type=str,
                        help="Name of the job without file extensions",
                        )
    
    args, _ = parser.parse_known_args()
    
    return vars(args)


def read_node_point_data(odb, steps, field_names, set_name):
    field_data = {}

    nset = odb.rootAssembly.instances["bearing race"].nodeSets[set_name]
    
    for field_name in field_names:
        data = []
        for step in steps:
            for frame in step.frames:
                fieldOut = frame.fieldOutputs[field_name].getSubset(region=nset)
                vals = fieldOut.values
                assert len(vals) == 1
                val = vals[0]
                
                try:
                    data.append(val.dataDouble)
                except:
                    data.append(val.data)

        field_data[field_name] = np.array(data)
    
    return field_data


def read_element_centroid_data(steps, field_names, coord):
    frame = steps[0].frames[0]
    fieldOut = frame.fieldOutputs["COORD"].getSubset(position=CENTROID)
    values = fieldOut.values
    
    try:
        centroid_coords = np.array([v.dataDouble for v in values])
    except:
        centroid_coords = np.array([v.data for v in values])
        
    coord = coord[0:centroid_coords.shape[1]]
    centroid_distances = np.linalg.norm(centroid_coords-coord, axis=1)
    imin = np.argmin(centroid_distances)
    
    enum = values[imin].elementLabel
    inum = values[imin].integrationPoint
    
    field_data = {}
    
    for field_name in field_names:
        data = []
        for step in steps:
            for frame in step.frames:
                fieldOut = frame.fieldOutputs[field_name]
                vals = [v for v in fieldOut.values if v.elementLabel == enum and v.integrationPoint == inum]
                assert len(vals) == 1
                val = vals[0]
                
                try:
                    data.append(val.dataDouble)
                except:
                    data.append(val.data)

        field_data[field_name] = np.array(data)
    
    return field_data


if __name__ == "__main__":
    args = parse_arguments()
        
    export_results(**args)
        
        