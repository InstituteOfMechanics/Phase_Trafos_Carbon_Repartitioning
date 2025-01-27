#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Export vtk files from the homogeneous state example.

"""

import argparse
import os

from paraqus.abaqus import OdbReader
from paraqus import BinaryWriter

FIELD_REQUESTS = [("U", "nodes", None),
                  ("NT11", "nodes", None),
                  ("S", "elements", None),
                  ("S", "elements", "mises"),
                  ("E", "elements", None),
                  ("SDV1", "elements", None),
                  ("SDV2", "elements", None),
                  ("SDV3", "elements", None),
                  ("SDV4", "elements", None),
                  ("SDV5", "elements", None),
                  ("SDV6", "elements", None),
                  ]

def export_results(job_name, output_dir):
    
    odb_path = job_name + ".odb"
    output_path = os.path.join(os.path.pardir, output_dir, job_name)

    reader = OdbReader(odb_path, job_name)
    
    for field_name, field_position, invariant in FIELD_REQUESTS:
        reader.add_field_export_request(field_name,
                                        field_position=field_position,
                                        invariant=invariant)
        
    step_name = "equalizing"
    
    models = list(reader.read_instances(step_name, -1))
    
    assert len(models) == 1
    
    model = models[0]
    
    writer = BinaryWriter(output_path, clear_output_dir=True)
    
    writer.write(model)
    
        
def parse_arguments():
    """Get a dict of parsed argument values."""
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--job_name",
                        type=str,
                        help="Name of the job without file extensions",
                        )
    
    parser.add_argument("--output_dir",
                        type=str,
                        help="Base directory for vtk output",
                        )
    
    args, _ = parser.parse_known_args()
    
    return vars(args)


if __name__ == "__main__":
    args = parse_arguments()
        
    export_results(**args)
        
        
