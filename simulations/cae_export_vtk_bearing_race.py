"""
Export field results from the bearing race example as vtk files.

The variables specified in FIELD_NAMES are exported either at the
node points or element centroids, as specified.

The open source software Paraqus is used to perform the actual export.

This script file is supposed to be executed in Abaqus Python.

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
    """
    Perform the actual export.
    
    Parameters
    ----------
    job_name : str
		Identifier of the Abaqus job.
	output_dir : str
		Directory for the vtk files.
		
	Returns
	-------
	None
	
    """
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
        
        
