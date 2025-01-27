"""
Functions that are used in both examples.

"""

from abaqusConstants import THERMOMECHANICAL

NUMBER_SDVS = 10

def create_material(model,
                    name,
                    density,
                    material_flag,
                    disable_martensite, 
                    disable_bainite):
    """
    Create the material calling the user subroutines UMAT and UMATHT.

    Parameters
    ----------
    model
        The Abaqus mode the material is created in.
    name : str
        The name of the new material.
    density : float
        Mass density for the material (must be set in Abaqus for
        thermomechanically coupled problems).
    material_flag : int
        0 for 100Cr6, 1 for 100CrMnSi6-4.
    disable_martensite : bool
        Whether to disable the evolution of martensite.
    disable_bainite : bool
        Whether to disable the evolution of bainite.

    Returns
    -------
    The new material.

    """
    # flags are 0 when evolution is active, 1 when evolution is suppressed
    martensite_flag = int(disable_martensite)
    bainite_flag = int(disable_bainite)
    
    print("Material: {mat}".format(mat=material_flag))
    
    # create the model in a fresh database
    
    
    
    
    # material definition
    material_steel = model.Material(name=name)
    
    # Density must be specified for thermal user material
    material_steel.Density(table=((density, ), ))
    
    material_steel.Depvar(n=NUMBER_SDVS)
    
    
    material_props = (density,
                      material_flag, 
                      int(martensite_flag),
                      int(bainite_flag),
                      )
    
    material_steel.UserMaterial(
        type=THERMOMECHANICAL,
        mechanicalConstants=material_props,
        thermalConstants=material_props,
        )
    
    return material_steel


def format_name(name):
        """Add quotation marks if name contains whitespace."""
        if ' ' in name:
            return '"{name}"'.format(name=name)
        else:
            return name
        

def add_keywords_intial_condition_sdvs(model, instance_name, set_name, step_name, sdvs):
    """
    Add a keyword block to specify homogeneous initial values for sdvs.

    Parameters
    ----------
    model
        The Abaqus model the keyword block is updated in.
    instance_name : str
        Name of the instance for which the initial condition is set.
    set_name : str
        Name of a set in the instance for which the initial conditions are set.
    step_name : str
        Set the initial condition for this step.
    sdvs : Sequence[float]
        Values for the initial sdvs.

    Raises
    ------
    Exception
        If the step is not found.

    Returns
    -------
    None.

    """
    modelkwb = model.keywordBlock
    modelkwb.synchVersions(storeNodesAndElements=False)
    
    # construct the keyword block for the sdvs
    kwds = '*INITIAL CONDITIONS,TYPE=SOLUTION \n{iname}.{sname}, '.format(
        iname=format_name(instance_name),
        sname=format_name(set_name)
        )
                                                                            
    for i, sdvi in enumerate(sdvs):
        kwds += "{val}, ".format(val=sdvi)
        if (i+1) % 7 == 0:
            kwds += "\n"
            
    # delete the last comma and whitespace
    kwds = kwds[0:-2]
            
    # find the correct position for the keyword block
    line_num = 0
    for n, line in enumerate(modelkwb.sieBlocks):
        if "STEP: {step}".format(step=step_name) in line:
            line_num = n
            break
        
    # insert the keyword block
    if line_num: 
        modelkwb.insert(position=line_num-1, text=kwds)
    else:
        e = ("Error: step... was not found",
              "in the Model KeywordBlock.")
        raise Exception(" ".join(e))
        
