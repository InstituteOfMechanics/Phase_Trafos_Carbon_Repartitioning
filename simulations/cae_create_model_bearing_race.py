# -*- coding: utf-8 -*-
"""
Create the example model for the inner race of a tapered bearing.

TODO: LICENSE INFO ETC

The quenching process is simulated in two steps. This choice is made to reduce
data size, since a reduced output frequency is prescribed in the second step,
where the evolution of the microstructure is slowed down considerably.

"""
# standard library imports
import os
import argparse

# abaqus imports
from abaqus import mdb, Mdb
from caeModules import mesh
from abaqusConstants import ( AVERAGED_AT_NODES, AXISYMMETRIC, CAX3T, CAX4T,
                              CENTROIDAL, CONSTANT_THROUGH_THICKNESS,
                              DEFORMABLE_BODY, FIXED, FROM_SECTION, FULL,
                              MEDIAL_AXIS, MIDDLE_SURFACE, ON, OFF,
                              PROPERTY_REF, QUAD, STANDARD, UNIFORM, UNSET,
                             )

# imports from utility module
from common import create_material, add_keywords_intial_condition_sdvs

# model constants - these are not changed between simulations, but can be used
# to rename things or slightly change the model for other applications
MODEL_NAME 	= 'bearing race quenching'

MATERIAL_NAME = 'steel'
SECTION_NAME = 'steel'

PART_NAME   = 'bearing race'
INSTANCE_NAME = 'bearing race'

STEP1_NAME = 'quenching'
STEP1_DURATION = 120e3 # ms
STEP1_MAX_TEMPERATURE_INC = 5
STEP1_MIN_TIME_INC = 0.5
STEP1_MAX_TIME_INC = 1000
STEP1_INIT_TIME_INC = 500
STEP1_MAX_NUM_INC = int(STEP1_DURATION/STEP1_MIN_TIME_INC)+1
STEP1_FIELDOUT_INTERVAL = 200.0

STEP2_NAME = 'equalizing'
STEP2_DURATION = 540e3 # ms
STEP2_MAX_TEMPERATURE_INC = 5
STEP2_TIME_INC = 1000
STEP2_MAX_NUM_INC = int(STEP2_DURATION/STEP2_TIME_INC)+1
STEP2_FIELDOUT_INTERVAL = 4000.0

INITIAL_TEMPERATURE = 850.0

BOLTZMANN_CONSTANT = 5.67E-14
ZERO_TEMP = -273.15

MASS_DENSITY = 0.00785

# the final job name will be constructed from the prefix, the material, and the
# quench temperature
JOB_NAME_PREFIX = "bearing_race_quenching"

# path to the Fortran UMAT main file
UMAT_PATH = os.path.join(os.getcwd(), "UMAT", "umat_main.f")

# initial values for state variables
INITIAL_SDVS = [0.0, # martensite volume fraction
                0.0, # bainite volume fraction
                1.0, # austenite volume fraction
                0.01, # austenite phase carbon content
                0.0, # thermal strain
                0.0, # transformation strain
                0.0, # derivative martensite volume fraction w.r.t. temperature
                0.0, # derivative bainite volume fraction w.r.t. temperature
                0.0, # derivative austenite phase carbon fraction w.r.t. temperature
                0.0, # factor for incomplete bainite trafo
                ]

def create_bearing_part(model,
                        r1=150.,
                        r2=160,
                        r3=195.,
                        r4=205.,
                        l1=25.,
                        l2=125.,
                        l3=150.,
                        mesh_size=3.
                        ):
    """
    Create the part object for the bearing race.

    Sketch of the parametrized race geometry as shown in the paper:

        
              r4
    |----------------->|

            r3
    |-------------->|
    
          r2
    |----------->|    
    
        r1
    |------->|
              ___       
    |        |   |           ^    ^    ^
             |   |           | l1 |    |
    |        |   \           v    | l2 |
             |    \               |    |
    |        |     \              |    | l3
             |      \__           v    |
    |        |         |               |
             |         |               |
    |        |_________|               v
    
    
    Unit for all lengths must be consistent with the simulation units.
    
    Parameters
    ----------
    model : Model
        The Abaqus model the part will be created in.
    r1 : float, optional
        Inner radius. Default: 150.
    r2 : float, optional
        Smallest outer radius. Default: 160.
    r3 : float, optional
        Outer radius at the end of the angled section. Default: 195.
    r4 : float, optional
        Largest outer radius. Default: 205.
    l1 : float, optional
        Length from top end to begin of angles section. Default: 25.
    l2 : float, optional
        Length from top end to end of angled section. Default: 125.
    l3 : float, optional
        Total length. Default: 150.
    mesh_size : float, optional
        Prescribed element size. Default: 3.
        
    Returns
    -------
    Part
        The new part object

    """
    sketch = model.ConstrainedSketch(name='__profile__', sheetSize=200.0)
    
    sketch.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    
    sketch.Line(point1=(r1,   0), point2=(r2,   0))
    sketch.Line(point1=(r2,   0), point2=(r2, -l1))
    sketch.Line(point1=(r2, -l1), point2=(r3, -l2))
    sketch.Line(point1=(r3, -l2), point2=(r4, -l2))
    sketch.Line(point1=(r4, -l2), point2=(r4, -l3))
    sketch.Line(point1=(r4, -l3), point2=(r1, -l3))
    sketch.Line(point1=(r1, -l3), point2=(r1,   0))
    
    bearing_part = model.Part(name=PART_NAME,
                              dimensionality=AXISYMMETRIC, 
                              type=DEFORMABLE_BODY)
    shell = bearing_part.BaseShell(sketch=sketch)
    
    
    # partitioning    
    sketch = model.ConstrainedSketch(name='__profile__', sheetSize=1.)

    # horizontal lines 
    sketch.Line(point1=(0, -l1/2.), point2=(r4, -l1/2.))
    sketch.Line(point1=(0, -l1), point2=(r4, -l1))
    sketch.Line(point1=(0, -l1-(l2-l1)/2.), point2=(r4, -l1-(l2-l1)/2.))
    sketch.Line(point1=(0, -l2), point2=(r4, -l2))
    sketch.Line(point1=(0, -l2-(l3-l2)/2.), point2=(r4, -l2-(l3-l2)/2.))
    
    # vertical lines
    sketch.Line(point1=((r1+r2)/2., 0), point2=((r1+r2)/2., -l1))
    sketch.Line(point1=((r1+r2)/2., -l1), point2=((r1+r3)/2., -l2))
    sketch.Line(point1=((r1+r4)/2., -l2), point2=((r1+r4)/2., -l3))
    
    sketch.Line(point1=(r3, 0), point2=(r3, -l3))
    
    bearing_part.PartitionFaceBySketch(faces=bearing_part.faces, sketch=sketch)

    del sketch
    
    
    # mesh creation
    bearing_part.seedPart(size=mesh_size,
                          deviationFactor=0.1,
                          minSizeFactor=0.1)
    
    
    
    
    bearing_part.setMeshControls(regions=bearing_part.faces,
                                 elemShape=QUAD,
                                 algorithm=MEDIAL_AXIS)
    elem_type_1 = mesh.ElemType(elemCode=CAX4T, elemLibrary=STANDARD)
    elem_type_2 = mesh.ElemType(elemCode=CAX3T, elemLibrary=STANDARD)
    
    bearing_part.setElementType(regions=(bearing_part.faces, ),
                                elemTypes=(elem_type_1, 
                                           elem_type_2)
                                )
    bearing_part.generateMesh()
    
    # sets and surfaces
    
    bearing_part.Set("all faces", faces=bearing_part.faces)
    outer_edges = bearing_part.getFeatureEdges(shell.name)
    bearing_part.Surface("whole surface", side1Edges=outer_edges)
    
    lower_left_vertex = bearing_part.vertices.findAt(((r1, -l3, 0),))
    bearing_part.Set("lower left vertex", vertices=lower_left_vertex)
    
    bearing_part.Set("all elements", elements=bearing_part.elements)


    p1 = ((r1+r2)/2., -l1/2., 0)
    p2 = ((2*r1+r2+r3)/4., -l1 - (l2-l1)/2., 0)
    p3 = ((r1+r4)/2., -l2 - (l3-l2)/2., 0)
    
    
    p1_vertex = bearing_part.vertices.findAt(p1)
    bearing_part.Set("P1", nodes=p1_vertex.getNodes())
    
    p2_vertex = bearing_part.vertices.findAt(p2)
    bearing_part.Set("P2", nodes=p2_vertex.getNodes())
    
    p3_vertex = bearing_part.vertices.findAt(p3)
    bearing_part.Set("P3", nodes=p3_vertex.getNodes())
    

    return bearing_part


def create_model(material_flag,
                 quench_temperature,
                 disable_martensite=False,
                 disable_bainite=False,
                 mesh_size=3):
    """
    Create the model for the quenching example and write job data.

    Parameters
    ----------
    material_flag : int
        0 for 100Cr6, 1 for 100CrMnSi6-4.
    quench_temperature : int
        Temperature of the quenching bath in degree celsius. Valid choices
        are 20 or 60.
    disable_martensite : bool, optional
        Whether to disable martensite evolution. Default: False.
    disable_bainite : bool, optional
        Whether to disable bainite evolution. Default: False.
    mesh_size : int, optional
        Typical element size.

    Returns
    -------
    job
        The Abaqus job for the requested model.

    """    
    # create the model in a fresh database
    Mdb()
    mdb.models.changeKey(fromName='Model-1', toName=MODEL_NAME)
    model = mdb.models[MODEL_NAME]
    
    model.setValues(stefanBoltzmann=BOLTZMANN_CONSTANT,
                    absoluteZero=ZERO_TEMP)
    
    bearing_part = create_bearing_part(model,
									   mesh_size=mesh_size,
									   )
    
    # material definition
    material_steel = create_material(model,
                                     MATERIAL_NAME,
                                     MASS_DENSITY,
                                     material_flag,
                                     disable_martensite,
                                     disable_bainite)
    
    # section definition and assignment
    model.HomogeneousSolidSection(name=SECTION_NAME,
                                  material=material_steel.name,
                                  thickness=None)
    
    bearing_part.SectionAssignment(bearing_part.sets['all faces'],
                                   sectionName=SECTION_NAME,
                                   offset=0.0, 
                                   offsetType=MIDDLE_SURFACE, 
                                   offsetField='', 
                                   thicknessAssignment=FROM_SECTION)
    
    # assembly
    bearing_instance = model.rootAssembly.Instance(name=INSTANCE_NAME,
                                                 part=bearing_part,
                                                 dependent=ON)
    
    # step 
    model.CoupledTempDisplacementStep(name=STEP1_NAME,
                                      previous='Initial',
                                      timePeriod=STEP1_DURATION,
                                      initialInc=STEP1_INIT_TIME_INC,
                                      minInc=STEP1_MIN_TIME_INC,
                                      maxInc=STEP1_MAX_TIME_INC,
                                      deltmx=STEP1_MAX_TEMPERATURE_INC,
                                      maxNumInc=STEP1_MAX_NUM_INC)
    
    model.CoupledTempDisplacementStep(name=STEP2_NAME,
                                      previous=STEP1_NAME,
                                      timePeriod=STEP2_DURATION,
                                      timeIncrementationMethod=FIXED,
                                      initialInc=STEP2_TIME_INC,
                                      maxNumInc=STEP2_MAX_NUM_INC)
    
    # boundary conditions
    
    # or fix only 3 sides to allow free expansion
    model.DisplacementBC(name='y displacement fixed',
                         createStepName='Initial',
                         region=bearing_instance.sets["lower left vertex"],
                         u1=UNSET, u2=0.0, u3=UNSET, 
                         amplitude=UNSET,
                         distributionType=UNIFORM,
                         fieldName='',
                         localCsys=None)
        
    model.Temperature(name='initial temperature', 
                      createStepName='Initial',
                      region=bearing_instance.sets['all faces'],
                      distributionType=UNIFORM, 
                      crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
                      magnitudes=(INITIAL_TEMPERATURE, ))
        
    
    if quench_temperature == 20:
        convection_property = model.FilmConditionProp(
            name='quenching water 20deg',
            temperatureDependency=ON,
            dependencies=0,
            property=((0.0043500,   0.0),
                      (0.0082071, 200.0),
                      (0.0119617, 400.0),
                      (0.0134917, 430.0),
                      (0.0125000, 500.0), 
                      (0.0102062, 560.0),
                      (0.0077930, 600.0),
                      (0.0025070, 700.0),
                      (0.0004371, 800.0), 
                      (0.0001353, 900.0))
            )
    
    elif quench_temperature == 60:
        convection_property = model.FilmConditionProp(
            name='quenching water 60deg',
            temperatureDependency=ON,
            dependencies=0,
            property=((0.0001353,   0.0), 
                      (0.0020292, 200.0),
                      (0.0028409, 400.0),
                      (0.0032912, 445.0),
                      (0.0034220, 500.0),
                      (0.0026109, 570.0),
                      (0.0021570, 600.0),
                      (0.0002303, 800.0), 
                      (0.0001352, 900.0))
            )
    else:
        raise NotImplementedError("Quench temperature must be 20 or 60.")
   
    model.FilmCondition(name='Convection', 
                        createStepName=STEP1_NAME,
                        surface=bearing_instance.surfaces['whole surface'],
                        definition=PROPERTY_REF, 
                        interactionProperty=convection_property.name,
                        sinkTemperature=quench_temperature, 
                        sinkDistributionType=UNIFORM, 
                        sinkFieldName='')
    
 
    
    del model.historyOutputRequests['H-Output-1']
    del model.fieldOutputRequests['F-Output-1']
    
    
    fo_nodes = model.FieldOutputRequest(name='Node Output 1', 
        createStepName=STEP1_NAME,
        variables=('U', 'NT'), 
        timeInterval=STEP1_FIELDOUT_INTERVAL,
        timeMarks=OFF)
    
    fo_nodes.setValuesInStep(stepName='equalizing',
                             timeInterval=STEP2_FIELDOUT_INTERVAL)
    
    fo_centroids = model.FieldOutputRequest(name='Element Center Output 1', 
        createStepName=STEP1_NAME,
        variables=('S', 'E', 'SDV', 'TEMP', 'COORD'),
        timeInterval=STEP1_FIELDOUT_INTERVAL,
        timeMarks=OFF,
        position=CENTROIDAL)
    
    fo_centroids.setValuesInStep(stepName='equalizing',
                             timeInterval=STEP2_FIELDOUT_INTERVAL)
    
    fo_plots = model.FieldOutputRequest(name='Plot Output 1',
        createStepName=STEP1_NAME,
        variables=('S', 'E', 'TEMP', 'SDV'), 
        timeInterval=STEP1_FIELDOUT_INTERVAL,
        timeMarks=OFF, 
        position=AVERAGED_AT_NODES)
    
    fo_plots.setValuesInStep(stepName='equalizing',
                             timeInterval=STEP2_FIELDOUT_INTERVAL)

    
    add_keywords_intial_condition_sdvs(model,
                                       INSTANCE_NAME,
                                       "all faces",
                                       STEP1_NAME,
                                       INITIAL_SDVS)
        
   
    
    
    # job definition
    job_name = get_job_name(JOB_NAME_PREFIX,
                            material_flag,
                            quench_temperature,
                            disable_martensite,
                            disable_bainite)
    
    job = mdb.Job(name=job_name,
                  model=MODEL_NAME,
                  userSubroutine=UMAT_PATH,
                  nodalOutputPrecision=FULL,
                  )    
        
    job.writeInput(consistencyChecking=OFF)

    return job


def get_job_name(prefix,
                 material_flag, 
                 quench_temperature,
                 disable_martensite,
                 disable_bainite):
    """
    Construct the name of the job based on all switches.

    Parameters
    ----------
    prefix : str
        Prefix for the job name.
    material_flag : int
        0 for 100Cr6, 1 for 100CrMnSi6-4.
    quench_temperature : int
        Temperature of the quenching bath in degree celsius. Valid choices
        are 20 or 60.
    disable_martensite : bool
        Whether to disable martensite transformation.
    disable_bainite : bool
        Whether to disable bainite transformation.

    Returns
    -------
    str
        The job name.

    """
    job_name = prefix
    job_name += '_100CrMnSi64' if material_flag==1 else '_100Cr6'
    
    job_name += '_Tinf_{Tinf}'.format(Tinf=quench_temperature)
    
    if disable_martensite:
        job_name += "_no_martensite"
        
    if disable_bainite:
        job_name += "_no_bainite"
        
    return job_name


def parse_arguments():
    """Get a dict of parsed argument values."""
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--material_flag",
                        type=int,
                        choices=[0, 1],
                        help="0: 100Cr6, 1: 100CrMnSi6-4",
                        default=1,
                        )
    
    parser.add_argument("--disable_martensite",
                        action="store_true",
                        default=False,
                        help="Diable martensite evolution",
                        )
    
    parser.add_argument("--disable_bainite",
                        action="store_true",
                        default=False,
                        help="Diable bainite evolution",
                        )
    
    parser.add_argument("--dir",
                        type=str,
                        help="(Relative) subdirectory execute the script in",
                        default=".",
                        )
    
    parser.add_argument("--quench_temperature",
                        type=int,
                        choices=[20,60],
                        default=60,
                        help="Temperature of the quenching bath. Must be 20 or 60."
                        )
    
    parser.add_argument("--mesh_size",
                        type=int,
                        default=3,
                        help="Typical mesh size."
                        )
    
    args, _ = parser.parse_known_args()
    
    args = vars(args)
    
    simdir = args.pop("dir")
    
    return simdir, args


if __name__ == "__main__":
    simdir, args = parse_arguments()

    if not os.path.isdir(simdir):
        os.mkdir(simdir)
        
    os.chdir(simdir)
    
    create_model(**args)
