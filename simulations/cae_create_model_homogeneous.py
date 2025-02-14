"""
Create the model for all homogeneous simulations.

Different cases and conditions can be selected through the parsed
arguments, e.g. isothermal or cooling etc.

This script file is supposed to be executed in Abaqus Python.

"""
# standard library imports
import os
import argparse

# abaqus imports
from abaqus import mdb, Mdb
from caeModules import mesh
from abaqusConstants import ( C3D8T, C3D6T, CENTROIDAL,
                              CONSTANT_THROUGH_THICKNESS, DEFORMABLE_BODY,
                              EMBEDDED_COEFF, FROM_SECTION, FULL, HEX,
                              MIDDLE_SURFACE, ON,OFF, STANDARD, STEP,
                              THREE_D, UNIFORM, UNKNOWN_TET, UNSET,
                             )

# imports from utility module
from common import create_material, add_keywords_intial_condition_sdvs

# model constants - these are not changed between simulations, but can be used
# to rename things or slightly change the model for other applications
MODEL_NAME  = 'homogeneous problem'

MATERIAL_NAME = 'steel'
SECTION_NAME = 'steel'

PART_NAME   = 'block'
INSTANCE_NAME = 'block'

STEP_NAME = 'quenching'

STEP_MAX_TEMPERATURE_INC = 2
STEP_MIN_TIME_INC = 0.05
STEP_INIT_TIME_INC = 100


END_TEMPERATURE = 60.0

BOLTZMANN_CONSTANT = 5.67E-14
ZERO_TEMP = -273.15

MASS_DENSITY = 0.00785

# the final job name will be constructed from the prefix, the material, and the
# quench temperature
JOB_NAME_PREFIX = "homogeneous"

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


def create_block_part(model,
                      size=10.0,
                      mesh_size=10.0,
                      ):
    """
    Create the part object for the cube specimen.

    Unit for all lengths must be consistent with the simulation units.

    Parameters
    ----------
    model : Model
        The Abaqus model the part will be created in.
    size : float
        Edge length of the cube.
    mesh_size : float, optional
        Prescribed element size. Default: 3.

    Returns
    -------
    Part
        The new part object

    """
    sketch = model.ConstrainedSketch(name='__profile__', sheetSize=200.0)

    sketch.rectangle(point1=(0, 0), point2=(size, size))

    block_part = model.Part(name=PART_NAME,
                            dimensionality=THREE_D,
                            type=DEFORMABLE_BODY)

    block_part.BaseSolidExtrude(sketch=sketch, depth=size)

    block_part.Set("all cells", cells=block_part.cells)
    block_part.Surface("whole surface", side1Faces=block_part.faces)

    left_face = block_part.faces.findAt(((0, size/2., size/2.),))
    block_part.Set("left face", faces=left_face)

    lower_face = block_part.faces.findAt(((size/2., 0., size/2.),))
    block_part.Set("lower face", faces=lower_face)

    back_face = block_part.faces.findAt(((size/2., size/2., 0),))
    block_part.Set("back face", faces=back_face)

    block_part.Set("all elements", elements=block_part.elements)

    # mesh creation
    block_part.seedPart(size=mesh_size,
                        deviationFactor=0.1,
                        minSizeFactor=0.1)

    block_part.setMeshControls(regions=block_part.cells,
                               elemShape=HEX)


    elemType1 = mesh.ElemType(elemCode=C3D8T, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=C3D6T, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=UNKNOWN_TET, elemLibrary=STANDARD)

    block_part.setElementType(regions=block_part.sets["all cells"],
                              elemTypes=(elemType1,
                                         elemType2,
                                         elemType3)
                              )

    block_part.generateMesh()

    return block_part


def create_model(material_flag,
                 initial_temperature,
                 htc,
                 step_duration,
                 max_timestep,
                 isothermal=False,
                 disable_martensite=False,
                 disable_bainite=False,
                 fix_temperature=False,
                 fix_displacements=False,
                 ):
    """
    Create the model for the homogeneous example and write job data.

    Parameters
    ----------
    material_flag : int
        0 for 100Cr6, 1 for 100CrMnSi6-4.
    initial_temperature : float
        Initial temperature.
    htc : float
        Heat transfer coefficient.
    step_duration : float
        Length of the step in seconds.
    max_timestep : float
        Maximum time step in seconds.
    isothermal : bool, optional
        If isothermal is True, the temperature is prescribed to be constant at
        initial_temperature. Default: False.
    disable_martensite : bool, optional
        Whether to disable martensite evolution. Default: False.
    disable_bainite : bool, optional
        Whether to disable bainite evolution. Default: False.
    fix_temperature : bool, optional
        Prescribe the temperature path by dirichlet condition.
        Default: False.
    fix_displacmenets : bool, optional
        Prescribe all displacements as 0 by dirichlet conditions.
        Default: False.

    Returns
    -------
    job
        The Abaqus job for the requested model.

    """
    Mdb()
    mdb.models.changeKey(fromName='Model-1', toName=MODEL_NAME)
    model = mdb.models[MODEL_NAME]

    model.setValues(stefanBoltzmann=BOLTZMANN_CONSTANT,
                    absoluteZero=ZERO_TEMP)

    block_part = create_block_part(model, size=1.0, mesh_size=1.0)

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

    block_part.SectionAssignment(block_part.sets['all cells'],
                                 sectionName=SECTION_NAME,
                                 offset=0.0,
                                 offsetType=MIDDLE_SURFACE,
                                 offsetField='',
                                 thicknessAssignment=FROM_SECTION)

    # assembly
    block_instance = model.rootAssembly.Instance(name=INSTANCE_NAME,
                                                 part=block_part,
                                                 dependent=ON)

    # convert step duration to ms
    step_duration *= 1000
    max_timestep *= 1000

    step_max_num_inc = int(step_duration / STEP_MIN_TIME_INC)

    model.CoupledTempDisplacementStep(name=STEP_NAME,
                                      previous='Initial',
                                      timePeriod=step_duration,
                                      initialInc=STEP_INIT_TIME_INC,
                                      minInc=STEP_MIN_TIME_INC,
                                      maxInc=max_timestep,
                                      deltmx=STEP_MAX_TEMPERATURE_INC,
                                      maxNumInc=step_max_num_inc)

    # boundary conditions

    if fix_displacements:
        # fix all displacements to disregard tangents
        model.DisplacementBC(name='all displacements fixed',
                             createStepName='Initial',
                             region=block_instance.sets["all cells"],
                             u1=0.0, u2=0.0, u3=0.0,
                             amplitude=UNSET,
                             distributionType=UNIFORM,
                             fieldName='',
                             localCsys=None)
    else:
        # or fix only 3 sides to allow free expansion
        model.DisplacementBC(name='x displacements fixed',
                             createStepName='Initial',
                             region=block_instance.sets["left face"],
                             u1=0.0, u2=UNSET, u3=UNSET,
                             amplitude=UNSET,
                             distributionType=UNIFORM,
                             fieldName='',
                             localCsys=None)

        model.DisplacementBC(name='y displacements fixed',
                             createStepName='Initial',
                             region=block_instance.sets["lower face"],
                             u1=UNSET, u2=0.0, u3=UNSET,
                             amplitude=UNSET,
                             distributionType=UNIFORM,
                             fieldName='',
                             localCsys=None)

        model.DisplacementBC(name='z displacements fixed',
                             createStepName='Initial',
                             region=block_instance.sets["back face"],
                             u1=UNSET, u2=UNSET, u3=0.0,
                             amplitude=UNSET,
                             distributionType=UNIFORM,
                             fieldName='',
                             localCsys=None)

    model.Temperature(name='initial temperature',
                      createStepName='Initial',
                      region=block_instance.sets['all cells'],
                      distributionType=UNIFORM,
                      crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
                      magnitudes=(initial_temperature, ))

    if fix_temperature:
        # prescribe uniform temperature

        if isothermal:
            # constant temperature
            model.TemperatureBC(name='prescribed temperature',
                                createStepName=STEP_NAME,
                                region=block_instance.sets["all cells"],
                                distributionType=UNIFORM,
                                fieldName='',
                                magnitude=initial_temperature)
        else:
            # exponential temperature decay
            area = 6.0
            volume = 1.0
            rho = MASS_DENSITY
            cp = 400.0

            a = htc*area/volume/rho/cp

            amp = model.DecayAmplitude(name='exponential temperature decay',
                                       timeSpan=STEP,
                                       initial=END_TEMPERATURE,
                                       maximum=initial_temperature-END_TEMPERATURE,
                                       start=0.0,
                                       decayTime=1/a)

            model.TemperatureBC(name='prescribed temperature',
                                createStepName=STEP_NAME,
                                region=block_instance.sets["all cells"],
                                fixed=OFF,
                                distributionType=UNIFORM,
                                fieldName='',
                                magnitude=1.0,
                                amplitude=amp.name)

    else:
        # prescribe heat flux to constant ambient temperature

        if isothermal:
            sink_temperature = initial_temperature
        else:
            sink_temperature = END_TEMPERATURE

        model.FilmCondition(name='Constant ambient temperature',
                            createStepName='quenching',
                            surface=block_instance.surfaces['whole surface'],
                            definition=EMBEDDED_COEFF,
                            filmCoeff=htc,
                            filmCoeffAmplitude='',
                            sinkTemperature=sink_temperature,
                            sinkAmplitude='',
                            sinkDistributionType=UNIFORM,
                            sinkFieldName='')



    # Field output request
    model.fieldOutputRequests['F-Output-1'].setValues(variables=(
        'S', 'E', 'U', 'RF', 'NT', 'HFL', 'RFL','SDV', 'TEMP'),
        timeInterval=500.0
        )

    model.FieldOutputRequest(name='F-Output-2',
        createStepName=STEP_NAME,
        variables=('S', 'E', 'SDV', 'TEMP'),
        numIntervals=100,
        position=CENTROIDAL)


    # insert the keyword for initial sdvs from user routine
    add_keywords_intial_condition_sdvs(model,
                                       INSTANCE_NAME,
                                       "all cells",
                                       STEP_NAME,
                                       INITIAL_SDVS)

    job_name = get_job_name(JOB_NAME_PREFIX,
                            material_flag,
                            initial_temperature,
                            htc,
                            isothermal,
                            disable_martensite,
                            disable_bainite,
                            fix_temperature)

    print("Job name is " + job_name)

    job = mdb.Job(name=job_name,
                  model=MODEL_NAME,
                  userSubroutine=UMAT_PATH,
                  nodalOutputPrecision=FULL,
                  )

    job.writeInput(consistencyChecking=OFF)

    return job


def get_job_name(prefix,
                 material_flag,
                 initial_temperature,
                 htc,
                 isothermal,
                 disable_martensite,
                 disable_bainite,
                 fix_temperature):
    """
    Construct the name of the job based on all switches.

    Parameters
    ----------
    prefix : str
        Prefix for the job name.
    material_flag : int
        0 for 100Cr6, 1 for 100CrMnSi6-4.
    initial_temperature : float
        Initial temperature.
    htc : float
        Heat transfer coefficient.
    isothermal : bool, optional
        If isothermal is True, the temperature is prescribed to be constant at
        initial_temperature. Default: False.
    disable_martensite : bool
        Whether to disable martensite transformation.
    disable_bainite : bool
        Whether to disable bainite transformation.
    fix_temperature
        Whether the temperature is prescribed as a Dirichlet condition.

    Returns
    -------
    str
        The job name.

    """
    job_name = prefix
    job_name += '_100CrMnSi64' if material_flag==1 else '_100Cr6'

    job_name += "_T0_{T0:03.0f}".format(T0=initial_temperature)

    job_name += "_htc_{htc:0.4f}".format(htc=1e5*htc).replace(".", "-")

    if isothermal:
        job_name += "_isothermal"

    if fix_temperature:
        job_name += "_fixed_temp"

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

    parser.add_argument("--fix_temperature",
                        action="store_true",
                        default=False,
                        help="Prescribe fixed temperature path",
                        )

    parser.add_argument("--fix_displacements",
                        action="store_true",
                        help="Prescribe all displacements as 0",
                        )

    parser.add_argument("--dir",
                        type=str,
                        help="(Relative) subdirectory execute the script in",
                        default=".",
                        )

    parser.add_argument("--initial_temperature",
                        type=float,
                        default=300,
                        help="Initial temperature",
                        )

    parser.add_argument("--htc",
                        type=float,
                        default=1e-5,
                        help="Heat transfer coefficient",
                        )

    parser.add_argument("--isothermal",
                        action="store_true",
                        help="Fix the temperature at the initial temperature"
                        )

    parser.add_argument("--step_duration",
                        type=float,
                        help="Step duration in seconds",
                        default=500.0,
                        )

    parser.add_argument("--max_timestep",
                        type=float,
                        help="Maximum time step in seconds",
                        default=5,
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
