class Options(object):
    def __init__(self, **kwargs):
        for key,value in kwargs.iteritems():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                raise KeyError('Unrecognized configuration option: %s' % key)


class Paths(Options):
    inputDir = "input"
    outputDir = "run/test1"
    logFile = None # write output to stdout


class General(Options):
    fixedBurnedVal = True
    fixedLeftLocation = False
    curvedFlame = False
    twinFlame = False

    unburnedLeft = True # applies to premixed flames
    fuelLeft = True # applies to diffusion flames

    chemistryIntegrator = "qss"
    errorStopCount = 100


class Chemistry(Options):
    mechanismFile = "chem.xml"
    phaseID = "gas"
    transportModel = "Mix"


class Adapchem(Options):
    mechanism = "chem.bin"
    infile = "adapchem.in"
    models = "modelsV2.chem"
    defaultModel = "default.model"
    donemodels = "donemodels"
    restart = "full.rstrt"


class TransportElimination(Options):
    diffusion = False
    convection = False
    atol = 1e-3
    stepInterval = 1


class Grid(Options):
    # High accuracy: vtol = 0.08, dvtol = 0.12, dampConst = 7
    # Medium accuracy: vtol = 0.1, dvtol = 0.25
    # Minimal accuracy: vtol = 0.2, dvtol = 0.4
    vtol = 0.12
    dvtol = 0.2
    vtolCont = 0.1
    dvtolCont = 0.2

    rmTol = 0.6
    dampConst = 7
    gridMin = 5e-7
    gridMax = 2e-4
    uniformityTol = 2.5
    absvtol = 1e-8

    boundaryTol = 5e-5
    boundaryTolRm = 1e-5
    addPointCount = 3

    centerGridMin = 1e-4;


class InitialCondition(Options):
    # Read Initial profile from a file.
    # file = "restart-steady3-split.h5"
    restartFile = None

    flameType = "premixed"

    # These properties will override parameters read from the file
    # Set them to None to use values from the restart file
    Tu = 300 # for premixed flames
    Tfuel = 300 # for diffusion flames
    Toxidizer = 300 # for diffusion flames
    fuel = "CH4:1.0"
    oxidizer = "N2:3.76, O2:1.0"
    equivalenceRatio = 0.75 # for premixed flames

    pressure = 101325

    # Generate a grid without a restart file
    nPoints = 100
    xLeft = -0.002
    xRight = 0.002

    # Sane defaults for premixed flames
    centerWidth = 0.000
    slopeWidth = 0.0005
    smoothCount = 4

    # Good choices for diffusion flames:
    #centerWidth = 0.002
    #slopeWidth = 0.001
    #smoothCount = 4


class StrainParameters(Options):
    initial = 400
    final = 400
    tStart = 0.000
    dt = 0.002
    rates = None
    # rates = [9216, 7680, 6144, 4608, 3840, 3072, 2304, 1920, 1536,
    #          1152, 960, 768, 576, 480, 384, 288, 240, 192, 144,
    #          120, 96, 72, 60, 48, 36, 30, 24, 18, 15, 12]


class PositionControl(Options):
    xInitial = 0.0025
    xFinal = 0.0025
    dt = 0.01
    tStart = 0
    proportionalGain = 10
    integralGain = 800


class Times(Options):
    tStart = 0
    globalTimestep = 5e-6
    diffusionTimestep = 1e-7
    maxTimestep = 1e-4

    # Intervals for regridding and grid adaptation
    regridTimeInterval = 100
    regridStepInterval = 20

    # Output intervals for integral flame data
    outputStepInterval = 1
    outputTimeInterval = 1e-5

    # Output intervals for flame profiles
    profileStepInterval = 5
    profileTimeInterval = 1e-4

    currentStateStepInterval = 10
    terminateStepInterval = 10


class CvodeTolerances(Options):
    relativeTolerance = 1e-6
    continuityAbsTol = 1e-10
    momentumAbsTol = 1e-7
    energyAbsTol = 1e-8
    speciesAbsTol = 1e-13
    minimumTimestep = 1e-18


class QssTolerances(Options):
    epsmax = 1e1
    epsmin = 2e-2
    dtmin = 1e-16
    iterationCount = 1
    abstol = 1e-11
    minval = 1e-60
    stabilityCheck = False


class Debug(Options):
    adaptation = False
    regridding = True
    sundials = True
    jacobian = False
    calcIC = False
    timesteps = True
    solverStats = False
    performanceStats = True
    veryVerbose = False
    flameRadiusControl = False
    sourcePoint = -1
    sourceTime = 0.0


class OutputFiles(Options):
    heatReleaseRate = True
    auxiliaryVariables = True
    timeDerivatives = True
    residualComponents = True
    firstFileNumber = 0
    saveProfiles = True
    debugIntegratorStages = False
    splitHeatReleaseRate = False


class TerminationCondition(Options):
    # Integrate until a specific time:
    tEnd = 0.010

    # Integrate until <measurement> is steady to within <tolerance> for <time>
    measurement = "Q"
    tolerance = 2e-4
    abstol = 0.5
    steadyPeriod = 0.005
    timeMax = 0.8


class Config(object):
    def __init__(self, *args):
        opts = {}
        for arg in args:
            assert isinstance(arg, Options)
            opts[arg.__class__.__name__] = arg

        get = lambda cls: opts.get(cls.__name__) or cls()

        self.paths = get(Paths)
        self.general = get(General)
        self.chemistry = get(Chemistry)
        self.adapchem = opts.get('Adapchem')
        self.transportElimination = get(TransportElimination)
        self.grid = get(Grid)
        self.initialCondition = get(InitialCondition)
        self.strainParameters = get(StrainParameters)
        self.positionControl = get(PositionControl)
        self.times = get(Times)
        self.cvodeTolerances = get(CvodeTolerances)
        self.qssTolerances = get(QssTolerances)
        self.debug = get(Debug)
        self.outputFiles = get(OutputFiles)
        self.terminationCondition = get(TerminationCondition)
