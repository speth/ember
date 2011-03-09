import numbers
import os
import Cantera

class Options(object):
    def __init__(self, **kwargs):
        for key,value in kwargs.iteritems():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                raise KeyError('Unrecognized configuration option: %s' % key)

    def _stringify(self, indent=0):
        ans = []
        spaces = None
        for attr in dir(self):
            if not attr.startswith('_'):
                value = getattr(self, attr)
                if not isinstance(value, numbers.Number):
                    value = '%r' % value

                if not spaces:
                    header = ' '*indent + self.__class__.__name__ + '('
                    spaces = ' '*len(header)

                else:
                    header = spaces

                ans.append('%s%s=%s,' % (header, attr, value))

        ans[-1] = ans[-1][:-1] + ')'

        return ans


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
    mechanismFile = "ucsd-methane.xml"
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
    relativeRestartPath = True # restart file path is relative to inputDir

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
    # Integrate until either of the following conditions is satisfied.
    # To disable steady-state check, set 'measurement=None', and give
    # a sane value for tEnd.

    # Integrate until a specific time:
    tEnd = 1000

    # Integrate until <measurement> is steady to within <tolerance> for <time>
    measurement = "Q" # integral heat release rate
    tolerance = 1e-4 # allowable relative RMS variation
    abstol = 0.5 # allowable absolute RMS variation
    steadyPeriod = 0.002 # time period over which heat release rate must be steady
    timeMax = 0.8 # give up at this time


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
        self.positionControl = opts.get('PositionControl')
        self.times = get(Times)
        self.cvodeTolerances = get(CvodeTolerances)
        self.qssTolerances = get(QssTolerances)
        self.debug = get(Debug)
        self.outputFiles = get(OutputFiles)
        self.terminationCondition = get(TerminationCondition)

    def stringify(self):
        ans = []
        for item in self.__dict__.itervalues():
            if isinstance(item, Options):
                ans.append('\n'.join(item._stringify(4)))
        ans[-1] = ans[-1]+')'
        return 'conf = Config(\n' + ',\n'.join(ans[:-1]) + ')\n'

    def validate(self):
        error = False

        # Check valid choices for string options
        error = self.checkString(self.terminationCondition,
                                 'measurement', ('Q', None)) or error
        error = self.checkString(self.initialCondition,
                                 'flameType', ('premixed', 'diffusion')) or error
        error = self.checkString(self.general,
                                 'chemistryIntegrator', ('cvode', 'qss')) or error
        erro = self.checkString(self.chemistry,
                                'transportModel', ('Multi', 'Mix')) or error

        # Make sure the mechanism file is in the correct place
        mech = os.path.join(self.paths.inputDir, self.chemistry.mechanismFile)
        if not os.path.exists(mech):
            error = True
            print "Error: Couldn't find mechanism file %r.\n" % mech

        # Make sure that the mechanism file actually works and contains the
        # specified fuel and oxidizer species
        gas = Cantera.IdealGasMix(mech, id=self.chemistry.phaseID)
        gas.set(X=self.initialCondition.fuel)
        gas.set(X=self.initialCondition.oxidizer)

        # Make sure the restart file is in the correct place (if specified)
        if self.initialCondition.restartFile:
            if self.initialCondition.relativeRestartPath:
                restar = os.path.join(self.paths.inputDir,
                                      self.initialCondition.restartFile)
            else:
                restart = self.initialCondition.restartFile

            if not os.path.exists(restart):
                error = True
                print "Error: Couldn't find restart file %r.\n" % restart

        # QSS Solver is incompatible with adapchem
        if self.adapchem and general.chemistryIntegrator == 'qss':
            error = True
            print 'Error: QSS integrator is incompatible with Adapchem.\n'

        if not error:
            print 'Validation completed successfully.'

    def checkString(self, options, attrname, choices):
        value = getattr(options, attrname)
        if value not in choices:
            print 'Error: Invalid option for %s.%s: %r' % (options.__class__.__name__,
                                                           attrname,
                                                           value)
            print '       Valid options are: %r\n' % list(choices)
            return True
        else:
            return False