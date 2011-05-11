""" \file input.py
A set of classes for specifiying input parameters for the Pyro solver.

Create a Config object to be passed to "run" or "multirun". Sane
defaults are given for most input parameters. Create a customized
configuration by passing Options objects to the constructor for Config:

\code
conf = Config(
    Paths(inputDir="somewhere",
          outputDir="somewhereElse"),
    InitialCondition(equivalenceRatio=0.9))
\endcode
"""

import numbers
import os
import Cantera
import numpy as np

class Options(object):
    """ Base class for elements of Config """
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
    """ Directories for input and output files """

    ## Relative path to the directory where input files (such as the
    ## Cantera mechanism file) are located.
    inputDir = "input"

    ## Relative path to the directory where output files
    ## (outNNNNNN.h5, profNNNNNN.h5) will be stored.
    outputDir = "run/test1"

    ## File to use for log messages. If *None*, write output to stdout
    logFile = None

class General(Options):
    """ High-level configuration options """

    ## True if the temperature and mass fractions on the burned gas
    ## side of the flame should be held constant. Applicable only to
    ## premixed flames.
    fixedBurnedVal = True

    ## True if the position of the leftmost grid point should be held
    ## constant. Should usually be True for Twin and Curved flame
    ## configurations.
    fixedLeftLocation = False

    ## True if solving a radially propagating flame in a cylindrical
    ## coordinate system.
    curvedFlame = False

    ## True if solving a planar flame that is symmetric about the x = 0 plane.
    twinFlame = False

    ## True if the unburned fuel/air mixture should be used as the
    ## left boundary condition. Applicable only to premixed flames.
    unburnedLeft = True # applies to premixed flames

    ## True if the fuel mixture should be used as the left boundary
    ## condition. Applicable only to diffusion flames.
    fuelLeft = True

    ## Integrator to use for the chemical source terms. Choices are
    ## "qss" (explicit, quasi-steady state) and "cvode" (implicit,
    ## variable-order BDF).
    chemistryIntegrator = "qss"

    ## Method to use for splitting the convection / diffusion / reaction
    ## terms. Options are "strang" and "balanced".
    splittingMethod = "balanced"

    ## Number of integration failures to tolerate in the chemistry
    ## integrator before aborting.
    errorStopCount = 100


class Chemistry(Options):
    """ Settings pertaining to the Cantera mechanism file """

    ## Path to the Cantera mechanism file in XML format
    mechanismFile = "ucsd-methane.xml"

    ## ID of the phase to use in the mechanism file.
    ## Found on a line that looks like:
    ## \code <phase dim="3" id="gas"> \endcode
    ## in the mechanism file. This is always "gas" for mechanisms
    ## converted using ck2cti and cti2ctml.
    phaseID = "gas"

    ## Transport model to use. Valid options are "Mix", "Multi", and "Approx"
    transportModel = "Mix"

    ## Mole fraction threshold for including species with transportModel = "Approx"
    threshold = 1e-5


class Adapchem(Options):
    """ Locations of files needed when using the AdapChem module """

    ## Chemkin-format mechanism. Must match the mechanism specified in
    ## Chemistry.mechanismFile
    mechanism = "chem.bin"

    ## Tolerance settings for AdapChem
    infile = "adapchem.in"

    ## A file listing a number of reduced models. For each model, this
    ## file contains a bit vector indicating which species and
    ## reactions are included, and a range of temperatures and species
    ## concentrations over which the model is valid.
    models = "modelsV2.chem"

    ## ??
    defaultModel = "default.model"

    ## ??
    donemodels = "donemodels"

    ## ??
    restart = "full.rstrt"

class TransportElimination(Options):
    """
    Options for integrating each of the species transport
    equations on a reduced set of grid points.
    """

    ## Integrate the diffusion equation on the reduced domain
    diffusion = False

    ## Integrate the convection equation on the reduced domain
    convection = False

    ## Absolute threshold on the time derivative of the mass fraction
    ## for each species at a point.
    atol = 1e-3

    ## Number of timesteps to take before re-evaluating the transport
    ## domain for each species.
    stepInterval = 1


class Grid(Options):
    """ Parameters controlling the adaptive grid """
    ## Maximum relative scalar variation of each state vector
    ## component between consecutive grid points. For high accuracy,
    ## vtol = 0.08; For minimal accuracy, vtol = 0.20.
    vtol = 0.12

    ## Maximum relative variation of the gradient of each state vector
    ## component between consecutive grid points. For high accuracy,
    ## dvtol = 0.12; For minimal accuracy, dvtol = 0.4.
    dvtol = 0.2

    ## Relative tolerance (compared to vtol and dvtol) for grid point
    ## removal.
    rmTol = 0.6

    ## Parameter to limit numerical diffusion in regions with high
    ## convective velocities.
    dampConst = 7

    ## Minimum grid spacing [m]
    gridMin = 5e-7

    ## Maximum grid spacing [m]
    gridMax = 2e-4

    ## Maximum ratio of the distances between adjacent pairs of grid
    ## points.
    ## \f[ \frac{1}{\tt uniformityTol} <
    ##     \frac{x_{j+1}-x_j}{x_j-x_{j-1}} <
    ##     {\tt uniformityTol} \f]
    uniformityTol = 2.5

    ## State vector components smaller than this value are not
    ## considered whether to add or remove a grid point.
    absvtol = 1e-8

    ## Tolerance for each state vector component for extending the
    ## domain to satisfy zero-gradient conditions at the left and
    ## right boundaries.
    boundaryTol = 5e-5

    ## Tolerance for removing points at the boundary. Must be smaller
    ## than boundaryTol.
    boundaryTolRm = 1e-5

    ## For unstrained flames, number of flame thicknesses (based on
    ## reaction zone width) downstream of the flame to keep the right
    ## edge of the domain.
    unstrainedDownstreamWidth = 5

    ## Number of points to add when extending a boundary to satisfy
    ## boundaryTol.
    addPointCount = 3

    ## For curved or twin flames, the minimum position of the first
    ## grid point past x = 0.
    centerGridMin = 1e-4;


class InitialCondition(Options):
    """
    Settings controlling the initial condition for the integrator. If
    no restartFile is specified, an initial profile is created based
    on the specified fuel and oxidizer compositions.

    If an input file is specified, then setting fuel and oxidizer
    compositions will cause new values to be used only at the
    boundaries of the domain.

    The grid parameters are only used if no restart file is specified.
    """

    ## Read initial profiles from the specified file, or if 'None',
    ## create a new initial profile.
    restartFile = None

    ## True if the restart file path is relative to Paths.inputDir.
    relativeRestartPath = True

    ## "premixed" or "diffusion"
    flameType = "premixed"

    ## Temperature of the unburned fuel/air mixture for premixed flames [K].
    Tu = 300

    ## Temperature of the fuel mixture for diffusion flames [K].
    Tfuel = 300

    ## Temperature of the oxidizer mixture for diffusion flames [K].
    Toxidizer = 300

    ## Molar composition of the fuel mixture.
    fuel = "CH4:1.0"

    ## Molar composition of the oxidizer mixture.
    oxidizer = "N2:3.76, O2:1.0"

    ## Equivalence ratio of the fuel/air mixture for premixed flames.
    equivalenceRatio = 0.75

    ## Thermodynamic pressure [Pa]
    pressure = 101325

    ## Number of points in the initial uniform grid.
    nPoints = 100

    ## Position of the leftmost point of the initial grid.
    xLeft = -0.002

    ## Position of the rightmost point of the initial grid.
    xRight = 0.002

    ## The width of the central plateau in the initial profile
    ## [m]. For premixed flames, this mixture is composed of oxidizer
    ## at the unburned gas temperature. Recommended value: 0.0. For
    ## diffusion flames, this mixture is composed of a stoichiometric
    ## fuel/air brought to equilibrium at constant enthalpy and
    ## pressure. Recommended value: 0.002.
    centerWidth = 0.000

    ## The width of the slope away from the central plateau in the
    ## initial profile [m]. Recommended value for premixed flames:
    ## 5e-4. Recommended value for diffusion flames: 1e-3.
    slopeWidth = 0.0005

    # Number of times to run the generated profile through a low pass
    # filter before starting the simulation.
    smoothCount = 4


class WallFlux(Options):
    ## Reference temperature for the wall heat flux
    Tinf = 300

    ## Conductance of the wall [W/m^2-K]
    Kwall = 100


class Ignition(Options):
    """
    Parameters for an artificial heat release rate function which can
    be used to simulate ignition. The heat release rate is a step
    function in time with a Gaussian spatial distribution.
    """

    ## Beginning of the external heat release rate pulse [s].
    tStart = 0

    ## Duration of the external heat release rate pulse [s].
    duration = 1e-3

    ## Integral amplitude of the pulse [W/m^2].
    energy = 0

    # Location of the center of the pulse [m].
    center = 0 # [m]

    # Characteristic width (standard deviation) of the pulse [m].
    stddev = 1e-4


class StrainParameters(Options):
    """
    Parameters defining the strain rate as a function of time.

    The strain rate changes linearly from #initial to #final over a
    period of #dt seconds, starting at #tStart.
    """
    initial = 400
    final = 400
    tStart = 0.000
    dt = 0.002

    ## A list of strain rates to use for a series of sequential
    ## integration periods (see #utils.multirun), with steady-state profiles
    ## generated for each strain rate before proceeding to the next. A
    ## typical list of strain rates to use:

    ## \code
    ## rates = [9216, 7680, 6144, 4608, 3840, 3072, 2304, 1920, 1536,
    ##          1152, 960, 768, 576, 480, 384, 288, 240, 192, 144, 120,
    ##          96, 72, 60, 48, 36, 30, 24, 18, 15, 12]
    ## \endcode
    rates = None


class PositionControl(Options):
    """
    Parameters defining the position of the flame as a function of
    time (for twin / curved flames).

    These parameters are used to adjust the mass flux at r = 0 to move
    the flame toward the desired location. The flame moves from
    #xInitial to #xFinal over #dt seconds, starting at #tStart.

    The feedback controller which determines the mass flux uses the
    distance between the current flame location and the desired flame
    location with the gains specified by #proportionalGain and
    #integralGain.
    """

    xInitial = 0.0025
    xFinal = 0.0025
    dt = 0.01
    tStart = 0

    proportionalGain = 10
    integralGain = 800


class Times(Options):
    """
    Paremeters controlling integrator timesteps and frequency of
    output profiles.
    """

    ## Integrator start time.
    tStart = 0

    ## Timestep used for operator splitting.
    globalTimestep = 5e-6

    ## Timestep used by the diffusion integrator.
    diffusionTimestep = 1e-7

    ## Maximum amount of time before regridding / adaptation.
    regridTimeInterval = 100

    ## Maximum number of timesteps before regridding / adaptation.
    regridStepInterval = 20

    ## Maximum number of steps between storing integral flame
    ## properties.
    outputStepInterval = 1

    ## Maximum time between storing integral flame properties.
    outputTimeInterval = 1e-5

    # Maximum number of timesteps before writing flame profiles.
    profileStepInterval = 5

    # Maximum time between writing flame profiles.
    profileTimeInterval = 1e-4

    # Number of timesteps between writing profNow.h5
    currentStateStepInterval = 10

    # Number of timesteps between checks of the steady-state
    # termination conditions.
    terminateStepInterval = 10


class CvodeTolerances(Options):
    """ Tolerances for the CVODE chemistry integrator """

    ## Relative tolerance for each state variable
    relativeTolerance = 1e-6

    ## Absolute tolerance for U velocity
    momentumAbsTol = 1e-7

    ## Absolute tolerance for T
    energyAbsTol = 1e-8

    ## Absolute tolerance for species mass fractions.
    speciesAbsTol = 1e-13

    ## Minimum interal timestep
    minimumTimestep = 1e-18


class QssTolerances(Options):
    """ Tolerances for the QSS chemistry integrator """

    ## Accuracy parameter for determining the next timestep.
    epsmin = 2e-2

    ## Accuracy parameter for repeating timesteps
    epsmax = 1e1

    ## Minimum internal timestep
    dtmin = 1e-16

    ## Maximum internal timestep
    dtmax = 1e-6

    ## Number of corrector iterations per timestep.
    iterationCount = 1

    ## Absolute threshold for including each component in the accuracy
    ## tests.
    abstol = 1e-11

    ## Lower limit on the value of each state vector component.
    minval = 1e-60

    ## Enable convergence-based stability check on timestep. Not
    ## enabled unless #iterationCount >= 3.
    stabilityCheck = False


class Debug(Options):
    """ Control of verbose debugging output """

    ## Addition / removal of internal grid points.
    adaptation = False

    ## Addition / removal of boundary grid points.
    regridding = True

    ## Print current time after each global timestep.
    timesteps = True

    ## Enable extensive timestep debugging output. Automatically
    ## enabled for debug builds.
    veryVerbose = False

    ## Print information about the flame radius feedback controller.
    flameRadiusControl = False

    ## Grid point to print debugging information about at #sourceTime.
    sourcePoint = -1

    ## Time at which to print extensive debugging information about
    ## the source term at j = #sourcePoint, then terminate.
    sourceTime = 0.0


class OutputFiles(Options):
    """ Control the contents of the periodic output files """

    ## Include the heat release rate as a function of space
    heatReleaseRate = True

    ## Include the reaction / diffusion / convection contributions to
    ## the net time derivative of each state variable
    timeDerivatives = True

    ## Include variables such as transport properties and grid
    ## parameters that can be recomputed from the state variables.
    extraVariables = False

    ## Include other miscellaneous variables
    auxiliaryVariables = False

    ## Used to generate a continuous sequence of output files after
    ## restarting the code.
    firstFileNumber = 0

    ## Generate profNNNNNN.h5 files
    saveProfiles = True

    ## Write profiles after each stage of the split integrator.
    debugIntegratorStages = False


class TerminationCondition(Options):
    """
    Integrate until either:

    (1) #tEnd is reached
    (2) The heat release rate (#measurement) reaches a steady-state
        value to within #tolerance (RMS) over a time period of
        #steadyPeriod, or the mean heat release rate over #steadyPeriod
        is less than #abstol.

    To disable steady-state check, set #measurement = None, and give a
    sane value for #tEnd.
    """

    tEnd = 0.8

    measurement = "Q"
    tolerance = 1e-4
    abstol = 0.5
    steadyPeriod = 0.002


class Config(object):
    """
    An object consisting of a set of Options objects which define a
    complete set of configuration options needed to run the flame
    solver.
    """
    def __init__(self, *args):
        opts = {}
        for arg in args:
            if not isinstance(arg, Options):
                raise TypeError('%r is not an instance of class Options' % arg)
            name = arg.__class__.__name__
            if name in opts:
                raise ValueError('Multiple instances of class %r encountered' % name)
            opts[name] = arg

        get = lambda cls: opts.get(cls.__name__) or cls()

        self.paths = get(Paths)
        self.general = get(General)
        self.chemistry = get(Chemistry)
        self.adapchem = opts.get('Adapchem')
        self.transportElimination = get(TransportElimination)
        self.grid = get(Grid)
        self.initialCondition = get(InitialCondition)
        self.wallFlux = opts.get('WallFlux')
        self.ignition = get(Ignition)
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
        error = self.checkString(self.chemistry,
                                'transportModel', ('Multi', 'Mix', 'Approx')) or error
        error = self.checkString(self.general,
                                 'splittingMethod', ('strang', 'balanced')) or error

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

        # Make sure that the mechanism file has sane rate coefficients
        if self.initialCondition.flameType == 'premixed':
            Tcheck = self.initialCondition.Tu
        else:
            Tcheck = min(self.initialCondition.Tfuel, self.initialCondition.Toxidizer)
        error |= self.checkRateConstants(gas, Tcheck)

        # Make sure the restart file is in the correct place (if specified)
        if self.initialCondition.restartFile:
            if self.initialCondition.relativeRestartPath:
                restart = os.path.join(self.paths.inputDir,
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

        if error:
            print 'Validation failed.'
        else:
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

    def checkRateConstants(self, gas, T):
        """
        A function for finding reactions with suspiciously high
        rate constants at low temperatures.
        """
        gas.set(T=300, P=101325, Y=np.ones(gas.nSpecies()))
        Rf = gas.fwdRateConstants()
        Rr = gas.revRateConstants()
        error = False
        for i in range(len(Rf)):
            if Rf[i] > 1e30:
                error = True
                print ('WARNING: Excessively high forward rate constant'
                       ' for reaction %i at T = %6.2f K' % (i+1,T))
                print '    Reaction equation: %s' % gas.reactionEqn(i)
                print '    Forward rate constant: %e' % Rf[i]

            if Rr[i] > 1e30:
                error = True
                print ('WARNING: Excessively high reverse rate constant'
                       ' for reaction %i at T = %6.2f K' % (i+1,T))
                print '    Reaction equation: %s' % gas.reactionEqn(i)
                print '    Reverse rate constant: %e' % Rr[i]

        return error
