"""
A set of classes for specifying input parameters for the Pyro solver.

Create a :class:`.Config` object to be passed to :func:`~pyro.utils.run` or
:func:`~pyro.utils.multirun`.
Sane defaults are given for most input parameters. Create a customized
configuration by passing :class:`.Options` objects to the constructor for
:class:`.Config`::

    conf = Config(
        Paths(inputDir="somewhere",
              outputDir="somewhereElse"),
        InitialCondition(equivalenceRatio=0.9))
"""

import numbers
import os
import types
import Cantera
import numpy as np
import utils

class Option(object):
    """
    Instances of this class are used as class members of descendants of class
    :class:`Options` to represent a single configurable value. When a
    user-specified value for an option is specified (as a keyword argument to
    the constructor of a class derived from :class:`Options`), that value is
    stored in this object and validated to make sure it satisfies any
    applicable constraints.

    :param default:
        The default value for this option
    :param choices:
        A sequence of valid values for this option, e.g. ``['Mix', 'Multi']``.
        The *default* choice is automatically included in *choices*.
    :param min:
        The minimum valid value for this option
    :param max:
        The maximum valid value for this option
    :param nullable:
        Set to *True* if *None* is a valid value for this option, regardless
        of any other restrictions. Automatically set to *True* if *default* is
        *None*.

    Requiring values of a particular type is done by using one of the derived
    classes: :class:`StringOption`, :class:`BoolOption`, :class:`IntegerOption`,
    :class:`FloatOption`.
    """
    def __init__(self, default, choices=None, min=None, max=None, nullable=False):
        self.value = default
        self.default = default
        if choices:
            self.choices = set(choices)
            self.choices.add(default)
        else:
            self.choices = None

        self.min = min
        self.max = max
        self.isSet = False
        self.nullable = nullable or self.default is None

    def validate(self):
        if self.choices and self.value not in self.choices:
            return '%r not in %r' % (self.value, list(self.choices))

        if self.min is not None and self.value < self.min:
            return 'Value (%s) must be greater than or equal to %s' % (self.value, self.min)

        if self.max is not None and self.value > self.max:
            return 'Value (%s) must be less than or equal to %s' % (self.value, self.max)

    def __repr__(self):
        return repr(self.value)

    def __nonzero__(self):
        return bool(self.value)

    def __eq__(self, other):
        try:
            return self.value == other.value
        except AttributeError:
            return self.value == other


class StringOption(Option):
    """ An option whose value must be a string. """
    def validate(self):
        if not (isinstance(self.value, types.StringTypes) or
                (self.value is None and self.nullable)):
            return 'Value must be a string. Got %r' % self.value
        return Option.validate(self)


class BoolOption(Option):
    """ An option whose value must be a boolean value. """
    def validate(self):
        if not (self.value in (True, False, 0, 1) or
                (self.value is None and self.nullable)):
            return 'Value must be a boolean. Got %r' % self.value
        return Option.validate(self)


class IntegerOption(Option):
    """ An option whose value must be a integer. """
    def validate(self):
        if not (isinstance(self.value, numbers.Integral) or
                (self.value is None and self.nullable)):
            return 'Value must be an integer. Got %r' % self.value
        return Option.validate(self)


class FloatOption(Option):
    """ An option whose value must be a floating point number. """
    def validate(self):
        if not (isinstance(self.value, numbers.Number) or
                (self.value is None and self.nullable)):
            return 'Value must be a number. Got %r' % self.value
        return Option.validate(self)


class Options(object):
    """ Base class for elements of :class:`.Config` """
    def __init__(self, **kwargs):
        for key,value in kwargs.iteritems():
            if hasattr(self, key):
                opt = getattr(self, key)
                if isinstance(opt, Option):
                    opt.value = value
                    opt.isSet = True
                    message = opt.validate()
                    if message:
                        raise ValueError('\nInvalid option specified for %s.%s:\n%s' %
                                         (self.__class__.__name__, key, message))
                else:
                    setattr(self, key, value)
            else:
                raise KeyError('Unrecognized configuration option: %s' % key)

    def _stringify(self, indent=0):
        ans = []
        spaces = None
        for attr in dir(self):
            if not attr.startswith('_'):
                value = getattr(self, attr)
                if isinstance(value, Option):
                    if value.value == value.default:
                        continue
                    else:
                        value = value.value

                if not isinstance(value, numbers.Number):
                    value = repr(value)

                if not spaces:
                    header = ' '*indent + self.__class__.__name__ + '('
                    spaces = ' '*len(header)

                else:
                    header = spaces

                ans.append('%s%s=%s,' % (header, attr, value))

        if ans:
            ans[-1] = ans[-1][:-1] + ')'
        else:
            ans = ''

        return ans

    def __iter__(self):
        for name,opt in self.__dict__.iteritems():
            if isinstance(opt, Option):
                yield name,opt

        for name,opt in self.__class__.__dict__.iteritems():
            if isinstance(opt, Option):
                yield name,opt


class Paths(Options):
    """ Directories for input and output files """

    #: Relative path to the directory where input files (such as the
    #: Cantera mechanism file) are located.
    inputDir = StringOption("input")

    #: Relative path to the directory where output files
    #: (outNNNNNN.h5, profNNNNNN.h5) will be stored.
    outputDir = StringOption("run/test1")

    #: File to use for log messages. If *None*, write output to stdout
    logFile = StringOption(None)

class General(Options):
    """ High-level configuration options """

    #: True if the temperature and mass fractions on the burned gas
    #: side of the flame should be held constant. Applicable only to
    #: premixed flames.
    fixedBurnedVal = BoolOption(True)

    #: True if the position of the leftmost grid point should be held
    #: constant. Should usually be True for Twin and Curved flame
    #: configurations.
    fixedLeftLocation = BoolOption(False)

    #: True if solving a radially propagating flame in a cylindrical
    #: coordinate system.
    curvedFlame = BoolOption(False)

    #: True if solving a planar flame that is symmetric about the x = 0 plane.
    twinFlame = BoolOption(False)

    #: Input file (HDF5 format) containing the interpolation data needed for
    #: the quasi2d mode. Contains:
    #:
    #: - vector *r* (length *N*)
    #: - vector *z* (length *M*)
    #: - Temperature array *T* (size *N* x *M*)
    #: - radial velocity array *vr* (size *N* x *M*)
    #: - axial velocity array *vz* (size *N* x *M*)
    interpFile = StringOption(None)

    #: *True* if the unburned fuel/air mixture should be used as the
    #: left boundary condition. Applicable only to premixed flames.
    unburnedLeft = BoolOption(True) # applies to premixed flames

    #: *True* if the fuel mixture should be used as the left boundary
    #: condition. Applicable only to diffusion flames.
    fuelLeft = BoolOption(True)

    #: Method for setting the boundary condition for the continuity equation
    #: Valid options are: ``fixedLeft``, ``fixedRight``, ``fixedQdot``,
    #: ``fixedTemperature``, and ``stagnationPoint``.
    continuityBC = StringOption("fixedLeft",
                                ("fixedRight", "fixedQdot", "fixedTemperature",
                                 "stagnationPoint"))

    #: Integrator to use for the chemical source terms. Choices are
    #: ``qss`` (explicit, quasi-steady state) and ``cvode`` (implicit,
    #: variable-order BDF).
    chemistryIntegrator = StringOption("qss", ("cvode",))

    #: Method to use for splitting the convection / diffusion / reaction
    #: terms. Options are ``strang`` and ``balanced``.
    splittingMethod = StringOption("balanced", ("strang",))

    #: Number of integration failures to tolerate in the chemistry
    #: integrator before aborting.
    errorStopCount = IntegerOption(100)

    #: Number of threads to use for evaluating terms in parallel
    nThreads = IntegerOption(1, min=1)


class Chemistry(Options):
    """ Settings pertaining to the Cantera mechanism file """

    #: Path to the Cantera mechanism file in XML format
    mechanismFile = StringOption("ucsd-methane.xml")

    #: ID of the phase to use in the mechanism file.
    #: Found on a line that looks like::
    #:
    #:     <phase dim="3" id="gas">
    #:
    #: in the mechanism file. This is always "gas" for mechanisms
    #: converted using ck2cti and cti2ctml.
    phaseID = StringOption("gas")

    #: Transport model to use. Valid options are ``Mix``, ``Multi``, and ``Approx``
    transportModel = StringOption("Mix", ("Multi", "Approx"))

    #: Mole fraction threshold for including species with ``transportModel = "Approx"``
    threshold = FloatOption(1e-5)


class Grid(Options):
    """ Parameters controlling the adaptive grid """
    #: Maximum relative scalar variation of each state vector
    #: component between consecutive grid points. For high accuracy,
    #: ``vtol = 0.08``; For minimal accuracy, ``vtol = 0.20``.
    vtol = FloatOption(0.12)

    #: Maximum relative variation of the gradient of each state vector
    #: component between consecutive grid points. For high accuracy,
    #: ``dvtol = 0.12``; For minimal accuracy, ``dvtol = 0.4``.
    dvtol = FloatOption(0.2)

    #: Relative tolerance (compared to vtol and dvtol) for grid point removal.
    rmTol = FloatOption(0.6)

    #: Parameter to limit numerical diffusion in regions with high
    #: convective velocities.
    dampConst = FloatOption(7)

    #: Minimum grid spacing [m]
    gridMin = FloatOption(5e-7)

    #: Maximum grid spacing [m]
    gridMax = FloatOption(2e-4)

    #: Maximum ratio of the distances between adjacent pairs of grid
    #: points.
    #:
    #: .. math:: \frac{1}{\tt uniformityTol} < \frac{x_{j+1}-x_j}{x_j-x_{j-1}}
    #:           < {\tt uniformityTol}
    uniformityTol = FloatOption(2.5)

    #: State vector components smaller than this value are not
    #: considered whether to add or remove a grid point.
    absvtol = FloatOption(1e-8)

    #: Tolerance for each state vector component for extending the
    #: domain to satisfy zero-gradient conditions at the left and
    #: right boundaries.
    boundaryTol = FloatOption(5e-5)

    #: Tolerance for removing points at the boundary. Must be smaller
    #: than boundaryTol.
    boundaryTolRm = FloatOption(1e-5)

    #: For unstrained flames, number of flame thicknesses (based on
    #: reaction zone width) downstream of the flame to keep the right
    #: edge of the domain.
    unstrainedDownstreamWidth = FloatOption(5)

    #: Number of points to add when extending a boundary to satisfy
    #: boundaryTol.
    addPointCount = IntegerOption(3, min=0)

    #: For curved or twin flames, the minimum position of the first
    #: grid point past x = 0.
    centerGridMin = FloatOption(1e-4, min=0)


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

    #: Read initial profiles from the specified file, or if 'None',
    #: create a new initial profile.
    restartFile = StringOption(None)

    #: True if the restart file path is relative to Paths.inputDir.
    relativeRestartPath = BoolOption(True)

    #: "premixed" or "diffusion"
    flameType = StringOption('premixed', ('diffusion', 'quasi2d'))

    #: Temperature of the unburned fuel/air mixture for premixed flames [K].
    Tu = FloatOption(300)

    #: Temperature of the fuel mixture for diffusion flames [K].
    Tfuel = FloatOption(300)

    #: Temperature of the oxidizer mixture for diffusion flames [K].
    Toxidizer = FloatOption(300)

    #: Molar composition of the fuel mixture.
    fuel = StringOption("CH4:1.0")

    #: Molar composition of the oxidizer mixture.
    oxidizer = StringOption("N2:3.76, O2:1.0")

    #: Equivalence ratio of the fuel/air mixture for premixed flames.
    equivalenceRatio = FloatOption(0.75, min=0)

    #: Thermodynamic pressure [Pa]
    pressure = FloatOption(101325, min=0)

    #: Number of points in the initial uniform grid.
    nPoints = IntegerOption(100)

    #: Position of the leftmost point of the initial grid.
    xLeft = FloatOption(-0.002)

    #: Position of the rightmost point of the initial grid.
    xRight = FloatOption(0.002)

    #: The width of the central plateau in the initial profile
    #: [m]. For premixed flames, this mixture is composed of oxidizer
    #: at the unburned gas temperature. Recommended value: 0.0. For
    #: diffusion flames, this mixture is composed of a stoichiometric
    #: fuel/air brought to equilibrium at constant enthalpy and
    #: pressure. Recommended value: 0.002.
    centerWidth = FloatOption(0.000)

    #: The width of the slope away from the central plateau in the
    #: initial profile [m]. Recommended value for premixed flames:
    #: 5e-4. Recommended value for diffusion flames: 1e-3.
    slopeWidth = FloatOption(0.0005)

    #: Number of times to run the generated profile through a low pass
    #: filter before starting the simulation.
    smoothCount = IntegerOption(4)

    #: True if initial profiles for x,T,U and Y (as well as rVzero) are given
    haveProfiles = BoolOption(False)

    x = Option(None)
    T = Option(None)
    U = Option(None)
    Y = Option(None)
    rVzero = FloatOption(None)


class WallFlux(Options):
    #: Reference temperature for the wall heat flux
    Tinf = FloatOption(300)

    #: Conductance of the wall [W/m^2-K]
    Kwall = FloatOption(100)


class Ignition(Options):
    """
    Parameters for an artificial heat release rate function which can
    be used to simulate ignition. The heat release rate is a step
    function in time with a Gaussian spatial distribution.
    """

    #: Beginning of the external heat release rate pulse [s].
    tStart = FloatOption(0)

    #: Duration of the external heat release rate pulse [s].
    duration = FloatOption(1e-3)

    #: Integral amplitude of the pulse [W/m^2].
    energy = FloatOption(0)

    #: Location of the center of the pulse [m].
    center = FloatOption(0)

    #: Characteristic width (standard deviation) of the pulse [m].
    stddev = FloatOption(1e-4)


class StrainParameters(Options):
    """
    Parameters defining the strain rate as a function of time.

    The strain rate changes linearly from *initial* to *final* over a
    period of *dt* seconds, starting at *tStart*.
    """
    initial = FloatOption(400) #: Initial strain rate [1/s]
    final = FloatOption(400) #: final strain rate [1/s]
    tStart = FloatOption(0.000) #: time at which strain rate starts to change [s]
    dt = FloatOption(0.002) #: time period over which strain rate changes [s]

    #: A list of strain rates to use for a series of sequential
    #: integration periods (see :func:`~pyro.utils.multirun`), with steady-state
    #: profiles generated for each strain rate before proceeding to the next.
    #: A typical list of strain rates to use::
    #:
    #:     rates = [9216, 7680, 6144, 4608, 3840, 3072, 2304, 1920, 1536,
    #:              1152, 960, 768, 576, 480, 384, 288, 240, 192, 144, 120,
    #:              96, 72, 60, 48, 36, 30, 24, 18, 15, 12]
    rates = Option(None)


class PositionControl(Options):
    """
    Parameters defining the position of the flame as a function of
    time (for twin / curved flames).

    These parameters are used to adjust the mass flux at r = 0 to move
    the flame toward the desired location. The flame moves from
    *xInitial* to *xFinal* over *dt* seconds, starting at *tStart*.

    The feedback controller which determines the mass flux uses the
    distance between the current flame location and the desired flame
    location with the gains specified by *proportionalGain* and
    *integralGain*.
    """

    xInitial = FloatOption(0.0025) #:
    xFinal = FloatOption(0.0025) #:
    dt = FloatOption(0.01) #:
    tStart = FloatOption(0) #:

    proportionalGain = FloatOption(10) #:
    integralGain = FloatOption(800) #:


class Times(Options):
    """
    Paremeters controlling integrator timesteps and frequency of
    output profiles.
    """

    #: Integrator start time.
    tStart = FloatOption(0)

    #: Timestep used for operator splitting.
    globalTimestep = FloatOption(5e-6)

    #: Control for timestep used by the diffusion integrator.
    #: Actual timestep will be this multiplier times the stability
    #: limit an explicit integrator.
    diffusionTimestepMultiplier = FloatOption(10)

    #: Maximum amount of time before regridding / adaptation.
    regridTimeInterval = FloatOption(100)

    #: Maximum number of timesteps before regridding / adaptation.
    regridStepInterval = IntegerOption(20)

    #: Maximum number of steps between storing integral flame
    #: properties.
    outputStepInterval = IntegerOption(1)

    #: Maximum time between storing integral flame properties.
    outputTimeInterval = FloatOption(1e-5)

    #: Maximum number of timesteps before writing flame profiles.
    profileStepInterval = IntegerOption(5)

    #: Maximum time between writing flame profiles.
    profileTimeInterval = FloatOption(1e-4)

    #: Number of timesteps between writing profNow.h5
    currentStateStepInterval = IntegerOption(10)

    #: Number of timesteps between checks of the steady-state
    #: termination conditions.
    terminateStepInterval = IntegerOption(10)


class CvodeTolerances(Options):
    """ Tolerances for the CVODE chemistry integrator """

    #: Relative tolerance for each state variable
    relativeTolerance = FloatOption(1e-6)

    #: Absolute tolerance for U velocity
    momentumAbsTol = FloatOption(1e-7)

    #: Absolute tolerance for T
    energyAbsTol = FloatOption(1e-8)

    #: Absolute tolerance for species mass fractions.
    speciesAbsTol = FloatOption(1e-13)

    #: Minimum internal timestep
    minimumTimestep = FloatOption(1e-18)


class QssTolerances(Options):
    """ Tolerances for the QSS chemistry integrator """

    #: Accuracy parameter for determining the next timestep.
    epsmin = FloatOption(2e-2)

    #: Accuracy parameter for repeating timesteps
    epsmax = FloatOption(1e1)

    #: Minimum internal timestep
    dtmin = FloatOption(1e-16)

    #: Maximum internal timestep
    dtmax = FloatOption(1e-6)

    #: Number of corrector iterations per timestep.
    iterationCount = IntegerOption(1)

    #: Absolute threshold for including each component in the accuracy
    #: tests.
    abstol = FloatOption(1e-11)

    #: Lower limit on the value of each state vector component.
    minval = FloatOption(1e-60)

    #: Enable convergence-based stability check on timestep. Not
    #: enabled unless *iterationCount* >= 3.
    stabilityCheck = BoolOption(False)


class Debug(Options):
    """ Control of verbose debugging output """

    #: Addition / removal of internal grid points.
    adaptation = BoolOption(False)

    #: Addition / removal of boundary grid points.
    regridding = BoolOption(True)

    #: Print current time after each global timestep.
    timesteps = BoolOption(True)

    #: Enable extensive timestep debugging output. Automatically
    #: enabled for debug builds.
    veryVerbose = BoolOption(False)

    #: Print information about the flame radius feedback controller.
    flameRadiusControl = BoolOption(False)

    #: Grid point to print debugging information about at *sourceTime*.
    sourcePoint = IntegerOption(-1)

    #: Time at which to print extensive debugging information about
    #: the source term at j = *sourcePoint*, then terminate.
    sourceTime = FloatOption(0.0)


class OutputFiles(Options):
    """ Control the contents of the periodic output files """

    #: Include the heat release rate as a function of space
    heatReleaseRate = BoolOption(True)

    #: Include the reaction / diffusion / convection contributions to
    #: the net time derivative of each state variable
    timeDerivatives = BoolOption(True)

    #: Include variables such as transport properties and grid
    #: parameters that can be recomputed from the state variables.
    extraVariables = BoolOption(False)

    #: Include other miscellaneous variables
    auxiliaryVariables = BoolOption(False)

    #: Used to generate a continuous sequence of output files after
    #: restarting the code.
    firstFileNumber = IntegerOption(0)

    #: Generate ``profNNNNNN.h5`` files
    saveProfiles = BoolOption(True)

    #: Write profiles after each stage of the split integrator.
    debugIntegratorStages = BoolOption(False)


class TerminationCondition(Options):
    """
    Integrate until either:

    (1) *tEnd* is reached
    (2) The heat release rate ( *measurement* ) reaches a steady-state
        value to within *tolerance* (RMS) over a time period of
        *steadyPeriod*, or the mean heat release rate over *steadyPeriod*
        is less than *abstol*.

    To disable steady-state check, set ``measurement = None``, and give a
    sane value for *tEnd*.
    """

    tEnd = FloatOption(0.8) #:

    measurement = Option("Q", (None,)) #:
    tolerance = FloatOption(1e-4) #:
    abstol = FloatOption(0.5, min=0) #:
    steadyPeriod = FloatOption(0.002, min=0) #:


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

    def evaluate(self):
        """
        Replace all of the :class:`Option` objects with their values so that
        they can be used by the C++ extension, in
        :meth:`Config.generateInitialCondition`, etc.
        """
        for opts in self.__dict__.itervalues():
            if not isinstance(opts, Options):
                continue
            for name in dir(opts):
                opt = getattr(opts, name)
                if isinstance(opt, Option):
                    setattr(opts, name, opt.value)

    def setup(self):
        """
        Perform any steps that need to be done before this object can be used
        to create a FlameSolver, such as generating initial reactant profiles.
        """
        self.evaluate()
        if (not self.initialCondition.restartFile and
            not self.initialCondition.haveProfiles):
            self.generateInitialCondition()

    def __iter__(self):
        for item in self.__dict__.itervalues():
            if isinstance(item, Options):
                yield item

    def stringify(self):
        ans = []
        for item in self:
            text = '\n'.join(item._stringify(4))
            if text:
                ans.append(text)

        return 'conf = Config(\n' + ',\n'.join(ans) + ')\n'

    def validate(self):
        error = False

        # Position control can only be used with "twin" or "curved" flames
        if (self.positionControl is not None and
            not self.general.twinFlame and
            not self.general.curvedFlame):
            error = True
            print ("Error: PositionControl can only be used when either 'twinFlame'\n"
                   "       or 'curvedFlame' is set to True.")

        # twinFlame and curvedFlame are mutually exclusive:
        if self.general.curvedFlame and self.general.twinFlame:
            error = True
            print "Error: 'twinFlame' and 'curvedFlame' are mutually exclusive."

        # the "fuelLeft" option only makes sense for diffusion flames
        if (self.initialCondition.flameType == 'premixed' and
            self.general.fuelLeft.isSet):
            error = True
            print "Error: 'general.fuelLeft' should not be specified for premixed flames."

        # the "unburnedLeft" option only makes sense for premixed flames
        if (self.initialCondition.flameType == 'diffusion' and
            self.general.unburnedLeft.isSet):
            error = True
            print "Error: 'general.unburnedLeft' should not be specified for diffusion flames."

        # Make sure the mechanism file is in the correct place
        mech = os.path.join(self.paths.inputDir.value,
                            self.chemistry.mechanismFile.value)
        if not os.path.exists(mech):
            error = True
            print "Error: Couldn't find mechanism file %r.\n" % mech

        # Make sure that the mechanism file actually works and contains the
        # specified fuel and oxidizer species
        gas = Cantera.IdealGasMix(mech, id=self.chemistry.phaseID.value)
        gas.set(X=self.initialCondition.fuel.value)
        gas.set(X=self.initialCondition.oxidizer.value)

        # Make sure that the mechanism file has sane rate coefficients
        if self.initialCondition.flameType == 'premixed':
            Tcheck = self.initialCondition.Tu.value
        else:
            Tcheck = min(self.initialCondition.Tfuel.value,
                         self.initialCondition.Toxidizer.value)
        error |= self.checkRateConstants(gas, Tcheck)

        # Make sure the restart file is in the correct place (if specified)
        if self.initialCondition.restartFile:
            if self.initialCondition.relativeRestartPath:
                restart = os.path.join(self.paths.inputDir.value,
                                       self.initialCondition.restartFile.value)
            else:
                restart = self.initialCondition.restartFile.value

            if not os.path.exists(restart):
                error = True
                print "Error: Couldn't find restart file %r.\n" % restart

        if error:
            print 'Validation failed.'
        else:
            print 'Validation completed successfully.'


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

    def generateInitialCondition(self):
        """
        Generate initial profiles for temperature, species mass fractions, and
        velocity using the specified fuel and oxidizer compositions and flame
        configuration parameters.
        """
        gas = Cantera.IdealGasMix(self.chemistry.mechanismFile,
                                  self.chemistry.phaseID)

        IC = self.initialCondition
        N = IC.nPoints

        xLeft = (0.0 if self.general.twinFlame or self.general.curvedFlame
                 else IC.xLeft)

        x = np.linspace(xLeft, IC.xRight, N)
        T = np.zeros(N)
        Y = np.zeros((gas.nSpecies(), N))

        jm = (N-1) // 2

        # make sure the initial profile fits comfortably in the domain
        scale = 0.8 * (x[-1] - x[0]) / (IC.centerWidth + 2 * IC.slopeWidth)
        if scale < 1.0:
            IC.slopeWidth *= scale
            IC.centerWidth *= scale

        # Determine the grid indices defining each profile segment
        dx = x[1]-x[0]
        centerPointCount = int(0.5 + 0.5 * IC.centerWidth / dx)
        slopePointCount = int(0.5 + IC.slopeWidth / dx)
        jl2 = jm - centerPointCount
        jl1 = jl2 - slopePointCount
        jr1 = jm + centerPointCount
        jr2 = jr1 + slopePointCount

        if IC.flameType == 'premixed':
            # Reactants
            reactants = utils.calculateReactantMixture(gas, IC.fuel, IC.oxidizer,
                                                       IC.equivalenceRatio)
            gas.set(X=reactants, T=IC.Tu, P=IC.pressure)
            rhou = gas.density()
            Yu = gas.massFractions()

            # Products
            gas.equilibrate('HP')
            Tb = gas.temperature()
            Yb = gas.massFractions()

            # Diluent in the middle
            gas.set(X=IC.oxidizer, T=IC.Tu, P=IC.pressure)
            Y[:,jm] = gas.massFractions()

            if self.general.unburnedLeft:
                Tleft = IC.Tu
                Yleft = Yu
                Tright = Tb
                Yright = Yb
            else:
                Tleft = Tb
                Yleft = Yb
                Tright = IC.Tu
                Yright = Yu

            T[0] = Tleft
            T[-1] = Tright
            T[jm] = IC.Tu

        elif IC.flameType == 'diffusion':
            # Stoichiometric mixture at the center
            IC.equivalenceRatio = 1.0
            products = utils.calculateReactantMixture(gas, IC.fuel, IC.oxidizer,
                                                      IC.equivalenceRatio)
            gas.set(X=products, T=0.5*(IC.Tfuel+IC.Toxidizer), P=IC.pressure)
            gas.equilibrate('HP')
            Tb = gas.temperature()
            Yb = gas.massFractions()
            Y[:,jm] = Yb

            # Fuel
            gas.set(X=IC.fuel, T=IC.Tfuel, P=IC.pressure)
            Yfuel = gas.massFractions()

            # Oxidizer
            gas.set(X=IC.oxidizer, T=IC.Toxidizer, P=IC.pressure)
            Yoxidizer = gas.massFractions()

            if self.general.fuelLeft:
                Tleft = IC.Tfuel
                Yleft = Yfuel
                Tright = IC.Toxidizer
                Yright = Yoxidizer
            else:
                Tleft = IC.Toxidizer
                Yleft = Yoxidizer
                Tright = IC.Tfuel
                Yright = Yfuel

            T[0] = Tleft
            T[-1] = Tright
            T[jm] = Tb

            gas.set(Y=Yleft, T=Tleft, P=IC.pressure)
            rhou = gas.density() # arbitrary for diffusion flame

        Y[:,0] = Yleft
        Y[:,-1] = Yright

        newaxis = np.newaxis
        Y[:,1:jl1] = Y[:,0,newaxis]
        T[1:jl1] = T[0]

        ramp = np.linspace(0, 1, jl2-jl1)
        Y[:,jl1:jl2] = Y[:,0,newaxis] + (Y[:,jm]-Y[:,0])[:,newaxis]*ramp
        T[jl1:jl2] = T[0] + (T[jm]-T[0]) * ramp

        Y[:,jl2:jr1] = Y[:,jm,newaxis]
        T[jl2:jr1] = T[jm]

        ramp = np.linspace(0, 1, jr2-jr1)
        Y[:,jr1:jr2] = Y[:,jm,newaxis] + (Y[:,-1]-Y[:,jm])[:,newaxis]*ramp
        T[jr1:jr2] = T[jm] + (T[-1]-T[jm]) * ramp

        Y[:,jr2:] = Y[:,-1,newaxis]
        T[jr2:] = T[-1]

        YT = Y.T
        for _ in range(IC.smoothCount):
            utils.smooth(YT)
            utils.smooth(T)
        Y = YT.T

        U = np.zeros(N)
        for j in range(N):
            gas.set(Y=Y[:,j], T=T[j], P=IC.pressure)
            rho = gas.density()
            U[j] = self.strainParameters.initial * np.sqrt(rhou/rho)

        for _ in range(2):
            utils.smooth(U)

        IC.x = x
        IC.Y = Y
        IC.T = T
        IC.U = U
        IC.rVzero = 0.0 # @todo: Better estimate? Is this currently used?
        IC.haveProfiles = True
