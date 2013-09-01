"""
A set of classes for specifying input parameters for the Ember solver.

Create a :class:`.Config` object to be passed to :func:`~ember.utils.run` or
:func:`~ember.utils.multirun`.
Sane defaults are given for most input parameters. Create a customized
configuration by passing :class:`.Options` objects to the constructor for
:class:`.Config`::

    conf = Config(
        Paths(outputDir="somewhere"),
        InitialCondition(equivalenceRatio=0.9))
"""

import numbers
import os
import sys
import types
import Cantera
import numpy as np
import utils
import copy
import time

import _ember

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
    :param label:
        A human readable label to be used in the GUI in place of the attribute
        name to which this Option is assigned.
    :param level:
        A number from 0-3 indicating the obscurity level of this option. Used
        to selectively hide advanced options in the GUI. 0 is always shown.
        3 is never shown.
    :param filter:
        A function that takes a :class:`Config` object as an argument and
        returns *True* if this option should be enabled

    Requiring values of a particular type is done by using one of the derived
    classes: :class:`StringOption`, :class:`BoolOption`, :class:`IntegerOption`,
    :class:`FloatOption`.
    """
    counter = [0]  # used to preserve options in the order they are defined

    def __init__(self, default, choices=None, min=None, max=None,
                 nullable=False, label=None, level=0, filter=None):
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
        self.label = label
        self.level = level
        self.filter = filter
        self.sortValue = self.counter[0]
        self.counter[0] += 1

    def validate(self):
        if self.choices and self.value not in self.choices:
            return '%r not in %r' % (self.value, list(self.choices))

        if self.min is not None and self.value < self.min:
            return ('Value (%s) must be greater than or equal to %s' %
                    (self.value, self.min))

        if self.max is not None and self.value > self.max:
            return ('Value (%s) must be less than or equal to %s' %
                    (self.value, self.max))

    def __repr__(self):
        return repr(self.value)

    def __nonzero__(self):
        return bool(self.value)

    def __eq__(self, other):
        try:
            return self.value == other.value
        except AttributeError:
            return self.value == other

    def shouldBeEnabled(self, conf):
        if self.filter:
            return bool(self.filter(conf))
        else:
            return True


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
        # Copy the defaults from the class's dictionary
        for name,value in self.__class__.__dict__.iteritems():
            if isinstance(value, Option):
                setattr(self, name, copy.deepcopy(value))

        # Apply the options specified in kwargs
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
        allOpts = [item for item in self.__dict__.items()
                   if isinstance(item[1], Option)]
        allOpts.sort(key=lambda item: item[1].sortValue)
        return allOpts.__iter__()


# frequently-used filter predicates
def _isPremixed(conf):
    return conf.initialCondition.flameType == 'premixed'


def _isDiffusion(conf):
    return conf.initialCondition.flameType == 'diffusion'


def _isSymmetric(conf):
    return conf.general.curvedFlame or conf.general.twinFlame


def _usingCvode(conf):
    return conf.general.chemistryIntegrator == 'cvode'


def _usingQss(conf):
    return conf.general.chemistryIntegrator == 'qss'


class Paths(Options):
    """ Directories for input and output files """

    #: Relative path to the directory where output files (outNNNNNN.h5,
    #: profNNNNNN.h5) will be stored. Automatically created if it doesn't
    #: already exist.
    outputDir = StringOption("run/test1", label='Output Directory')

    #: File to use for log messages. If *None*, write output to stdout
    logFile = StringOption(None, label='Log file', level=2)


class General(Options):
    """ High-level configuration options """

    #: True if the temperature and mass fractions on the burned gas
    #: side of the flame should be held constant. Applicable only to
    #: premixed flames.
    fixedBurnedVal = BoolOption(True, level=1, filter=_isPremixed)

    #: True if the position of the leftmost grid point should be held
    #: constant. Should usually be True for Twin and Curved flame
    #: configurations.
    fixedLeftLocation = BoolOption(False, level=1, filter=_isSymmetric)

    #: True if solving a radially propagating flame in a cylindrical
    #: coordinate system.
    curvedFlame = BoolOption(False,
        filter=lambda conf: not conf.general.twinFlame)

    #: True if solving a planar flame that is symmetric about the x = 0 plane.
    twinFlame = BoolOption(False,
        filter=lambda conf: not conf.general.curvedFlame)

    #: Input file (HDF5 format) containing the interpolation data needed for
    #: the quasi2d mode. Contains:
    #:
    #: - vector *r* (length *N*)
    #: - vector *z* (length *M*)
    #: - Temperature array *T* (size *N* x *M*)
    #: - radial velocity array *vr* (size *N* x *M*)
    #: - axial velocity array *vz* (size *N* x *M*)
    interpFile = StringOption(None, level=2,
        filter=lambda conf: conf.initialCondition.flameType == 'quasi2d')

    #: *True* if the unburned fuel/air mixture should be used as the
    #: left boundary condition. Applicable only to premixed flames.
    unburnedLeft = BoolOption(True, level=1, filter=_isPremixed)

    #: *True* if the fuel mixture should be used as the left boundary
    #: condition. Applicable only to diffusion flames.
    fuelLeft = BoolOption(True, level=1, filter=_isDiffusion)

    #: Method for setting the boundary condition for the continuity equation
    #: Valid options are: ``fixedLeft``, ``fixedRight``, ``fixedQdot``,
    #: ``fixedTemperature``, and ``stagnationPoint``.
    continuityBC = StringOption("fixedLeft",
                                ("fixedRight", "fixedQdot", "fixedTemperature",
                                 "stagnationPoint"),
                                level=2)

    #: Integrator to use for the chemical source terms. Choices are
    #: ``qss`` (explicit, quasi-steady state) and ``cvode`` (implicit,
    #: variable-order BDF).
    chemistryIntegrator = StringOption("qss", ("cvode",), level=1)

    #: Method to use for splitting the convection / diffusion / reaction
    #: terms. Options are ``strang`` and ``balanced``.
    splittingMethod = StringOption("balanced", ("strang",), level=2)

    #: Number of integration failures to tolerate in the chemistry
    #: integrator before aborting.
    errorStopCount = IntegerOption(100, level=2)

    #: Number of threads to use for evaluating terms in parallel
    nThreads = IntegerOption(1, min=1, level=1,
                             label='Number of Processors')

    #: Order of Chebyshev polynomial to use in approximating input
    #: functions over a single global time step.
    chebyshevOrder = IntegerOption(5, min=2, level=3)


class Chemistry(Options):
    """ Settings pertaining to the Cantera mechanism file """

    #: Path to the Cantera mechanism file in XML format
    mechanismFile = StringOption("ucsd-methane.xml")

    #: ID of the phase to use in the mechanism file.
    #: Found on a line that looks like::
    #:
    #:     <phase dim="3" id="gas">
    #:
    #: in the mechanism file. This is always "gas" for mechanisms converted
    #: using ck2cti and cti2ctml. This option only needs to be specified if
    #: the desired phase is not the first phase defined in the input file.
    phaseID = StringOption("", level=1)

    #: Transport model to use. Valid options are ``Mix``, ``Multi``, and ``Approx``
    transportModel = StringOption("Mix", ("Multi", "Approx"), level=1)

    #: Kinetics model to use. Valid options are ``standard`` and ``interp``.
    kineticsModel = StringOption("standard", ("interp",), level=2)

    #: Mole fraction threshold for including species with ``transportModel = "Approx"``
    threshold = FloatOption(1e-5, level=2, label='Approx. transport threshold',
        filter=lambda conf: conf.chemistry.transportModel == 'Approx')

    #: Set a scalar multiplier for the reaction rate term as a function time
    rateMultiplierFunction = Option(None, level=3)


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
    rmTol = FloatOption(0.6, level=1)

    #: Parameter to limit numerical diffusion in regions with high
    #: convective velocities.
    dampConst = FloatOption(7, level=2)

    #: Minimum grid spacing [m]
    gridMin = FloatOption(5e-7, level=1)

    #: Maximum grid spacing [m]
    gridMax = FloatOption(2e-4, level=1)

    #: Maximum ratio of the distances between adjacent pairs of grid
    #: points.
    #:
    #: .. math:: \frac{1}{\tt uniformityTol} < \frac{x_{j+1}-x_j}{x_j-x_{j-1}}
    #:           < {\tt uniformityTol}
    uniformityTol = FloatOption(2.5, level=2)

    #: State vector components smaller than this value are not
    #: considered whether to add or remove a grid point.
    absvtol = FloatOption(1e-8, level=2)

    #: Tolerance for each state vector component for extending the
    #: domain to satisfy zero-gradient conditions at the left and
    #: right boundaries.
    boundaryTol = FloatOption(5e-5, level=2)

    #: Tolerance for removing points at the boundary. Must be smaller
    #: than boundaryTol.
    boundaryTolRm = FloatOption(1e-5, level=2)

    #: For unstrained flames, number of flame thicknesses (based on
    #: reaction zone width) downstream of the flame to keep the right
    #: edge of the domain.
    unstrainedDownstreamWidth = FloatOption(5, level=2)

    #: Number of points to add when extending a boundary to satisfy
    #: boundaryTol.
    addPointCount = IntegerOption(3, min=0, level=2)

    #: For curved or twin flames, the minimum position of the first
    #: grid point past x = 0.
    centerGridMin = FloatOption(1e-4, min=0, level=2, filter=_isSymmetric)


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
    restartFile = StringOption(None, level=1)

    #: "premixed", "diffusion", or "quasi2d"
    flameType = StringOption('premixed', ('diffusion', 'quasi2d'))

    #: Temperature of the unburned fuel/air mixture for premixed flames [K].
    Tu = FloatOption(300, label='Reactant Temperature', filter=_isPremixed)

    #: Temperature of the fuel mixture for diffusion flames [K].
    Tfuel = FloatOption(300, label='Fuel Temperature', filter=_isDiffusion)

    #: Temperature of the oxidizer mixture for diffusion flames [K].
    Toxidizer = FloatOption(300, label='Oxidizer Temperature', filter=_isDiffusion)

    #: Molar composition of the fuel mixture.
    fuel = Option("CH4:1.0", label='Molar Fuel Composition')

    #: Molar composition of the oxidizer mixture.
    oxidizer = Option("N2:3.76, O2:1.0", label='Molar Oxidizer Composition')

    #: Equivalence ratio of the fuel/air mixture for premixed flames.
    equivalenceRatio = FloatOption(0.75, min=0, label='Equivalence Ratio',
        filter=_isPremixed)

    #: Molar composition of the flow opposite the premixed reactant
    #: stream, if different from the equilibrium composition
    counterflow = StringOption(None, label='Composition of counterflow',
                               filter=_isPremixed, level=1)

    #: Temperature of the flow opposite the premixed reactant stream, if
    #: different from the equilibrium temperature
    Tcounterflow = FloatOption(None, label='Temperature of counteflow',
                               filter=_isPremixed, level=1)

    #: Adjust the composition of the counterflow stream so that the components
    #: are at equilibrium. This option specifies the property pair to hold
    #: constant during equilibration, or *False* to skip equilibration. The
    #: boundary condition is not consistent if this mixture has reactions that
    #: are proceeding at finite rates. For diffusion flames, this option is applied
    #: to the state oxidizer stream.
    equilibrateCounterflow = Option('TP', ('HP','UV','SV','TV','SP',False),
                                    label='Equilibrate specified counterflow mixture',
                                    level=1)

    #: Thermodynamic pressure [Pa]
    pressure = FloatOption(101325, min=0)

    #: Number of points in the initial uniform grid.
    nPoints = IntegerOption(100, level=1)

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
    centerWidth = FloatOption(0.000, level=1)

    #: The width of the slope away from the central plateau in the
    #: initial profile [m]. Recommended value for premixed flames:
    #: 5e-4. Recommended value for diffusion flames: 1e-3.
    slopeWidth = FloatOption(0.0005, level=1)

    #: Number of times to run the generated profile through a low pass
    #: filter before starting the simulation.
    smoothCount = IntegerOption(4, level=2)

    #: True if initial profiles for x,T,U,V and Y are given
    haveProfiles = BoolOption(False, level=3)

    x = Option(None, level=3)
    T = Option(None, level=3)
    U = Option(None, level=3)
    Y = Option(None, level=3)
    V = Option(None, level=3)


class WallFlux(Options):
    #: Reference temperature for the wall heat flux
    Tinf = FloatOption(300, level=2)

    #: Conductance of the wall [W/m^2-K]
    Kwall = FloatOption(100, level=2)


class Ignition(Options):
    """
    Parameters for an artificial heat release rate function which can
    be used to simulate ignition. The heat release rate is a step
    function in time with a Gaussian spatial distribution.
    """

    #: Beginning of the external heat release rate pulse [s].
    tStart = FloatOption(0, level=2)

    #: Duration of the external heat release rate pulse [s].
    duration = FloatOption(1e-3, level=2)

    #: Integral amplitude of the pulse [W/m^2].
    energy = FloatOption(0, level=2)

    #: Location of the center of the pulse [m].
    center = FloatOption(0, level=2)

    #: Characteristic width (standard deviation) of the pulse [m].
    stddev = FloatOption(1e-4, level=2)


class StrainParameters(Options):
    """
    Parameters defining the strain rate as a function of time.

    The strain rate changes linearly from *initial* to *final* over a
    period of *dt* seconds, starting at *tStart*.
    """
    initial = FloatOption(400)  #: Initial strain rate [1/s]
    final = FloatOption(400)  #: final strain rate [1/s]
    tStart = FloatOption(0.000)  #: time at which strain rate starts to change [s]
    dt = FloatOption(0.002)  #: time period over which strain rate changes [s]

    #: A list of strain rates to use for a series of sequential
    #: integration periods (see :func:`~ember.utils.multirun`), with steady-state
    #: profiles generated for each strain rate before proceeding to the next.
    #: A typical list of strain rates to use::
    #:
    #:     rates = [9216, 7680, 6144, 4608, 3840, 3072, 2304, 1920, 1536,
    #:              1152, 960, 768, 576, 480, 384, 288, 240, 192, 144, 120,
    #:              96, 72, 60, 48, 36, 30, 24, 18, 15, 12]
    rates = Option(None, level=3)

    #: Specify the strain rate as a function of time, using any callable
    #! Python object.
    function = Option(None, level=2)


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

    xInitial = FloatOption(0.0025, filter=_isSymmetric)  #:
    xFinal = FloatOption(0.0025, filter=_isSymmetric)  #:
    dt = FloatOption(0.01, filter=_isSymmetric)  #:
    tStart = FloatOption(0, filter=_isSymmetric)  #:

    proportionalGain = FloatOption(10, filter=_isSymmetric)  #:
    integralGain = FloatOption(800, filter=_isSymmetric)  #:


class Times(Options):
    """
    Paremeters controlling integrator timesteps and frequency of
    output profiles.
    """

    #: Integrator start time.
    tStart = FloatOption(0, level=1)

    #: Timestep used for operator splitting.
    globalTimestep = FloatOption(2e-5, level=1)

    #: Control for timestep used by the diffusion integrator.
    #: Actual timestep will be this multiplier times the stability
    #: limit an explicit integrator.
    diffusionTimestepMultiplier = FloatOption(10, level=2)

    #: Maximum amount of time before regridding / adaptation.
    regridTimeInterval = FloatOption(100, level=2)

    #: Maximum number of timesteps before regridding / adaptation.
    regridStepInterval = IntegerOption(20, level=2)

    #: Maximum number of steps between storing integral flame properties.
    outputStepInterval = IntegerOption(1, level=2)

    #: Maximum time between storing integral flame properties.
    outputTimeInterval = FloatOption(1e-5, level=2)

    #: Maximum number of timesteps before writing flame profiles.
    profileStepInterval = IntegerOption(1000, level=2)

    #: Maximum time between writing flame profiles.
    profileTimeInterval = FloatOption(1e-3, level=1)

    #: Number of timesteps between writing profNow.h5
    currentStateStepInterval = IntegerOption(20, level=2)

    #: Number of timesteps between checks of the steady-state
    #: termination conditions.
    terminateStepInterval = IntegerOption(10, level=2)


class CvodeTolerances(Options):
    """ Tolerances for the CVODE chemistry integrator """

    #: Relative tolerance for each state variable
    relativeTolerance = FloatOption(1e-6, level=2, filter=_usingCvode)

    #: Absolute tolerance for U velocity
    momentumAbsTol = FloatOption(1e-7, level=2, filter=_usingCvode)

    #: Absolute tolerance for T
    energyAbsTol = FloatOption(1e-8, level=2, filter=_usingCvode)

    #: Absolute tolerance for species mass fractions.
    speciesAbsTol = FloatOption(1e-13, level=2, filter=_usingCvode)

    #: Minimum internal timestep
    minimumTimestep = FloatOption(1e-18, level=2, filter=_usingCvode)


class QssTolerances(Options):
    """ Tolerances for the QSS chemistry integrator """

    #: Accuracy parameter for determining the next timestep.
    epsmin = FloatOption(2e-2, level=2, filter=_usingQss)

    #: Accuracy parameter for repeating timesteps
    epsmax = FloatOption(1e1, level=2, filter=_usingQss)

    #: Minimum internal timestep
    dtmin = FloatOption(1e-16, level=2, filter=_usingQss)

    #: Maximum internal timestep
    dtmax = FloatOption(1e-6, level=2, filter=_usingQss)

    #: Number of corrector iterations per timestep.
    iterationCount = IntegerOption(1, level=2, filter=_usingQss)

    #: Absolute threshold for including each component in the accuracy
    #: tests.
    abstol = FloatOption(1e-11, level=2, filter=_usingQss)

    #: Lower limit on the value of each state vector component.
    minval = FloatOption(1e-60, level=2, filter=_usingQss)

    #: Enable convergence-based stability check on timestep. Not
    #: enabled unless *iterationCount* >= 3.
    stabilityCheck = BoolOption(False, level=2, filter=_usingQss)


class Debug(Options):
    """ Control of verbose debugging output """

    #: Addition / removal of internal grid points.
    adaptation = BoolOption(False, level=2)

    #: Addition / removal of boundary grid points.
    regridding = BoolOption(True, level=2)

    #: Print current time after each global timestep.
    timesteps = BoolOption(True, level=2)

    #: Enable extensive timestep debugging output. Automatically
    #: enabled for debug builds.
    veryVerbose = BoolOption(False, level=2)

    #: Print information about the flame radius feedback controller.
    flameRadiusControl = BoolOption(False, level=2)

    #: Grid point to print debugging information about at *sourceTime*.
    sourcePoint = IntegerOption(-1, level=3)

    #: Time at which to print extensive debugging information about
    #: the source term at j = *sourcePoint*, then terminate.
    sourceTime = FloatOption(0.0, level=3)

    #: Time at which to start saving intermediate integrator profiles when
    #: OutputFiles.debugIntegratorStages is True
    startTime = FloatOption(0.0, level=3)

    #: Time at which to stop saving intermediate integrator profiles when
    #: OutputFiles.debugIntegratorStages is True
    stopTime = FloatOption(100.0, level=3)

class OutputFiles(Options):
    """ Control the contents of the periodic output files """

    #: Include the heat release rate as a function of space
    heatReleaseRate = BoolOption(True, level=2)

    #: Include the reaction / diffusion / convection contributions to
    #: the net time derivative of each state variable
    timeDerivatives = BoolOption(True, level=2)

    #: Include variables such as transport properties and grid
    #: parameters that can be recomputed from the state variables.
    extraVariables = BoolOption(False, level=2)

    #: Include other miscellaneous variables
    auxiliaryVariables = BoolOption(False, level=2)

    #: Used to generate a continuous sequence of output files after
    #: restarting the code.
    firstFileNumber = IntegerOption(0, level=2)

    #: Generate ``profNNNNNN.h5`` files
    saveProfiles = BoolOption(True, level=1)

    #: Write profiles after each stage of the split integrator.
    debugIntegratorStages = BoolOption(False, level=2)


class TerminationCondition(Options):
    r"""
    Integrate until either *tEnd* is reached or a criterion based on
    *measurement* is satisfied. The *measurement*-based check is not enabled
    until *tMin* is reached.

    - If `measurement == None`, integration will proceed to *tEnd*.

    - If `measurement == 'Q'`, integration will terminate when the the heat
      release rate reaches a steady-state value to within *tolerance* (RMS)
      over a time period of *steadyPeriod*, or the mean heat release rate over
      *steadyPeriod* is less than *abstol*.

    - If `measurement == 'dTdt', integration will terminate when
      `||1/T * dT/dt|| / sqrt(nPoints)`  is less than `dTdt_tol`.
    """

    tEnd = FloatOption(0.8)  #:

    measurement = Option("Q", (None,'dTdt'))  #:
    tolerance = FloatOption(1e-4, level=2)  #:
    abstol = FloatOption(0.5, min=0, level=2)  #:
    steadyPeriod = FloatOption(0.002, min=0, level=1)  #:
    tMin = FloatOption(0.0, level=1)
    dTdtTol = FloatOption(10.0)


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
        return ConcreteConfig(self)

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

        # Make sure that the mechanism file actually works and contains the
        # specified fuel and oxidizer species
        gas = Cantera.IdealGasMix(self.chemistry.mechanismFile.value,
                                  id=self.chemistry.phaseID.value)
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

    def run(self):
        """
        Run the simulation using the parameters set in this Config.

        If a list strain rates is provided by the field
        :attr:`strainParameters.rates`, a sequence of flame simulations at the
        given strain rates will be run. Otherwise, a single simulation will be
        run.

        If the script which calls this function is passed the argument
        *validate*, then the configuration will be checked for errors and the
        script will exit without running the simulation.
        """
        if len(sys.argv) > 1 and sys.argv[1].lower() == 'validate':
            # Validate the configuration and exit
            self.validate()
            return

        concrete = self.evaluate()
        if self.strainParameters.rates:
            concrete.multirun()
        else:
            concrete.run()


class ConcreteConfig(_ember.ConfigOptions):
    """
    Same structure as class Config, but all the Option objects are replaced
    with their actual values, and these values are propagated to an underlying
    C++ object as necessary.
    """
    def __init__(self, config):
        super(ConcreteConfig, self).__init__()
        self.original = config

        for name, opts in config.__dict__.iteritems():
            if isinstance(opts, Options):
                group = opts.__class__()
                setattr(self, name, group)

                for name in dir(opts):
                    opt = getattr(opts, name)
                    if isinstance(opt, Option):
                        setattr(group, name, opt.value)
            elif opts is None:
                setattr(self, name, None)

        self.gas = Cantera.IdealGasMix(self.chemistry.mechanismFile,
                                       self.chemistry.phaseID)

        if self.initialCondition.restartFile:
            self.readInitialCondition(self.initialCondition.restartFile)
        elif not self.initialCondition.haveProfiles:
            self.generateInitialCondition()
        self.apply_options()

    def readInitialCondition(self, restartFile):
        """
        Read the initial profiles for temperature, species mass fractions, and
        velocity from the specified input file.
        """
        IC = self.initialCondition
        data = utils.HDFStruct(restartFile)
        IC.x = data.x
        IC.Y = data.Y
        IC.T = data.T
        IC.U = data.U
        IC.V = data.V

        IC.haveProfiles = True
        if (IC.fuel or IC.oxidizer or IC.Tfuel or IC.Toxidizer or IC.equivalenceRatio
            or IC.Tcounterflow or IC.counterflow):
            self.setBoundaryValues(IC.T, IC.Y)

    def setBoundaryValues(self, T, Y):
        IC = self.initialCondition
        jm = (IC.nPoints-1) // 2
        gas = self.gas

        if IC.flameType == 'premixed':
            # Reactants
            reactants = utils.calculateReactantMixture(gas, IC.fuel, IC.oxidizer,
                                                       IC.equivalenceRatio)
            gas.set(X=reactants, T=IC.Tu, P=IC.pressure)
            rhou = gas.density()
            Yu = gas.massFractions()

            # Products
            gas.equilibrate('HP')
            if IC.Tcounterflow is None:
                Tb = gas.temperature()
            else:
                Tb = IC.Tcounterflow
            if IC.counterflow is None:
                Yb = gas.massFractions()
            else:
                gas.set(X=IC.counterflow, T=Tb, P=IC.pressure)
                Yb = gas.massFractions()

            if IC.equilibrateCounterflow:
                gas.set(Y=Yb, T=Tb, P=IC.pressure)
                gas.equilibrate(IC.equilibrateCounterflow)
                Yb = gas.massFractions()
                Tb = gas.temperature()

            if self.general.unburnedLeft:
                T[0] = IC.Tu
                Y[:,0] = Yu
                T[-1] = Tb
                Y[:,-1] = Yb
            else:
                T[0] = Tb
                Y[:,0] = Yb
                T[-1] = IC.Tu
                Y[:,-1] = Yu

        elif IC.flameType == 'diffusion':
            # Fuel
            gas.set(X=IC.fuel, T=IC.Tfuel, P=IC.pressure)
            Yfuel = gas.massFractions()

            # Oxidizer
            gas.set(X=IC.oxidizer, T=IC.Toxidizer, P=IC.pressure)
            if IC.equilibrateCounterflow:
                gas.equilibrate(IC.equilibrateCounterflow)

            Yoxidizer = gas.massFractions()
            Toxidizer = gas.temperature()

            if self.general.fuelLeft:
                T[0] = IC.Tfuel
                Y[:,0] = Yfuel
                T[-1] = Toxidizer
                Y[:,-1] = Yoxidizer
            else:
                T[0] = Toxidizer
                Y[:,0] = Yoxidizer
                T[-1] = IC.Tfuel
                Y[:,-1] = Yfuel

            gas.set(Y=Y[:,0], T=T[0], P=IC.pressure)
            rhou = gas.density()  # arbitrary for diffusion flame

        return rhou

    def generateInitialCondition(self):
        """
        Generate initial profiles for temperature, species mass fractions, and
        velocity using the specified fuel and oxidizer compositions and flame
        configuration parameters.
        """
        IC = self.initialCondition
        N = IC.nPoints
        gas = self.gas

        xLeft = (0.0 if self.general.twinFlame or self.general.curvedFlame
                 else IC.xLeft)

        x = np.linspace(xLeft, IC.xRight, N)
        T = np.zeros(N)
        Y = np.zeros((self.gas.nSpecies(), N))
        V = np.zeros(N)

        jm = (IC.nPoints-1) // 2

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

        rhou = self.setBoundaryValues(T, Y)

        if IC.flameType == 'premixed':
            T[jm] = IC.Tu
            # Diluent in the middle
            gas.set(X=IC.oxidizer, T=IC.Tu, P=IC.pressure)
            Y[:,jm] = gas.massFractions()

        elif IC.flameType == 'diffusion':
            # Stoichiometric mixture at the center
            IC.equivalenceRatio = 1.0
            products = utils.calculateReactantMixture(gas, IC.fuel, IC.oxidizer,
                                                      IC.equivalenceRatio)
            gas.set(X=products, T=0.5*(IC.Tfuel+IC.Toxidizer), P=IC.pressure)
            gas.equilibrate('HP')
            T[jm] = gas.temperature()
            Y[:,jm] = gas.massFractions()


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

        rho = np.zeros(N)
        U = np.zeros(N)
        if self.strainParameters.function:
            a0 = self.strainParameters.function(self.times.tStart)
        else:
            a0 = self.strainParameters.initial
        for j in range(N):
            gas.set(Y=Y[:,j], T=T[j], P=IC.pressure)
            rho[j] = gas.density()
            U[j] = a0 * np.sqrt(rhou/rho[j])

        for _ in range(2):
            utils.smooth(U)

        if self.general.twinFlame or self.general.curvedFlame:
            # Stagnation point at x = 0
            V[0] = 0
            for j in range(1, N):
                V[j] = V[j-1] - rho[j]*U[j]*(x[j] - x[j-1])
        else:
            # Put the stagnation point on the hot side of the flame
            if T[-1] > T[0]:
                jz = 3 * N // 4
            else:
                jz = N // 4
            V[jz] = 0
            for j in range(jz+1, N):
                V[j] = V[j-1] - rho[j]*U[j]*(x[j] - x[j-1])
            for j in range(jz-1, -1, -1):
                V[j] = V[j] + rho[j-1]*U[j-1]*(x[j] - x[j-1])

        IC.x = x
        IC.Y = Y
        IC.T = T
        IC.U = U
        IC.V = V
        IC.haveProfiles = True

    def run(self):
        """
        Run a single flame simulation using the parameters set in this Config.
        """
        confString = self.original.stringify()

        if not os.path.isdir(self.paths.outputDir):
            os.makedirs(self.paths.outputDir, 0755)
        confOutPath = os.path.join(self.paths.outputDir, 'config')
        if (os.path.exists(confOutPath)):
            os.unlink(confOutPath)
        confOut = file(confOutPath, 'w')
        confOut.write(confString)

        solver = _ember.FlameSolver(self)
        solver.initialize()
        done = 0
        while not done:
            done = solver.step()
        solver.finalize()
        return solver

    def multirun(self):
        """
        Run a sequence of flame simulations at different strain rates using
        the parameters set in this Config object. The list of strain rates is
        defined by the configuration field :attr:`strainParameters.rates`.
        """

        confString = self.original.stringify()
        strainRates = self.strainParameters.rates
        if not strainRates:
            print 'No strain rate list specified'
            return

        self.strainParameters.rates = None
        _logFile = file(self.paths.logFile, 'w')
        def log(message):
            _logFile.write(message)
            _logFile.write('\n')
            _logFile.flush()

        if not os.path.exists(self.paths.outputDir):
            os.mkdir(self.paths.outputDir, 0755)

        self.strainParameters.initial = strainRates[0]

        Q = []
        Sc = []
        xFlame = []

        for a in strainRates:
            restartFile = 'prof_eps%04i' % a
            historyFile = 'out_eps%04i' % a
            configFile = 'conf_eps%04i' % a

            restartPath = os.path.join(self.paths.outputDir, restartFile)
            historyPath = os.path.join(self.paths.outputDir, historyFile)
            configPath = os.path.join(self.paths.outputDir, configFile)

            if os.path.exists(restartPath) and os.path.exists(historyPath):
                # If the output files already exist, we simply retrieve the
                # integral flame properties from the existing profiles and
                # advance to the next strain rate.

                log('Skipping run at strain rate a = %g'
                    ' because the output file "%s" already exists.' % (a, restartFile))

                # Compute integral properties using points from the last half
                # of the termination-check period
                data = HDFStruct(historyPath)
                mask = data.t > data.t[-1] - 0.5*self.terminationCondition.steadyPeriod
                if not any(mask):
                    log('Warning: old data file did not contain data'
                        ' spanning the requested period.')
                    mask = data.t > 0.5*data.t[-1]

                Q.append(np.mean(data.Q[mask]))
                Sc.append(np.mean(data.Sc[mask]))
                xFlame.append(np.mean(data.xFlame[mask]))
                del data

            else:
                # Data is not already present, so run the flame solver for this strain rate

                log('Beginning run at strain rate a = %g s^-1' % a)
                confOut = file(configPath, 'w')
                confOut.write(confString)

                self.strainParameters.initial = a
                self.strainParameters.final = a
                self.paths.logFile = os.path.join(self.paths.outputDir, 'log-eps%04i.txt' % a)
                solver = _ember.FlameSolver(self)
                t1 = time.time()
                solver.initialize()
                done = 0
                while not done:
                    done = solver.step()
                solver.finalize()
                t2 = time.time()

                log('Completed run at strain rate a = %g s^-1' % a)
                log('Integration took %.1f seconds.' % (t2-t1))

                solver.writeStateFile(restartFile, False, False)
                solver.writeTimeseriesFile(historyFile)
                tRun = np.array(solver.timeVector)
                QRun = np.array(solver.heatReleaseRate)
                ScRun = np.array(solver.consumptionSpeed)
                xFlameRun = np.array(solver.flamePosition)

                # Compute integral properties using points from the last half
                # of the termination-check period
                mask = tRun > tRun[-1] - 0.5*self.terminationCondition.steadyPeriod
                Q.append(np.mean(QRun[mask]))
                Sc.append(np.mean(ScRun[mask]))
                xFlame.append(np.mean(xFlameRun[mask]))

            self.initialCondition.restartFile = restartPath

            # Sort by strain rate:
            strainRates, Q, Sc, xFlame = map(list, zip(*sorted(zip(strainRates, Q, Sc, xFlame))))

            integralFile = os.path.join(self.paths.outputDir, "integral.h5")
            if os.path.exists(integralFile):
                os.unlink(integralFile)
            data = h5py.File(integralFile)
            data['a'] = strainRates
            data['Q'] = Q
            data['Sc'] = Sc
            data['xFlame'] = xFlame
            data.close()
