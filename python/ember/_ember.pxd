from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool as cbool
from cpython import bool as pybool

import numpy as np
cimport numpy as np

cdef extern from "eigen_utils.h":
    cdef cppclass CxxEigenVec "dvec":
        pass

    cdef cppclass CxxEigenMatrix "dmatrix":
        pass

    # Technically, these functions return Map<...> objects, but it's better
    # if Cython doesn't know about that.
    cdef CxxEigenVec map_vector(double*, int, int)
    cdef CxxEigenMatrix map_matrix(double*, int, int, int, int)


cdef extern from "readConfig.h":
    cdef cppclass CxxConfigOptions "ConfigOptions":
        void setContinuityBC(string&)
        void setLogFile(string&)

        string outputDir, restartFile
        cbool overrideTu, overrideReactants, haveRestartFile
        cbool fixedBurnedVal, fixedLeftLoc, twinFlame
        int continuityBC
        cbool curvedFlame, unburnedLeft, fuelLeft
        string flameType
        int regridStepInterval, outputStepInterval, profileStepInterval
        int currentStateStepInterval, terminateStepInterval
        double regridTimeInterval, outputTimeInterval, profileTimeInterval
        double globalTimestep, diffusionTimestepMultiplier

        string splittingMethod, chemistryIntegrator
        double integratorRelTol, integratorMomentumAbsTol, integratorEnergyAbsTol
        double integratorSpeciesAbsTol, integratorMinTimestep

        double qss_epsmax, qss_epsmin, qss_dtmin, qss_dtmax
        int qss_iterationCount
        double qss_abstol, qss_minval
        cbool qss_stabilityCheck

        string gasMechanismFile, gasPhaseID, transportModel
        double transportThreshold

        string fuel, oxidizer
        double equivalenceRatio, pressure

        double Tu, Tfuel, Toxidizer
        double strainRateInitial, strainRateFinal, strainRateDt, strainRateT0

        cbool haveInitialProfiles
        CxxEigenVec x_initial, T_initial, U_initial
        CxxEigenMatrix Y_initial
        double rVzero_initial

        cbool quasi2d
        string interpFile

        cbool wallFlux
        double Tinf, Kwall

        double ignition_tStart, ignition_duration, ignition_energy
        double ignition_center, ignition_stddev

        double vtol, dvtol, rmTol, dampConst, gridMin, gridMax
        double uniformityTol, absvtol
        double boundaryTol, boundaryTolRm, unstrainedDownstreamWidth
        int addPointCount

        double tStart, tEnd
        cbool haveTStart

        int gridAlpha

        cbool outputAuxiliaryVariables, outputTimeDerivatives,
        cbool outputHeatReleaseRate, outputExtraVariables, outputProfiles

        int debugSourcePoint
        double debugSourceTime, debugStartTime, debugStopTime

        cbool terminateForSteadyQdot
        double terminationTolerance, terminationAbsTol, terminationPeriod

        int errorStopCount
        cbool stopIfError
        int nThreads

        int outputFileNumber
        cbool fileNumberOverride

        double centerGridMin

        cbool xFlameControl
        double xFlameInitial, xFlameFinal, xFlameDt, xFlameT0
        double xFlameIntegralGain, xFlameProportionalGain

        cbool xStagControl
        double xStag


cdef extern from "debugUtils.h":
    cdef cbool CxxSetDebugParameters "debugParameters::setParameters" (cbool, cbool, cbool, cbool, cbool)
    cdef cppclass CxxLogFile "LogFile":
        void open(string)
    cdef CxxLogFile CxxSingletonLogfile "logFile"

cdef extern from "flameSolver.h":
    cdef cppclass CxxFlameSolver "FlameSolver":
        void setOptions(CxxConfigOptions&)
        void initialize()
        void finalize()
        int step() except +


cdef class ConfigOptions:
    cdef CxxConfigOptions* opts

cdef class FlameSolver:
    cdef ConfigOptions options
    cdef CxxFlameSolver* solver
