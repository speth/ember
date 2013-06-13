from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool as cbool
from cpython import bool as pybool

import numpy as np
cimport numpy as np

cdef extern from "eigen_utils.h":
    cdef cppclass CxxEigenVecMap "VecMap":
        double& item "operator()" (int)
        int size()
        void resize(int)

    cdef cppclass CxxEigenVec "dvec":
        double& item "operator()" (int)
        int size()
        void resize(int)

    cdef cppclass CxxEigenMatrix "dmatrix":
        double& item "operator()" (int, int)
        int rows()
        int cols()
        void resize(int, int)

    cdef cppclass CxxEigenMatrixMap "MatrixMap":
        double& item "operator()" (int, int)
        int rows()
        int cols()
        void resize(int, int)


    # Technically, these functions return Map<...> objects, but it's better
    # if Cython doesn't know about that.
    cdef CxxEigenVec map_vector(double*, int, int)
    cdef CxxEigenMatrix map_matrix(double*, int, int, int, int)


cdef extern from "chemistry0d.h":
    cdef void CxxAddCanteraDirectory "Cantera::addDirectory" (string&)

cdef extern from "readConfig.h":
    cdef cppclass CxxConfigOptions "ConfigOptions":
        void setContinuityBC(string&)
        void setLogFile(string&)

        string outputDir
        cbool fixedBurnedVal, fixedLeftLoc, twinFlame
        int continuityBC
        cbool curvedFlame, unburnedLeft, fuelLeft
        string flameType
        int regridStepInterval, outputStepInterval, profileStepInterval
        int currentStateStepInterval, terminateStepInterval
        double regridTimeInterval, outputTimeInterval, profileTimeInterval
        double globalTimestep, diffusionTimestepMultiplier

        string splittingMethod, chemistryIntegrator, rateMultiplierFunctionType
        double integratorRelTol, integratorMomentumAbsTol, integratorEnergyAbsTol
        double integratorSpeciesAbsTol, integratorMinTimestep

        double qss_epsmax, qss_epsmin, qss_dtmin, qss_dtmax
        int qss_iterationCount
        double qss_abstol, qss_minval
        cbool qss_stabilityCheck

        string gasMechanismFile, gasPhaseID, transportModel, kineticsModel
        double transportThreshold

        double pressure

        string strainFunctionType
        double strainRateInitial, strainRateFinal, strainRateDt, strainRateT0

        CxxEigenVec x_initial, T_initial, U_initial, V_initial
        CxxEigenMatrix Y_initial

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

        double tStart, tEnd, tEndMin
        cbool haveTStart

        int gridAlpha

        cbool outputAuxiliaryVariables, outputTimeDerivatives,
        cbool outputHeatReleaseRate, outputExtraVariables, outputProfiles

        int debugSourcePoint
        double debugSourceTime, debugStartTime, debugStopTime
        cbool outputDebugIntegratorStages

        string terminationMeasurement
        double terminationTolerance, terminationAbsTol, terminationPeriod
        double termination_dTdtTol

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
        void write(string&)
    cdef CxxLogFile CxxSingletonLogfile "logFile"


cdef extern from "grid.h":
    cdef cppclass CxxOneDimGrid "OneDimGrid":
        CxxEigenVec x


cdef extern from "scalarFunction.h":
    cdef cppclass CxxScalarFunction "ScalarFunction":
        void setCoefficients(int, double*)
        double a(double)
        double dadt(double)


cdef extern from "flameSolver.h":
    cdef cppclass CxxFlameSolver "FlameSolver":
        void setOptions(CxxConfigOptions&) except +
        void initialize() except +
        void finalize() except +
        int step() nogil except +

        CxxScalarFunction* strainfunc
        CxxScalarFunction* rateMultiplierFunction

        double tNow, dt
        vector[double] timeVector
        vector[double] heatReleaseRate
        vector[double] consumptionSpeed
        vector[double] flamePosition
        double terminationCondition

        CxxOneDimGrid grid

        CxxEigenVecMap T
        CxxEigenVecMap U
        CxxEigenMatrixMap Y

        CxxEigenMatrix ddtDiff
        CxxEigenMatrix ddtConv
        CxxEigenMatrix ddtProd
        CxxEigenMatrix ddtCross
        CxxEigenVec drhodt
        CxxEigenVec rho
        CxxEigenVec jCorr
        CxxEigenVec sumcpj
        CxxEigenVec qDot
        CxxEigenMatrix wDot
        CxxEigenVec Wmx
        CxxEigenVec W
        CxxEigenVec mu
        CxxEigenVec k "lambda"
        CxxEigenVec cp
        CxxEigenMatrix cpSpec
        CxxEigenMatrix rhoD
        CxxEigenMatrix Dkt
        CxxEigenMatrix hk
        CxxEigenMatrix jFick
        CxxEigenMatrix jSoret


cdef class ConfigOptions:
    cdef CxxConfigOptions* opts

cdef class FlameSolver:
    cdef ConfigOptions options
    cdef CxxFlameSolver* solver
    cdef object strainFunction
    cdef object rateMultiplierFunction
    cdef object strainInterpPoints

    # usedby the GUI
    cdef public object lock
    cdef public object progress
