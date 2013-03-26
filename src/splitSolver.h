#pragma once

#include "mathUtils.h"
#include "perfTimer.h"
#include "readConfig.h"

class SplitSolver
{
public:
    SplitSolver() {}
    virtual ~SplitSolver() {}

    void setOptions(const ConfigOptions& _options);

    //! Set the size of all internal arrays when the problem size changes
    void resize(index_t nRows, index_t nCols);

    //! Compute estimated time derivatives for each split term based on deltas
    //! for each integrator phase.
    void calculateTimeDerivatives(double dt);

    int step(); //!< Take one global timestep

    virtual void setupStep() = 0;
    virtual int finishStep() = 0;
    virtual void prepareIntegrators() = 0;

    //! Take one global timestep by integrating the split terms and computing
    //! the deltas for each integrator phase.
    void integrateSplitTerms(double t, double dt);

    //! @name Balanced Strang-split Integration Terms
    //! These functions should be overloaded by derived classes to implement the
    //! integration of the split terms.
    //! @{
    virtual void integrateConvectionTerms(double tStart, double tEnd) = 0;
    virtual void integrateProductionTerms(double tStart, double tEnd) = 0;
    virtual void integrateDiffusionTerms(double tStart, double tEnd) = 0;
    //! @}

    //! Create prof.h5 file.
    //! @param fileName The name of the output file.
    //! @param errorFile Setting this to 'true' will generate a more verbose output file
    //! @param updateDerivatives Setting this to 'true' will compute time derivatives based
    //!        on the current state. This flag should only set if all of the solvers and auxillary
    //!        arrays are sized corresponding to the current state, e.g. after integration but
    //!        but before regridding.
    virtual void writeStateFile(const std::string& fileName="", bool errorFile=false,
                                bool updateDerivatives=true) {}

    // Current state
    dmatrix state;

    // State at the start of the current integrator stage
    dmatrix startState;

    // Changes in each state variable for each terms of the governing equations
    dmatrix deltaConv, deltaDiff, deltaProd;

    // Estimated time derivatives for each split term
    dmatrix ddtConv, ddtDiff, ddtProd, ddtCross;

    // Computed splitting constants for each term
    dmatrix splitConstConv, splitConstDiff, splitConstProd;

    PerfTimer splitTimer;
    ConfigOptions options; //!< Options read from the input file

    double tStart; //!< Integrator start time
    double tEnd; //!< Integrator termination time (upper bound)
    double tNow; //!< Current time reached by the integrator
    double t; //!< start time of the current global timestep
    double dt; //!< Global timestep used by the split solver

private:
    //! 'true' to use rebalanced splitting, 'false' to use regular Strang splitting
    bool useBalancedSplitting;

    //! @name Internal methods for each stages
    //! These methods do the neccessary actions before and after each split
    //! integration stage.
    //! @{
    void _integrateConvectionTerms(double tStart, double tEnd, int stage);
    void _integrateProductionTerms(double tStart, double tEnd, int stage);
    void _integrateDiffusionTerms(double tStart, double tEnd, int stage);
    //! @}
};
