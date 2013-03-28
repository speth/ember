#pragma once

#include "mathUtils.h"
#include "perfTimer.h"
#include "readConfig.h"

//! Implements operator split integration for reacting flow problems.
/*!
 *  Solves systems of equations with convection, reaction, and diffusion terms
 *  using either Strang splitting or rebalanced splitting.
 */
class SplitSolver
{
public:
    SplitSolver() {}
    virtual ~SplitSolver() {}

    //! Set configuration options needed by the solver.
    void setOptions(const ConfigOptions& _options);

    //! Set the size of all internal arrays when the problem size changes.
    void resize(index_t nRows, index_t nCols);

    //! Compute estimated time derivatives for each split term based on deltas
    //! for each integrator phase.
    void calculateTimeDerivatives(double dt);

    int step(); //!< Take one global timestep

    //! Take actions at the start of the time step. Used by derived classes to
    //! perform actions that need to be done at the start of a time step,
    //! before the splitting constants are computed.
    virtual void setupStep() = 0;

    //! Take actions at the end of a time step. Used by derived classes to
    //! perform any actions that need to be done at the end of the time step,
    //! after the time derivatives have been calculated based on the results
    //! of integrating the split equations.
    virtual int finishStep() = 0;

    //! Prepare split term integrators. Used by derived classes to assign
    //! split constants to the integrators for the individual terms and
    //! perform any other actions that should take place
    virtual void prepareIntegrators() = 0;

    //! Take one global timestep by integrating the split terms and computing
    //! the deltas for each integrator phase.
    void integrateSplitTerms(double t, double dt);

    //! @name Balanced Strang-split Integration Terms
    //! These functions should be overloaded by derived classes to implement the
    //! integration of the split terms.
    //! @{
    virtual void integrateConvectionTerms() = 0;
    virtual void integrateProductionTerms() = 0;
    virtual void integrateDiffusionTerms() = 0;
    //! @}

    //! Create prof.h5 file.
    //! @param fileName The name of the output file.
    //! @param errorFile Setting this to 'true' will generate a more verbose
    //!     output file
    //! @param updateDerivatives Setting this to 'true' will compute time
    //!        derivatives based on the current state. This flag should only
    //!        set if all of the solvers and auxillary arrays are sized
    //!        corresponding to the current state, e.g. after integration but
    //!        but before regridding.
    virtual void writeStateFile(const std::string& fileName="",
                                bool errorFile=false,
                                bool updateDerivatives=true) {}

    //! Current state
    dmatrix state;

    //! State at the start of the current integrator stage
    dmatrix startState;

    //! Change in state due to convection terms during the current time step
    dmatrix deltaConv;

    //! Change in state due to diffusion terms during the current time step
    dmatrix deltaDiff;

    //! Change in state due to production terms during the current time step
    dmatrix deltaProd;

    //! Estimated time derivatives of convection terms, based on the
    //! #deltaConv of the previous time step
    dmatrix ddtConv;

    //! Estimated time derivatives of the diffusion terms, based on the
    //! #deltaDiff of the previous time step
    dmatrix ddtDiff;

    //! Estimated time derivatives of production terms, based on the
    //! #deltaProd of the previous time step
    dmatrix ddtProd;

    //! Time derivative associated with terms taken as constant over a time
    //! step. These terms are held constant because they depend on state
    //! variables that cross the split dimensions, e.g. Soret mass fluxes.
    dmatrix ddtCross;

    //! Splitting constants for the convection terms
    dmatrix splitConstConv;

    //! Splitting constants for the diffusion terms
    dmatrix  splitConstDiff;

    //! Splitting constants for the production terms
    dmatrix splitConstProd;

    PerfTimer splitTimer; //!< Timer for amount of time used by splitting procedure.
    ConfigOptions options; //!< Options read from the input file.

    double tStart; //!< Integrator start time
    double tEnd; //!< Integrator termination time (upper bound)
    double tNow; //!< Current time reached by the integrator
    double tStageStart; //!< Start time for the active integrator stage
    double tStageEnd; //!< End time for the active integrator stage
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
