#include "splitSolver.h"
#include "readConfig.h"

void SplitSolver::setOptions(const ConfigOptions& options_)
{
    options = options_;
    useBalancedSplitting = (options.splittingMethod == "balanced");
    dt = options.globalTimestep;
}

void SplitSolver::resize(index_t nRows, index_t nCols)
{
    state.resize(nRows, nCols);
    startState.resize(nRows, nCols);
    deltaConv.resize(nRows, nCols);
    deltaDiff.resize(nRows, nCols);
    deltaProd.resize(nRows, nCols);
    ddtConv.resize(nRows, nCols);
    ddtDiff.resize(nRows, nCols);
    ddtProd.resize(nRows, nCols);
    ddtCross.resize(nRows, nCols);
    splitConstConv.setZero(nRows, nCols);
    splitConstDiff.setZero(nRows, nCols);
    splitConstProd.setZero(nRows, nCols);
}

void SplitSolver::calculateTimeDerivatives(double dt)
{
    double split = useBalancedSplitting;

    if (useBalancedSplitting) {
        deltaConv -= 0.25 * dt * (ddtProd + ddtDiff + ddtCross);
        deltaDiff -= 0.25 * dt * (ddtProd + ddtConv + ddtCross);
        deltaProd -= 0.5 * dt * (ddtConv + ddtDiff + ddtCross);
    }

    ddtConv = deltaConv / dt + split * 0.75 * ddtConv;
    ddtDiff = deltaDiff / dt + split * 0.75 * ddtDiff;
    ddtProd = deltaProd / dt + split * 0.5 * ddtProd;
}

int SplitSolver::step()
{
    setupStep();

    if (useBalancedSplitting) {
        splitConstDiff = 0.25 * (ddtProd + ddtConv + ddtCross - 3 * ddtDiff);
        splitConstConv = 0.25 * (ddtProd + ddtDiff + ddtCross - 3 * ddtConv);
        splitConstProd = 0.5 * (ddtConv + ddtDiff + ddtCross - ddtProd);
    } else {
        splitConstDiff = ddtCross;
    }

    if (t == tStart && options.outputProfiles) {
        writeStateFile("", false, false);
    }

    if (options.debugIntegratorStages(tNow)) {
        writeStateFile((format("start_t%.6f") % t).str(), true, false);
    }

    prepareIntegrators();
    integrateSplitTerms(t, dt);
    calculateTimeDerivatives(dt);
    return finishStep();
}

void SplitSolver::integrateSplitTerms(double t, double dt)
{
    splitTimer.start();
    deltaConv.setZero();
    deltaProd.setZero();
    deltaDiff.setZero();
    splitTimer.stop();

    _integrateDiffusionTerms(t, t + 0.25*dt, 1); // step 1/4
    _integrateConvectionTerms(t, t + 0.5*dt, 1); // step 1/2
    _integrateDiffusionTerms(t + 0.25*dt, t + 0.5*dt, 2); // step 2/4
    _integrateProductionTerms(t, t + dt, 1); // full step
    _integrateDiffusionTerms(t + 0.5*dt, t + 0.75*dt, 3); // step 3/4
    _integrateConvectionTerms(t + 0.5*dt, t + dt, 2); // step 2/2
    _integrateDiffusionTerms(t + 0.75*dt, t + dt, 4); // step 4/4
}

void SplitSolver::_integrateDiffusionTerms(double tStart, double tEnd, int stage)
{
    tStageStart = tStart;
    tStageEnd = tEnd;
    assert(mathUtils::notnan(state));
    logFile.verboseWrite(format("diffusion terms %i/4...") % stage, false);
    startState = state;
    integrateDiffusionTerms();
    assert(mathUtils::notnan(state));
    deltaDiff += state - startState;
    if (stage && options.debugIntegratorStages(tNow)) {
        writeStateFile((format("diff%i_t%.6f") % stage % tNow).str(), true, false);
    }
}

void SplitSolver::_integrateProductionTerms(double tStart, double tEnd, int stage)
{
    tStageStart = tStart;
    tStageEnd = tEnd;
    logFile.verboseWrite("Source term...", false);
    assert(mathUtils::notnan(state));
    startState = state;
    integrateProductionTerms();
    assert(mathUtils::notnan(state));
    deltaProd += state - startState;
    if (options.debugIntegratorStages(tNow)) {
        writeStateFile((format("prod_t%.6f") % tNow).str(), true, false);
    }
}

void SplitSolver::_integrateConvectionTerms(double tStart, double tEnd, int stage)
{
    tStageStart = tStart;
    tStageEnd = tEnd;
    assert(stage == 1 || stage == 2);
    assert(mathUtils::notnan(state));
    logFile.verboseWrite(format("convection term %i/2...") % stage, false);
    startState = state;
    integrateConvectionTerms();
    assert(mathUtils::notnan(state));
    deltaConv += state - startState;
    if (options.debugIntegratorStages(tNow)) {
        writeStateFile((format("conv%i_t%.6f") % stage % tNow).str(), true, false);
    }
}
