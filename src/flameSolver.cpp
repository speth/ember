#include "flameSolver.h"
#include "dataFile.h"
#include "scalarFunction.h"
#include "tbb_tools.h"

#include <boost/foreach.hpp>

#define foreach BOOST_FOREACH

FlameSolver::FlameSolver()
    : U(0, 0, Stride1X(1 ,1))
    , T(0, 0, Stride1X(1 ,1))
    , Y(0, 0, 0, StrideXX(1 ,1))
    , jCorrSolver(jCorrSystem)
    , strainfunc(NULL)
    , tbbTaskSched(tbb::task_scheduler_init::deferred)
{
}

FlameSolver::~FlameSolver()
{
    delete strainfunc;
}

void FlameSolver::setOptions(const ConfigOptions& _options)
{
    SplitSolver::setOptions(_options);
    tStart = options.tStart;
    tEnd = options.tEnd;

    gas.setOptions(_options);
    grid.setOptions(_options);
}

void FlameSolver::initialize(void)
{
    tbbTaskSched.initialize (options.nThreads);
    delete strainfunc;
    strainfunc = newStrainFunction(options);

    flamePosIntegralError = 0;
    terminationCondition = 1e10;

    // Cantera initialization
    gas.initialize();
    nSpec = gas.nSpec;
    nVars = nSpec + 2;
    W.resize(nSpec);
    gas.getMolecularWeights(W);

    // Get Initial Conditions
    loadProfile();

    // Interpolation data for quasi-2d problem
    if (options.quasi2d) {
        TInterp.reset(new BilinearInterpolator);
        vrInterp.reset(new BilinearInterpolator);
        vzInterp.reset(new BilinearInterpolator);
        TInterp->open(options.interpFile, "T", "r", "z");
        vrInterp->open(options.interpFile, "vr", "r", "z");
        vzInterp->open(options.interpFile, "vz", "r", "z");
        TInterp->get(0.05, 0.0);
    }

    grid.setSize(x.size());
    convectionSystem.setGas(gas);
    convectionSystem.setLeftBC(Tleft, Yleft);
    convectionSystem.setTolerances(options);

    for (size_t k=0; k<nVars; k++) {
        DiffusionSystem* term = new DiffusionSystem();
        TridiagonalIntegrator* integrator = new TridiagonalIntegrator(*term);
        integrator->resize(nPoints);
        diffusionTerms.push_back(term);
        diffusionSolvers.push_back(integrator);
    }
    if (options.wallFlux) {
        diffusionTerms[kEnergy].yInf = options.Tinf;
        diffusionTerms[kEnergy].wallConst = options.Kwall;
    }

    resizeAuxiliary();

    ddtConv.setZero();
    ddtDiff.setZero();
    ddtProd.setZero();

    updateChemicalProperties();
    calculateQdot();
    convectionSystem.utwSystem.updateContinuityBoundaryCondition(qDot, options.continuityBC);

    t = tStart;
    tOutput = t;
    tRegrid = t + options.regridTimeInterval;
    tProfile = t + options.profileTimeInterval;
    nTotal = 0;
    nRegrid = 0;
    nOutput = 0;
    nProfile = 0;
    nTerminate = 0;
    nCurrentState = 0;

    grid.updateValues();
    resizeAuxiliary();

    tFlamePrev = t;
    tNow = t;

    totalTimer.start();
}

void FlameSolver::setupStep()
{
    setupTimer.start();

    // Debug sanity check
    #ifndef NDEBUG
        bool error = false;
        for (size_t j=0; j<nPoints; j++) {
            if (T(j) < 295 || T(j) > 3000) {
                logFile.write(format(
                    "WARNING: Unexpected Temperature: T = %f at j = %i") % T(j) % j);
                error = true;
            }
        }
        if (error) {
            writeStateFile();
        }
    #endif

    // Reset boundary conditions to prevent numerical drift
    if (grid.leftBC == BoundaryCondition::FixedValue) {
        T(0) = Tleft;
        Y.col(0) = Yleft;
    }

    if (grid.rightBC == BoundaryCondition::FixedValue) {
        T(jj) = Tright;
        Y.col(jj) = Yright;
    }

    updateChemicalProperties();

    updateBC();
    if (options.xFlameControl) {
        update_xStag(t, true); // calculate the value of rVzero
    }
    convectionSystem.set_rVzero(rVzero);
    setupTimer.stop();

    // Set up solvers for split integration
    updateCrossTerms();
}

void FlameSolver::prepareIntegrators()
{
    splitTimer.resume();
    // Diffusion terms
    if (!options.quasi2d) {
        // Diffusion solvers: Energy and momentum
        diffusionTerms[kMomentum].B = rho.inverse();
        diffusionTerms[kEnergy].B = (rho * cp).inverse();

        diffusionTerms[kMomentum].D = mu;
        diffusionTerms[kEnergy].D = lambda;

        // Diffusion solvers: Species
        for (size_t k=0; k<nSpec; k++) {
            diffusionTerms[kSpecies+k].B =rho.inverse();
            diffusionTerms[kSpecies+k].D = rhoD.row(k);
        }
    } else {
        // Diffusion solvers: Energy and momentum
        diffusionTerms[kMomentum].B.setZero(nPoints);
        diffusionTerms[kEnergy].B.setZero(nPoints);

        diffusionTerms[kMomentum].D.setZero(nPoints);
        diffusionTerms[kEnergy].D.setZero(nPoints);

        // Diffusion solvers: Species
        for (size_t k=0; k<nSpec; k++) {
            DiffusionSystem& sys = diffusionTerms[kSpecies+k];
            sys.D = rhoD.row(k);
            for (size_t j = 0; j <= jj; j++) {
                sys.B[j] = 1 / (rho[j] * vzInterp->get(x[j], tNow));
            }
        }
    }

    setDiffusionSolverState(tNow);
    for (size_t i=0; i<nVars; i++) {
        diffusionTerms[i].splitConst = splitConstDiff.row(i);
    }

    // Production terms
    setProductionSolverState(tNow);
    for (size_t j=0; j<nPoints; j++) {
        sourceTerms[j].splitConst = splitConstProd.col(j);
    }

    // Convection terms
    setConvectionSolverState(tNow);
    dmatrix ddt = ddtConv + ddtDiff + ddtProd;
    if (options.splittingMethod == "balanced") {
        ddt += ddtCross;
    }

    dvec tmp = (W.inverse().matrix().transpose() * ddt.bottomRows(nSpec).matrix()).array();
    drhodt = - rho * (ddt.row(kEnergy).transpose() / T + tmp * Wmx);

    assert(mathUtils::notnan(drhodt));
    convectionSystem.setDensityDerivative(drhodt);
    convectionSystem.setSplitConstants(splitConstConv);
    convectionSystem.utwSystem.updateContinuityBoundaryCondition(qDot, options.continuityBC);
    splitTimer.stop();
}

int FlameSolver::finishStep()
{
    logFile.verboseWrite("done!");

    // *** End of Strang-split integration step ***
    correctMassFractions();

    t = tNow + dt;
    tNow += dt;

    nOutput++;
    nRegrid++;
    nProfile++;
    nTerminate++;
    nCurrentState++;

    if (debugParameters::debugTimesteps) {
        int nSteps = convectionSystem.getNumSteps();
        logFile.write(format("t = %8.6f (dt = %9.3e) [C: %i]") % t % dt % nSteps);
    }

    setupTimer.resume();
    if (t + 0.5 * dt > tOutput || nOutput >= options.outputStepInterval) {
        calculateQdot();
        timeVector.push_back(t);
        timestepVector.push_back(dt);
        heatReleaseRate.push_back(getHeatReleaseRate());
        consumptionSpeed.push_back(getConsumptionSpeed());
        flamePosition.push_back(getFlamePosition());
        strainRateVector.push_back(strainfunc->a(t));
        dStrainRateVector.push_back(strainfunc->dadt(t));

        tOutput = t + options.outputTimeInterval;
        nOutput = 0;
    }

    // Periodic check for terminating the integration (based on steady heat
    // release rate, etc.) Quit now to skip grid adaptation on the last step
    if (t >= tEnd) {
        return 1;
    } else if (nTerminate >= options.terminateStepInterval) {
        nTerminate = 0;
        if (checkTerminationCondition()) {
            return 1;
        }
    }

    // Save the current integral and profile data in files that are
    // automatically overwritten, and save the time-series data (out.h5)
    if (nCurrentState >= options.currentStateStepInterval) {
        calculateQdot();
        nCurrentState = 0;
        writeTimeseriesFile("out");
        writeStateFile("profNow");
    }

    // *** Save flame profiles
    if (t + 0.5 * dt > tProfile || nProfile >= options.profileStepInterval) {
        if (options.outputProfiles) {
            writeStateFile();
        }
        tProfile = t + options.profileTimeInterval;
        nProfile = 0;
    }
    setupTimer.stop();

    if (t > tRegrid || nRegrid >= options.regridStepInterval) {
        regridTimer.start();
        tRegrid = t + options.regridTimeInterval;
        nRegrid = 0;

        // If the left grid point moves, a new boundary value for rVzero
        // needs to be calculated from the mass flux V on the current grid points
        dvec x_prev = grid.x;
        convectionSystem.evaluate();
        dvec V_prev = convectionSystem.V;

        // dampVal sets a limit on the maximum grid size
        grid.dampVal.resize(grid.x.rows());
        for (size_t j=0; j<nPoints; j++) {
            double num = std::min(mu[j],lambda[j]/cp[j]);
            for (size_t k=0; k<nSpec; k++) {
                if (rhoD(k,j) > 0) {
                    num = std::min(num, rhoD(k,j));
                }
            }
            double den = std::max(std::abs(rho[j]*strainfunc->a(t)), 1e-100);
            grid.dampVal[j] = sqrt(num/den);
        }
        dvec dampVal_prev = grid.dampVal;

        vector<dvector> currentSolution;
        rollVectorVector(currentSolution, state);
        rollVectorVector(currentSolution, ddtConv);
        rollVectorVector(currentSolution, ddtDiff);
        rollVectorVector(currentSolution, ddtProd);

        grid.nAdapt = nVars;
        if (options.quasi2d) {
            // do not change grid extents in this case
        } else if (strainfunc->a(tNow) == 0) {
            calculateQdot();
            grid.regridUnstrained(currentSolution, qDot);
        } else {
            grid.regrid(currentSolution);
        }

        // Interpolate dampVal onto the modified grid
        grid.dampVal = mathUtils::interp1(x_prev, dampVal_prev, grid.x);

        grid.adapt(currentSolution);

        // Perform updates that are necessary if the grid has changed
        if (grid.updated) {
            logFile.write(format("Grid size: %i points.") % nPoints);
            resizeMappedArrays();
            unrollVectorVector(currentSolution, state, 0);
            unrollVectorVector(currentSolution, ddtConv, 1);
            unrollVectorVector(currentSolution, ddtDiff, 2);
            unrollVectorVector(currentSolution, ddtProd, 3);
            correctMassFractions();

            // Update the mass flux (including the left boundary value)
            rVzero = mathUtils::interp1(x_prev, V_prev, grid.x[0]);
            convectionSystem.utwSystem.V = mathUtils::interp1(x_prev, V_prev, grid.x);

            // Allocate the solvers and arrays for auxiliary variables
            resizeAuxiliary();
            calculateQdot();

            if (debugParameters::debugAdapt || debugParameters::debugRegrid) {
                writeStateFile("postAdapt", false, false);
            }
            grid.updated = false;
        }
        regridTimer.stop();
    }

    if (nTotal % 10 == 0) {
        printPerformanceStats();
    }
    nTotal++;
    return 0;
}

void FlameSolver::finalize()
{
    calculateQdot();
    writeTimeseriesFile("out");

    // *** Integration has reached the termination condition
    if (options.outputProfiles) {
        writeStateFile();
    }

    totalTimer.stop();
    printPerformanceStats();
    logFile.write(format("Runtime: %f seconds.") % totalTimer.getTime());
}

bool FlameSolver::checkTerminationCondition(void)
{
    if (options.terminateForSteadyQdot) {
        int j1 = mathUtils::findLast(timeVector < (tNow - options.terminationPeriod));
        if (j1 == -1) {
            logFile.write(format(
                    "Continuing integration: t (%8.6f) < terminationPeriod (%8.6f)") %
                    (tNow-timeVector[0]) % options.terminationPeriod);
            return false;
        }

        int j2 = timeVector.size()-1;
        double qMean = mathUtils::mean(heatReleaseRate,j1,j2);
        double hrrError = 0;
        for (int j=j1; j<=j2; j++) {
            hrrError += pow(heatReleaseRate[j]-qMean, 2);
        }
        hrrError = sqrt(hrrError) / (j2-j1+1);
        terminationCondition = hrrError / abs(qMean);

        logFile.write(format(
            "Heat release rate RMS error = %6.3f%%. absolute error: %9.4e") %
            (hrrError/qMean*100) % hrrError);

        if (hrrError/abs(qMean) < options.terminationTolerance) {
            logFile.write("Terminating integration: "
                    "Heat release RMS variation less than relative tolerance.");
            return true;
        } else if (hrrError < options.terminationAbsTol) {
            logFile.write("Terminating integration: "
                    "Heat release rate RMS variation less than absolute tolerance.");
            return true;
        } else {
            logFile.write(format(
                "Continuing integration. t = %8.6f") % (tNow-timeVector[0]));
        }
    }
    return false;
}

void FlameSolver::writeStateFile
(const std::string& fileNameStr, bool errorFile, bool updateDerivatives)
{
    std::string filename;
    bool incrementFileNumber = false;

    if (fileNameStr.length() == 0) {
        // Determine the name of the output file (profXXXXXX.h5)
        incrementFileNumber = true;
        std::string prefix = (errorFile) ? "error" : "prof";
        filename = (format("%s/%s%06i.h5") %
            options.outputDir % prefix % options.outputFileNumber).str();
    } else {
        filename = (format("%s/%s.h5") % options.outputDir % fileNameStr).str();
    }

    if (errorFile) {
        logFile.write(format("Writing error output file: %s") % filename);
    } else {
        logFile.write(format("Writing output file: %s") % filename);
    }

    if (updateDerivatives) {
        updateChemicalProperties();
        convectionSystem.evaluate();
    }

    // Write the state data to the output file:
    DataFile outFile(filename, std::ios_base::trunc);
    outFile.writeScalar("t", tNow);
    outFile.writeVec("x", x);
    outFile.writeVec("T", T);
    outFile.writeVec("U", U);
    outFile.writeArray2D("Y", Y, true);
    outFile.writeVec("V", convectionSystem.V);
    outFile.writeScalar("P", options.pressure);
    outFile.writeScalar("gridAlpha", grid.alpha);
    outFile.writeScalar("a", strainfunc->a(tNow));
    outFile.writeScalar("dadt", strainfunc->dadt(tNow));
    outFile.writeScalar("fileNumber", options.outputFileNumber);

    if (options.outputHeatReleaseRate || errorFile) {
        outFile.writeVec("q", qDot);
        outFile.writeVec("rho", rho);
    }

    if (options.outputTimeDerivatives || errorFile) {
        outFile.writeVec("dUdtDiff", ddtDiff.row(kMomentum));
        outFile.writeVec("dUdtConv", ddtConv.row(kMomentum));
        outFile.writeVec("dUdtProd", ddtProd.row(kMomentum));

        outFile.writeVec("dTdtDiff", ddtDiff.row(kEnergy));
        outFile.writeVec("dTdtConv", ddtConv.row(kEnergy));
        outFile.writeVec("dTdtProd", ddtProd.row(kEnergy));
        outFile.writeVec("dTdtCross", ddtCross.row(kEnergy));

        outFile.writeArray2D("dYdtConv", ddtConv.bottomRows(nSpec), true);
        outFile.writeArray2D("dYdtDiff", ddtDiff.bottomRows(nSpec), true);
        outFile.writeArray2D("dYdtProd", ddtProd.bottomRows(nSpec), true);
        outFile.writeArray2D("dYdtCross", ddtCross.bottomRows(nSpec), true);

        outFile.writeVec("dWdt", convectionSystem.dWdt);
        outFile.writeVec("drhodt", drhodt);
    }

    if (options.outputAuxiliaryVariables || errorFile) {
        outFile.writeVec("sumcpj", sumcpj);
        outFile.writeScalar("Tleft", Tleft);
        outFile.writeVec("Yleft", Yleft);

        outFile.writeVec("dWdx", convectionSystem.utwSystem.dWdx);
        outFile.writeVec("dTdx", convectionSystem.utwSystem.dTdx);
    }

    if (options.outputExtraVariables || errorFile) {
        // These variables can be recomputed from the state variables
        outFile.writeArray2D("wdot", wDot, true);
        outFile.writeArray2D("rhoD", rhoD, true);
        outFile.writeVec("lambda", lambda);
        outFile.writeVec("cp", cp);
        outFile.writeVec("mu", mu);
        outFile.writeVec("Wmx", Wmx);
        outFile.writeVec("W", W);
        outFile.writeVec("cfp", grid.cfp);
        outFile.writeVec("cf", grid.cf);
        outFile.writeVec("cfm", grid.cfm);
        outFile.writeVec("hh", hh);
        outFile.writeVec("rphalf", grid.rphalf);
        outFile.writeArray2D("jFick", jFick, true);
        outFile.writeArray2D("jSoret", jSoret, true);
        outFile.writeVec("jCorr", jCorr);
    }

    if (incrementFileNumber) {
        options.outputFileNumber++;
    }
}

void FlameSolver::writeTimeseriesFile(const std::string& filename)
{
    DataFile outFile(options.outputDir+"/"+filename+".h5", std::ios_base::trunc);
    outFile.writeVector("t", timeVector);
    outFile.writeVector("dt", timestepVector);
    outFile.writeVector("Q", heatReleaseRate);
    outFile.writeVector("Sc", consumptionSpeed);
    outFile.writeVector("xFlame", flamePosition);
    outFile.writeVector("a", strainRateVector);
    outFile.writeVector("dadt", dStrainRateVector);
}

void FlameSolver::resizeAuxiliary()
{
    size_t nPointsOld = rho.size();
    grid.setSize(T.size());

    if (nPoints == nPointsOld && !grid.updated) {
        return; // nothing to do
    }

    resizeMappedArrays();

    ddtCross.topRows(2).setZero();
    ddtCross.bottomRows(nSpec) *= NaN;

    rho.setZero(nPoints);
    drhodt.setZero(nPoints);
    Wmx.resize(nPoints);
    mu.resize(nPoints);
    lambda.setZero(nPoints);
    cp.setZero(nPoints);
    jCorr.resize(nPoints);
    sumcpj.setZero(nPoints);
    qDot.resize(nPoints);
    cpSpec.resize(nSpec, nPoints);
    rhoD.resize(nSpec, nPoints);
    Dkt.resize(nSpec, nPoints);
    wDot.resize(nSpec, nPoints);
    hk.resize(nSpec, nPoints);
    jFick.setZero(nSpec, nPoints);
    jSoret.setZero(nSpec, nPoints);

    grid.jj = nPoints-1;
    grid.updateBoundaryIndices();

    if (nPoints > nPointsOld) {
        for (size_t j=nPointsOld; j<nPoints; j++) {
            useCVODE.push_back(options.chemistryIntegrator == "cvode");

            SourceSystem* system;
            if (useCVODE[j]) {
                system = new SourceSystemCVODE();
            } else {
                system = new SourceSystemQSS();
            }
            // initialize the new SourceSystem
            system->setGas(&gas);
            system->initialize(nSpec);
            system->setOptions(options);
            system->setTimers(&reactionRatesTimer, &thermoTimer, &jacobianTimer);
            system->setStrainFunction(strainfunc);
            system->setRhou(rhou);
            system->setPosition(j, x[j]);
            if (options.quasi2d) {
                system->setupQuasi2d(vzInterp, TInterp);
            }

            // Store the solver and system
            sourceTerms.push_back(system);
        }

    } else {
        // Delete the unwanted solvers and systems
        sourceTerms.erase(sourceTerms.begin()+nPoints, sourceTerms.end());
        useCVODE.erase(useCVODE.begin()+nPoints, useCVODE.end());
    }

    // Resize solution vector for diffusion systems / solvers
    for (size_t k=0; k<nVars; k++) {
        diffusionSolvers[k].resize(nPoints);
        diffusionTerms[k].setGrid(grid);
    }

    convectionSystem.setGrid(grid);
    convectionSystem.resize(nPoints, nSpec, state);
    convectionSystem.setLeftBC(Tleft, Yleft);

    if (options.quasi2d) {
        convectionSystem.setupQuasi2D(vzInterp, vrInterp);
    }

    // Resize the jCorr stabilizer
    jCorrSolver.resize(nPoints);
    jCorrSystem.setGrid(grid);

    // Set the grid position for each of the source solvers
    for (size_t j=0; j<nPoints; j++) {
        sourceTerms[j].setPosition(j, x[j]);
    }
}

void FlameSolver::resizeMappedArrays()
{
    resize(nVars, nPoints);
    remap(state, T, nPoints, kEnergy);
    remap(state, U, nPoints, kMomentum);
    remap(state, Y, nSpec, nPoints, kSpecies);
}

void FlameSolver::updateCrossTerms()
{
    assert(mathUtils::notnan(state));
    assert(mathUtils::notnan(rhoD));

    for (size_t j=0; j<jj; j++) {
        jCorr[j] = 0;
        for (size_t k=0; k<nSpec; k++) {
            jFick(k,j) = -0.5*(rhoD(k,j)+rhoD(k,j+1)) * ((Y(k,j+1)-Y(k,j))/hh[j]);
            jSoret(k,j) = -0.5*(Dkt(k,j)/T(j) + Dkt(k,j+1)/T(j+1))
                * (T(j+1)-T(j))/hh[j];
            jCorr[j] -= jFick(k,j) + jSoret(k,j);
        }
    }
    jCorr[jj] = 0;

    assert(mathUtils::notnan(jFick));
    assert(mathUtils::notnan(jSoret));
    assert(mathUtils::notnan(lambda));
    assert(mathUtils::notnan(rho));
    assert(mathUtils::notnan(cp));

    // Add a bit of artificial diffusion to jCorr to improve stability
    jCorrSolver.y = jCorr;
    jCorrSystem.B = lambda / (rho * cp); // treat as Le = 1
    jCorrSystem.D.setConstant(nPoints, 1.0);
    jCorrSystem.splitConst.setZero(nPoints);

    double dt = options.globalTimestep;
    for (size_t j=1; j<jj; j++) {
        dt = std::min(dt, options.diffusionTimestepMultiplier*dlj[j]*dlj[j]/(jCorrSystem.B[j]));
    }

    jCorrSolver.initialize(0, dt);
    jCorrSolver.integrateToTime(options.globalTimestep);
    assert(mathUtils::notnan(jCorrSolver.y));

    jCorr = jCorrSolver.y;

    // dYdt due to gradients in other species and temperature
    Eigen::Block<dmatrix> dYdtCross = ddtCross.bottomRows(nSpec);

    // dTdt due to gradients in species composition
    Eigen::Block<dmatrix, 1> dTdtCross = ddtCross.row(kEnergy);

    dYdtCross.col(0).setZero();
    dYdtCross.col(jj).setZero();
    for (size_t j=1; j<jj; j++) {
        sumcpj[j] = 0;
        for (size_t k=0; k<nSpec; k++) {
            dYdtCross(k,j) = -0.5 / (r[j] * rho[j] * dlj[j]) *
                (rphalf[j] * (Y(k,j) + Y(k,j+1)) * jCorr[j] -
                 rphalf[j-1] * (Y(k,j-1) + Y(k,j)) * jCorr[j-1]);
            dYdtCross(k,j) -= 1 / (r[j] * rho[j] * dlj[j]) *
                (rphalf[j] * jSoret(k,j) - rphalf[j-1] * jSoret(k,j-1));
            sumcpj[j] += 0.5*(cpSpec(k,j) + cpSpec(k,j+1)) / W[k] *
                (jFick(k,j) + jSoret(k,j) + 0.5 * (Y(k,j) + Y(k,j+1)) * jCorr[j]);
        }
        double dTdx = cfm[j] * T(j-1) + cf[j] * T(j) + cfp[j] * T(j+1);
        if (!options.quasi2d) {
            dTdtCross[j] = - 0.5 * (sumcpj[j] + sumcpj[j-1]) * dTdx / (cp[j] * rho[j]);
        }
    }

    assert(mathUtils::notnan(ddtCross));
    assert(mathUtils::notnan(sumcpj));
}

void FlameSolver::updateBC()
{
    BoundaryCondition::BC leftPrev = grid.leftBC;
    BoundaryCondition::BC rightPrev = grid.rightBC;

    if (options.wallFlux && x[0] >= 0.0 && x[0] <= options.centerGridMin) {
        grid.leftBC = BoundaryCondition::WallFlux;
    } else if ((options.twinFlame || options.curvedFlame) &&
                x[0] >= 0.0 && x[0] <= options.centerGridMin) {
        grid.leftBC = BoundaryCondition::ControlVolume;
    } else if (grid.ju == 0 || (grid.jb == 0 && grid.fixedBurnedVal)) {
        grid.leftBC = BoundaryCondition::FixedValue;
    } else {
        grid.leftBC = BoundaryCondition::ZeroGradient;
    }

    if (options.flameType == "premixed" && grid.jb == jj && !grid.fixedBurnedVal) {
        grid.rightBC = BoundaryCondition::Floating;
    } else {
        grid.rightBC = BoundaryCondition::FixedValue;
    }

    if (leftPrev != grid.leftBC) {
        logFile.write(format("updateBC: Left BC changed from %i to %i.") %
                      leftPrev % grid.leftBC);
    }

    if (rightPrev != grid.rightBC) {
        logFile.write(format("updateBC: Right BC changed from %i to %i.") %
                      rightPrev % grid.rightBC);
    }
}

void FlameSolver::updateChemicalProperties()
{
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nPoints,1),
                      TbbWrapper<FlameSolver>(&FlameSolver::updateChemicalProperties, this));
}

void FlameSolver::updateChemicalProperties(size_t j1, size_t j2)
{
    CanteraGas& gas = gases.local();
    if (!gas.initialized()) {
        tbb::mutex::scoped_lock lock(gasInitMutex);
        gas.setOptions(options);
        gas.initialize();
    }

    // Calculate auxiliary data
    for (size_t j=j1; j<j2; j++) {
        // Thermodynamic properties
        thermoTimer.start();
        gas.setStateMass(&Y(0,j), T(j));
        rho[j] = gas.getDensity();
        Wmx[j] = gas.getMixtureMolecularWeight();
        cp[j] = gas.getSpecificHeatCapacity();
        gas.getSpecificHeatCapacities(&cpSpec(0,j));
        gas.getEnthalpies(&hk(0,j));
        thermoTimer.stop();

        // Transport Properties
        transportTimer.start();

        conductivityTimer.start();
        lambda[j] = gas.getThermalConductivity();
        conductivityTimer.stop();

        viscosityTimer.start();
        mu[j] = gas.getViscosity();
        viscosityTimer.stop();

        diffusivityTimer.start();
        gas.getWeightedDiffusionCoefficientsMass(&rhoD(0,j));
        gas.getThermalDiffusionCoefficients(&Dkt(0,j));
        diffusivityTimer.stop();
        transportTimer.stop();
    }
}


void FlameSolver::setDiffusionSolverState(double tInitial)
{
    splitTimer.resume();
    size_t k = 0;
    foreach (TridiagonalIntegrator& integrator, diffusionSolvers) {
        double dt = (dlj.square().tail(nPoints - 2) /
            (diffusionTerms[k].B * diffusionTerms[k].D).segment(1, nPoints-2)).minCoeff();
        dt = std::min(options.diffusionTimestepMultiplier * dt, options.globalTimestep);
        integrator.initialize(tInitial, dt);
        k++;
    }

    for (size_t i=0; i<nVars; i++) {
        diffusionSolvers[i].y = state.row(i);
    }
    splitTimer.stop();
}

void FlameSolver::setConvectionSolverState(double tInitial)
{
    splitTimer.resume();
    convectionSystem.setState(tInitial);
    splitTimer.stop();
}

void FlameSolver::setProductionSolverState(double tInitial)
{
    splitTimer.resume();
    for (size_t j=0; j<nPoints; j++) {
        sourceTerms[j].setState(tInitial, U(j), T(j), Y.col(j));
    }
    splitTimer.stop();
}

void FlameSolver::integrateConvectionTerms()
{
    setConvectionSolverState(tStageStart);
    convectionTimer.start();
    try {
        convectionSystem.integrateToTime(tStageEnd);
    } catch (DebugException& e) {
        logFile.write(e.errorString);
        writeStateFile("err_convectionIntegration", true, false);
        throw;
    }
    convectionTimer.stop();

    splitTimer.resume();
    convectionSystem.unroll_y();
    splitTimer.stop();
}

void FlameSolver::integrateProductionTerms()
{
    setProductionSolverState(tStageStart);
    reactionTimer.start();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nPoints,1),
                     TbbWrapper<FlameSolver>(&FlameSolver::integrateProductionTerms, this));
    reactionTimer.stop();
}

void FlameSolver::integrateProductionTerms(size_t j1, size_t j2)
{
    CanteraGas& gas = gases.local();
    if (!gas.initialized()) {
        tbb::mutex::scoped_lock lock(gasInitMutex);
        gas.setOptions(options);
        gas.initialize();
    }

    int err = 0;
    for (size_t j=j1; j<j2; j++) {
        SourceSystem& system = sourceTerms[j];
        logFile.verboseWrite(format("%i") % j, false);
        system.setGas(&gas);
        if (int(j) == options.debugSourcePoint &&
            tStageEnd >= options.debugSourceTime) {
            system.setDebug(true);
            std::ofstream steps("sourceTermSteps.py");
            system.writeState(steps, true);

            while (system.time() < tStageEnd - tStageStart && err >= 0) {
                err = system.integrateOneStep(tStageEnd - tStageStart);
                system.writeState(steps, false);
            }

            system.writeJacobian(steps);
            steps.close();
            std::terminate();

        } else {
            err = system.integrateToTime(tStageEnd - tStageStart);
        }
        if (err < 0) {
            logFile.write(format("Error at j = %i") % j);
            logFile.write(format("T = %s") % system.T);
            logFile.write(format("U = %s") % system.U);
            logFile.write("Y = ", false);
            logFile.write(system.Y);
            writeStateFile((format("prod%i_error_t%.6f_j%03i") %
                    1 % tStageEnd % j).str(), true, false);
        }
        logFile.verboseWrite(format(" [%s]...") % system.getStats(), false);
        sourceTerms[j].unroll_y();
        U(j) = sourceTerms[j].U;
        T(j) = sourceTerms[j].T;
        Y.col(j) = sourceTerms[j].Y;
    }
}

void FlameSolver::integrateDiffusionTerms()
{
    setDiffusionSolverState(tStageStart);
    diffusionTimer.start();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nVars, 1),
                      TbbWrapper<FlameSolver>(&FlameSolver::integrateDiffusionTerms, this));
    diffusionTimer.stop();
}

void FlameSolver::integrateDiffusionTerms(size_t k1, size_t k2)
{
    for (size_t k=k1; k<k2; k++) {
        diffusionSolvers[k].integrateToTime(tStageEnd);
        assert(mathUtils::almostEqual(diffusionSolvers[k].t, tStageEnd));
        state.row(k) = diffusionSolvers[k].y;
    }
}

void FlameSolver::rollVectorVector
(vector<dvector>& vv, const dmatrix& M) const
{
    size_t N = vv.size();
    vv.resize(N + nVars, dvector(nPoints));

    for (size_t k=0; k<nVars; k++) {
        Eigen::Map<dvec>(&vv[N+k][0], nPoints) = M.row(k);
    }
}


void FlameSolver::unrollVectorVector
(vector<dvector>& vv, dmatrix& M, size_t i) const
{
    for (size_t k=0; k<nVars; k++) {
        M.row(k) = Eigen::Map<dvec>(&vv[i*nVars+k][0], nPoints);
    }
}


void FlameSolver::update_xStag(const double t, const bool updateIntError)
{
    calculateQdot();
    xFlameActual = getFlamePosition();
    xFlameTarget = targetFlamePosition(t);
    if (updateIntError) {
        flamePosIntegralError += (xFlameTarget-xFlameActual)*(t-tFlamePrev);
        tFlamePrev = t;
    }

    // controlSignal is approximately a*xStag
    double controlSignal = options.xFlameProportionalGain *
        ((xFlameTarget - xFlameActual) +
        (flamePosIntegralError + (xFlameTarget - xFlameActual) * (t - tFlamePrev)) *
        options.xFlameIntegralGain);

    if (debugParameters::debugFlameRadiusControl) {
        logFile.write(format(
            "rFlameControl: rF=%g;  control=%g;  P=%g;  I=%g;  dt=%g") %
            xFlameActual %
            controlSignal %
            (options.xFlameProportionalGain * (xFlameTarget - xFlameActual)) %
            (options.xFlameProportionalGain * flamePosIntegralError *
            options.xFlameIntegralGain) %
            (t - tFlamePrev));
    }

    double a = strainfunc->a(t);
    if (alpha == 1) {
        rVzero = 0.5*rhoLeft*(controlSignal*abs(controlSignal)-a*x[0]*x[0]);
    } else {
        rVzero = rhoLeft*(controlSignal-a*x[0]);
    }
}


double FlameSolver::targetFlamePosition(double t)
{
    if (t <= options.xFlameT0) {
        return options.xFlameInitial;
    } else if (t >= options.xFlameT0 + options.xFlameDt) {
        return options.xFlameFinal;
    } else {
        return options.xFlameInitial + (options.xFlameFinal - options.xFlameInitial) *
            (t - options.xFlameT0) / options.xFlameDt;
    }
}


void FlameSolver::calculateQdot()
{
    reactionRatesTimer.start();
    for (size_t j=0; j<nPoints; j++) {
        gas.setStateMass(&Y(0,j), T(j));
        gas.getEnthalpies(&hk(0,j));
        gas.getReactionRates(&wDot(0,j));
        qDot[j] = - (wDot.col(j) * hk.col(j)).sum();
    }
    reactionRatesTimer.stop();
}


void FlameSolver::correctMassFractions() {
    setupTimer.resume();
    for (size_t j=0; j<nPoints; j++) {
        gas.setStateMass(&Y(0,j), T(j));
        gas.getMassFractions(&Y(0,j));
    }
    setupTimer.stop();
}

double FlameSolver::getHeatReleaseRate(void)
{
    return mathUtils::integrate(x, qDot);
}

double FlameSolver::getConsumptionSpeed(void)
{
    double QoverCp = mathUtils::integrate(x, qDot / cp);
    double rhouDeltaT = rhou*(T(grid.jb)-T(grid.ju));
    return QoverCp/rhouDeltaT;
}

double FlameSolver::getFlamePosition(void)
{
    return mathUtils::trapz(x, (x*qDot).eval())/mathUtils::trapz(x,qDot);
}

void FlameSolver::loadProfile(void)
{
    grid.alpha = (options.curvedFlame) ? 1 : 0;
    grid.unburnedLeft = options.unburnedLeft;

    // Read initial condition specified in the configuration file
    x = options.x_initial;
    nPoints = x.size();
    resizeMappedArrays();

    U = options.U_initial;
    T = options.T_initial;
    Y = options.Y_initial;
    if (Y.rows() == static_cast<dmatrix::Index>(x.size())) {
        Y = Y.transpose().eval();
    }
    convectionSystem.V = options.V_initial;
    convectionSystem.utwSystem.V = options.V_initial;
    rVzero = convectionSystem.utwSystem.V[0];

    grid.setSize(x.size());
    grid.updateValues();
    grid.updateBoundaryIndices();

    tNow = (options.haveTStart) ? options.tStart : 0.0;

    if (options.flameType == "premixed") {
        gas.setStateMass(&Y(0,grid.ju), T(grid.ju));
        rhou = gas.getDensity();
        gas.setStateMass(&Y(0,grid.jb), T(grid.ju));
        rhob = gas.getDensity();

        if (options.unburnedLeft) {
            rhoLeft = rhou;
            Tleft = T(grid.ju);
            Yleft = Y.col(grid.ju);
            rhoRight = rhob;
            Tright = T(grid.jb);
            Yright = Y.col(grid.jb);
        } else {
            rhoLeft = rhob;
            Tleft = T(grid.jb);
            Yleft = Y.col(grid.jb);
            rhoRight = rhou;
            Tright = T(grid.ju);
            Yright = Y.col(grid.ju);
        }

    } else if (options.flameType == "diffusion") {
        // Fuel composition
        size_t jFuel = (options.fuelLeft) ? 0 : jj;
        gas.thermo.setState_TPY(T(jFuel), options.pressure, &Y(0,jFuel));
        double rhoFuel = gas.getDensity();
        dvec Yfuel(nSpec);
        gas.getMassFractions(Yfuel);
        double Tfuel = gas.thermo.temperature();

        // Oxidizer composition
        size_t jOxidizer = (options.fuelLeft) ? jj : 0;
        gas.thermo.setState_TPY(T(jOxidizer),
                                options.pressure,
                                &Y(0,jOxidizer));
        double rhoOxidizer = gas.getDensity();
        dvec Yoxidizer(nSpec);
        gas.getMassFractions(Yoxidizer);
        double Toxidizer = gas.thermo.temperature();

        if (options.fuelLeft) {
            rhoLeft = rhoFuel;
            Tleft = Tfuel;
            Yleft = Yfuel;
            rhoRight = rhoOxidizer;
            Tright = Toxidizer;
            Yright = Yoxidizer;
        } else {
            rhoLeft = rhoOxidizer;
            Tleft = Toxidizer;
            Yleft = Yoxidizer;
            rhoRight = rhoFuel;
            Tright = Tfuel;
            Yright = Yfuel;
        }

        rhou = rhoLeft;
    } else if (options.flameType == "quasi2d") {
        gas.thermo.setState_TPY(T(0), options.pressure, &Y(0,0));
        rhoLeft = gas.thermo.density();
        Tleft = T(0);
        Yleft.resize(nSpec);
        gas.getMassFractions(Yleft);
        gas.thermo.setState_TPY(T(jj), options.pressure, &Y(0,jj));
        rhoRight = gas.thermo.density();
        Tright = T(jj);
        Yright.resize(nSpec);
        gas.getMassFractions(Yright);

        rhou = rhoRight;

    } else {
        throw DebugException("Invalid flameType: " + options.flameType);
    }

    updateBC();

    if (grid.leftBC == BoundaryCondition::ControlVolume &&
        options.xFlameControl)
    {
        double controlSignal;
        if (alpha == 0) {
            controlSignal = rVcenter/rhoLeft;
        } else {
            double tmp = pow(x[0],2) + 2*rVcenter/rhoLeft;
            controlSignal = mathUtils::sign(tmp)*sqrt(abs(tmp));
        }
        flamePosIntegralError = controlSignal /
            (options.xFlameProportionalGain * options.xFlameIntegralGain);
    }
}

void FlameSolver::printPerformanceStats(void)
{
    std::string filename = options.outputDir + "/stats";
    std::ofstream stats(filename.c_str(), std::ios::trunc | std::ios::out);
    totalTimer.stop();
    totalTimer.resume();
    stats << "\n   *** Performance Stats ***       time   ( call count )\n";
    printPerfString(stats, "                General Setup: ", setupTimer);
    printPerfString(stats, "             Split Term Setup: ", splitTimer);
    printPerfString(stats, "              Grid Adaptation: ", regridTimer);
    printPerfString(stats, "    Reaction Term Integration: ", reactionTimer);
    printPerfString(stats, "   Diffusion Term Integration: ", diffusionTimer);
    printPerfString(stats, "  Convection Term Integration: ", convectionTimer);
    printPerfString(stats, "                        Total: ", totalTimer);
    stats << "\n Subcomponents:\n";
    printPerfString(stats, "               Reaction Rates: ", reactionRatesTimer);
    printPerfString(stats, "         Transport Properties: ", transportTimer);
    printPerfString(stats, "          - thermal cond.    : ", conductivityTimer);
    printPerfString(stats, "          - viscosity        : ", viscosityTimer);
    printPerfString(stats, "          - diffusion coeff. : ", diffusivityTimer);
    printPerfString(stats, "     Thermodynamic Properties: ", thermoTimer);
    printPerfString(stats, "   Source Jacobian Evaluation: ", jacobianTimer);
    printPerfString(stats, "   UTW Convection Integration: ", convectionSystem.utwTimer);
    printPerfString(stats, "    Yk Convection Integration: ", convectionSystem.speciesTimer);
}

void FlameSolver::printPerfString(std::ostream& stats, const std::string& label,
                                  const PerfTimer& T)
{
    if (T.getTime() > 0.0) {
        stats << format("%s %9.3f (%12i)\n") % label % T.getTime() % T.getCallCount();
    }
}
