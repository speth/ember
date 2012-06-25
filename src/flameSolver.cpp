#include "flameSolver.h"
#include "debugUtils.h"
#include "perfTimer.h"
#include "dataFile.h"
#include "mathUtils.h"

#include "tbb/task_scheduler_init.h"

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/python.hpp>

#define foreach BOOST_FOREACH

FlameSolver::FlameSolver()
    : jCorrSolver(jCorrSystem)
{
}

FlameSolver::FlameSolver(const boost::python::api::object& config)
    : jCorrSolver(jCorrSystem)
{
    setOptions(configOptions(config));
}

void FlameSolver::setOptions(const configOptions& _options)
{
    options = _options;
    tStart = options.tStart;
    tEnd = options.tEnd;

    gas.setOptions(_options);
    grid.setOptions(_options);
}

void FlameSolver::initialize(void)
{
    try {
        tbb::task_scheduler_init tbbTaskSched(options.nThreads);
        strainfunc.setOptions(options);

        flamePosIntegralError = 0;

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

        dUdtConv.setZero(nPoints);
        dUdtDiff.setZero(nPoints);
        dUdtProd.setZero(nPoints);

        dTdtConv.setZero(nPoints);
        dTdtDiff.setZero(nPoints);
        dTdtProd.setZero(nPoints);

        dYdtConv.setZero(nSpec, nPoints);
        dYdtDiff.setZero(nSpec, nPoints);
        dYdtProd.setZero(nSpec, nPoints);

        updateChemicalProperties();
        // Determine initial condition for V
        dvec& V(convectionSystem.utwSystem.V);
        V.resize(nPoints);
        if ((options.twinFlame || options.curvedFlame) && !options.xFlameControl) {
            rVzero = 0;
            V[0] = 0;
            for (size_t j=1; j<nPoints; j++) {
                V[j] = V[j-1] - rho[j]*U[j]*(x[j]-x[j-1]);
            }
        } else {
            // Put the stagnation point on the burned side of the flame
            size_t jz = (grid.ju + 3 * grid.jb) / 4;
            convectionSystem.utwSystem.xVzero = x[jz];
            V[jz] = 0;

            for (size_t j=jz+1; j<nPoints; j++) {
                V[j] = V[j-1] - rho[j]*U[j]*(x[j]-x[j-1]);
            }

            if (jz != 0) {
                for (size_t j=jz; j>0; j--) {
                    V[j-1] = V[j] + rho[j-1]*U[j-1]*(x[j]-x[j-1]);
                }
            }
            rVzero = V[0];
        }
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
    catch (Cantera::CanteraError& err) {
        std::cout << err.what() << std::endl;
        throw;
    } catch (debugException& e) {
        std::cout << e.errorString << std::endl;
        logFile.write(e.errorString);
        throw;
    } catch (std::exception& e) {
        std::string message(e.what());
        std::cout << message << std::endl;
        logFile.write(message);
        throw;
    } catch (...) {
        std::string message = "I have no idea what went wrong!";
        std::cout << message << std::endl;
        logFile.write(message);
        throw;
    }
}

int FlameSolver::step(void)
{
    try {
        return step_internal();
    }
    catch (Cantera::CanteraError& err) {
        std::cout << err.what() << std::endl;
        throw;
    } catch (debugException& e) {
        std::cout << e.errorString << std::endl;
        logFile.write(e.errorString);
        throw;
    } catch (tbb::tbb_exception& e) {
        std::string message("\nUnhandled TBB exception:\nname:");
        message.append(e.name());
        message.append("\nwhat:");
        message.append(e.what());
        std::cout << message << std::endl;
        logFile.write(message);
        throw;
    } catch (std::exception& e) {
        std::string message(e.what());
        std::cout << message << std::endl;
        logFile.write(message);
    } catch (...) {
        std::string message = "I have no idea what went wrong!";
        std::cout << message << std::endl;
        logFile.write(message);
        throw;
    }
    return -1;
}

int FlameSolver::step_internal()
{
    setupTimer.start();

    // Debug sanity check
    #ifndef NDEBUG
        bool error = false;
        for (size_t j=0; j<nPoints; j++) {
            if (T[j] < 295 || T[j] > 3000) {
                logFile.write(format(
                    "WARNING: Unexpected Temperature: T = %f at j = %i") % T[j] % j);
                error = true;
            }
        }
        if (error) {
            writeStateFile();
        }
    #endif

    // Reset boundary conditions to prevent numerical drift
    if (grid.leftBC == BoundaryCondition::FixedValue) {
        T[0] = Tleft;
        Y.col(0) = Yleft;
    }

    if (grid.rightBC == BoundaryCondition::FixedValue) {
        T[jj] = Tright;
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
    prepareDiffusionTerms();
    prepareProductionTerms();
    prepareConvectionTerms();

    double dt = options.globalTimestep;
    double tNext = tNow + dt;

    if (t == tStart && options.outputProfiles) {
        writeStateFile("", false, false);
    }

    if (options.outputDebugIntegratorStages) {
        writeStateFile((format("start_t%.6f") % t).str(), true, false);
    }

    // *** Balanced Strang-split Integration ***

    integrateDiffusionTerms(t, t + 0.25*dt, 1); // step 1/4
    integrateConvectionTerms(t, t + 0.5*dt, 1); // step 1/2
    integrateDiffusionTerms(t + 0.25*dt, t + 0.5*dt, 2); // step 2/4
    integrateProductionTerms(t, t + dt, 0); // full step
    integrateDiffusionTerms(t + 0.5*dt, t + 0.75*dt, 3); // step 3/4
    integrateConvectionTerms(t + 0.5*dt, t + dt, 2); // step 2/2
    integrateDiffusionTerms(t + 0.75*dt, t + dt, 4); // step 4/4

    if (debugParameters::veryVerbose) {
        logFile.write("done!");
    }

    // *** End of Strang-split integration step ***
    correctMassFractions();
    calculateTimeDerivatives(dt);

    t = tNext;
    tNow = tNext;

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
    if (t > tOutput || nOutput >= options.outputStepInterval) {
        calculateQdot();
        timeVector.push_back(t);
        timestepVector.push_back(dt);
        heatReleaseRate.push_back(getHeatReleaseRate());
        consumptionSpeed.push_back(getConsumptionSpeed());
        flamePosition.push_back(getFlamePosition());

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
    if (t > tProfile || nProfile >= options.profileStepInterval) {
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
            num = std::min(num, rhoD.col(j).minCoeff());
            grid.dampVal[j] = sqrt(num/(rho[j]*strainfunc.a(t)));
        }
        dvec dampVal_prev = grid.dampVal;

        vector<dvector> currentSolution;
        rollVectorVector(currentSolution, U, T, Y);
        rollVectorVector(currentSolution, dUdtConv, dTdtConv, dYdtConv);
        rollVectorVector(currentSolution, dUdtDiff, dTdtDiff, dYdtDiff);
        rollVectorVector(currentSolution, dUdtProd, dTdtProd, dYdtProd);

        grid.nAdapt = nVars;
        if (options.quasi2d) {
            // do not change grid extents in this case
        } else if (strainfunc.a(tNow) == 0) {
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

            unrollVectorVector(currentSolution, U, T, Y, 0);
            unrollVectorVector(currentSolution, dUdtConv, dTdtConv, dYdtConv, 1);
            unrollVectorVector(currentSolution, dUdtDiff, dTdtDiff, dYdtDiff, 2);
            unrollVectorVector(currentSolution, dUdtProd, dTdtProd, dYdtProd, 3);
            correctMassFractions();

            // Update the mass flux (including the left boundary value)
            rVzero = mathUtils::interp1(x_prev, V_prev, grid.x[0]);
            convectionSystem.utwSystem.V = mathUtils::interp1(x_prev, V_prev, grid.x);

            // Allocate the solvers and arrays for auxiliary variables
            resizeAuxiliary();

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

    // Erase the existing file and create a new one
    if (boost::filesystem::exists(filename)) {
        boost::filesystem::remove(filename);
    }
    DataFile outFile(filename);

    updateChemicalProperties();
    if (updateDerivatives) {
        convectionSystem.evaluate();
    }

    // Write the state data to the output file:
    outFile.writeScalar("t", tNow);
    outFile.writeVec("x", x);
    outFile.writeVec("T", T);
    outFile.writeVec("U", U);
    outFile.writeArray2D("Y", Y, true);
    outFile.writeVec("V", convectionSystem.V);
    outFile.writeScalar("P", options.pressure);
    outFile.writeScalar("gridAlpha", grid.alpha);
    outFile.writeScalar("a", strainfunc.a(tNow));
    outFile.writeScalar("dadt", strainfunc.dadt(tNow));
    outFile.writeScalar("fileNumber", options.outputFileNumber);

    if (options.outputHeatReleaseRate || errorFile) {
        outFile.writeVec("q", qDot);
        outFile.writeVec("rho", rho);
    }

    if (options.outputTimeDerivatives || errorFile) {
        outFile.writeVec("dUdtDiff", dUdtDiff);
        outFile.writeVec("dUdtConv", dUdtConv);
        outFile.writeVec("dUdtProd", dUdtProd);

        outFile.writeVec("dTdtDiff", dTdtDiff);
        outFile.writeVec("dTdtConv", dTdtConv);
        outFile.writeVec("dTdtProd", dTdtProd);
        outFile.writeVec("dTdtCross", dTdtCross);

        outFile.writeArray2D("dYdtConv", dYdtConv, true);
        outFile.writeArray2D("dYdtDiff", dYdtDiff, true);
        outFile.writeArray2D("dYdtProd", dYdtProd, true);
        outFile.writeArray2D("dYdtCross", dYdtCross, true);

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

    outFile.close();
    if (incrementFileNumber) {
        options.outputFileNumber++;
    }
}

void FlameSolver::writeTimeseriesFile(const std::string& filename)
{
    std::string fullfilename = options.outputDir+"/"+filename+".h5";
    if (boost::filesystem::exists(fullfilename)) {
        boost::filesystem::remove(fullfilename);
    }
    DataFile outFile(fullfilename);
    outFile.writeVector("t", timeVector);
    outFile.writeVector("dt", timestepVector);
    outFile.writeVector("Q", heatReleaseRate);
    outFile.writeVector("Sc", consumptionSpeed);
    outFile.writeVector("xFlame", flamePosition);
    outFile.close();
}

void FlameSolver::resizeAuxiliary()
{
    size_t nPointsOld = rho.size();
    grid.setSize(T.size());

    if (nPoints == nPointsOld && !grid.updated) {
        return; // nothing to do
    }

    nSpec = gas.nSpec;
    nVars = 2+nSpec;
    N = nVars*nPoints;

    dYdtCross.resize(nSpec, nPoints);

    deltaYconv.resize(nSpec, nPoints);
    deltaYdiff.resize(nSpec, nPoints);
    deltaYprod.resize(nSpec, nPoints);

    dYdtCross *= NaN;

    deltaYconv *= NaN;
    deltaYdiff *= NaN;
    deltaYprod *= NaN;

    dTdtCross.setZero(nPoints);
    dTdtDiff.resize(nPoints);
    dTdtProd.resize(nPoints);
    dTdtConv.resize(nPoints);

    dUdtDiff.resize(nPoints);
    dUdtProd.resize(nPoints);
    dUdtConv.resize(nPoints);

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
    convectionSystem.resize(nPoints, nSpec);
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

void FlameSolver::updateCrossTerms()
{
    assert(mathUtils::notnan(Y));
    assert(mathUtils::notnan(T));
    assert(mathUtils::notnan(rhoD));

    for (size_t j=0; j<jj; j++) {
        jCorr[j] = 0;
        for (size_t k=0; k<nSpec; k++) {
            jFick(k,j) = -0.5*(rhoD(k,j)+rhoD(k,j+1)) * ((Y(k,j+1)-Y(k,j))/hh[j]);
            jSoret(k,j) = -0.5*(Dkt(k,j)/T[j] + Dkt(k,j+1)/T[j+1])
                * (T[j+1]-T[j])/hh[j];
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
        double dTdx = cfm[j] * T[j-1] + cf[j] * T[j] + cfp[j] * T[j+1];
        if (!options.quasi2d) {
            dTdtCross[j] = - 0.5 * (sumcpj[j] + sumcpj[j-1]) * dTdx / (cp[j] * rho[j]);
        }
    }

    assert(mathUtils::notnan(dYdtCross));
    assert(mathUtils::notnan(sumcpj));
    assert(mathUtils::notnan(dTdtCross));
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
    // Calculate auxiliary data
    for (size_t j=0; j<nPoints; j++) {
        // Thermodynamic properties
        thermoTimer.start();
        gas.setStateMass(&Y(0,j), T[j]);
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

void FlameSolver::prepareDiffusionTerms()
{
    splitTimer.resume();

    deltaUdiff.setZero(nPoints);
    deltaTdiff.setZero(nPoints);
    deltaYdiff.setZero(nSpec, nPoints);

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

    if (options.splittingMethod == "balanced") {
        diffusionTerms[kEnergy].splitConst = -0.75 * dTdtDiff +
                0.25 * (dTdtProd + dTdtConv + dTdtCross);
        diffusionTerms[kMomentum].splitConst = -0.75 * dUdtDiff +
                0.25 * (dUdtProd + dUdtConv);

        for (size_t k=0; k<nSpec; k++) {
            diffusionTerms[kSpecies+k].splitConst = - 0.75 * dYdtDiff.row(k) +
                0.25 * (dYdtProd + dYdtConv + dYdtCross).row(k);
        }
    } else { // options.splittingMethod == "strang"
        diffusionTerms[kEnergy].splitConst = dTdtCross;
        for (size_t k=0; k<nSpec; k++) {
            diffusionTerms[kSpecies+k].splitConst = dYdtCross.row(k);
        }
    }
    splitTimer.stop();
}

void FlameSolver::prepareProductionTerms()
{
    deltaUprod.setZero(nPoints);
    deltaTprod.setZero(nPoints);
    deltaYprod.setZero(nSpec, nPoints);

    setProductionSolverState(tNow);
    if (options.splittingMethod == "strang") {
        for (size_t j=0; j<nPoints; j++) {
            sourceTerms[j].resetSplitConstants();
        }
    } else { // options.splittingMethod == "balanced"
        for (size_t j=0; j<nPoints; j++) {
            sourceTerms[j].splitConst[kEnergy] = -0.5 * dTdtProd[j] +
                0.5 * (dTdtConv[j] + dTdtDiff[j] + dTdtCross[j]);
            sourceTerms[j].splitConst[kMomentum] = -0.5 * dUdtProd[j] +
                0.5 * (dUdtConv[j] + dUdtDiff[j]);
            sourceTerms[j].splitConst.tail(nSpec) = -0.5 * dYdtProd.col(j) +
                0.5 * (dYdtConv + dYdtDiff + dYdtCross).col(j);
        }
    }
}

void FlameSolver::prepareConvectionTerms()
{
    deltaUconv.setZero(nPoints);
    deltaTconv.setZero(nPoints);
    deltaYconv.setZero(nSpec, nPoints);

    setConvectionSolverState(tNow);
    dvec dTdt = dTdtConv + dTdtDiff + dTdtProd;
    Eigen::MatrixXd dYdt = dYdtConv + dYdtDiff + dYdtProd;
    if (options.splittingMethod == "balanced") {
        dTdt += dTdtCross;
        dYdt += dYdtCross.matrix();
    }

    dvec tmp = (W.inverse().matrix().transpose() * dYdt).array();
    drhodt = - rho * (dTdt / T + tmp * Wmx);

    assert(mathUtils::notnan(drhodt));
    convectionSystem.setDensityDerivative(drhodt);

    dvec splitConstT = dvec::Zero(nPoints);
    dvec splitConstU = dvec::Zero(nPoints);
    dmatrix splitConstY = dmatrix::Zero(nSpec, nPoints);

    if (options.splittingMethod == "balanced") {
        splitConstT = 0.25 * (dTdtProd + dTdtDiff + dTdtCross - 3 * dTdtConv);
        splitConstU = 0.25 * (dUdtProd + dUdtDiff - 3 * dUdtConv);
        splitConstY = 0.25 * (dYdtProd + dYdtDiff + dYdtCross - 3 * dYdtConv);
    }

    convectionSystem.setSplitConstants(splitConstU, splitConstT, splitConstY);
    convectionSystem.utwSystem.updateContinuityBoundaryCondition(qDot, options.continuityBC);
}

void FlameSolver::setDiffusionSolverState(double tInitial)
{
    splitTimer.resume();
    Ustart = U;
    Tstart = T;
    Ystart = Y;

    size_t k = 0;
    foreach (TridiagonalIntegrator& integrator, diffusionSolvers) {
        double dt = (dlj.square().tail(nPoints - 2) /
            (diffusionTerms[k].B * diffusionTerms[k].D).segment(1, nPoints-2)).minCoeff();
        dt = std::min(options.diffusionTimestepMultiplier * dt, options.globalTimestep);
        integrator.initialize(tInitial, dt);
        k++;
    }

    diffusionSolvers[kMomentum].y = U;
    diffusionSolvers[kEnergy].y = T;
    for (size_t k=0; k<nSpec; k++) {
        diffusionSolvers[kSpecies+k].y = Y.row(k);
    }
    splitTimer.stop();
}

void FlameSolver::setConvectionSolverState(double tInitial)
{
    splitTimer.resume();
    Ustart = U;
    Tstart = T;
    Ystart = Y;

    convectionSystem.setState(U, T, Y, tInitial);
    splitTimer.stop();
}

void FlameSolver::setProductionSolverState(double tInitial)
{
    splitTimer.resume();
    Ustart = U;
    Tstart = T;
    Ystart = Y;

    for (size_t j=0; j<nPoints; j++) {
        sourceTerms[j].setState(tInitial, U[j], T[j], Y.col(j));
    }
    splitTimer.stop();
}

void FlameSolver::integrateConvectionTerms(double tStart, double tEnd, int stage)
{
    assert(stage == 1 || stage == 2);
    if (debugParameters::veryVerbose) {
        logFile.write(format("convection term %i/2...") % stage, false);
    }
    setConvectionSolverState(tStart);

    convectionTimer.start();
    try {
        convectionSystem.integrateToTime(tEnd);
    } catch (debugException& e) {
        logFile.write(e.errorString);
        writeStateFile("err_convectionIntegration", true, false);
        throw;
    }
    convectionTimer.stop();

    splitTimer.resume();
    convectionSystem.unroll_y();
    U = convectionSystem.U;
    T = convectionSystem.T;
    Y = convectionSystem.Y;

    assert(mathUtils::notnan(T));
    assert(mathUtils::notnan(U));
    assert(mathUtils::notnan(Y));

    deltaUconv += U - Ustart;
    deltaTconv += T - Tstart;
    deltaYconv += Y - Ystart;
    splitTimer.stop();

    if (options.outputDebugIntegratorStages) {
        writeStateFile((format("conv%i_t%.6f") % stage % tNow).str(), true, false);
    }
}

void FlameSolver::integrateProductionTerms(double tStart, double tEnd, int stage)
{
    setProductionSolverState(tStart);

    reactionTimer.start();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nPoints,1),
                      SourceTermWrapper(this, tEnd, stage));
    reactionTimer.stop();

    splitTimer.resume();
    for (size_t j=0; j<nPoints; j++) {
        sourceTerms[j].unroll_y();
        U[j] = sourceTerms[j].U;
        T[j] = sourceTerms[j].T;
        Y.col(j) = sourceTerms[j].Y;
    }
    assert(mathUtils::notnan(T));
    assert(mathUtils::notnan(U));
    assert(mathUtils::notnan(Y));

    deltaUprod += U - Ustart;
    deltaTprod += T - Tstart;
    deltaYprod += Y - Ystart;
    splitTimer.stop();

    if (options.outputDebugIntegratorStages) {
        writeStateFile((format("prod_t%.6f") % tNow).str(), true, false);
    }
}


void FlameSolver::integrateDiffusionTerms(double tStart, double tEnd, int stage)
{
    assert(stage >= 0 && stage <= 4);
    if (debugParameters::veryVerbose) {
        logFile.write(format("diffusion terms %i/4...") % stage, false);
    }

    setDiffusionSolverState(tStart);
    diffusionTimer.start();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nVars, 1),
                       DiffusionTermWrapper(this, tEnd));
    diffusionTimer.stop();

    splitTimer.resume();
    U = diffusionSolvers[kMomentum].y;
    T = diffusionSolvers[kEnergy].y;
    for (size_t k=0; k<nSpec; k++) {
        Y.row(k) = diffusionSolvers[kSpecies+k].y;
    }
    assert(mathUtils::notnan(T));
    assert(mathUtils::notnan(U));
    assert(mathUtils::notnan(Y));

    deltaUdiff += U - Ustart;
    deltaTdiff += T - Tstart;
    deltaYdiff += Y - Ystart;
    splitTimer.stop();

    if (stage && options.outputDebugIntegratorStages) {
        writeStateFile((format("diff%i_t%.6f") % stage % tNow).str(), true, false);
    }
}


void FlameSolver::rollVectorVector
(vector<dvector>& vv, const dvec& u, const dvec& t, const dmatrix& y) const
{
    size_t N = vv.size();
    vv.resize(N + nSpec + 2, dvector(nPoints));

    Eigen::Map<dvec>(&vv[N][0], nPoints) = u;
    Eigen::Map<dvec>(&vv[N+1][0], nPoints) = t;
    for (size_t k=0; k<nSpec; k++) {
        Eigen::Map<dvec>(&vv[N+k+2][0], nPoints) = y.row(k);
    }
}


void FlameSolver::unrollVectorVector
(vector<dvector>& vv, dvec& u, dvec& t, dmatrix& y, size_t i) const
{
    u = Eigen::Map<dvec>(&vv[i*nVars+kMomentum][0], nPoints);
    t = Eigen::Map<dvec>(&vv[i*nVars+kEnergy][0], nPoints);
    y.resize(nSpec, nPoints);
    for (size_t k=0; k<nSpec; k++) {
        y.row(k) = Eigen::Map<dvec>(&vv[i*nVars+kSpecies+k][0], nPoints);
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

    double a = strainfunc.a(t);
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
        gas.setStateMass(&Y(0,j), T[j]);
        gas.getEnthalpies(&hk(0,j));
        gas.getReactionRates(&wDot(0,j));
        qDot[j] = - (wDot.col(j) * hk.col(j)).sum();
    }
    reactionRatesTimer.stop();
}


void FlameSolver::correctMassFractions() {
    setupTimer.resume();
    for (size_t j=0; j<nPoints; j++) {
        gas.setStateMass(&Y(0,j), T[j]);
        gas.getMassFractions(&Y(0,j));
    }
    setupTimer.stop();
}

void FlameSolver::calculateTimeDerivatives(double dt)
{
    double split(options.splittingMethod == "balanced");

    if (options.splittingMethod == "balanced") {
        deltaUconv -= 0.25 * dt * (dUdtProd + dUdtDiff);
        deltaUdiff -= 0.25 * dt * (dUdtProd + dUdtConv);
        deltaUprod -= 0.5 * dt * (dUdtConv + dUdtDiff);

        deltaTconv -= 0.25 * dt * (dTdtProd + dTdtDiff + dTdtCross);
        deltaTdiff -= 0.25 * dt * (dTdtProd + dTdtConv + dTdtCross);
        deltaTprod -= 0.5 * dt * (dTdtConv + dTdtDiff + dTdtCross);

        deltaYconv -= 0.25 * dt * (dYdtProd + dYdtDiff + dYdtCross);
        deltaYdiff -= 0.25 * dt * (dYdtProd + dYdtConv + dYdtCross);
        deltaYprod -= 0.5 * dt * (dYdtConv + dYdtDiff + dYdtCross);
    }

    dUdtConv = deltaUconv / dt + split * 0.75 * dUdtConv;
    dUdtDiff = deltaUdiff / dt + split * 0.75 * dUdtDiff;
    dUdtProd = deltaUprod / dt + split * 0.5 * dUdtProd;

    dTdtConv = deltaTconv / dt + split * 0.75 * dTdtConv;
    dTdtDiff = deltaTdiff / dt + split * 0.75 * dTdtDiff;
    dTdtProd = deltaTprod / dt + split * 0.5 * dTdtProd;

    dYdtConv = deltaYconv / dt + split * 0.75 * dYdtConv;
    dYdtDiff = deltaYdiff / dt + split * 0.75 * dYdtDiff;
    dYdtProd = deltaYprod / dt + split * 0.5 * dYdtProd;
}

double FlameSolver::getHeatReleaseRate(void)
{
    return mathUtils::integrate(x, qDot);
}

double FlameSolver::getConsumptionSpeed(void)
{
    double QoverCp = mathUtils::integrate(x, qDot / cp);
    double rhouDeltaT = rhou*(Tb-Tu);
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

    if (options.haveRestartFile) {
        std::string inputFilename;
        if (options.useRelativeRestartPath) {
            inputFilename = options.inputDir + "/" + options.restartFile;
        } else {
            inputFilename = options.restartFile;
        }

        logFile.write(format("Reading initial condition from %s") % inputFilename);
        DataFile infile(inputFilename);
        x = infile.readVec("x");

        grid.setSize(x.size());
        grid.updateValues();
        grid.updateBoundaryIndices();

        U = infile.readVec("U");
        T = infile.readVec("T");
        Y = infile.readArray2D("Y", true);
        tNow = infile.readScalar("t");
        if (!options.haveTStart) {
            // If tStart is not in the input file, use the time from the restart file.
            tStart = tNow;
        }

        dvec V = infile.readVec("V");
        rVzero = V[0];

        if (!options.fileNumberOverride) {
            options.outputFileNumber = (int) infile.readScalar("fileNumber");
        }

        infile.close();
    } else if (options.haveInitialProfiles) {
        // Read initial condition specified in the configuration file
        logFile.write("Reading initial condition from configuration file.");
        x = options.x_initial;
        U = options.U_initial;
        T = options.T_initial;
        Y = options.Y_initial;
        if (Y.rows() == static_cast<dmatrix::Index>(x.size())) {
            Y = Y.transpose().eval();
        }
        rVzero = options.rVzero_initial;

        grid.setSize(x.size());
        grid.updateValues();
        grid.updateBoundaryIndices();

        tNow = (options.haveTStart) ? options.tStart : 0.0;
    } else {
        throw debugException("Initial profile data required but not provided.");
    }

    Tu = (options.overrideTu) ? options.Tu : T[grid.ju];

    CanteraGas gas;
    gas.setOptions(options);
    gas.initialize();
    size_t nSpec = gas.nSpec;

    if (options.flameType == "premixed") {
        // Save the burned gas properties for the case where the burned
        // values are not fixed.
        double TbSave = T[grid.jb];
        dvec YbSave = Y.col(grid.jb);

        if (options.overrideReactants) {
            dvec reactants = gas.calculateReactantMixture(
                options.fuel, options.oxidizer, options.equivalenceRatio);
            gas.thermo.setState_TPX(Tu,gas.pressure, &reactants[0]);
            gas.thermo.getMassFractions(&Y(0,grid.ju));
            gas.thermo.setState_TPX(Tu,gas.pressure, &reactants[0]);
            Cantera::equilibrate(gas.thermo,"HP");
            gas.thermo.getMassFractions(&Y(0,grid.jb));
            T[grid.jb] = gas.thermo.temperature();
        }
        Tb = T[grid.jb];

        gas.setStateMass(&Y(0,grid.ju), T[grid.ju]);
        rhou = gas.getDensity();
        gas.setStateMass(&Y(0,grid.jb), T[grid.ju]);
        rhob = gas.getDensity();
        Yu = Y.col(grid.ju);
        Yb = Y.col(grid.jb);

        if (options.unburnedLeft) {
            rhoLeft = rhou;
            Tleft = Tu;
            Yleft = Yu;
            rhoRight = rhob;
            Tright = Tb;
            Yright = Yb;
        } else {
            rhoLeft = rhob;
            Tleft = Tb;
            Yleft = Yb;
            rhoRight = rhou;
            Tright = Tu;
            Yright = Yu;
        }

        if (!options.fixedBurnedVal) {
            T[grid.jb] = TbSave;
            Y.col(grid.jb) = YbSave;
        }

        if (options.overrideTu) {
            T[grid.ju] = Tu;
        }

    } else if (options.flameType == "diffusion") {
        // Fuel composition
        size_t jFuel = (options.fuelLeft) ? 0 : jj;
        if (options.overrideReactants) {
            gas.thermo.setState_TPX(options.Tfuel, options.pressure, options.fuel);
        } else {
            gas.thermo.setState_TPY(T[jFuel], options.pressure, &Y(0,jFuel));
        }
        double rhoFuel = gas.getDensity();
        dvec Yfuel(nSpec);
        gas.getMassFractions(Yfuel);
        double Tfuel = gas.thermo.temperature();

        // Oxidizer composition
        size_t jOxidizer = (options.fuelLeft) ? jj : 0;
        if (options.overrideReactants) {
            gas.thermo.setState_TPX(options.Toxidizer,
                                    options.pressure,
                                    options.oxidizer);
        } else {
            gas.thermo.setState_TPY(T[jOxidizer],
                                    options.pressure,
                                    &Y(0,jOxidizer));
        }
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
        gas.thermo.setState_TPY(T[0], options.pressure, &Y(0,0));
        rhoLeft = gas.thermo.density();
        Tleft = T[0];
        Yleft.resize(nSpec);
        gas.getMassFractions(Yleft);
        gas.thermo.setState_TPY(T[jj], options.pressure, &Y(0,jj));
        rhoRight = gas.thermo.density();
        Tright = T[jj];
        Yright.resize(nSpec);
        gas.getMassFractions(Yright);

        rhou = rhoRight;

    } else {
        throw debugException("Invalid flameType: " + options.flameType);
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
    if (boost::filesystem::exists(filename)) {
        boost::filesystem::remove(filename);
    }
    statsFile.open(filename.c_str(), std::ios::trunc | std::ios::out);
    statsFile << "\n   *** Performance Stats ***       time   ( call count )\n";
    printPerfString("                General Setup: ", setupTimer);
    printPerfString("             Split Term Setup: ", splitTimer);
    printPerfString("    Reaction Term Integration: ", reactionTimer);
    printPerfString("   Diffusion Term Integration: ", diffusionTimer);
    printPerfString("  Convection Term Integration: ", convectionTimer);
    statsFile << "\n Subcomponents:\n";
    printPerfString("               Reaction Rates: ", reactionRatesTimer);
    printPerfString("         Transport Properties: ", transportTimer);
    printPerfString("          - thermal cond.    : ", conductivityTimer);
    printPerfString("          - viscosity        : ", viscosityTimer);
    printPerfString("          - diffusion coeff. : ", diffusivityTimer);
    printPerfString("     Thermodynamic Properties: ", thermoTimer);
    printPerfString("   Source Jacobian Evaluation: ", jacobianTimer);
    printPerfString("   UTW Convection Integration: ", convectionSystem.utwTimer);
    printPerfString("    Yk Convection Integration: ", convectionSystem.speciesTimer);
    statsFile.close();
}

void FlameSolver::printPerfString(const std::string& label, const PerfTimer& T)
{
    statsFile << format("%s %9.3f (%12i)\n") % label % T.getTime() % T.getCallCount();
}

void SourceTermWrapper::operator()(const tbb::blocked_range<size_t>& r) const
{
    if (debugParameters::veryVerbose) {
        if (stage_) {
            logFile.write(format("Source term %i/2...") % stage_, false);
        } else {
            logFile.write("Source term...", false);
        }
    }
    CanteraGas& gas = parent_->gases.local();
    if (!gas.initialized()) {
        tbb::mutex::scoped_lock lock(parent_->gasInitMutex);
        gas.setOptions(parent_->options);
        gas.initialize();
    }

    int err = 0;
    for (size_t j=r.begin(); j<r.end(); j++) {
        SourceSystem& system = parent_->sourceTerms[j];
        if (debugParameters::veryVerbose) {
            logFile.write(format("%i") % j, false);
        }
        system.setGas(&gas);
        if (int(j) == parent_->options.debugSourcePoint &&
            t_ >= parent_->options.debugSourceTime) {
            system.setDebug(true);
            std::ofstream steps("sourceTermSteps.py");
            system.writeState(steps, true);

            while (system.time() < t_ - parent_->tNow && err >= 0) {
                err = system.integrateOneStep(t_ - parent_->tNow);
                system.writeState(steps, false);
            }

            system.writeJacobian(steps);
            steps.close();
            std::terminate();

        } else {
            err = system.integrateToTime(t_ - parent_->tNow);
        }
        if (err < 0) {
            logFile.write(format("Error at j = %i") % j);
            logFile.write(format("T = %s") % system.T);
            logFile.write(format("U = %s") % system.U);
            logFile.write("Y = ", false);
            logFile.write(system.Y);
            parent_->writeStateFile((format("prod%i_error_t%.6f_j%03i") %
                    stage_ % t_ % j).str(), true, false);
        }

        if (debugParameters::veryVerbose) {
            logFile.write(format(" [%s]...") % system.getStats(), false);
        }
    }
}

void DiffusionTermWrapper::operator()(const tbb::blocked_range<size_t>& r) const
{
    for (size_t k=r.begin(); k<r.end(); k++) {
        parent_->diffusionSolvers[k].integrateToTime(t_);
    }
}
