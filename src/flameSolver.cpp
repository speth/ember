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

    gases.resize(options.nThreads);

    // Cantera has problems with the multithreaded initialization of objects.
    std::vector<CanteraGas**> ptrs(options.nThreads);
    for (int i=0; i<options.nThreads; i++) {
        logFile.write(format("Intializing Cantera object %i/%i") %
                      (i+1) % options.nThreads);
        gases[i]->setOptions(options);
        gases[i]->initialize();
    }

    gas.setOptions(_options);
    grid.setOptions(_options);
}

void FlameSolver::initialize(void)
{
    try {
        strainfunc.setOptions(options);

        flamePosIntegralError = 0;

        // Cantera initialization
        gas.initialize();
        nSpec = gas.nSpec;
        nVars = nSpec + 2;
        W.resize(nSpec);
        gas.getMolecularWeights(W);

        // Initial Conditions
        if (options.haveRestartFile || options.haveInitialProfiles) {
            loadProfile();
        } else {
            generateProfile();
        }

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

void FlameSolver::tryrun(void)
{
    try {
        run();
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
}

void FlameSolver::run(void)
{
    totalTimer.start();
    tbb::task_scheduler_init tbbTaskSched(options.nThreads);

    double t = tStart;

    long int nTotal = 0; // total number of timesteps taken
    int nRegrid = 0; // number of time steps since regridding/adaptation
    int nOutput = 0; // number of time steps since storing integral flame parameters
    int nProfile = 0; // number of time steps since saving flame profiles
    int nTerminate = 1; // number of steps since last checking termination condition
    int nCurrentState = 0; // number of time steps since profNow.h5 and out.h5 were written

    double tOutput = t; // time of next integral flame parameters output (this step)
    double tRegrid = t + options.regridTimeInterval; // time of next regridding
    double tProfile = t + options.profileTimeInterval; // time of next profile output

    grid.updateValues();
    resizeAuxiliary();

    tFlamePrev = t;
    tNow = t;

    while (true) {
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

        splitTimer.start();
        if (!options.quasi2d) {
            // Diffusion solvers: Energy and momentum
            diffusionTerms[kMomentum].B = rho.inverse();
            diffusionTerms[kEnergy].B = (rho * cp).inverse();

            diffusionTerms[kMomentum].D = mu;
            diffusionTerms[kEnergy].D = lambda;

            // Diffusion solvers: Species
            for (size_t k=0; k<nSpec; k++) {
                DiffusionSystem& sys = diffusionTerms[kSpecies+k];
                sys.B = rho.inverse();
                sys.D = rhoD.row(k);
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
        splitTimer.stop();

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

        // Diffusion Integration: first quarter-step
        setDiffusionSolverState(t);
        integrateDiffusionTerms(t + 0.25*dt, 1);
        extractDiffusionState(0.25*dt, 1);

        // Convection Integration: first half-step
        setConvectionSolverState(t, 1);
        integrateConvectionTerms(t + 0.5*dt, 1);
        extractConvectionState(0.5*dt, 1);

        // Diffusion Integration: second quarter-step
        setDiffusionSolverState(t + 0.25*dt);
        integrateDiffusionTerms(t + 0.5*dt, 2);
        extractDiffusionState(0.25*dt, 2);

        // Reaction Integration: full step
        setProductionSolverState(t);
        integrateProductionTerms(t + dt, 0);
        extractProductionState(dt);

        // Diffusion Integration: third quarter-step
        setDiffusionSolverState(t + 0.5*dt);
        integrateDiffusionTerms(t + 0.75*dt, 3);
        extractDiffusionState(0.25*dt, 3);

        // Convection Integration: second half-step
        setConvectionSolverState(t + 0.5*dt, 2);
        integrateConvectionTerms(t + dt, 2);
        extractConvectionState(0.5*dt, 2);

        // Diffusion Integration: fourth quarter-step
        setDiffusionSolverState(t + 0.75*dt);
        integrateDiffusionTerms(t + dt, 4);
        extractDiffusionState(0.25*dt, 4);

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

        // *** Periodic check for terminating the integration
        //     (based on steady heat release rate, etc.)
        //     Quit now to skip grid adaptation on the last step
        if (t >= tEnd) {
            break;
        } else if (nTerminate >= options.terminateStepInterval) {
            nTerminate = 0;
            if (checkTerminationCondition()) {
                break;
            }
        }

        // *** Save the current integral and profile data
        //     in files that are automatically overwritten,
        //     and save the time-series data (out.h5)
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
            dvector x_prev = grid.x;
            convectionSystem.evaluate();
            dvector V_prev(nPoints);

            // dampVal sets a limit on the maximum grid size
            grid.dampVal.resize(grid.x.size());
            for (size_t j=0; j<nPoints; j++) {
                double num = std::min(mu[j],lambda[j]/cp[j]);
                num = std::min(num, rhoD.col(j).minCoeff());
                grid.dampVal[j] = sqrt(num/(rho[j]*strainfunc.a(t)));
                V_prev[j] = convectionSystem.V[j];
            }
            dvector dampVal_prev = grid.dampVal;

            // "rollVectorVector"
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
                dvector V = mathUtils::interp1(x_prev, V_prev, grid.x);
                convectionSystem.utwSystem.V = Eigen::Map<dvec>(&V[0], V.size());

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
    }

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
    outFile.writeVector("x", x);
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
        outFile.writeVector("q", qDot);
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

        // Number of timesteps in the chemistry solver in the last global timestep
        dvector chemSteps(nPoints, 0);
        for (size_t j=0; j<nPoints; j++) {
            if (sourceSolvers[j].initialized()) {
                chemSteps[j] = sourceSolvers[j].getNumSteps();
            }
        }
        outFile.writeVector("chemSteps", chemSteps);
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
        outFile.writeVector("cfp", grid.cfp);
        outFile.writeVector("cf", grid.cf);
        outFile.writeVector("cfm", grid.cfm);
        outFile.writeVector("hh", hh);
        outFile.writeVector("rphalf", grid.rphalf);
        outFile.writeArray2D("jFick", jFick, true);
        outFile.writeArray2D("jSoret", jSoret, true);
        outFile.writeVec("jCorr", jCorr);
    }

    outFile.close();
    if (incrementFileNumber) {
        options.outputFileNumber++;
    }

//    if (errorFile && options.stopIfError) {
//        logFile.write(format("Error outputs remaining until termination: %i") %
//                options.errorStopCount);
//      if (options.errorStopCount-- <= 0) {
//        throw debugException("Too many integration failures.");
//      }
//    }
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
            // Create and initialize the new SourceSystem
            SourceSystem* system = new SourceSystem();
            system->resize(nSpec);
            system->gas = &gas;
            system->options = &options;
            system->thermoTimer = &thermoTimer;
            system->reactionRatesTimer = &reactionRatesTimer;
            system->jacobianTimer = &jacobianTimer;
            system->strainFunction = strainfunc;
            system->rhou = rhou;
            system->W = W;
            system->j = j;
            system->x = x[j];
            if (options.quasi2d) {
                system->setupQuasi2d(vzInterp, TInterp);
            }

            // Create and initialize the new Sundials solver
            sundialsCVODE* solver = new sundialsCVODE(nVars);
            solver->setODE(system);
            solver->abstol[kMomentum] = options.integratorMomentumAbsTol;
            solver->abstol[kEnergy] = options.integratorEnergyAbsTol;
            for (size_t k=0; k<nSpec; k++) {
                solver->abstol[kSpecies+k] = options.integratorSpeciesAbsTol;
            }
            solver->reltol = options.integratorRelTol;
            solver->linearMultistepMethod = CV_BDF;
            solver->nonlinearSolverMethod = CV_NEWTON;
            solver->maxNumSteps = 1000000;
            solver->minStep = options.integratorMinTimestep;

            // Store the solver and system
            sourceTerms.push_back(system);
            sourceSolvers.push_back(solver);

            // Create the new SourceSystemQSS
            SourceSystemQSS* qssSolver = new SourceSystemQSS();
            qssSolver->initialize(nSpec);
            qssSolver->setOptions(options);
            qssSolver->gas = &gas;
            qssSolver->thermoTimer = &thermoTimer;
            qssSolver->reactionRatesTimer = &reactionRatesTimer;
            qssSolver->strainFunction = strainfunc;
            qssSolver->rhou = rhou;
            qssSolver->W = W;
            qssSolver->j = j;
            qssSolver->x = x[j];
            qssSolver->ymin.setConstant(nVars, options.qss_minval);
            qssSolver->ymin[kMomentum] = -1e4;
            if (options.quasi2d) {
                qssSolver->setupQuasi2d(vzInterp, TInterp);
            }

            sourceTermsQSS.push_back(qssSolver);
            useCVODE.push_back(options.chemistryIntegrator == "cvode");
        }

    } else {
        // Delete the unwanted solvers and systems
        sourceTerms.erase(sourceTerms.begin()+nPoints, sourceTerms.end());
        sourceSolvers.erase(sourceSolvers.begin()+nPoints, sourceSolvers.end());
        sourceTermsQSS.erase(sourceTermsQSS.begin()+nPoints, sourceTermsQSS.end());
        useCVODE.erase(useCVODE.begin()+nPoints, useCVODE.end());
    }

    // Resize solution vector for diffusion systems / solvers
    for (size_t k=0; k<nVars; k++) {
        diffusionSolvers[k].resize(nPoints);
        diffusionTerms[k].setGrid(grid);
    }

    // Diffusion solvers: Energy and Momentum
    diffusionSolvers[kMomentum].resize(nPoints);
    diffusionSolvers[kEnergy].resize(nPoints);

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
        sourceTerms[j].x = x[j];
        sourceTermsQSS[j].x = x[j];
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

    setDiffusionSolverState(tNow);

    if (options.splittingMethod == "balanced") {
        diffusionTerms[kEnergy].splitConst = -0.75 * dTdtDiff +
                0.25 * (dTdtProd + dTdtConv + dTdtCross);
        diffusionTerms[kMomentum].splitConst = -0.75 * dUdtDiff +
                0.25 * (dUdtProd + dUdtConv);

        for (size_t k=0; k<nSpec; k++) {
            dvec& splitConstY = diffusionTerms[kSpecies+k].splitConst;
            for (size_t j = 0; j <= jj; j++)
            {
                splitConstY[j] = -0.75 * dYdtDiff(k,j) +
                    0.25 * (dYdtProd(k,j) + dYdtConv(k,j) + dYdtCross(k,j));
            }
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
            if (useCVODE[j]) {
                sourceTerms[j].resetSplitConstants();
            } else {
                sourceTermsQSS[j].resetSplitConstants();
            }
        }
    } else { // options.splittingMethod == "balanced"
        for (size_t j=0; j<nPoints; j++) {
            if (useCVODE[j]) {
                sourceTerms[j].splitConst[kEnergy] = -0.5 * dTdtProd[j] +
                    0.5 * (dTdtConv[j] + dTdtDiff[j] + dTdtCross[j]);
                sourceTerms[j].splitConst[kMomentum] = -0.5 * dUdtProd[j] +
                    0.5 * (dUdtConv[j] + dUdtDiff[j]);
                for (size_t k=0; k<nSpec; k++) {
                    sourceTerms[j].splitConst[kSpecies+k] = -0.5 * dYdtProd(k,j) +
                        0.5 * (dYdtConv(k,j) + dYdtDiff(k,j) + dYdtCross(k,j));
                }
            } else {
                sourceTermsQSS[j].splitConstT = -0.5 * dTdtProd[j] +
                    0.5 * (dTdtConv[j] + dTdtDiff[j] + dTdtCross[j]);
                sourceTermsQSS[j].splitConstU = -0.5 * dUdtProd[j] +
                    0.5 * (dUdtConv[j] + dUdtDiff[j]);
                for (size_t k=0; k<nSpec; k++) {
                    sourceTermsQSS[j].splitConstY[k] = -0.5 * dYdtProd(k,j) +
                        0.5 * (dYdtConv(k,j) + dYdtDiff(k,j) + dYdtCross(k,j));
                }
            }
        }
    }
}

void FlameSolver::prepareConvectionTerms()
{
    deltaUconv.setZero(nPoints);
    deltaTconv.setZero(nPoints);
    deltaYconv.setZero(nSpec, nPoints);

    setConvectionSolverState(tNow, 0);
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
        double dt = options.globalTimestep;
        for (size_t j=1; j<jj; j++) {
            dt = std::min(dt, options.diffusionTimestepMultiplier*dlj[j]*dlj[j] /
                              (diffusionTerms[k].B[j]*diffusionTerms[k].D[j]));
        }
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

void FlameSolver::setConvectionSolverState(double tInitial, int stage)
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
        if (useCVODE[j]) {
            sourceSolvers[j].t0 = tInitial;
            sdVector& ySource = sourceSolvers[j].y;
            ySource[kMomentum] = U[j];
            ySource[kEnergy] = T[j];
            for (size_t k=0; k<nSpec; k++) {
                ySource[kSpecies+k] = Y(k,j);
            }
            sourceSolvers[j].initialize();

        } else {
            dvec ySource(nVars);
            ySource[kMomentum] = U[j];
            ySource[kEnergy] = T[j];
            for (size_t k=0; k<nSpec; k++) {
                ySource[kSpecies+k] = Y(k,j);
            }
            sourceTermsQSS[j].setState(ySource, tInitial);
        }
    }
    splitTimer.stop();
}

void FlameSolver::extractConvectionState(double dt, int stage)
{
    assert(stage == 1 || stage == 2);

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

void FlameSolver::extractDiffusionState(double dt, int stage)
{
    assert(stage >= 0 && stage <= 4);

    splitTimer.resume();
    U = diffusionSolvers[kMomentum].y;
    T = diffusionSolvers[kEnergy].y;
    for (size_t k=0; k<nSpec; k++) {
        const TridiagonalIntegrator& solver = diffusionSolvers[kSpecies+k];
        Y.row(k) = solver.y;
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

void FlameSolver::extractProductionState(double dt)
{
    splitTimer.resume();
    for (size_t j=0; j<nPoints; j++) {
        if (useCVODE[j]) {
            sourceTerms[j].unroll_y(sourceSolvers[j].y, tNow);
            U[j] = sourceTerms[j].U;
            T[j] = sourceTerms[j].T;
            for (size_t k=0; k<nSpec; k++) {
                Y(k,j) = sourceTerms[j].Y[k];
            }
        } else {
            sourceTermsQSS[j].unroll_y(sourceTermsQSS[j].y);
            U[j] = sourceTermsQSS[j].U;
            T[j] = sourceTermsQSS[j].T;
            for (size_t k=0; k<nSpec; k++) {
                Y(k,j) = sourceTermsQSS[j].Y[k];
            }
        }
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

void FlameSolver::integrateConvectionTerms(double t, int stage)
{
    if (debugParameters::veryVerbose) {
        logFile.write(format("convection term %i/2...") % stage, false);
    }

    convectionTimer.start();
    try {
        convectionSystem.integrateToTime(t);
    } catch (debugException& e) {
        logFile.write(e.errorString);
        writeStateFile("err_convectionIntegration", true, false);
        throw;
    }
    convectionTimer.stop();
}

void FlameSolver::integrateProductionTerms(double t, int stage)
{
    reactionTimer.start();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nPoints,1),
                      SourceTermWrapper(this, t, stage));
    reactionTimer.stop();
}


void FlameSolver::integrateDiffusionTerms(double t, int stage)
{
    if (debugParameters::veryVerbose) {
        logFile.write(format("diffusion terms %i/2...") % stage, false);
    }

    diffusionTimer.start();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nVars, 1),
                       DiffusionTermWrapper(this, t));
    diffusionTimer.stop();
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
        //
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
    dvector q_cp(nPoints);
    for (size_t j=0; j<nPoints; j++) {
        q_cp[j] = qDot[j]/cp[j];
    }
    double QoverCp = mathUtils::integrate(x,q_cp);
    double rhouDeltaT = rhou*(Tb-Tu);
    return QoverCp/rhouDeltaT;
}

double FlameSolver::getFlamePosition(void)
{
    return mathUtils::trapz(x,x*qDot)/mathUtils::trapz(x,qDot);
}

void FlameSolver::generateProfile(void)
{
    logFile.write("Generating initial profiles from given fuel and oxidizer compositions.");

    // Set up a CanteraGas object to use for calculating the initial profiles
    nSpec = gas.nSpec;
    grid.setSize(options.nPoints);

    // Create a uniform initial grid
    x.resize(nPoints);
    double xRight, xLeft;
    xRight = options.xRight;
    if (options.twinFlame || options.curvedFlame) {
        x[0] = 0;
        xLeft = options.centerGridMin;
    } else {
        xLeft = options.xLeft;
        x[0] = xLeft;
    }
    // Uniform initial grid
    double dx = (xRight-xLeft)/((double) nPoints-1);
    for (size_t j=1; j<nPoints; j++) {
        x[j] = xLeft + j * dx;
    }

    U.resize(nPoints);
    T.resize(nPoints);
    Yb.resize(nSpec);
    Yu.resize(nSpec);
    Y.resize(nSpec, nPoints);

    double a = options.strainRateInitial;

    Tu = options.Tu;
    gas.pressure = options.pressure;

    grid.alpha = (options.curvedFlame) ? 1 : 0;
    grid.unburnedLeft = options.unburnedLeft;
    grid.updateValues();
    grid.updateBoundaryIndices();

    int jm = (grid.ju+grid.jb)/2; // midpoint of the profiles.

    // Make sure the initial profile fits in the domain
    double scale = 0.8*(x[jj]-x[0]) /
            (options.initialCenterWidth + 2 * options.initialSlopeWidth);
    if (scale < 1.0) {
        options.initialCenterWidth *= scale;
        options.initialSlopeWidth *= scale;
    }

    // Determine the grid indices that define each segment of the initial profile
    int centerPointCount = floor(0.5 + 0.5 * options.initialCenterWidth / dx);
    int slopePointCount = floor(0.5 + options.initialSlopeWidth / dx);
    int jl2 = jm - centerPointCount;
    int jl1 = jl2 - slopePointCount;
    int jr1 = jm + centerPointCount;
    int jr2 = jr1 + slopePointCount;
    int nSmooth = options.initialSmoothCount;

    if (options.flameType == "premixed") {
        // Reactants
        dvector reactants = gas.calculateReactantMixture(
            options.fuel, options.oxidizer, options.equivalenceRatio);
        gas.setStateMole(reactants, Tu);
        rhou = gas.getDensity();
        gas.getMassFractions(Yu);

        // Products
        Cantera::equilibrate(gas.thermo,"HP");
        Tb = gas.thermo.temperature();
        rhob = gas.getDensity();
        gas.thermo.getMassFractions(&Yb[0]);

        // Diluent in the middle
        gas.thermo.setState_TPY(Tu, gas.pressure, options.oxidizer);
        gas.thermo.getMassFractions(&Y(0,jm));

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

        T[0] = Tleft;
        T[grid.jj] = Tright;
        T[jm] = T[grid.ju];

    } else if (options.flameType == "diffusion") {
        // Stoichiometric mixture at the center
        options.equivalenceRatio = 1.0;
        dvector products = gas.calculateReactantMixture(
            options.fuel, options.oxidizer, options.equivalenceRatio);
        gas.setStateMole(products, 0.5*(options.Tfuel+options.Toxidizer));
        Cantera::equilibrate(gas.thermo,"HP");

        Tb = gas.thermo.temperature();
        rhob = gas.getDensity();
        gas.thermo.getMassFractions(&Yb[0]);
        gas.thermo.getMassFractions(&Y(0,jm));

        gas.thermo.setState_TPX(options.Tfuel, options.pressure, options.fuel);
        double rhoFuel = gas.getDensity();
        dvec Yfuel(nSpec);
        gas.getMassFractions(Yfuel);

        gas.thermo.setState_TPX(options.Toxidizer, options.pressure, options.oxidizer);
        double rhoOxidizer = gas.getDensity();
        dvec Yoxidizer(nSpec);
        gas.getMassFractions(Yoxidizer);

        if (options.fuelLeft) {
            rhoLeft = rhoFuel;
            Tleft = options.Tfuel;
            Yleft = Yfuel;
            rhoRight = rhoOxidizer;
            Tright = options.Toxidizer;
            Yright = Yoxidizer;
        } else {
            rhoLeft = rhoOxidizer;
            Tleft = options.Toxidizer;
            Yleft = Yoxidizer;
            rhoRight = rhoFuel;
            Tright = options.Tfuel;
            Yright = Yfuel;
        }

        rhou = rhoLeft; // arbitrary

        T[0] = Tleft;
        T[grid.jj] = Tright;
        T[jm] = Tb;

    } else {
        throw debugException("Invalid flameType: " + options.flameType);
    }

    for (size_t k=0; k<nSpec; k++) {
        Y(k, 0) = Yleft[k];
        Y(k, jj) = Yright[k];
    }

    for (int j=1; j<jl1; j++) {
        for (size_t k=0; k<nSpec; k++) {
            Y(k,j) = Y(k,0);
        }
        T[j] = T[0];
    }

    for (int j=jl1; j<jl2; j++) {
        double delta_x = x[jl2]-x[jl1];
        for (size_t k=0; k<nSpec; k++) {
            Y(k,j) = Y(k,0) + (Y(k,jm)-Y(k,0))*(x[j]-x[jl1])/delta_x;
        }
        T[j] = T[0] + (T[jm]-T[0])*(x[j]-x[jl1])/delta_x;
    }

    for (int j=jl2; j<jr1; j++) {
        for (size_t k=0; k<nSpec; k++) {
            Y(k,j) = Y(k,jm);
        }
        T[j] = T[jm];
    }

    for (int j=jr1; j<jr2; j++) {
        double delta_x = x[jr2]-x[jr1];
        for (size_t k=0; k<nSpec; k++) {
            Y(k,j) = Y(k,jm) + (Y(k,grid.jj)-Y(k,jm))*(x[j]-x[jr1])/delta_x;
        }
        T[j] = T[jm] + (T[grid.jj]-T[jm])*(x[j]-x[jr1])/delta_x;
    }

    for (size_t j=jr2; j<nPoints; j++) {
        for (size_t k=0; k<nSpec; k++) {
            Y(k,j) = Y(k,grid.jj);
        }
        T[j] = T[grid.jj];
    }

    dvector yTemp(nPoints);
    for (size_t k=0; k<nSpec; k++) {
        for (size_t j=0; j<nPoints; j++) {
            yTemp[j] = Y(k,j);
        }

        for (int i=0; i<nSmooth; i++) {
            mathUtils::smooth(yTemp);
        }

        for (size_t j=0; j<nPoints; j++) {
            Y(k,j) = yTemp[j];
        }
    }

    for (int i=0; i<nSmooth; i++) {
        mathUtils::smooth(T);
    }

    // Grid and initial profiles of T, U and V
    dvector rho(nPoints);
    for (size_t j=0; j<nPoints; j++) {
        gas.setStateMass(&Y(0,j), T[j]);
        rho[j] = gas.getDensity();
        U[j] = a*sqrt(rhou/rho[j]);
    }

    for (int i=0; i<2; i++) {
        mathUtils::smooth(U);
    }

    if (options.fixedLeftLoc)
    {
        jm = 0;
    }
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
        x = infile.readVector("x");

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
        dvector YbSave(nSpec);
        for (size_t k=0; k<nSpec; k++) {
            YbSave[k] = Y(k,grid.jb);
        }

        if (options.overrideReactants) {
            dvector reactants = gas.calculateReactantMixture(
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
        Yu.resize(nSpec); Yb.resize(nSpec);
        for (size_t k=0; k<nSpec; k++) {
            Yu[k] = Y(k,grid.ju);
            Yb[k] = Y(k,grid.jb);
        }

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
            for (size_t k=0; k<nSpec; k++) {
                Y(k,grid.jb) = YbSave[k];
            }
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
    configOptions& options = parent_->options;
    boost::ptr_vector<SourceSystem>& sourceTerms = parent_->sourceTerms;
    boost::ptr_vector<sundialsCVODE>& sourceSolvers = parent_->sourceSolvers;
    boost::ptr_vector<SourceSystemQSS>& sourceTermsQSS = parent_->sourceTermsQSS;

    CanteraGas** gasHandle = parent_->gases.acquire();
    CanteraGas* gas = *gasHandle;
    if (!gas->initialized()) {
        gas->setOptions(options);
        gas->initialize();
    }

    int err = 0;
    //std::cout << "pointRange: " << r.begin() << ", " << r.end() << std::endl;
    for (size_t j=r.begin(); j<r.end(); j++) {
        if (debugParameters::veryVerbose) {
            logFile.write(format("%i") % j, false);
        }
        if (parent_->useCVODE[j]) {
            sourceTerms[j].gas = gas;
            if (int(j) == options.debugSourcePoint && t_ >= options.debugSourceTime) {
                std::ofstream steps;
                steps.open("cvodeSteps.py");
                sourceTerms[j].writeState(sourceSolvers[j], steps, true);

                while (sourceSolvers[j].tInt < t_) {
                    err = sourceSolvers[j].integrateOneStep(t_);
                    sourceTerms[j].writeState(sourceSolvers[j], steps, false);
                    if (err != CV_SUCCESS) {
                        break;
                    }
                }

                sourceTerms[j].writeJacobian(sourceSolvers[j], steps);

                steps.close();
                std::terminate();

            } else {
                err = sourceSolvers[j].integrateToTime(t_);
            }
            if (err && err != CV_TSTOP_RETURN) {
                logFile.write(format("Error at j = %i") % j);
                logFile.write(format("T = %s") % sourceTerms[j].T);
                logFile.write(format("U = %s") % sourceTerms[j].U);
                logFile.write("Y = ", false);
                logFile.write(sourceTerms[j].Y);
                parent_->writeStateFile((format("prod%i_error_t%.6f_j%03i") %
                        stage_ % t_ % j).str(), true, false);
            }

            if (debugParameters::veryVerbose) {
                logFile.write(format(
                    " [%i]...") % sourceSolvers[j].getNumSteps(), false);
            }

        } else {
            sourceTermsQSS[j].gas = gas;
            if (int(j) == options.debugSourcePoint &&
                t_ >= options.debugSourceTime)
            {
                sourceTermsQSS[j].debug = true;
                std::ofstream steps;
                steps.open("cvodeSteps.py");
                sourceTermsQSS[j].writeState(steps, true);

                while (sourceTermsQSS[j].tn < (t_-parent_->tNow)) {
                    err = sourceTermsQSS[j].integrateOneStep(t_-parent_->tNow);
                    sourceTermsQSS[j].writeState(steps, false);
                    if (err) {
                        break;
                    }
                }

                steps.close();
                std::terminate();

            } else {
                err = sourceTermsQSS[j].integrateToTime(t_-parent_->tNow);
            }
            if (err) {
                logFile.write(format("Error at j = %i") % j);
                logFile.write(format("T = %s") % sourceTermsQSS[j].T);
                logFile.write(format("U = %s") % sourceTermsQSS[j].U);
                logFile.write("Y = ", false);
                logFile.write(sourceTermsQSS[j].Y);
                parent_->writeStateFile((format("prod%i_error_t%.6f_j%03i") %
                        stage_ % t_ % j).str(), true, false);
            }

            if (debugParameters::veryVerbose) {
                logFile.write(format(" [%i/%i]...") % sourceTermsQSS[j].gcount
                        % sourceTermsQSS[j].rcount, false);
            }

        }
    }

    parent_->gases.release(gasHandle);
}

void DiffusionTermWrapper::operator()(const tbb::blocked_range<size_t>& r) const
{
    for (size_t k=r.begin(); k<r.end(); k++) {
        parent_->diffusionSolvers[k].integrateToTime(t_);
    }
}
