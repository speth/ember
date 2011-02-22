#include "flameSolver.h"
#include "debugUtils.h"
#include "perfTimer.h"
#include "dataFile.h"
#include "mathUtils.h"

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#define foreach BOOST_FOREACH

using boost::format;

FlameSolver::FlameSolver()
    : jCorrSolver(jCorrSystem)
    , diffusionTestSolver(diffusionTestTerm)
{
    convectionSystem.setSpeciesDomains(convectionStartIndices, convectionStopIndices);
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
    //theSys.copyOptions(); TODO: copy to self?
    strainfunc.setOptions(options);

    // Cantera initialization
    gas.initialize();
    nSpec = gas.nSpec;
    nVars = nSpec + 2;
    W.resize(nSpec);
    gas.getMolecularWeights(W);

    // Chemkin & Adapchem Initialization
    if (options.usingAdapChem) {
        ckGas.reset(new AdapChem(options.inputDir+"/"+options.chemkinMechanismFile,
                                 true,
                                 options.inputDir+"/"+options.adapchemInputFile,
                                 options.inputDir+"/"+options.adapchemModelsFile,
                                 options.inputDir+"/"+options.adapchemDefaultModelFile,
                                 options.outputDir+"/"+options.adapchemDonemodelsFile,
                                 options.outputDir+"/"+options.adapchemRestartFile));
        ckGas->setPressure(options.pressure);
    }

    // Initial Conditions
    if (options.haveRestartFile) {
        loadProfile();
    } else {
        generateProfile();
    }

    grid.setSize(x.size());
    convectionSystem.setLeftBC(Tleft, Yleft);
    convectionSystem.setGas(gas);
    convectionSystem.setTolerances(options);

    for (size_t k=0; k<nVars; k++) {
        DiffusionSystem* term = new DiffusionSystem();
        BDFIntegrator* integrator = new BDFIntegrator(*term);
        integrator->resize(nPoints, 1, 1);
        diffusionTerms.push_back(term);
        diffusionSolvers.push_back(integrator);
    }

    resizeAuxiliary();
}

void FlameSolver::run(void)
{
    totalTimer.start();

    double t = tStart;
    double dt = options.globalTimestep;

    long int nTotal = 0; // total number of timesteps taken
    int nRegrid = 0; // number of time steps since regridding/adaptation
    int nOutput = 0; // number of time steps since storing integral flame parameters
    int nProfile = 0; // number of time steps since saving flame profiles
    int nIntegrate = 0; // number of time steps since restarting integrator
    int nTerminate = 1; // number of steps since last checking termination condition
    int nCurrentState = 0; // number of time steps since profNow.h5 and outNow.h5 were written

    // number of time steps since performing transport species elimination
    int nEliminate = options.transportEliminationStepInterval;

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
                    cout << format("WARNING: Unexpected Temperature: T = %f at j = %i")
                            % T[j] % j << endl;
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
            for (size_t k=0; k<nSpec; k++) {
                Y(k,0) = Yleft[k];
            }
        }
        T[jj] = Tright;
        for (size_t k=0; k<nSpec; k++) {
            Y(k,jj) = Yright[k];
        }

        updateChemicalProperties();

        if (options.usingAdapChem) {
            // TODO: Determine whether this should be called before or after initializeStep
            ckGas->incrementStep();

            // Because AdapChem only assigns values for species in the
            // reduced mechanisms, we need to zero the reaction rates
            // whenever the mechansim at each point could change.
            wDot.data().assign(wDot.data().size(), 0);
        }

        updateLeftBC();
        if (options.xFlameControl) {
            update_xStag(t, true); // calculate the value of rVzero
        }
        convectionSystem.set_rVzero(rVzero);

        setupTimer.stop();

        splitTimer.start();
        if (nEliminate >= options.transportEliminationStepInterval) {
            updateTransportDomain();
            nEliminate = 0;
        } else {
            nEliminate++;
        }

        for (size_t j=0; j<nPoints; j++) {
            diffusionTerms[kMomentum].B[j] = 1/rho[j];
            diffusionTerms[kEnergy].B[j] = 1/(rho[j]*cp[j]);

            diffusionTerms[kMomentum].D[j] = mu[j];
            diffusionTerms[kEnergy].D[j] = lambda[j];
        }

        // Diffusion solvers: Species
        for (size_t k=0; k<nSpec; k++) {
            DiffusionSystem& sys = diffusionTerms[kSpecies+k];
            size_t i = 0;
            sys.jLeft = diffusionStartIndices[k];
            sys.jRight = diffusionStopIndices[k];

            for (size_t j = diffusionStartIndices[k];
                 j <= diffusionStopIndices[k];
                 j++)
            {
                sys.B[i] = 1/rho[j];
                sys.D[i] = rhoD(k,j);
                i++;
            }
        }

        splitTimer.stop();

        if (t == tStart && options.outputProfiles) {
            writeStateFile();
        }

        if (options.outputDebugIntegratorStages) {
            writeStateFile((format("start_t%.6f") % t).str(), true, false);
        }

        double tNext = tNow + dt;

        // *** Strang-split Integration ***

        // Convection Integration: first half-step
        setConvectionSolverState(t, 1);
        integrateConvectionTerms(t + 0.5*dt, 1);
        extractConvectionState(1);

        // Diffusion Integration: first half-step
        setDiffusionSolverState(t);
        integrateDiffusionTerms(t + 0.5*dt, 1);
        extractDiffusionState(1);

        // Reaction Integration
        setProductionSolverState(t);

        if (options.outputSplitHeatReleaseRate) {
            // First half-step
            integrateProductionTerms(t + 0.5*dt, 1);
            extractProductionState(1);

            // Second half-step
            integrateProductionTerms(tNext, 2);
            extractProductionState(2);

        } else {
            // Full step
            integrateProductionTerms(tNext, 0);
            extractProductionState(2);
        }

        // Diffusion: second half-step
        setDiffusionSolverState(t + 0.5*dt); // assign initial conditions
        integrateDiffusionTerms(tNext, 2); // integrate
        extractDiffusionState(2); // extract solution components

        // Convection Integration: second half-step
        setConvectionSolverState(t + 0.5*dt, 2); // assign initial conditions
        integrateConvectionTerms(tNext, 2); // integrate
        extractConvectionState(2); // extract solution components

        if (debugParameters::veryVerbose) {
            cout << "done!" << endl;
        }

        // *** End of Strang-split integration step ***
        correctMassFractions();

        t = tNext;
        tNow = tNext;

        nOutput++;
        nRegrid++;
        nProfile++;
        nIntegrate++;
        nTerminate++;
        nCurrentState++;

        if (debugParameters::debugTimesteps) {
            int nSteps = convectionSystem.getNumSteps();
            cout << format("t = %8.6f (dt = %9.3e) [C: %i]") %
                    t % dt % nSteps << endl;
        }

        setupTimer.resume();
        if (t > tOutput || nOutput >= options.outputStepInterval) {
            calculateQdot();
            timeVector.push_back(t);
            timestepVector.push_back(dt);
            heatReleaseRate.push_back(getHeatReleaseRate());
            consumptionSpeed.push_back(getConsumptionSpeed());
            flamePosition.push_back(getFlamePosition());

            if (options.outputSplitHeatReleaseRate) {
                qDotProd1Vec.push_back(qDotProd1);
                qDotProd2Vec.push_back(qDotProd2);
                qDotDiff1Vec.push_back(qDotDiff1);
                qDotDiff2Vec.push_back(qDotDiff2);
                qDotConv1Vec.push_back(qDotConv1);
                qDotConv2Vec.push_back(qDotConv2);
            }

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

        // call ckGas::adapBox

        // *** Save the current integral and profile data
        //     in files that are automatically overwritten,
        //     and save the time-series data (out.h5)
        if (nCurrentState >= options.currentStateStepInterval) {
            nCurrentState = 0;
            writeTimeseriesFile("outNow");
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
            dvector& V_prev = convectionSystem.V;

            // "rollVectorVector"
            vector<dvector> currentSolution;
            currentSolution.push_back(U);
            currentSolution.push_back(T);
            for (size_t k=0; k<nSpec; k++) {
                dvector tmp(nPoints);
                for (size_t j=0; j<nPoints; j++) {
                    tmp[j] = Y(k,j);
                }
                currentSolution.push_back(tmp);
            }

            grid.regrid(currentSolution);
            grid.dampVal.resize(grid.x.size());

            // dampVal sets a limit on the maximum grid size
            for (size_t j=0; j<nPoints; j++) {
                double num = min(mu[j],lambda[j]/cp[j]);
                for (size_t k=0; k<nSpec; k++) {
                    num = min(num, rhoD(k,j));
                }
                grid.dampVal[j] = sqrt(num/(rho[j]*strainfunc.a(t)));
            }

            grid.adapt(currentSolution);

            // Perform updates that are necessary if the grid has changed
            if (grid.updated) {
                nIntegrate = 0;

                // Update species elimination at the start of the next time step
                nEliminate = options.transportEliminationStepInterval;

                cout << format("Grid size: %i points.") % nPoints << endl;

                // "unrollVectorVector"
                U.resize(nPoints);
                T.resize(nPoints);
                Y.resize(nSpec, nPoints);

                for (size_t j=0; j<nPoints; j++) {
                    U[j] = currentSolution[kMomentum][j];
                    T[j] = currentSolution[kEnergy][j];
                    for (size_t k=0; k<nSpec; k++) {
                        Y(k,j) = currentSolution[kSpecies+k][j];
                    }
                }

                correctMassFractions();

                // Update the mass flux at the left boundary
                rVzero = mathUtils::interp1(x_prev, V_prev, grid.x[0]);

                // Allocate the solvers and arrays for auxiliary variables
                resizeAuxiliary();

                if (debugParameters::debugAdapt || debugParameters::debugRegrid) {
                    writeStateFile("postAdapt", false, false);
                }
                grid.updated = false;
            }
            regridTimer.stop();
        }

        if (debugParameters::debugPerformanceStats && (nTotal % 10 == 0)) {
            printPerformanceStats();
        }
        nTotal++;
    }

    // *** Integration has reached the termination condition
    if (options.outputProfiles) {
        writeStateFile();
    }

    totalTimer.stop();
    printPerformanceStats();
    cout << format("Runtime: %f seconds.") % totalTimer.getTime() << endl;
}

bool FlameSolver::checkTerminationCondition(void)
{

    if (options.terminateForSteadyQdot) {
        int j1 = mathUtils::findLast(timeVector < (tNow - options.terminationPeriod));

        if (j1 == -1)
        {
            cout << format("Continuing integration: t (%8.6f) < terminationPeriod (%8.6f)")
                    % (tNow-timeVector[0]) % options.terminationPeriod << endl;
            return false;
        }

        int j2 = timeVector.size()-1;
        double qMean = mathUtils::mean(heatReleaseRate,j1,j2);
        double hrrError = 0;
        for (int j=j1; j<=j2; j++) {
            hrrError += abs(heatReleaseRate[j]-qMean);
        }
        hrrError /= (j2-j1+1);

        cout << format("Heat release rate deviation = %6.3f%%. hrrError = %9.4e")
                % (hrrError/qMean*100) % hrrError << endl;

        if (hrrError/abs(qMean) < options.terminationTolerance) {
            cout << "Terminating integration: ";
            cout << "Heat release deviation less than relative tolerance." << endl;
            return true;
        } else if (hrrError < options.terminationAbsTol) {
            cout << "Terminating integration: ";
            cout << "Heat release rate deviation less than absolute tolerance." << endl;
            return true;
        } else if (tNow-tStart > options.terminationMaxTime ) {
          cout << "Terminating integration: Maximum integration time reached." << endl;
          return true;
        } else {
            cout << format("Continuing integration. t = %8.6f") % (tNow-timeVector[0]) << endl;
        }

    }
    return false;
}

void FlameSolver::writeStateFile(const std::string fileNameStr, bool errorFile, bool updateDerivatives)
{
    std::ostringstream fileName(ostringstream::out);
    bool incrementFileNumber = false;

    if (fileNameStr.length() == 0) {
        // Determine the name of the output file (outXXXXXX.h5)
        incrementFileNumber = true;
        if (errorFile) {
            fileName << options.outputDir << "/error";
        } else {
            fileName << options.outputDir << "/prof";
        }
        fileName.flags(ios_base::right);
        fileName.fill('0');
        fileName.width(6);
        fileName << options.outputFileNumber << ".h5";
    } else {
        fileName << options.outputDir << "/" << fileNameStr << ".h5";
    }
    if (errorFile) {
        cout << "Writing error output file: " << fileName.str() << endl;
    } else {
        cout << "Writing output file: " << fileName.str() << endl;
    }

    // Erase the existing file and create a new one
    if (boost::filesystem::exists(fileName.str())) {
        boost::filesystem::remove(fileName.str());
    }
    DataFile outFile(fileName.str());

    updateChemicalProperties();
    convectionSystem.evaluate();

    // Write the state data to the output file:
    outFile.writeScalar("t", tNow);
    outFile.writeVector("x", x);
    outFile.writeVector("T", T);
    outFile.writeVector("U", U);
    outFile.writeArray2D("Y", Y);
    outFile.writeScalar("a", strainfunc.a(tNow));
    outFile.writeScalar("dadt", strainfunc.dadt(tNow));
    outFile.writeScalar("fileNumber", options.outputFileNumber);

    if (updateDerivatives && (options.outputTimeDerivatives || options.outputAuxiliaryVariables)) {
        calculateTimeDerivatives();
    }

    if (options.outputHeatReleaseRate || errorFile) {
        outFile.writeVector("q", qDot);
        outFile.writeVector("rho", rho);
    }

    if (options.outputTimeDerivatives || errorFile) {
        outFile.writeVector("dUdt", dUdt);
        outFile.writeVector("dTdt", dTdt);
        outFile.writeArray2D("dYdt", dYdt);
    }

    if (options.outputAuxiliaryVariables || errorFile) {
        outFile.writeVector("V", convectionSystem.V);
        outFile.writeArray2D("wdot", wDot);
        outFile.writeArray2D("rhoD", rhoD);
        outFile.writeVector("lambda", lambda);

        outFile.writeVector("cp", cp);
        outFile.writeVector("mu", mu);
        outFile.writeVector("Wmx", Wmx);
        outFile.writeVector("W", W);
        outFile.writeVector("sumcpj", sumcpj);
        outFile.writeArray2D("jFick", jFick);
        outFile.writeArray2D("jSoret", jSoret);
        outFile.writeVector("jCorr", jCorr);
        outFile.writeVector("cfp", grid.cfp);
        outFile.writeVector("cf", grid.cf);
        outFile.writeVector("cfm", grid.cfm);
        outFile.writeVector("hh", hh);
        outFile.writeVector("rphalf", grid.rphalf);
        outFile.writeScalar("Tleft", Tleft);
        outFile.writeVector("Yleft", Yleft);

        if (options.usingAdapChem) {
            dvector nSpecReduced(nPoints);
            for (size_t j=0; j<nPoints; j++) {
                nSpecReduced[j] = ckGas->getNumSpec(j);
            }
            outFile.writeVector("nSpecReduced", nSpecReduced);
        }

        // Number of timesteps in the chemistry solver in the last global timestep
        dvector chemSteps(nPoints, 0);
        for (size_t j=0; j<nPoints; j++) {
            if (sourceSolvers[j].initialized()) {
                chemSteps[j] = sourceSolvers[j].getNumSteps();
            }
        }
        outFile.writeVector("chemSteps", chemSteps);

        if (options.transportEliminationDiffusion) {
            dvector jStart, jStop;
            for (size_t k=0; k<nSpec; k++) {
                jStart.push_back(diffusionStartIndices[k]);
                jStop.push_back(diffusionStopIndices[k]);
            }
            outFile.writeVector("diffusionStartIndices", jStart);
            outFile.writeVector("diffusionStopIndices", jStop);
        }

        if (options.transportEliminationConvection) {
            dvector jStart, jStop;
            for (size_t k=0; k<nSpec; k++) {
                jStart.push_back(convectionStartIndices[k]);
                jStop.push_back(convectionStopIndices[k]);
            }
            outFile.writeVector("convectionStartIndices", jStart);
            outFile.writeVector("convectionStopIndices", jStop);
        }
    }

    if (options.outputResidualComponents || errorFile) {

        outFile.writeVector("dUdtDiff", diffusionSolvers[kMomentum].get_ydot());
        outFile.writeVector("dUdtConv", convectionSystem.dUdt);
        outFile.writeVector("dUdtProd", dUdtProd);

        outFile.writeVector("dTdtDiff", diffusionSolvers[kEnergy].get_ydot());
        outFile.writeVector("dTdtConv", convectionSystem.dTdt);
        outFile.writeVector("dTdtProd", dTdtProd);
        outFile.writeVector("dTdtCross", dTdtCross);

        outFile.writeArray2D("dYdtConv", convectionSystem.dYdt);
        outFile.writeArray2D("dYdtDiff", dYdtDiff);
        outFile.writeArray2D("dYdtProd", dYdtProd);
        outFile.writeArray2D("dYdtCross", dYdtCross);

        outFile.writeVector("dWdt", convectionSystem.dWdt);
        outFile.writeVector("dWdx", convectionSystem.utwSystem.dWdx);
        outFile.writeVector("dTdx", convectionSystem.utwSystem.dTdx);
        outFile.writeVector("dWdtSplit", convectionSystem.utwSystem.dWdtSplit);
        outFile.writeVector("dTdtSplit", convectionSystem.utwSystem.dTdtSplit);
    }

    outFile.close();
    if (incrementFileNumber) {
        options.outputFileNumber++;
    }

    if (errorFile && options.stopIfError) {
      cout << "Error outputs remaining until termination: " << options.errorStopCount << endl;
      if (options.errorStopCount-- <= 0) {
        throw debugException("Too many integration failures.");
      }
    }
}

void FlameSolver::writeTimeseriesFile(const std::string filename)
{
    DataFile outFile(options.outputDir+"/"+filename+".h5");
    outFile.writeVector("t", timeVector);
    outFile.writeVector("dt", timestepVector);
    outFile.writeVector("Q", heatReleaseRate);
    outFile.writeVector("Sc", consumptionSpeed);
    outFile.writeVector("xFlame", flamePosition);

    if (options.outputSplitHeatReleaseRate) {
        outFile.writeVector("qDotProd1", qDotProd1Vec);
        outFile.writeVector("qDotProd2", qDotProd2Vec);
        outFile.writeVector("qDotDiff1", qDotDiff1Vec);
        outFile.writeVector("qDotDiff2", qDotDiff2Vec);
        outFile.writeVector("qDotConv1", qDotConv1Vec);
        outFile.writeVector("qDotConv2", qDotConv2Vec);
    }

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

    dUdt.resize(nPoints,0);
    dTdt.resize(nPoints,0);
    dYdt.resize(nSpec,nPoints);

    dTdtCross.resize(nPoints,0);
    dYdtCross.resize(nSpec, nPoints, 0);

    rho.resize(nPoints);
    Wmx.resize(nPoints);
    mu.resize(nPoints);
    lambda.resize(nPoints);
    cp.resize(nPoints);
    jCorr.resize(nPoints);
    sumcpj.resize(nPoints);
    qDot.resize(nPoints);
    cpSpec.resize(nSpec, nPoints);
    rhoD.resize(nSpec, nPoints);
    Dkt.resize(nSpec, nPoints);
    wDot.resize(nSpec, nPoints);
    hk.resize(nSpec, nPoints);
    jFick.resize(nSpec, nPoints);
    jSoret.resize(nSpec, nPoints);

    grid.jj = nPoints-1;
    grid.updateBoundaryIndices();
    diffusionStartIndices.assign(nSpec, 0);
    diffusionStopIndices.assign(nSpec, jj);
    nPointsDiffusion.assign(nSpec, nPoints);

    convectionStartIndices.assign(nSpec, 0);
    convectionStopIndices.assign(nSpec, jj);
    nPointsConvection.assign(nSpec, nPoints);

    if (options.usingAdapChem) {
        ckGas->setGridSize(nPoints);
    }

    if (nPoints > nPointsOld) {
        for (size_t j=nPointsOld; j<nPoints; j++) {
            // Create and initialize the new SourceSystem
            SourceSystem* system = new SourceSystem();
            system->resize(nSpec);
            system->gas = &gas;
            system->ckGas = ckGas;
            system->usingAdapChem = options.usingAdapChem;
            system->thermoTimer = &thermoTimer;
            system->reactionRatesTimer = &reactionRatesTimer;
            system->jacobianTimer = &jacobianTimer;
            system->strainFunction = strainfunc;
            system->rhou = rhou;
            system->W = W;
            system->j = j;

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
            qssSolver->resize(nSpec);
            qssSolver->gas = &gas;
            qssSolver->ckGas = ckGas;
            qssSolver->usingAdapChem = options.usingAdapChem;
            qssSolver->thermoTimer = &thermoTimer;
            qssSolver->reactionRatesTimer = &reactionRatesTimer;
            qssSolver->strainFunction = strainfunc;
            qssSolver->rhou = rhou;
            qssSolver->W = W;
            qssSolver->j = j;

            qssSolver->ymin.assign(nVars, 1e-24);
            qssSolver->ymin[kMomentum] = -1e4;
            qssSolver->epsmin = options.integratorRelTol;
            qssSolver->dtmin = options.integratorMinTimestep;

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
    // With transport species elimination enabled, this is done later
    // to handle the varying number of grid points for each species
    if (options.transportEliminationDiffusion || options.transportEliminationConvection) {
        diffusionTestSolver.resize(nPoints, 1, 1);
        diffusionTestTerm.setGrid(grid);
    } else {
        for (size_t k=0; k<nVars; k++) {
            diffusionSolvers[k].resize(nPoints, 1, 1);
        }
    }

    for (size_t k=0; k<nVars; k++) {
        diffusionTerms[k].setGrid(grid);
    }

    // Diffusion solvers: Energy and Momentum
    diffusionSolvers[kMomentum].resize(nPoints, 1, 1);
    diffusionSolvers[kEnergy].resize(nPoints, 1, 1);

    convectionSystem.resize(nPoints, nPointsConvection, nSpec);
    convectionSystem.setGrid(grid);

    // Resize the jCorr stabilizer
    jCorrSolver.resize(nPoints, 1, 1);
    jCorrSystem.setGrid(grid);
}

void FlameSolver::updateCrossTerms()
{
    assert(mathUtils::notnan(Y.data()));
    assert(mathUtils::notnan(T));
    assert(mathUtils::notnan(rhoD.data()));

    for (size_t j=0; j<jj; j++) {
        jCorr[j] = 0;
        for (size_t k=0; k<nSpec; k++) {
            jFick(k,j) = -0.5*(rhoD(k,j)+rhoD(k,j+1)) * ((Y(k,j+1)-Y(k,j))/hh[j]);
            jSoret(k,j) = -0.5*(Dkt(k,j)/T[j] + Dkt(k,j+1)/T[j+1])
                * (T[j+1]-T[j])/hh[j];
            jCorr[j] -= jFick(k,j) + jSoret(k,j);
        }
    }
    assert(mathUtils::notnan(jFick.data()));
    assert(mathUtils::notnan(jSoret.data()));
    assert(mathUtils::notnan(lambda));
    assert(mathUtils::notnan(rho));
    assert(mathUtils::notnan(cp));

    // Add a bit of artificial diffusion to jCorr to improve stability
    jCorrSolver.y = jCorr;
    for (size_t j=0; j<jj; j++) {
        jCorrSystem.B[j] = lambda[j]/(rho[j]*cp[j]); // treat as Le = 1
        jCorrSystem.D[j] = 1.0;
    }

    jCorrSolver.initialize(0, options.diffusionTimestep);
    jCorrSolver.integrateToTime(options.globalTimestep);
    jCorr = jCorrSolver.y;

    assert(mathUtils::notnan(jCorr));

    for (size_t j=1; j<jj; j++) {
        sumcpj[j] = 0;
        for (size_t k=0; k<nSpec; k++) {
            dYdtCross(k,j) = -0.5/(r[j]*rho[j]*dlj[j]) *
                (rphalf[j]*(Y(k,j)+Y(k,j+1))*jCorr[j] - rphalf[j-1]*(Y(k,j-1)+Y(k,j))*jCorr[j-1]);
            dYdtCross(k,j) -= 1/(r[j]*rho[j]*dlj[j]) *
                (rphalf[j]*jSoret(k,j) - rphalf[j-1]*jSoret(k,j-1));
            sumcpj[j] += 0.5*(cpSpec(k,j)+cpSpec(k,j+1))/W[k]*(jFick(k,j)
                    + jSoret(k,j) + 0.5*(Y(k,j)+Y(k,j+1))*jCorr[j]);
        }
        double dTdx = cfm[j]*T[j-1] + cf[j]*T[j] + cfp[j]*T[j+1];
        dTdtCross[j] = - 0.5*(sumcpj[j]+sumcpj[j-1]) * dTdx / (cp[j]*rho[j]);
    }

    assert(mathUtils::notnan(dYdtCross.data()));
    assert(mathUtils::notnan(sumcpj));
    assert(mathUtils::notnan(dTdtCross));
}

void FlameSolver::updateLeftBC()
{
    BoundaryCondition::BC prev = grid.leftBC;

    if ((options.twinFlame || options.curvedFlame) && x[0] >= 0.0 && x[0] <= options.centerGridMin) {
        grid.leftBC = BoundaryCondition::ControlVolume;
    } else if (grid.ju == 0 || (grid.jb == 0 && grid.fixedBurnedVal)) {
        grid.leftBC = BoundaryCondition::FixedValue;
    } else {
        grid.leftBC = BoundaryCondition::ZeroGradient;
    }

    if (prev != grid.leftBC) {
        cout << "updateLeftBC: BC changed from " << prev << " to " << grid.leftBC << "." << endl;
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

        // Reaction rates
        reactionRatesTimer.start();
        if (options.usingAdapChem) {
            ckGas->initializeStep(&Y(0,j), T[j], j);
            ckGas->getReactionRates(&wDot(0,j));
        } else {
            gas.getReactionRates(&wDot(0,j));
        }
        qDot[j] = 0;
        for (size_t k=0; k<nSpec; k++) {
            qDot[j] -= wDot(k,j)*hk(k,j);
        }
        reactionRatesTimer.stop();
    }
}

void FlameSolver::setDiffusionSolverState(double tInitial)
{
    splitTimer.resume();
    for (size_t k=0; k<nVars; k++) {
        // TODO: Use timestep that is based on each component's diffusivity
        diffusionSolvers[k].initialize(tInitial, options.diffusionTimestep);
    }

    updateCrossTerms();
    diffusionSolvers[kMomentum].y = U;
    diffusionSolvers[kEnergy].y = T;
    diffusionTerms[kEnergy].splitConst = dTdtCross;
    for (size_t k=0; k<nSpec; k++) {
        dvector& yDiff_Y = diffusionSolvers[kSpecies+k].y;
        dvector& splitConstY = diffusionTerms[kSpecies+k].splitConst;
        size_t i = 0;
        for (size_t j = diffusionStartIndices[k];
             j <= diffusionStopIndices[k];
             j++)
        {
            yDiff_Y[i] = Y(k,j);
            splitConstY[i] = dYdtCross(k,j);
            i++;
        }
    }
    splitTimer.stop();
}

void FlameSolver::setConvectionSolverState(double tInitial, int stage)
{
    assert(stage == 1 || stage == 2);

    splitTimer.resume();
    convectionSystem.setState(U, T, Y);
    convectionSystem.initialize(tInitial);
    splitTimer.stop();
    setDiffusionSolverState(tInitial);

    if (stage == 1) {
        calculateSplitDerivatives(tInitial);
    }
}

void FlameSolver::setProductionSolverState(double tInitial)
{
    splitTimer.resume();
    for (size_t j=0; j<nPoints; j++) {
        if (useCVODE[j]) {
            sdVector& ySource = sourceSolvers[j].y;
            ySource[kMomentum] = U[j];
            ySource[kEnergy] = T[j];
            for (size_t k=0; k<nSpec; k++) {
                ySource[kSpecies+k] = Y(k,j);
            }
            sourceSolvers[j].t0 = tInitial;
            sourceSolvers[j].initialize();
            sourceTerms[j].wDot.assign(nSpec, 0);
            sourceTerms[j].strainFunction.pin(tInitial);
        } else {
            dvector ySource(nVars);
            ySource[kMomentum] = U[j];
            ySource[kEnergy] = T[j];
            for (size_t k=0; k<nSpec; k++) {
                ySource[kSpecies+k] = Y(k,j);
            }
            sourceTermsQSS[j].initialize(ySource, tInitial);
            sourceTermsQSS[j].wDotQ.assign(nSpec, 0);
            sourceTermsQSS[j].wDotD.assign(nSpec, 0);
            sourceTermsQSS[j].strainFunction.pin(tInitial);
        }
    }
    splitTimer.stop();
}

void FlameSolver::calculateSplitDerivatives(double t)
{
    Array2D dYdtSplit(nSpec, nPoints, 0);

    // Evaluate derivatives from the diffusion terms
    dvector dTdtSplit = diffusionSolvers[kEnergy].get_ydot();
    for (size_t k=0; k<nSpec; k++) {
        const dvector& ydotk = diffusionSolvers[kSpecies+k].get_ydot();
        size_t i = 0;
        for (size_t j = diffusionStartIndices[k];
                         j <= diffusionStopIndices[k];
                         j++)
        {
            dYdtSplit(k,j) = ydotk[i];
            i++;
        }
    }

    // Evaluate the derivatives from the source terms
    for (size_t j=0; j<nPoints; j++) {
        if (useCVODE[j]) {
            sdVector& ySource = sourceSolvers[j].y;
            ySource[kMomentum] = U[j];
            ySource[kEnergy] = T[j];
            for (size_t k=0; k<nSpec; k++) {
                ySource[kSpecies+k] = Y(k,j);
            }
            sourceTerms[j].wDot.assign(nSpec, 0);
            sourceTerms[j].strainFunction.pin(t);

            sdVector ydotSource(nVars);
            sourceTerms[j].f(t, ySource, ydotSource);
            dTdtSplit[j] += ydotSource[kEnergy];
            for (size_t k=0; k<nSpec; k++) {
                dYdtSplit(k,j) += ydotSource[kSpecies+k];
            }
        } else {
            dvector ySource(nVars);
            ySource[kMomentum] = U[j];
            ySource[kEnergy] = T[j];
            for (size_t k=0; k<nSpec; k++) {
                ySource[kSpecies+k] = Y(k,j);
            }
            sourceTermsQSS[j].wDotQ.assign(nSpec, 0);
            sourceTermsQSS[j].wDotD.assign(nSpec, 0);
            sourceTermsQSS[j].strainFunction.pin(t);

            dvector q(nVars), d(nVars);
            sourceTermsQSS[j].odefun(t, ySource, q, d);
            dTdtSplit[j] += (q[kEnergy] - d[kEnergy]);
            for (size_t k=0; k<nSpec; k++) {
                dYdtSplit(k,j) += (q[kSpecies+k] - d[kSpecies+k]);
            }
        }
    }

    convectionSystem.setSplitDerivatives(dTdtSplit, dYdtSplit);
}

void FlameSolver::extractConvectionState(int stage)
{
    assert(stage == 1 || stage == 2);

    splitTimer.resume();
    convectionSystem.unroll_y();
    U = convectionSystem.U;
    T = convectionSystem.T;
    for (size_t k=0; k<nSpec; k++) {
        for (size_t j = convectionStartIndices[k];
             j <= convectionStopIndices[k];
             j++)
        {
            Y(k,j) = convectionSystem.Y(k,j);
        }
    }
    splitTimer.stop();

    if (options.outputSplitHeatReleaseRate) {
        calculateQdot();
        (stage == 1) ? qDotConv1 : qDotConv2 = getHeatReleaseRate();
    }

    if (options.outputDebugIntegratorStages) {
        writeStateFile((format("conv%i_t%.6f") % stage % tNow).str(), true, false);
    }
}

void FlameSolver::extractDiffusionState(int stage)
{
    assert(stage == 1 || stage == 2);

    splitTimer.resume();
    U = diffusionSolvers[kMomentum].y;
    T = diffusionSolvers[kEnergy].y;
    for (size_t k=0; k<nSpec; k++) {
        size_t i = 0;
        const BDFIntegrator& solver = diffusionSolvers[kSpecies+k];
        for (size_t j = diffusionStartIndices[k];
             j <= diffusionStopIndices[k];
             j++)
        {
            Y(k,j) = solver.y[i];
            i++;
        }
    }
    splitTimer.stop();

    if (options.outputSplitHeatReleaseRate) {
        calculateQdot();
        (stage == 1) ? qDotDiff1 : qDotDiff2 = getHeatReleaseRate();
    }

    if (options.outputDebugIntegratorStages) {
        writeStateFile((format("diff%i_t%.6f") % stage % tNow).str(), true, false);
    }
}

void FlameSolver::extractProductionState(int stage)
{
    assert(stage == 1 || stage == 2);

    splitTimer.resume();
    for (size_t j=0; j<nPoints; j++) {
        if (useCVODE[j]) {
            sourceTerms[j].unroll_y(sourceSolvers[j].y);
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
    splitTimer.stop();

    if (options.outputSplitHeatReleaseRate) {
        calculateQdot();
        (stage == 1) ? qDotProd1 : qDotProd2 = getHeatReleaseRate();
    }

    if (options.outputDebugIntegratorStages) {
        writeStateFile((format("prod%i_t%.6f") % stage % tNow).str(), true, false);
    }
}

void FlameSolver::integrateConvectionTerms(double t, int stage)
{
    if (debugParameters::veryVerbose) {
        cout << format("convection term %i/2...") % stage;
        cout.flush();
    }

    convectionTimer.start();
    convectionSystem.integrateToTime(t);
    convectionTimer.stop();
}

void FlameSolver::integrateProductionTerms(double t, int stage)
{
    if (debugParameters::veryVerbose) {
        if (stage) {
            cout << format("Source term %i/2...") % stage;
        } else {
            cout << format("Source term...");
        }
        cout.flush();
    }

    reactionTimer.start();
    int err = 0;
    for (size_t j=0; j<nPoints; j++) {
        if (debugParameters::veryVerbose) {
            cout << j;
            cout.flush();
        }
        if (useCVODE[j]) {
            if (j == options.debugSourcePoint && t >= options.debugSourceTime) {
                ofstream steps;
                steps.open("cvodeSteps.py");
                sourceTerms[j].writeState(sourceSolvers[j], steps, true);

                while (sourceSolvers[j].tInt < t) {
                    err = sourceSolvers[j].integrateOneStep(t);
                    sourceTerms[j].writeState(sourceSolvers[j], steps, false);
                    if (err != CV_SUCCESS) {
                        break;
                    }
                }

                sourceTerms[j].writeJacobian(sourceSolvers[j], steps);

                steps.close();
                terminate();

            } else {
                err = sourceSolvers[j].integrateToTime(t);
            }
            if (err && err != CV_TSTOP_RETURN) {
                cout << "Error at j = " << j << endl;
                cout << "T = " << sourceTerms[j].T << endl;
                cout << "U = " << sourceTerms[j].U << endl;
                cout << "Y = " << sourceTerms[j].Y << endl;
                writeStateFile((format("prod%i_error_t%.6f_j%03i") % stage % t % j).str(), true, false);
            }

            if (debugParameters::veryVerbose) {
                cout << " [" << sourceSolvers[j].getNumSteps() << "]...";
                cout.flush();
            }

        } else {
            if (j == options.debugSourcePoint && t >= options.debugSourceTime) {
                ofstream steps;
                steps.open("cvodeSteps.py");
                sourceTermsQSS[j].writeState(steps, true);

                while (sourceTermsQSS[j].tn < (t-tNow)) {
                    err = sourceTermsQSS[j].integrateOneStep(t-tNow);
                    sourceTermsQSS[j].writeState(steps, false);
                    if (err) {
                        break;
                    }
                }

                steps.close();
                terminate();

            } else {
                err = sourceTermsQSS[j].integrateToTime(t-tNow);
            }
            if (err) {
                cout << "Error at j = " << j << endl;
                cout << "T = " << sourceTermsQSS[j].T << endl;
                cout << "U = " << sourceTermsQSS[j].U << endl;
                cout << "Y = " << sourceTermsQSS[j].Y << endl;
                writeStateFile((format("prod%i_error_t%.6f_j%03i") % stage % t % j).str(), true, false);
            }

            if (debugParameters::veryVerbose) {
                cout << " [" << sourceTermsQSS[j].gcount << "]...";
                cout.flush();
            }

        }
    }
    reactionTimer.stop();
}

void FlameSolver::integrateDiffusionTerms(double t, int stage)
{
    if (debugParameters::veryVerbose) {
        cout << format("diffusion terms %i/2...") % stage;
        cout.flush();
    }

    diffusionTimer.start();
    for (size_t k=0; k<nVars; k++) {
        diffusionSolvers[k].integrateToTime(t);
    }
    diffusionTimer.stop();
}

void FlameSolver::update_xStag(const double t, const bool updateIntError)
{
    xFlameActual = getFlamePosition();
    xFlameTarget = targetFlamePosition(t);
    if (updateIntError) {
        flamePosIntegralError += (xFlameTarget-xFlameActual)*(t-tFlamePrev);
        tFlamePrev = t;
    }

    // controlSignal is approximately a*xStag
    double controlSignal = options.xFlameProportionalGain *
        ( (xFlameTarget-xFlameActual) + (flamePosIntegralError + (xFlameTarget-xFlameActual)*(t-tFlamePrev))*options.xFlameIntegralGain );

    if (debugParameters::debugFlameRadiusControl) {
        cout << "rFlameControl: " << "rF = " << xFlameActual << "   control = " << controlSignal;
        cout << "   P = " <<  options.xFlameProportionalGain*(xFlameTarget-xFlameActual);
        cout << "   I = " << options.xFlameProportionalGain*flamePosIntegralError*options.xFlameIntegralGain;
        cout << "  dt = " << t-tFlamePrev << endl;
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
    return (t <= options.xFlameT0) ? options.xFlameInitial
        :  (t >= options.xFlameT0+options.xFlameDt) ? options.xFlameFinal
        : options.xFlameInitial + (options.xFlameFinal-options.xFlameInitial)*(t-options.xFlameT0)/options.xFlameDt;
}

void FlameSolver::calculateQdot()
{
    for (size_t j=0; j<nPoints; j++) {
        gas.setStateMass(&Y(0,j), T[j]);
        gas.getEnthalpies(&hk(0,j));

        if (options.usingAdapChem) {
            ckGas->initializeStep(&Y(0,j), T[j], j);
            ckGas->getReactionRates(&wDot(0,j));
        } else {
            gas.getReactionRates(&wDot(0,j));
        }
        qDot[j] = 0;
        for (size_t k=0; k<nSpec; k++) {
            qDot[j] -= wDot(k,j)*hk(k,j);
        }
    }
}

void FlameSolver::calculateTimeDerivatives()
{
    dYdtDiff.resize(nSpec, nPoints);
    dYdtProd.resize(nSpec, nPoints);
    dTdtProd.resize(nPoints);
    dUdtProd.resize(nPoints);

    // Diffusion term contribution
    setDiffusionSolverState(tNow);
    dUdt = diffusionSolvers[kMomentum].get_ydot();
    dTdt = diffusionSolvers[kEnergy].get_ydot();
    dYdt.data().assign(dYdt.data().size(), 0);
    for (size_t k=0; k<nSpec; k++) {
        const dvector& ydot = diffusionSolvers[kSpecies+k].get_ydot();

        size_t i = 0;
        for (size_t j = diffusionStartIndices[k];
             j <= diffusionStopIndices[k];
             j++)
        {
            dYdt(k,j) = dYdtDiff(k,j) = ydot[i];
            i++;
        }
    }

    // Production term contribution
    setProductionSolverState(tNow);
    for (size_t j=0; j<nPoints; j++) {
        if (useCVODE[j]) {
            sdVector& ySource = sourceSolvers[j].y;
            sdVector ydotSource(nVars);
            sourceTerms[j].f(tNow, ySource, ydotSource);

            dUdt[j] += dUdtProd[j] = ydotSource[kMomentum];
            dTdt[j] += dTdtProd[j] = ydotSource[kEnergy];
            for (size_t k=0; k<nSpec; k++) {
                dYdt(k,j) += dYdtProd(k,j) = ydotSource[kSpecies+k];
            }
        } else {
            dvector& ySource = sourceTermsQSS[j].y;
            dvector q(nVars), d(nVars);
            sourceTermsQSS[j].odefun(tNow, ySource, q, d);

            dUdt[j] += dUdtProd[j] = (q[kMomentum] - d[kMomentum]);
            dTdt[j] += dTdtProd[j] = (q[kEnergy] - d[kEnergy]);
            for (size_t k=0; k<nSpec; k++) {
                dYdt(k,j) += dYdtProd(k,j) = (q[kSpecies+k] - d[kSpecies+k]);
            }
        }
    }

    // Convection term contribution
    setConvectionSolverState(tNow, 1);
    convectionSystem.evaluate();
    dUdt += convectionSystem.dUdt;
    dTdt += convectionSystem.dTdt;
    dYdt += convectionSystem.dYdt;
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

double FlameSolver::getHeatReleaseRate(void)
{
    return mathUtils::integrate(x, qDot);
}

double FlameSolver::getConsumptionSpeed(void)
{
    double QoverCp = mathUtils::integrate(x,qDot/cp);
    double rhouDeltaT = rhou*(Tb-Tu);
    return QoverCp/rhouDeltaT;
}

double FlameSolver::getFlamePosition(void)
{
    return mathUtils::trapz(x,x*qDot)/mathUtils::trapz(x,qDot);
}

void FlameSolver::generateProfile(void)
{
    cout << "Generating initial profiles from given fuel and oxidizer compositions." << endl;

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

    grid.alpha = (options.curvedFlame) ? 1 : 0;
    grid.unburnedLeft = options.unburnedLeft;
    grid.updateValues();
    grid.updateBoundaryIndices();

    int jm = (grid.ju+grid.jb)/2; // midpoint of the profiles.
    int jl = jm - 4;
    int jr = jm + 4;

    U.resize(nPoints);
    T.resize(nPoints);
    Yb.resize(nSpec);
    Yu.resize(nSpec);
    Y.resize(nSpec, nPoints);

    double a = options.strainRateInitial;

    Tu = options.Tu;
    gas.pressure = options.pressure;

    // Reactants
    dvector reactants = calculateReactantMixture();
    gas.setStateMole(reactants, Tu);
    rhou = gas.getDensity();
    gas.getMassFractions(Yu);
    gas.getMassFractions(&Y(0,grid.ju));

    // Products
    Cantera::equilibrate(gas.thermo,"HP");
    Tb = gas.thermo.temperature();
    rhob = gas.getDensity();
    gas.thermo.getMassFractions(&Yb[0]);
    gas.thermo.getMassFractions(&Y(0,grid.jb));

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

    T[0] = Tleft; T[grid.jj] = Tright;
    T[jm] = T[grid.ju];


    for (int j=1; j<jl; j++) {
        for (size_t k=0; k<nSpec; k++) {
            Y(k,j) = Y(k,0);
        }
        T[j] = T[0];
    }

    for (int j=jl; j<jm; j++) {
        for (size_t k=0; k<nSpec; k++) {
            Y(k,j) = Y(k,0) + (Y(k,jm)-Y(k,0))*(x[j]-x[jl])/(x[jm]-x[jl]);
        }
        T[j] = T[0] + (T[jm]-T[0])*(x[j]-x[jl])/(x[jm]-x[jl]);
    }

    for (int j=jm+1; j<jr; j++) {
        for (size_t k=0; k<nSpec; k++) {
            Y(k,j) = Y(k,jm) + (Y(k,grid.jj)-Y(k,jm))*(x[j]-x[jm])/(x[jr]-x[jm]);
        }
        T[j] = T[jm] + (T[grid.jj]-T[jm])*(x[j]-x[jm])/(x[jr]-x[jm]);
    }

    for (size_t j=jr; j<nPoints; j++) {
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

        for (size_t i=0; i<10; i++) {
            mathUtils::smooth(yTemp);
        }

        for (size_t j=0; j<nPoints; j++) {
            Y(k,j) = yTemp[j];
        }
    }

    for (size_t i=0; i<5; i++) {
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

    // TODO: Generate a better profile for the curved flame case
    dvector V(nPoints);
    V[jm] = 0;
    for (size_t j=jm+1; j<nPoints; j++) {
        V[j] = V[j-1] - rho[j]*U[j]*(x[j]-x[j-1]);
    }

    if (jm != 0) {
        for (size_t j=jm; j>0; j--) {
            V[j-1] = V[j] + rho[j-1]*U[j-1]*(x[j]-x[j-1]);
        }
    }

    rVzero = V[0];
}

void FlameSolver::loadProfile(void)
{
    std::string inputFilename;
    if (options.useRelativeRestartPath) {
        inputFilename = options.inputDir + "/" + options.restartFile;
    } else {
        inputFilename = options.restartFile;
    }

    cout << "Reading initial condition from " << inputFilename << endl;
    DataFile infile(inputFilename);
    x = infile.readVector("x");

    grid.setSize(x.size());
    grid.alpha = (options.curvedFlame) ? 1 : 0;
    grid.unburnedLeft = options.unburnedLeft;
    grid.updateValues();
    grid.updateBoundaryIndices();

    U = infile.readVector("U");
    T = infile.readVector("T");
    Y = infile.readArray2D("Y");
    tNow = infile.readScalar("t");
    if (!options.haveTStart) {
        // If tStart is not in the input file, use the time from the restart file.
        tStart = tNow;
    }

    dvector V = infile.readVector("V");
    rVzero = V[0];

    if (!options.fileNumberOverride) {
        options.outputFileNumber = (int) infile.readScalar("fileNumber");
    }

    infile.close();

    if (options.overrideTu) {
        T[grid.ju] = options.Tu;
    }
    Tu = T[grid.ju];

    CanteraGas gas;
    gas.setOptions(options);
    gas.initialize();
    size_t nSpec = gas.nSpec;

    // save the burned gas properties for the case where burned values are not fixed
    double TbSave = T[grid.jb];
    dvector YbSave(nSpec);
    for (size_t k=0; k<nSpec; k++) {
        YbSave[k] = Y(k,grid.jb);
    }

    if (options.overrideReactants) {
        dvector reactants = calculateReactantMixture();
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

    updateLeftBC();
    double controlSignal;
    if (grid.leftBC == BoundaryCondition::ControlVolume && options.xFlameControl) {
        if (alpha == 0) {
            controlSignal = rVcenter/rhoLeft;
        } else {
            double tmp = pow(x[0],2) + 2*rVcenter/rhoLeft;
            controlSignal = mathUtils::sign(tmp)*sqrt(abs(tmp));
        }
        flamePosIntegralError = controlSignal/(options.xFlameProportionalGain*options.xFlameIntegralGain);
    }
}

dvector FlameSolver::calculateReactantMixture(void)
{
    // Calculate the composition of the reactant mixture from compositions of
    // the fuel and oxidizer mixtures and the equivalence ratio.

    int mC = gas.thermo.elementIndex("C");
    int mO = gas.thermo.elementIndex("O");
    int mH = gas.thermo.elementIndex("H");

    double Cf(0), Hf(0), Of(0); // moles of C/H/O in fuel
    double Co(0), Ho(0), Oo(0); // moles of C/H/O in oxidizer

    dvector Xf(nSpec), Xo(nSpec), Xr(nSpec);

    gas.thermo.setState_TPX(options.Tu, options.pressure, options.fuel);
    gas.getMoleFractions(Xf);
    gas.thermo.setState_TPX(options.Tu, options.pressure, options.oxidizer);
    gas.getMoleFractions(Xo);

    dvector a(gas.thermo.nElements());
    for (size_t k=0; k<nSpec; k++) {
        gas.thermo.getAtoms(k,&a[0]);
        Cf += a[mC]*Xf[k];
        Co += a[mC]*Xo[k];
        Hf += a[mH]*Xf[k];
        Ho += a[mH]*Xo[k];
        Of += a[mO]*Xf[k];
        Oo += a[mO]*Xo[k];
    }
    double stoichAirFuelRatio = -(Of-2*Cf-Hf/2)/(Oo-2*Co-Ho/2);
    Xr = Xf*options.equivalenceRatio + stoichAirFuelRatio*Xo;
    Xr /= mathUtils::sum(Xr);

    return Xr;
}

void FlameSolver::printPerformanceStats(void)
{
    cout << endl << "   *** Performance Stats ***       time   ( call count )" << endl;
    printPerfString("                General Setup: ", setupTimer);
    printPerfString("             Split Term Setup: ", splitTimer);
    printPerfString("    Reaction Term Integration: ", reactionTimer);
    printPerfString("   Diffusion Term Integration: ", diffusionTimer);
    printPerfString("  Convection Term Integration: ", convectionTimer);
    cout << endl << " Subcomponents:" << endl;
    printPerfString("               Reaction Rates: ", reactionRatesTimer);
    printPerfString("         Transport Properties: ", transportTimer);
    printPerfString("          - thermal cond.    : ", conductivityTimer);
    printPerfString("          - viscosity        : ", viscosityTimer);
    printPerfString("          - diffusion coeff. : ", diffusivityTimer);
    printPerfString("     Thermodynamic Properties: ", thermoTimer);
    printPerfString("   Source Jacobian Evaluation: ", jacobianTimer);
    printPerfString("   UTW Convection Integration: ", convectionSystem.utwTimer);
    printPerfString("    Yk Convection Integration: ", convectionSystem.speciesTimer);
    cout << endl;
}

void FlameSolver::printPerfString(const std::string& label, const perfTimer& T) const
{
    cout << format("%s %9.3f (%12i)") % label % T.getTime() % T.getCallCount() << endl;
}

void FlameSolver::updateTransportDomain()
{
    if (!options.transportEliminationDiffusion && !options.transportEliminationConvection) {
        return;
    }

    Array2D dYdtTransport(nSpec, nPoints);
    Array2D dYdtProduction(nSpec, nPoints);

    // Evaluate the full diffusion term for each species
    diffusionTestSolver.initialize(0, options.diffusionTimestep);
    for (size_t k=0; k<nSpec; k++) {
        for (size_t j=0; j<nPoints; j++) {
            diffusionTestSolver.y[j] = Y(k,j);
            diffusionTestTerm.B[j] = 1/rho[j];
            diffusionTestTerm.D[j] = rhoD(k,j);
        }
        const dvector& dYkdt_diff = diffusionTestSolver.get_ydot();
        dYdt.setRow(k, const_cast<double*>(&dYkdt_diff[0]));
        dYdtTransport.setRow(k, const_cast<double*>(&dYkdt_diff[0]));
    }

    // Evaluate the reaction term (using the current reduced mechanism, if any)
    // for reach component
    for (size_t j=0; j<nPoints; j++) {
        if (useCVODE[j]) {
            sdVector& ySource = sourceSolvers[j].y;
            ySource[kMomentum] = U[j];
            ySource[kEnergy] = T[j];
            for (size_t k=0; k<nSpec; k++) {
                ySource[kSpecies+k] = Y(k,j);
            }
            sourceTerms[j].wDot.assign(nSpec, 0);
            sourceTerms[j].strainFunction.pin(tNow);

            sdVector ydotSource(nVars);
            sourceTerms[j].f(tNow, ySource, ydotSource);
            for (size_t k=0; k<nSpec; k++) {
                dYdt(k,j) += ydotSource[kSpecies+k];
                dYdtProduction(k,j) = ydotSource[kSpecies+k];
            }
        } else {
            dvector& ySource = sourceTermsQSS[j].y;
            ySource[kMomentum] = U[j];
            ySource[kEnergy] = T[j];
            for (size_t k=0; k<nSpec; k++) {
                ySource[kSpecies+k] = Y(k,j);
            }
            sourceTermsQSS[j].wDotD.assign(nSpec, 0);
            sourceTermsQSS[j].wDotQ.assign(nSpec, 0);
            sourceTermsQSS[j].strainFunction.pin(tNow);

            dvector q(nVars), d(nVars);
            sourceTermsQSS[j].odefun(tNow, ySource, q, d);
            for (size_t k=0; k<nSpec; k++) {
                dYdt(k,j) += (q[kSpecies+k] = d[kSpecies+k]);
                dYdtProduction(k,j) = (q[kSpecies+k] - d[kSpecies+k]);
            }

        }
    }

    // Evaluate the full convection term
    // Because the convection solver includes the continuity equation containing
    // drho/dt, we need to include the time derivatives from the other terms when
    // evaluating the derivatives of this term, then subtract them from the output
    convectionStartIndices.assign(nSpec, 0);
    convectionStopIndices.assign(nSpec, jj);
    nPointsConvection.assign(nSpec, nPoints);
    if (options.transportEliminationConvection) {
        convectionSystem.resize(nPoints, nPointsConvection, nSpec);
    }

    convectionSystem.setState(U, T, Y);
    convectionSystem.setSplitDerivatives(dTdt, dYdt);
    convectionSystem.evaluate();

    for (size_t j=0; j<nPoints; j++) {
        for (size_t k=0; k<nSpec; k++) {
            dYdt(k,j) += convectionSystem.dYdt(k,j);
            dYdtTransport(k,j) += convectionSystem.dYdt(k,j);
        }
    }

    for (size_t k=0; k<nSpec; k++) {
        // Find the left boundary for species k
        int jStart = jj;
        for (size_t j=0; j<nPoints; j++) {
            if (abs(dYdtTransport(k,j)) > options.adapchem_atol ||
                abs(dYdtProduction(k,j)) > options.adapchem_atol) {
                jStart = j;
                break;
            }
        }

        jStart = max(0, jStart-2);
        if (options.transportEliminationDiffusion) {
            diffusionStartIndices[k] = jStart;
        }
        if (options.transportEliminationConvection) {
            convectionStartIndices[k] = jStart;
        }

        // Find the right boundary for species k
        size_t jStop = jStart;
        for (int j=jj; j>jStart; j--) {
            if (abs(dYdtTransport(k,j)) > options.adapchem_atol ||
                abs(dYdtProduction(k,j)) > options.adapchem_atol) {
                jStop = j;
                break;
            }
        }

        jStop = min(jj, jStop+2);
        if (options.transportEliminationDiffusion) {
            diffusionStopIndices[k] = jStop;
            nPointsDiffusion[k] = jStop - jStart + 1;
        }
        if (options.transportEliminationConvection) {
            convectionStopIndices[k] = jStop;
            nPointsConvection[k] = jStop - jStart + 1;
        }
    }

    for (size_t k=0; k<nSpec; k++) {
        // size the Diffusion systems appropriately
        diffusionSolvers[kSpecies+k].resize(nPointsDiffusion[k], 1, 1);
    }

    if (options.transportEliminationConvection) {
        convectionSystem.resize(nPoints, nPointsConvection, nSpec);
    }

    // Eliminated species don't count for regridding
    grid.leftComponents.resize(nVars, true);
    grid.rightComponents.resize(nVars, true);
    for (size_t k=0; k<nSpec; k++) {
        grid.leftComponents[kSpecies+k] = (diffusionStartIndices[k] == 0 && convectionStartIndices[k] == 0);
        grid.rightComponents[kSpecies+k] = (diffusionStopIndices[k] == jj && convectionStopIndices[k] == jj);
    }
}
