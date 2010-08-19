#include "flameSolver.h"
#include "debugUtils.h"
#include "perfTimer.h"
#include "dataFile.h"
#include "mathUtils.h"

#include <boost/format.hpp>
#include <boost/filesystem.hpp>

using boost::format;

FlameSolver::FlameSolver()
    : convectionSolver(NULL)
{
}

FlameSolver::~FlameSolver()
{
    delete convectionSolver;
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

    // Initial Conditions
    if (options.haveRestartFile) {
        loadProfile();
    } else {
        generateProfile();
    }

    grid.setSize(x.size());
    convectionTerm.Yleft = Yleft;
    convectionTerm.Tleft = Tleft;
    convectionTerm.gas = &gas;

    for (size_t k=0; k<nVars; k++) {
        DiffusionSystem* term = new DiffusionSystem();
        BDFIntegrator* integrator = new BDFIntegrator(*term);
        diffusionTerms.push_back(term);
        diffusionSolvers.push_back(integrator);
    }

    resizeAuxiliary();
}

void FlameSolver::run(void)
{
    clock_t tIDA1, tIDA2;
    perfTimer runTime;
    runTime.start();

    double t = tStart;
    double dt;

    int nRegrid = 0; // number of time steps since regridding/adaptation
    int nOutput = 0; // number of time steps since storing integral flame parameters
    int nProfile = 0; // number of time steps since saving flame profiles
    int nIntegrate = 0; // number of time steps since restarting integrator
    int nTerminate = 1; // number of steps since last checking termination condition
    int nCurrentState = 0; // number of time steps since profNow.h5 and outNow.h5 were written

    double tOutput = t; // time of next integral flame parameters output (this step)
    double tRegrid = t + options.regridTimeInterval; // time of next regridding
    double tProfile = t + options.profileTimeInterval; // time of next profile output

    grid.updateValues();

    tFlamePrev = t;
    tNow = t;
    tPrev = t;
    aPrev = strainfunc.a(t);

    resizeAuxiliary();

    if (options.outputProfiles) {
        writeStateFile();
    }

    while (t < tEnd) {

        // Debug sanity check
        bool error = false;
        for (size_t j=0; j<nPoints; j++) {
            if (T[j] < 295 || T[j] > 3000) {
                cout << "WARNING: Unexpected Temperature: T = " << T[j] << "at j = " << j << endl;
                error = true;
            }
        }
        if (error) {
            writeStateFile();
        }

        tIDA1 = clock();

        // Calculate auxiliary data
        for (size_t j=0; j<nPoints; j++) {
            gas.setStateMass(&Y(0,j), T[j]);
            rho[j] = gas.getDensity();
            Wmx[j] = gas.getMixtureMolecularWeight();
            lambda[j] = gas.getThermalConductivity();
            cp[j] = gas.getSpecificHeatCapacity();
            mu[j] = gas.getViscosity();
            gas.getWeightedDiffusionCoefficients(&rhoD(0,j));
            gas.getThermalDiffusionCoefficients(&Dkt(0,j));
            gas.getSpecificHeatCapacities(&cpSpec(0,j));
            gas.getReactionRates(&wDot(0,j));
            gas.getEnthalpies(&hk(0,j));
            qDot[j] = 0;
            for (size_t k=0; k<nSpec; k++) {
                qDot[j] -= wDot(k,j)*hk(k,j);
            }
        }

        updateLeftBC();
        updateDiffusionFluxes();
        if (options.xFlameControl) {
            update_xStag(t, true); // calculate the value of rVzero
        }
        convectionTerm.rVzero = rVzero;
        dt = options.globalTimestep;

        // *** Set the initial conditions and auxiliary variables for each
        // solver/system, and set the value of the "constant" term(s) to zero

        // Convection solver
        sdVector& yConv = convectionSolver->y;
        for (size_t j=0; j<nPoints; j++) {
            yConv[nVars*j+kMomentum] = U[j];
            yConv[nVars*j+kEnergy] = T[j];
            for (size_t k=0; k<nSpec; k++) {
                yConv[nVars*j+kSpecies+k] = Y(k,j);
            }
        }
        convectionSolver->t0 = t;
        convectionSolver->maxNumSteps = 1000000;
        convectionSolver->minStep = 1e-16;
        convectionSolver->initialize();

        convectionTerm.splitConstU.assign(nPoints, 0);
        convectionTerm.splitLinearU.assign(nPoints, 0);
        convectionTerm.splitConstT.assign(nPoints, 0);
        convectionTerm.splitLinearT.assign(nPoints, 0);
        convectionTerm.splitConstY.data().assign(nPoints*nSpec, 0);
        convectionTerm.splitLinearY.data().assign(nPoints*nSpec, 0);

        // Source solvers
        for (size_t j=0; j<nPoints; j++) {
            sdVector& ySource = sourceSolvers[j].y;
            ySource[kMomentum] = U[j];
            ySource[kEnergy] = T[j];
            for (size_t k=0; k<nSpec; k++) {
                ySource[kSpecies+k] = Y(k,j);
            }
            sourceSolvers[j].t0 = t;
            sourceSolvers[j].initialize();
            sourceTerms[j].splitConst.assign(nVars, 0);
            sourceTerms[j].splitLinear.assign(nVars, 0);
            sourceTerms[j].strainFunction.pin(t);
        }

        // Diffusion solvers
        dvector& yDiff_U = diffusionSolvers[kMomentum].y;
        dvector& yDiff_T = diffusionSolvers[kEnergy].y;
        for (size_t j=0; j<nPoints; j++) {
            yDiff_U[j] = U[j];
            yDiff_T[j] = T[j];

            diffusionTerms[kMomentum].B[j] = 1/rho[j];
            diffusionTerms[kEnergy].B[j] = 1/(rho[j]*cp[j]);

            diffusionTerms[kMomentum].D[j] = mu[j];
            diffusionTerms[kEnergy].D[j] = lambda[j];
        }

        for (size_t k=0; k<nSpec; k++) {
            dvector& yDiff_Y = diffusionSolvers[kSpecies+k].y;
            DiffusionSystem& sys = diffusionTerms[kSpecies+k];
            for (size_t j=0; j<nPoints; j++) {
                yDiff_Y[j] = Y(k,j);
                sys.B[j] = 1/rho[j];
                sys.D[j] = rhoD(k,j);
            }
        }

        // TODO: Use timestep that is based on each component's diffusivity
        for (size_t k=0; k<nVars; k++) {
            diffusionSolvers[k].set_dt(options.diffusionTimestep);
        }

        for (size_t k=0; k<nVars; k++) {
            diffusionTerms[k].splitConst.assign(nPoints, 0);
            diffusionTerms[k].splitLinear.assign(nPoints, 0);
            diffusionSolvers[k].t = t;
            diffusionSolvers[k].initialize();
        }

        // *** Get the current value of each term's time derivative and
        // diagonalized Jacobian needed to compute the diagonalized
        // approximation used in the other terms, as well as the
        // "predicted" values Uextrap, Textrap, and Yextrap

        // Convection solver
        sdVector ydotConv(nVars*nPoints);
        convectionTerm.f(tNow, yConv, ydotConv);
        // TODO: Can this loop just be replaced by copying dUdt etc. from the system?
        for (size_t j=0; j<nPoints; j++) {
            constUconv[j] = ydotConv[nVars*j+kMomentum];
            constTconv[j] = ydotConv[nVars*j+kEnergy];
            for (size_t k=0; k<nSpec; k++) {
                constYconv(k,j) = ydotConv[nVars*j+kSpecies+k];
            }
        }
        convectionTerm.get_diagonal(tNow, linearUconv, linearTconv, linearYconv);

        // Source solvers
        sdVector ydotSource(nVars);
        sdMatrix Jtmp(nVars, nVars);
        for (size_t j=0; j<nPoints; j++) {
            sourceTerms[j].updateDiagonalJac = true;
            sourceTerms[j].f(tNow, sourceSolvers[j].y, ydotSource);
            sourceTerms[j].denseJacobian(tNow, sourceSolvers[j].y, ydotSource, Jtmp);
            constUprod[j] = ydotSource[kMomentum];
            constTprod[j] = ydotSource[kEnergy];
            for (size_t k=0; k<nSpec; k++) {
                constYprod(k,j) = ydotSource[kSpecies+k];
                linearYprod(k,j) = sourceTerms[j].diagonalJac[kSpecies+k];
            }
        }

        // Diffusion solvers
        constUdiff = diffusionSolvers[kMomentum].get_ydot();
        constTdiff = diffusionSolvers[kEnergy].get_ydot();
        linearUdiff = diffusionSolvers[kMomentum].get_diagonal();
        linearTdiff = diffusionSolvers[kEnergy].get_diagonal();

        dvector Ydiff_tmp(nVars);
        for (size_t k=0; k<nSpec; k++) {
            const dvector& constantTerm = diffusionSolvers[kSpecies+k].get_ydot();
            const dvector& linearTerm = diffusionSolvers[kSpecies+k].get_diagonal();
            // const_cast required because Array2D::setRow is missing a const qualifier
            constYdiff.setRow(k, const_cast<double*>(&constantTerm[0]));
            linearYdiff.setRow(k, const_cast<double*>(&linearTerm[0]));
        }

        // *** Use the time derivatives to calculate the values for the
        //     state variables based on the diagonalized approximation

        // Offset the value of the linear terms at the current solution
        constUprod -= linearUprod*U;
        constUconv -= linearUconv*U;
        constUdiff -= linearUdiff*U;

        constTprod -= linearTprod*T;
        constTconv -= linearTconv*T;
        constTdiff -= linearTdiff*T;

        for (size_t j=0; j<nPoints; j++) {
            for (size_t k=0; k<nSpec; k++) {
                constYprod(k,j) -= linearYprod(k,j)*Y(k,j);
                constYconv(k,j) -= linearYconv(k,j)*Y(k,j);
                constYdiff(k,j) -= linearYdiff(k,j)*Y(k,j);
            }
        }

        // This is the exact solution to the linearized problem
        for (size_t j=0; j<nPoints; j++) {
            double K = constUprod[j] + constUdiff[j] + constUconv[j];
            double L = linearUprod[j] + linearUdiff[j] + linearUconv[j];
            double L2; // L2 handles the limit as L->0 for the second term in the solution
            L2 = (abs(L*dt) < 1e-10) ? 1e-10/dt : L;

            Uextrap[j] = U[j]*exp(L*dt) + K/L2*(exp(L2*dt)-1);

            K = constTprod[j] + constTdiff[j] + constTconv[j];
            L = linearTprod[j] + linearTdiff[j] + linearTconv[j];
            L2 = (abs(L*dt) < 1e-10) ? 1e-10/dt : L;
            Textrap[j] = T[j]*exp(L*dt) + K/L2*(exp(L2*dt)-1);

            for (size_t k=0; k<nSpec; k++) {
                K = constYprod(k,j) + constYdiff(k,j) + constYconv(k,j);
                L = linearYconv(k,j) + linearYdiff(k,j) + linearYprod(k,j);
                L2 = (abs(L*dt) < 1e-10) ? 1e-10/dt : L;
                Yextrap(k,j) = Y(k,j)*exp(L*dt) + K/L2*(exp(L2*dt)-1);
            }
        }

        // *** Calculate the constant and linear terms for each system of
        //     equations.

        // Convection term
        convectionTerm.splitConstU = constUprod + constUdiff;
        convectionTerm.splitLinearU = linearUprod + linearUdiff;
        convectionTerm.splitConstT = constTprod + constTdiff;
        convectionTerm.splitLinearT = linearTprod + linearTdiff;
        for (size_t j=0; j<nPoints; j++) {
            for (size_t k=0; k<nSpec; k++) {
                convectionTerm.splitConstY(k,j) = constYprod(k,j) + constYdiff(k,j);
                convectionTerm.splitLinearY(k,j) = linearYprod(k,j) + linearYdiff(k,j);
            }
        }

        // Source terms
        for (size_t j=0; j<nPoints; j++) {
            SourceSystem& term = sourceTerms[j];
            term.splitConst[kMomentum] = constUconv[j] + constUdiff[j];
            term.splitLinear[kMomentum] = linearUconv[j] + linearUdiff[j];
            term.splitConst[kEnergy] = constTconv[j] + constTdiff[j];
            term.splitLinear[kEnergy] = linearTconv[j] + linearTdiff[j];
            for (size_t k=0; k<nSpec; k++) {
                term.splitConst[kSpecies+k] = constYconv(k,j) + constYdiff(k,j);
                term.splitLinear[kSpecies+k] = linearYconv(k,j) + linearYdiff(k,j);
            }
        }

        // Diffusion terms
        diffusionTerms[kMomentum].splitConst = constUconv + constUprod;
        diffusionTerms[kMomentum].splitLinear = linearUconv + linearUprod;
        diffusionTerms[kEnergy].splitConst = constTconv + constTprod;
        diffusionTerms[kEnergy].splitLinear = linearTconv + linearTprod;
        for (size_t k=0; k<nSpec; k++) {
            dvector& splitConst = diffusionTerms[kSpecies+k].splitConst;
            dvector& splitLinear = diffusionTerms[kSpecies+k].splitLinear;
            for (size_t j=0; j<nPoints; j++) {
                splitConst[j] = constYconv(k,j) + constYprod(k,j);
                splitLinear[j] = linearYconv(k,j) + linearYprod(k,j);
            }
        }

        // *** Take one global timestep
        double tNext = tNow + dt;

        int cvode_flag = 0;
        try {
            for (size_t j=0; j<nPoints; j++) {
                cvode_flag |= sourceSolvers[j].integrateToTime(tNext);
            }

            for (size_t k=0; k<nVars; k++) {
                diffusionSolvers[k].integrateToTime(tNext);
            }
            cvode_flag |= convectionSolver->integrateToTime(tNext);

        } catch (Cantera::CanteraError) {
            writeStateFile("errorOutput", true);
            cout << "Integration failed at t = " << tNow << std::endl;
        }

        t = tNext;
        tNow = tNext;
        aPrev = strainfunc.a(t);

        cvode_flag = CV_SUCCESS; // CHARGE!
        if (cvode_flag == CV_SUCCESS) {
            nOutput++;
            nRegrid++;
            nProfile++;
            nIntegrate++;
            nTerminate++;
            nCurrentState++;

            if (debugParameters::debugTimesteps) {
                cout << "t = " << format("%8.6f") % t;
                cout << "  (dt = " << format("%9.3e") % dt;
                cout << ")" << endl;
            }
        } else {
            cout << "CVODE Solver failed at time t = " << format("%8.6f") % t;
            cout << "  (dt = " << format("%9.3e") % dt << ")" << endl;
            writeStateFile("errorOutput",true);
            break;
        }

        // *** Combine the solutions from the split integrators and form
        //     the new state vector
        convectionTerm.unroll_y(convectionSolver->y);
        for (size_t j=0; j<nPoints; j++) {
            sourceTerms[j].unroll_y(sourceSolvers[j].y);
        }

        for (size_t j=0; j<nPoints; j++) {
            U[j] = convectionTerm.U[j] + sourceTerms[j].U + diffusionSolvers[kMomentum].y[j] - 2*Uextrap[j];
            T[j] = convectionTerm.T[j] + sourceTerms[j].T + diffusionSolvers[kEnergy].y[j] - 2*Textrap[j];
            for (size_t k=0; k<nSpec; k++) {
                Y(k,j) = convectionTerm.Y(k,j) + sourceTerms[j].Y[k] +
                         diffusionSolvers[kSpecies+k].y[j] - 2*Yextrap(k,j);
            }
        }

        if (t > tOutput || nOutput >= options.outputStepInterval) {
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
        if (nTerminate >= options.terminateStepInterval) {
            nTerminate = 0;
            if (checkTerminationCondition()) {
                tIDA2 = clock();
                //theSolver.printStats(tIDA2-tIDA1);
                if (options.outputProfiles) {
                    writeStateFile();
                }
                runTime.stop();
                cout << "Runtime: " << runTime.getTime() << " seconds." << endl;
                return;
            }
        }

        // *** Save the current integral and profile data
        //     in files that are automatically overwritten,
        //     and save the time-series data (out.h5)
        if (nCurrentState >= options.currentStateStepInterval) {
            nCurrentState = 0;
            DataFile outFile(options.outputDir+"/outNow.h5");
            outFile.writeVector("t",timeVector);
            outFile.writeVector("dt",timestepVector);
            outFile.writeVector("Q",heatReleaseRate);
            outFile.writeVector("Sc",consumptionSpeed);
            outFile.writeVector("xFlame",flamePosition);
            outFile.close();
            writeStateFile("profNow");
        }

        // *** Save flame profiles
        if (t > tProfile || nProfile >= options.profileStepInterval) {
            if (options.outputProfiles) {
//                sdVector resTemp(theSys.N);
//                theSys.f(t, theSolver.y, theSolver.ydot, resTemp);
                writeStateFile();
            }
            tProfile = t + options.profileTimeInterval;
            nProfile = 0;
        }

        if (t > tRegrid || nRegrid >= options.regridStepInterval) {
            tRegrid = t + options.regridTimeInterval;
            nRegrid = 0;

            // "rollVectorVector"
            vector<dvector> currentSolution;
            vector<dvector> currentSolutionDot; // TODO: Remove this, as we don't need ydot anymore
            currentSolution.push_back(U);
            currentSolution.push_back(T);
            currentSolutionDot.push_back(dUdt);
            currentSolutionDot.push_back(dTdt);
            for (size_t k=0; k<nSpec; k++) {
                dvector tmp(nPoints);
                dvector dtmp(nPoints);
                for (size_t j=0; j<nPoints; j++) {
                    tmp[j] = Y(k,j);
                    dtmp[j] = dYdt(k,j);
                }
                currentSolution.push_back(tmp);
                currentSolutionDot.push_back(dtmp);
            }

            grid.regrid(currentSolution, currentSolutionDot);

            grid.dampVal.resize(grid.x.size());
            cout << nPoints << endl;
            // dampVal sets a limit on the maximum grid size
            for (size_t j=0; j<nPoints; j++) {
                double num = min(mu[j],lambda[j]/cp[j]);
                for (size_t k=0; k<nSpec; k++) {
                    num = min(num, rhoD(k,j));
                }
                grid.dampVal[j] = sqrt(num/(rho[j]*strainfunc.a(t)));
            }

            grid.adapt(currentSolution, currentSolutionDot);

            // Perform updates that are necessary if the grid has changed
            if (grid.updated) {
                grid.updated = false;
                nIntegrate = 0;
                cout << "Grid size: " << nPoints << " points." << endl;

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
                    // Correct the drift of the total mass fractions
                    gas.setStateMass(&Y(0,j), T[j]);
                    gas.getMassFractions(&Y(0,j));
                }

                // Assign the new grid to the terms that need it
                for (size_t k=0; k<nVars; k++) {
                    diffusionTerms[k].grid = grid;
                }
                convectionTerm.grid = grid;

                // Allocate the solvers and arrays for auxiliary variables
                resizeAuxiliary();

                writeStateFile("postAdapt");
            }

        }

        tIDA2 = clock();
        // *** This is the end for the current instance of the IDA solver
        // theSolver.printStats(tIDA2-tIDA1);

        if (debugParameters::debugPerformanceStats) {
            printPerformanceStats();
        }
    }

    // *** Integration has reached the termination condition
    if (options.outputProfiles) {
        writeStateFile();
    }

    runTime.stop();
    cout << "Runtime: " << runTime.getTime() << " seconds." << endl;

}

bool FlameSolver::checkTerminationCondition(void)
{

    if (options.terminateForSteadyQdot) {
        int j1 = mathUtils::findLast(timeVector < (tNow - options.terminationPeriod));

        if (j1 == -1)
        {
            cout << "Continuing integration: t (" << format("%8.6f") % (tNow-timeVector[0]) <<
                ") < terminationPeriod (" << format("%8.6f") % options.terminationPeriod << ")" << endl;
            return false;
        }

        int j2 = timeVector.size()-1;
        double qMean = mathUtils::mean(heatReleaseRate,j1,j2);
        double hrrError = 0;
        for (int j=j1; j<=j2; j++) {
            hrrError += abs(heatReleaseRate[j]-qMean);
        }
        hrrError /= (j2-j1+1);

        cout << "Heat release rate deviation =  " << format("%6.3f") % (hrrError/qMean*100) << "%    ";
        cout << "hrrError = " << format("%9.4e") % hrrError << endl;

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
            cout << "Continuing integration. t = "<< format("%8.6f") % (tNow-timeVector[0]) << endl;
        }

    }
    return false;
}

void FlameSolver::writeStateFile(const std::string fileNameStr, bool errorFile)
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

    // Write the state data to the output file:
    outFile.writeScalar("t", tNow);
    outFile.writeVector("x", x);
    outFile.writeVector("T", T);
    outFile.writeVector("U", U);
    outFile.writeArray2D("Y", Y);
    outFile.writeScalar("a", strainfunc.a(tNow));
    outFile.writeScalar("dadt", strainfunc.dadt(tNow));
    outFile.writeScalar("fileNumber", options.outputFileNumber);

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
        outFile.writeVector("V", convectionTerm.V);
        outFile.writeArray2D("Yextrap", Yextrap);
        outFile.writeVector("Textrap", Textrap);
        outFile.writeVector("Uextrap", Uextrap);
        outFile.writeArray2D("wdot", wDot);
        outFile.writeArray2D("rhoD", rhoD);
        outFile.writeVector("lambda", lambda);

        outFile.writeVector("cp", cp);
        outFile.writeVector("mu", mu);
        outFile.writeVector("Wmx", Wmx);
        outFile.writeVector("W", W);
//        outFile.writeArray2D("jFick", jFick);
//        outFile.writeArray2D("jSoret", jSoret);
//        outFile.writeVector("qFourier",qFourier);
        outFile.writeVector("cfp", grid.cfp);
        outFile.writeVector("cf", grid.cf);
        outFile.writeVector("cfm", grid.cfm);
        outFile.writeVector("hh", hh);
        outFile.writeVector("rphalf", grid.rphalf);
        outFile.writeScalar("Tleft", Tleft);
        outFile.writeVector("Yleft", Yleft);
        outFile.writeVector("sumcpj", sumcpj);

        dvector Uprod(nPoints);
        dvector Tprod(nPoints);
        Array2D Ydiff(nSpec, nPoints), Yprod(nSpec, nPoints);

        for (size_t j=0; j<nPoints; j++) {
            Uprod[j] = sourceTerms[j].U;
            Tprod[j] = sourceTerms[j].T;
            for (size_t k=0; k<nSpec; k++) {
                Yprod(k,j) = sourceTerms[j].Y[k];
                Ydiff(k,j) = diffusionSolvers[kSpecies+k].y[j];
            }
        }

        outFile.writeVector("Uconv", convectionTerm.U);
        outFile.writeVector("Udiff", diffusionSolvers[kMomentum].y);
        outFile.writeVector("Uprod", Uprod);
        outFile.writeVector("Tconv", convectionTerm.T);
        outFile.writeVector("Tdiff", diffusionSolvers[kEnergy].y);
        outFile.writeVector("Tprod", Tprod);
        outFile.writeArray2D("Yconv", convectionTerm.Y);
        outFile.writeArray2D("Ydiff", Ydiff);
        outFile.writeArray2D("Yprod", Yprod);

        outFile.writeVector("Uextrap", Uextrap);
        outFile.writeVector("Textrap", Textrap);
        outFile.writeArray2D("Yextrap", Yextrap);
    }

    if (options.outputResidualComponents || errorFile) {
        outFile.writeVector("dTdt_diff", constTdiff);
        outFile.writeVector("dTdt_conv", constTconv);
        outFile.writeVector("dTdt_prod", constTprod);
        outFile.writeVector("dUdt_diff", constUdiff);
        outFile.writeVector("dUdt_conv", constUconv);
        outFile.writeVector("dUdt_prod", constUprod);
        outFile.writeArray2D("dYdt_diff", constYdiff);
        outFile.writeArray2D("dYdt_conv", constYconv);
        outFile.writeArray2D("dYdt_prod", constYprod);

        outFile.writeVector("linearTdiff", linearTdiff);
        outFile.writeVector("linearTconv", linearTconv);
        outFile.writeVector("linearTprod", linearTprod);
        outFile.writeVector("linearUdiff", linearUdiff);
        outFile.writeVector("linearUconv", linearUconv);
        outFile.writeVector("linearUprod", linearUprod);
        outFile.writeArray2D("linearYdiff", linearYdiff);
        outFile.writeArray2D("linearYconv", linearYconv);
        outFile.writeArray2D("linearYprod", linearYprod);
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

void FlameSolver::resizeAuxiliary()
{
    perfTimerResize.start();
    size_t nPointsOld = rho.size();
    grid.setSize(T.size());

    if (nPoints == nPointsOld) {
        return; // nothing to do
    }

    nSpec = gas.nSpec;
    nVars = 2+nSpec;
    N = nVars*nPoints;

    Uextrap.resize(nPoints, 0);
    Textrap.resize(nPoints, 0);
    Yextrap.resize(nSpec,nPoints, 0);

    dUdt.resize(nPoints,0);
    dTdt.resize(nPoints,0);
    dYdt.resize(nSpec,nPoints);

    constTdiff.resize(nPoints,0);
    constTconv.resize(nPoints,0);
    constTprod.resize(nPoints,0);
    linearTdiff.resize(nPoints,0);
    linearTconv.resize(nPoints,0);
    linearTprod.resize(nPoints,0);

    constUdiff.resize(nPoints,0);
    constUconv.resize(nPoints,0);
    constUprod.resize(nPoints,0);
    linearUdiff.resize(nPoints,0);
    linearUconv.resize(nPoints,0);
    linearUprod.resize(nPoints,0);

    constYdiff.resize(nSpec,nPoints,0);
    constYconv.resize(nSpec,nPoints,0);
    constYprod.resize(nSpec,nPoints,0);
    linearYdiff.resize(nSpec,nPoints,0);
    linearYconv.resize(nSpec,nPoints,0);
    linearYprod.resize(nSpec,nPoints,0);

    rho.resize(nPoints);
    Wmx.resize(nPoints);
    mu.resize(nPoints);
    lambda.resize(nPoints);
    cp.resize(nPoints);
    sumcpj.resize(nPoints);
    jCorr.resize(nPoints);
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

    if (nPoints > nPointsOld) {
        for (size_t i=nPointsOld; i<nPoints; i++) {
            // Create and initialize the new SourceSystem
            SourceSystem* system = new SourceSystem();
            system->resize(nSpec);
            system->gas = &gas;
            system->strainFunction = strainfunc;
            system->rhou = rhou;
            system->W = W;

            // Create and initialize the new Sundials solver
            sundialsCVODE* solver = new sundialsCVODE(nVars);
            solver->setODE(system);
            solver->abstol[kMomentum] = options.idaMomentumAbsTol;
            solver->abstol[kEnergy] = options.idaEnergyAbsTol;
            for (size_t k=0; k<nSpec; k++) {
                solver->abstol[kSpecies+k] = options.idaSpeciesAbsTol;
            }
            solver->reltol = options.idaRelTol;
            solver->linearMultistepMethod = CV_BDF;
            solver->nonlinearSolverMethod = CV_NEWTON;
            solver->maxNumSteps = 1000000;
            solver->minStep = 1e-16;

            // Store the solver and system
            sourceTerms.push_back(system);
            sourceSolvers.push_back(solver);
        }

    } else {
        // Delete the unwanted solvers and systems
        sourceTerms.erase(sourceTerms.begin()+nPoints, sourceTerms.end());
        sourceSolvers.erase(sourceSolvers.begin()+nPoints, sourceSolvers.end());
    }

    // Resize solution vector for diffusion systems / solvers
    for (size_t k=0; k<nVars; k++) {
        diffusionSolvers[k].set_size(nPoints, 1, 1);
        // TODO: Refactor this into class DiffusionSystem
        diffusionTerms[k].B.resize(nPoints);
        diffusionTerms[k].splitConst.resize(nPoints);
        diffusionTerms[k].splitLinear.resize(nPoints);
        diffusionTerms[k].D.resize(nPoints);
        diffusionTerms[k].grid = grid;
    }

    convectionTerm.resize(nSpec, nPoints);
    convectionTerm.grid = grid;

    delete convectionSolver;
    convectionSolver = new sundialsCVODE(N);
    convectionSolver->setODE(&convectionTerm);
    convectionSolver->reltol = options.idaRelTol;
    convectionSolver->linearMultistepMethod = CV_ADAMS;
    convectionSolver->nonlinearSolverMethod = CV_FUNCTIONAL;
    for (size_t j=0; j<nPoints; j++) {
        convectionSolver->abstol[nVars*j+kMomentum] = options.idaMomentumAbsTol;
        convectionSolver->abstol[nVars*j+kEnergy] = options.idaEnergyAbsTol;
        for (size_t k=0; k<nSpec; k++) {
            convectionSolver->abstol[nVars*j+kSpecies+k] = options.idaSpeciesAbsTol;
        }
    }

    perfTimerResize.stop();
}

void FlameSolver::updateDiffusionFluxes()
{
    for (size_t j=0; j<jj; j++) {
        sumcpj[j] = 0;
        for (size_t k=0; k<nSpec; k++) {
            jFick(k,j) = -0.5*(rhoD(k,j)+rhoD(k,j+1)) * ((Y(k,j+1)-Y(k,j))/hh[j]) -
                0.5*(rhoD(k,j)*Y(k,j)/Wmx[j]+Y(k,j+1)*rhoD(k,j+1)/Wmx[j+1])*(Wmx[j+1]-Wmx[j])/hh[j];
            jSoret(k,j) = -0.5*(Dkt(k,j)/T[j] + Dkt(k,j+1)/T[j+1])
                * (T[j+1]-T[j])/hh[j];
        }

        jCorr[j] = 0;
        for (size_t k=0; k<nSpec; k++) {
            jCorr[j] -= jFick(k,j);
        }
        for (size_t k=0; k<nSpec; k++) {
            jFick(k,j) += 0.5*(Y(k,j)+Y(k,j+1))*jCorr[j]; // correction to ensure that sum of mass fractions equals 1
            sumcpj[j] += 0.5*(cpSpec(k,j)+cpSpec(k,j+1))/W[k]*(jFick(k,j) + jSoret(k,j));
        }
    }
}

void FlameSolver::updateLeftBC()
{
    // Boundary condition for left edge of domain
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

    U = infile.readVector("U");
    T = infile.readVector("T");
    Y = infile.readArray2D("Y");
    tNow = infile.readScalar("t");

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
    cout << endl << " **Performance Stats**:  time  (call count)" << endl;
    printPerfString("         General Setup: ", perfTimerResize);
    printPerfString("  Preconditioner Setup: ", perfTimerPrecondSetup);
    printPerfString("  Factorizing Jacobian: ", perfTimerLU);
    printPerfString("  Preconditioner Solve: ", perfTimerPrecondSolve);
    printPerfString("   Residual Evaluation: ", perfTimerResFunc);
    printPerfString("        Reaction Rates: ", perfTimerRxnRates);
    printPerfString("  Transport Properties: ", perfTimerTransportProps);
    cout << endl;
}

void FlameSolver::printPerfString(const std::string& label, const perfTimer& T) const
{
    cout << label << mathUtils::stringify(T.getTime(),6) << " (" << T.getCallCount() << ")" << endl;
}
