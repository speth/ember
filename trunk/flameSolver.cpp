#include "flameSolver.h"
#include "debugUtils.h"
#include "perfTimer.h"
#include "matlabFile.h"

void flameSolver::setOptions(const configOptions& theOptions)
{
    options = theOptions;
}

void flameSolver::initialize(void)
{
    theSys.options = options;
    theSys.copyOptions();

    // Cantera initialization
    theSys.gas.initialize(options.usingMultiTransport);

    // Initial Conditions
    theSys.setup();
    if (options.haveRestartFile) {
        theSys.loadInitialProfiles();
    } else {
        theSys.generateInitialProfiles();
    }
}

void flameSolver::run(void)
{
    clock_t tIDA1, tIDA2;
    perfTimer runTime;
    runTime.start();

    double integratorTimestep = 0;
    double t = theSys.tStart;
    double dt;

    int nRegrid = 0; // number of time steps since regridding/adaptation
    int nOutput = 0; // number of time steps since storing integral flame parameters
    int nProfile = 0; // number of time steps since saving flame profiles
    int nIntegrate = 0; // number of time steps since restarting integrator
    int nTerminate = 1; // number of steps since last checking termination condition
    int nCurrentState = 0; // number of time steps since profNow.mat and outNow.mat were written

    double tOutput = t; // time of last integral flame parameters output
    double tRegrid = t; // time of last regridding
    double tProfile = t; // time of last profile output

    theSys.grid.updateValues();

    theSys.gas.setStateMass(theSys.Y,theSys.T);
    theSys.gas.getReactionRates(theSys.wDot);
    theSys.updateThermoProperties();

    for (int j=0; j<=theSys.nPoints-1; j++) {
        theSys.qDot[j] = 0;
        for (int k=0; k<theSys.nSpec; k++) {
            theSys.qDot[j] -= theSys.wDot(k,j)*theSys.hk(k,j);
        }
    }

    theSys.tFlamePrev = t;
    theSys.tPrev = t;
    theSys.aPrev = theSys.strainRate(t);

    if (options.outputProfiles) {
        theSys.writeStateMatFile();
    }

    while (t < theSys.tEnd) {

        theSys.setup();

        // **************************************
        // *** Set up the Sundials IDA solver ***
        // **************************************

        tIDA1 = clock();

        sundialsIDA theSolver(theSys.N);
        theSolver.reltol = options.idaRelTol;
        theSolver.nRoots = 0;
        theSolver.findRoots = false;
        theSolver.t0 = theSolver.tInt = t;

        if (options.enforceNonnegativeSpecies) {
            theSolver.imposeConstraints = true;
            theSys.updateConstraints(theSolver.constraints);
        }

        int N = theSys.nVars;

        theSys.rollY(theSolver.y);
        for (int j=0; j<theSys.nPoints; j++) {
            theSolver.abstol(N*j) = options.idaContinuityAbsTol;
            theSolver.abstol(N*j+1) = options.idaMomentumAbsTol;
            theSolver.abstol(N*j+2) = options.idaEnergyAbsTol;
            for (int k=0; k<theSys.nSpec; k++) {
                theSolver.abstol(N*j+k+3) = options.idaSpeciesAbsTol;
            }
        }

        theSys.reltol = theSolver.reltol;
        theSys.abstol = &theSolver.abstol;

        for (int j=0; j<theSys.N; j++) {
            theSolver.ydot(j) = 0;
        }

        theSys.updateLeftBC();
        theSolver.setDAE(&theSys);
        theSolver.calcIC = false;

        // ******************************************
        // *** Get a consistent initial condition ***
        // ******************************************

        theSys.updateAlgebraicComponents();
        if (options.xFlameControl) {
            theSys.gas.setStateMass(theSys.Y,theSys.T);
            theSys.gas.getReactionRates(theSys.wDot);
            theSys.updateThermoProperties();

            for (int j=0; j<=theSys.nPoints-1; j++) {
                theSys.qDot[j] = 0;
                for (int k=0; k<theSys.nSpec; k++) {
                    theSys.qDot[j] -= theSys.wDot(k,j)*theSys.hk(k,j);
                }
            }
            theSys.update_xStag(t, true);
        }
        int ICflag = -1;
        int ICcount = 0;

        while (ICflag!=0 && ICcount < 5) {
            ICcount++;

            // This corrects the drift of the total mass fractions
            theSys.unrollY(theSolver.y);
            theSys.gas.setStateMass(theSys.Y,theSys.T);
            theSys.gas.getMassFractions(theSys.Y);
            theSys.rollY(theSolver.y);
            ICflag = theSys.getInitialCondition(t, theSolver.y, theSolver.ydot, theSys.algebraic);
        }
        if (ICflag != 0) {
            theSys.debugFailedTimestep(theSolver.y);
        }
        if (ICflag == 100) {
            theSys.writeStateMatFile("",true);
            throw;
        }

        // *** Final preparation of IDA solver
        theSolver.initialize();
        theSolver.disableErrorOutput();
        theSolver.setMaxStepSize(options.maxTimestep);

        if (integratorTimestep == 0.0) {
            theSolver.setInitialStepSize(1e-9);
        }

        // ************************************************************
        // *** Integrate until the termination condition is reached ***
        // ************************************************************

        int IDAflag;

        while (t < theSys.tEnd) {
            // *** Take a time step
            try {
                IDAflag = theSolver.integrateOneStep();
            } catch (Cantera::CanteraError) {
                theSys.writeStateMatFile("errorOutput",true);
            }

            dt = integratorTimestep = theSolver.getStepSize();
            t = theSys.tPrev = theSolver.tInt;
            theSys.aPrev = theSys.strainRate(t);

            // *** See if it worked
            if (IDAflag == CV_SUCCESS) {
                nOutput++;
                nRegrid++;
                nProfile++;
                nIntegrate++;
                nTerminate++;
                nCurrentState++;

                if (debugParameters::debugTimesteps) {
                    int order = theSolver.getLastOrder();
                    cout << "t = " << t << "  (dt = " << dt << ") [" << order << "]" << endl;
                }
                if (options.xFlameControl) {
                    theSys.update_xStag(t, true);
                }

            } else {
                cout << "IDA Solver failed at time t = " << t << "  (dt = " << dt << ")" << endl;
                theSys.debugFailedTimestep(theSolver.y);
                theSys.writeStateMatFile("errorOutput",true);
                integratorTimestep = 0;
                break;
            }

            // *** Save the time-series data (out.mat)
            if (t > tOutput || nOutput >= options.outputStepInterval) {
                timeVector.push_back(t);
                timestepVector.push_back(dt);
                heatReleaseRate.push_back(theSys.getHeatReleaseRate());
                consumptionSpeed.push_back(theSys.getConsumptionSpeed());
                flamePosition.push_back(theSys.getFlamePosition());

                tOutput = t + options.outputTimeInterval;
                nOutput = 0;
            }

            // *** Periodic check for terminating the integration
            //     (based on steady heat release rate, etc.)
            if (nTerminate >= options.terminateStepInterval) {
                nTerminate = 0;
                if (checkTerminationCondition()) {
                    tIDA2 = clock();
                    theSolver.printStats(tIDA2-tIDA1);
                    if (options.outputProfiles) {
                        theSys.writeStateMatFile();
                    }
                    runTime.stop();
                    cout << "Runtime: " << runTime.getTime() << " seconds." << endl;
                    return;
                }
            }

            // *** Save the current integral and profile data
            //     in files that are automatically overwritten.
            if (nCurrentState >= options.currentStateStepInterval) {
                nCurrentState = 0;
                matlabFile outFile(options.outputDir+"/outNow.mat");
                outFile.writeVector("t",timeVector);
                outFile.writeVector("dt",timestepVector);
                outFile.writeVector("Q",heatReleaseRate);
                outFile.writeVector("Sc",consumptionSpeed);
                outFile.writeVector("xFlame",flamePosition);
                outFile.close();
                theSys.writeStateMatFile("profNow");
            }

            // *** Save flame profiles
            if (t > tProfile || nProfile >= options.profileStepInterval) {
                if (options.outputProfiles) {
                    sdVector resTemp(theSys.N);
                    theSys.f(t, theSolver.y, theSolver.ydot, resTemp);
                    theSys.writeStateMatFile();
                }

                tProfile = t + options.profileTimeInterval;
                nProfile = 0;
            }

            // *** Adapt the grid if necessary
            if (t > tRegrid || nRegrid >= options.regridStepInterval) {
                tRegrid = t + options.regridTimeInterval;
                nRegrid = 0;

                // dampVal sets a limit on the maximum grid size
                for (int j=0; j<theSys.nPoints; j++) {
                    double num = min(theSys.mu[j],theSys.lambda[j]/theSys.cp[j]);
                    for (int k=0; k<theSys.nSpec; k++) {
                        num = min(num,theSys.rhoD(k,j));
                    }
                    theSys.grid.dampVal[j] = sqrt(num/(theSys.rho[j]*theSys.strainRate(t)));
                }

                 vector<dvector> currentSolution, currentSolutionDot;
                theSys.rollVectorVector(theSolver.y, theSys.qDot, currentSolution);
                theSys.rollVectorVector(theSolver.ydot, theSys.qDot*0, currentSolutionDot);

                bool regridFlag = theSys.grid.regrid(currentSolution, currentSolutionDot);
                bool adaptFlag = theSys.grid.adapt(currentSolution, currentSolutionDot);

                // Perform updates that are necessary if the grid has changed
                if (adaptFlag || regridFlag) {
                    nIntegrate = 0;
                    theSys.nPoints = theSys.grid.jj+1;
                    cout << "Grid size: " << theSys.nPoints << " points." << endl;
                    theSys.setup();

                    theSys.unrollVectorVector(currentSolution);
                    theSys.unrollVectorVectorDot(currentSolutionDot);

                    // Correct the drift of the total mass fractions
                    theSys.gas.setStateMass(theSys.Y,theSys.T);
                    theSys.gas.getMassFractions(theSys.Y);

                    // exit the inner loop to reinitialize the integrator for the new problem size
                    break;
                }
            }

            if (nIntegrate > options.integratorRestartInterval) {
              nIntegrate = 0;

              // exit inner loop to reinitialize the integrator
              break;
            }
        }

        // *** This is the end for the current instance of the IDA solver
        tIDA2 = clock();
        theSolver.printStats(tIDA2-tIDA1);
        if (debugParameters::debugPerformanceStats) {
            theSys.printPerformanceStats();
        }

    }

    // *** Integration has reached the termination condition
    if (options.outputProfiles) {
        theSys.writeStateMatFile();
    }
    runTime.stop();
    cout << "Runtime: " << runTime.getTime() << " seconds." << endl;
}

void flameSolver::calculateReactantMixture(void)
{
    // Calculate the composition of the reactant mixture from compositions of
    // the fuel and oxidizer mixtures and the equivalence ratio.

    Cantera_CXX::IdealGasMix fuel(options.gasMechanismFile,options.gasPhaseID);
    Cantera_CXX::IdealGasMix oxidizer(options.gasMechanismFile,options.gasPhaseID);

    fuel.setState_TPX(options.Tu, options.pressure, options.fuel);
    oxidizer.setState_TPX(options.Tu, options.pressure, options.oxidizer);

    double Cf(0), Hf(0), Of(0); // moles of C/H/O in fuel
    double Co(0), Ho(0), Oo(0); // moles of C/H/O in oxidizer

    int nSpec = fuel.nSpecies();
    int mC = fuel.elementIndex("C");
    int mO = fuel.elementIndex("O");
    int mH = fuel.elementIndex("H");

    dvector Xf(nSpec), Xo(nSpec), Xr(nSpec);
    fuel.getMoleFractions(&Xf[0]);
    oxidizer.getMoleFractions(&Xo[0]);
    dvector a(fuel.nElements());
    for (int k=0; k<nSpec; k++) {
        fuel.getAtoms(k,&a[0]);
        Cf += a[mC]*Xf[k];
        Co += a[mC]*Xo[k];
        Hf += a[mH]*Xf[k];
        Ho += a[mH]*Xo[k];
        Of += a[mO]*Xf[k];
        Oo += a[mO]*Xo[k];
    }
    double stoichAirFuelRatio = -(Of-2*Cf-Hf/2)/(Oo-2*Co-Ho/2);
    options.reactants = Xf*options.equivalenceRatio + stoichAirFuelRatio*Xo;
    options.reactants /= mathUtils::sum(options.reactants);
}

bool flameSolver::checkTerminationCondition(void)
{

    if (options.terminateForSteadyQdot) {
        int j1 = mathUtils::findLast(timeVector < (theSys.tNow - options.terminationPeriod));

        if (j1 == -1)
        {
            cout << "Continuing integration: t (" << theSys.tNow-timeVector[0] <<
                ") < terminationPeriod (" << options.terminationPeriod << ")" << endl;
            return false;
        }

        int j2 = timeVector.size()-1;
        double qMean = mathUtils::mean(heatReleaseRate,j1,j2);
        double hrrError = 0;
        for (int j=j1; j<=j2; j++) {
            hrrError += abs(heatReleaseRate[j]-qMean);
        }
        hrrError /= (j2-j1+1);

        cout << "Heat release rate deviation =  " << hrrError/qMean*100 << "%" << endl;
        cout << "hrrError = " << hrrError << endl;

        if (hrrError/abs(qMean) < options.terminationTolerance) {
            cout << "Terminating integration: ";
            cout << "Heat release deviation less than relative tolerance." << endl;
            return true;
        } else if (hrrError < options.terminationAbsTol) {
            cout << "Terminating integration: ";
            cout << "Heat release rate deviation less than absolute tolerance." << endl;
            return true;
        } else if (theSys.tNow-theSys.tStart > options.terminationMaxTime ) {
          cout << "Terminating integration: Maximum integration time reached." << endl;
          return true;
        } else {
            cout << "Continuing integration. t = "<< theSys.tNow-timeVector[0] << endl;
        }

    }
    return false;
}
