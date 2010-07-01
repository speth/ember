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
}

void FlameSolver::initialize(void)
{
    //theSys.copyOptions(); TODO: copy to self?
    strainfunc.setOptions(options);

    // Cantera initialization
    gas.initialize();
    nSpec = gas.nSpec;
    W.resize(nSpec);
    gas.getMolecularWeights(W);

    // Initial Conditions
    if (options.haveRestartFile) {
        loadProfile();
    } else {
        generateProfile();
    }

    for (size_t i=0; i<nSpec; i++) {
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

    double integratorTimestep = 0;
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
    tPrev = t;
    aPrev = strainfunc.a(t);

    if (options.outputProfiles) {
        writeStateFile();
    }

    while (t < tEnd) {

        // Allocate the solvers and arrays for auxiliary variables
        resizeAuxiliary();

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

        updateDiffusionFluxes();

        // Set the initial conditions for each solver
        tIDA1 = clock();

        sundialsIDA theSolver(theSys.N);
        theSolver.reltol = options.idaRelTol;
        theSolver.nRoots = 0;
        theSolver.findRoots = false;
        theSolver.t0 = theSolver.tInt = t;

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
            ICflag = theSys.getInitialCondition(t, theSolver.y, theSolver.ydot);
        }
        if (ICflag != 0) {
            theSys.debugFailedTimestep(theSolver.y);
        }
        if (ICflag == 100) {
            theSys.writeStateFile("",true);
            throw debugException("Initial condition calculation failed repeatedly.");
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

        int IDAflag(0);

        while (t < theSys.tEnd) {
            // *** Take a time step
            try {
                IDAflag = theSolver.integrateOneStep();
            } catch (Cantera::CanteraError) {
                theSys.writeStateFile("errorOutput",true);
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
                    cout << "t = " << format("%8.6f") % t;
                    cout << "  (dt = " << format("%9.3e") % dt;
                    cout << ")  [" << order << "]" << endl;
                }
                if (options.xFlameControl) {
                    theSys.update_xStag(t, true);
                }

            } else {
                cout << "IDA Solver failed at time t = " << format("%8.6f") % t;
                cout << "  (dt = " << format("%9.3e") % dt << ")" << endl;
                theSys.debugFailedTimestep(theSolver.y);
                theSys.writeStateFile("errorOutput",true);
                integratorTimestep = 0;
                break;
            }

            // *** Save the time-series data (out.h5)
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
                        theSys.writeStateFile();
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
                DataFile outFile(options.outputDir+"/outNow.h5");
                outFile.writeVector("t",timeVector);
                outFile.writeVector("dt",timestepVector);
                outFile.writeVector("Q",heatReleaseRate);
                outFile.writeVector("Sc",consumptionSpeed);
                outFile.writeVector("xFlame",flamePosition);
                outFile.close();
                theSys.writeStateFile("profNow");
            }

            // *** Save flame profiles
            if (t > tProfile || nProfile >= options.profileStepInterval) {
                if (options.outputProfiles) {
                    sdVector resTemp(theSys.N);
                    theSys.f(t, theSolver.y, theSolver.ydot, resTemp);
                    theSys.writeStateFile();
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
    }

    if (options.outputResidualComponents || errorFile) {
        outFile.writeVector("resEnergyDiff", energyDiff);
        outFile.writeVector("resEnergyConv", energyConv);
        outFile.writeVector("resEnergyProd", energyProd);
        outFile.writeVector("resMomentumDiff", momentumDiff);
        outFile.writeVector("resMomentumConv", momentumConv);
        outFile.writeVector("resMomentumProd", momentumProd);
        outFile.writeArray2D("resSpeciesDiff", speciesDiff);
        outFile.writeArray2D("resSpeciesConv", speciesConv);
        outFile.writeArray2D("resSpeciesProd", speciesProd);
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
    nPoints = T.size();

    if (nPoints == nPointsOld) {
        return; // nothing to do
    }

    nSpec = gas.nSpec;
    nVars = 2+nSpec;
    N = nVars*nPoints;

//    U.resize(nPoints);
//    T.resize(nPoints);
//    Y.resize(nSpec,nPoints);

    dUdt.resize(nPoints,0);
    dTdt.resize(nPoints,0);
    dYdt.resize(nSpec,nPoints);

    energyDiff.resize(nPoints,0);
    energyConv.resize(nPoints,0);
    energyProd.resize(nPoints,0);

    momentumDiff.resize(nPoints,0);
    momentumConv.resize(nPoints,0);
    momentumProd.resize(nPoints,0);

    speciesDiff.resize(nSpec,nPoints,0);
    speciesConv.resize(nSpec,nPoints,0);
    speciesProd.resize(nSpec,nPoints,0);

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

            // Create and initialize the new Sundials solver
            sundialsCVODE* solver = new sundialsCVODE(nVars);
            solver->setODE(system);
            solver->abstol[kMomentum] = options.idaMomentumAbsTol;
            solver->abstol[kEnergy] = options.idaEnergyAbsTol;
            for (size_t k=0; k<nSpec; k++) {
                solver->abstol[kSpecies+k] = options.idaSpeciesAbsTol;
            }
            solver->reltol = options.idaRelTol;

            // Store the solver and system
            sourceTerms.push_back(system);
            sourceSolvers.push_back(solver);
        }

    } else {
        // Delete the unwanted solvers and systems
        sourceTerms.erase(sourceTerms.begin()+nPoints, sourceTerms.end());
        sourceSolvers.erase(sourceSolvers.begin()+nPoints, sourceSolvers.end());
    }

    convectionTerm.resize(nSpec, nPoints);
    delete convectionSolver;
    convectionSolver = new sundialsCVODE(N);
    convectionSolver->setODE(&convectionTerm);
    for (size_t j=0; j<nPoints; j++) {
        convectionSolver->abstol[N*j+kMomentum] = options.idaMomentumAbsTol;
        convectionSolver->abstol[N*j+kEnergy] = options.idaEnergyAbsTol;
        for (size_t k=0; k<nSpec; k++) {
            convectionSolver->abstol[N*j+kSpecies+k] = options.idaSpeciesAbsTol;
        }
    }
    convectionSolver->reltol = options.idaRelTol;

    perfTimerResize.stop();
}

void FlameSolver::updateDiffusionFluxes()
{
    for (int j=0; j<jj; j++) {
        sumcpj[j] = 0;
        for (int k=0; k<nSpec; k++) {
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
    CanteraGas gas;
    gas.setOptions(options);
    gas.initialize();
    size_t nSpec = gas.nSpec;

    // Create a uniform initial grid
    x.resize(options.nPoints);
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

    grid.updateBoundaryIndices();
    grid.updateValues();

    int jm = (grid.ju+grid.jb)/2; // midpoint of the profiles.
    int jl = jm - 4;
    int jr = jm + 4;

    Yb.resize(nSpec);
    Yu.resize(nSpec);

    double a = options.strainRateInitial;

    Tu = options.Tu;

    // Reactants
    dvector reactants = calculateReactantMixture();
    gas.thermo.setState_TPX(Tu, gas.pressure, &reactants[0]);
    rhou = gas.thermo.density();
    gas.thermo.getMassFractions(&Yu[0]);
    gas.thermo.getMassFractions(&Y(0,grid.ju));

    // Products
    gas.thermo.setState_TPX(Tu, gas.pressure, &reactants[0]);
    Cantera::equilibrate(gas.thermo,"HP");
    Tb = gas.thermo.temperature();
    rhob = gas.thermo.density();
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

    dvector V(nPoints);
    V[jm] = 0;
    for (size_t j=jm+1; j<nPoints; j++) {
        V[j] = V[j-1] - rho[j]*U[j]*(x[j]-x[j-1]);
    }

    for (size_t j=jm-1; j>=0; j--) {
        V[j] = V[j+1] + rho[j]*U[j]*(x[j+1]-x[j]);
    }
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

    nPoints = x.size();

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
    Xr = Xf*options.equivalenceRatio + stoichAirFuelRatio*Xo;
    Xr /= mathUtils::sum(Xr);

    return Xr;
}
