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

#ifdef NDEBUG
    static const bool VERY_VERBOSE = false;
#else
    static const bool VERY_VERBOSE = true;
#endif

FlameSolver::FlameSolver()
    : jCorrSolver(jCorrSystem)
    , diffusionTestSolver(diffusionTestTerm)
{
    convectionSystem.setSpeciesDomains(transportStartIndices, transportStopIndices);
}

FlameSolver::~FlameSolver()
{
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
        diffusionTerms.push_back(term);
        diffusionSolvers.push_back(integrator);
    }

    resizeAuxiliary();
}

void FlameSolver::run(void)
{
    totalTimer.start();

    double t = tStart;
    double dt;

    long int nTotal = 0; // total number of timesteps taken
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

    while (true) {

        setupTimer.start();
        // Debug sanity check
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

        if (options.usingAdapChem) {
            ckGas->incrementStep();

            // Because AdapChem only assigns values for species in the
            // reduced mechanisms, we need to zero the reaction rates
            // whenever the mechansim at each point could change.
            wDot.data().assign(wDot.data().size(), 0);
        }

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

        // Calculate auxiliary data
        for (size_t j=0; j<nPoints; j++) {
            thermoTimer.start();
            gas.setStateMass(&Y(0,j), T[j]);
            rho[j] = gas.getDensity();
            Wmx[j] = gas.getMixtureMolecularWeight();
            cp[j] = gas.getSpecificHeatCapacity();
            gas.getSpecificHeatCapacities(&cpSpec(0,j));
            gas.getEnthalpies(&hk(0,j));
            thermoTimer.stop();

            transportTimer.start();
            lambda[j] = gas.getThermalConductivity();
            mu[j] = gas.getViscosity();
            gas.getWeightedDiffusionCoefficientsMass(&rhoD(0,j));
            gas.getThermalDiffusionCoefficients(&Dkt(0,j));
            transportTimer.stop();

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

        updateLeftBC();
        updateCrossTerms();
        if (options.xFlameControl) {
            update_xStag(t, true); // calculate the value of rVzero
        }
        convectionSystem.set_rVzero(rVzero);
        dt = options.globalTimestep;

        setupTimer.stop();
        splitTimer.start();

        // *** Set the initial conditions and auxiliary variables for each
        // solver/system, and set the value of the "constant" term(s) to zero.
        // Get the current value of each term's time derivative and
        // diagonalized Jacobian needed to compute the diagonalized
        // approximation used in the other terms, as well as the "predicted"
        // values Uextrap, Textrap, and Yextrap.

        // Source solvers
        // Here, we only compute the constant term, which is needed to calculate
        // the convection solver's constant term. The linear term is calculated
        // only while integrating the source systems.

        // Source solvers initialization
        for (size_t j=0; j<nPoints; j++) {
            sdVector& ySource = sourceSolvers[j].y;
            ySource[kMomentum] = U[j];
            ySource[kEnergy] = T[j];
            for (size_t k=0; k<nSpec; k++) {
                ySource[kSpecies+k] = Y(k,j);
            }
            sourceSolvers[j].t0 = t;
            sourceSolvers[j].initialize();
            sourceTerms[j].wDot.assign(nSpec, 0);
            sourceTerms[j].strainFunction.pin(t);
        }

        sdVector ydotSource(nVars);
        for (size_t j=0; j<nPoints; j++) {
            sourceTerms[j].f(tNow, sourceSolvers[j].y, ydotSource);
            constUprod[j] = ydotSource[kMomentum];
            constTprod[j] = ydotSource[kEnergy];
            for (size_t k=0; k<nSpec; k++) {
                constYprod(k,j) = ydotSource[kSpecies+k];
            }
        }

        // constYprod is needed to determine the region where the diffusion term
        // needs to be integrated for each species
        updateTransportDomain();

        // Convection solver: Initialization
        convectionSystem.setState(U, T, Y);
        convectionSystem.initialize(t);

        // Diffusion solvers: Energy and Momentum
        diffusionSolvers[kMomentum].resize(nPoints, 1, 1);
        diffusionSolvers[kEnergy].resize(nPoints, 1, 1);

        for (size_t k=0; k<nVars; k++) {
            // TODO: Use timestep that is based on each component's diffusivity
            diffusionSolvers[k].initialize(t, options.diffusionTimestep);
        }

        diffusionSolvers[kMomentum].y = U;
        diffusionSolvers[kEnergy].y = T;
        for (size_t j=0; j<nPoints; j++) {
            diffusionTerms[kMomentum].B[j] = 1/rho[j];
            diffusionTerms[kEnergy].B[j] = 1/(rho[j]*cp[j]);

            diffusionTerms[kMomentum].D[j] = mu[j];
            diffusionTerms[kEnergy].D[j] = lambda[j];
        }

        constUdiff = diffusionSolvers[kMomentum].get_ydot();
        constTdiff = diffusionSolvers[kEnergy].get_ydot();
        linearUdiff = diffusionSolvers[kMomentum].get_diagonal();
        linearTdiff = diffusionSolvers[kEnergy].get_diagonal();

        // Diffusion solvers: Species
        for (size_t k=0; k<nSpec; k++) {
            dvector& yDiff_Y = diffusionSolvers[kSpecies+k].y;
            DiffusionSystem& sys = diffusionTerms[kSpecies+k];
            size_t i = 0;
            sys.jLeft = transportStartIndices[k];
            sys.jRight = transportStopIndices[k];

            for (size_t j = transportStartIndices[k];
                 j <= transportStopIndices[k];
                 j++)
            {
                yDiff_Y[i] = Y(k,j);
                sys.B[i] = 1/rho[j];
                sys.D[i] = rhoD(k,j);
                i++;
            }
        }

        for (size_t k=0; k<nSpec; k++) {
            const dvector& constantTerm = diffusionSolvers[kSpecies+k].get_ydot();
            const dvector& linearTerm = diffusionSolvers[kSpecies+k].get_diagonal();

            if (options.transportElimination) {
                // left excluded region
                for (size_t j=0; j<transportStartIndices[k]; j++) {
                    constYdiff(k,j) = 0;
                    linearYdiff(k,j) = 0;
                }

                // transport solution region
                size_t i = 0;
                for (size_t j = transportStartIndices[k];
                     j <= transportStopIndices[k];
                     j++)
                {
                    constYdiff(k,j) = constantTerm[i];
                    linearYdiff(k,j) = linearTerm[i];
                    i++;
                }

                // right excluded region
                for (size_t j=transportStopIndices[k]+1; j<nPoints; j++) {
                    constYdiff(k,j) = 0;
                    linearYdiff(k,j) = 0;
                }

            } else {
                // const_cast required because Array2D::setRow is missing a const qualifier
                constYdiff.setRow(k, const_cast<double*>(&constantTerm[0]));
                linearYdiff.setRow(k, const_cast<double*>(&linearTerm[0]));
            }
        }


        // Convection solver: constant determination
        // Because the convection solver includes the continuity equation containing
        // drho/dt, we need to include the time derivatives from the other terms when
        // evaluating the derivatives of this term, then subtract them from the output
        dvector Utmp = constUprod + constUdiff;
        dvector Ttmp = constTprod + constTdiff + constTcross;
        Array2D Ytmp(nSpec, nPoints);
        for (size_t j=0; j<nPoints; j++) {
            for (size_t k=0; k<nSpec; k++) {
                Ytmp(k,j) = constYprod(k,j) + constYdiff(k,j) + constYcross(k,j);
            }
        }

        convectionSystem.setSplitConst(Utmp, Ttmp, Ytmp, false);

        convectionSystem.evaluate();
        constUconv = convectionSystem.dUdt - Utmp;
        constTconv = convectionSystem.dTdt - Ttmp;

        for (size_t j=0; j<nPoints; j++) {
            for (size_t k=0; k<nSpec; k++) {
                constYconv(k,j) = convectionSystem.dYdt(k,j) - Ytmp(k,j);
            }
        }
        convectionSystem.get_diagonal(tNow, linearUconv, linearTconv, linearYconv);

        // Store the time derivatives for output files
        // TODO: Only do this if we're actually about to write an output file
        dUdtdiff = constUdiff;
        dUdtconv = constUconv;
        dTdtconv = constTconv;
        dTdtdiff = constTdiff;
        dYdtconv = constYconv;
        dYdtdiff = constYdiff;
        dTdtcross = constTcross;
        dYdtcross = constYcross;

        // *** Use the time derivatives to calculate the values for the
        //     state variables based on the diagonalized approximation

        // Offset the value of the linear terms at the current solution
        constUconv -= linearUconv*U;
        constUdiff -= linearUdiff*U;

        constTconv -= linearTconv*T;
        constTdiff -= linearTdiff*T;
        constTcross -= linearTcross*T;

        for (size_t j=0; j<nPoints; j++) {
            for (size_t k=0; k<nSpec; k++) {
                constYconv(k,j) -= linearYconv(k,j)*Y(k,j);
                constYdiff(k,j) -= linearYdiff(k,j)*Y(k,j);
                constYcross(k,j) -= linearYcross(k,j)*Y(k,j);
            }
        }

        // Source terms: calculate constant and linear terms
        for (size_t j=0; j<nPoints; j++) {
            SourceSystem& term = sourceTerms[j];
            term.splitConst[kMomentum] = constUconv[j] + constUdiff[j];
            term.splitLinear[kMomentum] = linearUconv[j] + linearUdiff[j];
            term.splitConst[kEnergy] = constTconv[j] + constTdiff[j] + constTcross[j];
            term.splitLinear[kEnergy] = linearTconv[j] + linearTdiff[j] + linearTcross[j];
            for (size_t k=0; k<nSpec; k++) {
                term.splitConst[kSpecies+k] = constYconv(k,j) + constYdiff(k,j) + constYcross(k,j);
                term.splitLinear[kSpecies+k] = linearYconv(k,j) + linearYdiff(k,j) + linearYcross(k,j);
            }
        }

        splitTimer.stop();

        if (VERY_VERBOSE) {
            cout << "Starting Integration: source terms...";
            cout.flush();
        }

        double tNext = tNow + dt;
        // Source solvers: Integrate, and extract the linear terms
        // needed by the diffusion/convection solvers
        reactionTimer.start();
        for (size_t j=0; j<nPoints; j++) {
            sourceTerms[j].updateDiagonalJac = true;

            int err = sourceSolvers[j].integrateToTime(tNext);
            if (err) {
                cout << "Error at j = " << j << endl;
                cout << "T = " << sourceTerms[j].T << endl;
                cout << "U = " << sourceTerms[j].U << endl;
                cout << "Y = " << sourceTerms[j].Y << endl;
            }
            linearUprod[j] = sourceTerms[j].diagonalJac[kMomentum] -
                    sourceTerms[j].splitLinear[kMomentum];
            linearTprod[j] = sourceTerms[j].diagonalJac[kEnergy] -
                    sourceTerms[j].splitLinear[kEnergy];
            for (size_t k=0; k<nSpec; k++) {
                linearYprod(k,j) = sourceTerms[j].diagonalJac[kSpecies+k] -
                        sourceTerms[j].splitLinear[kSpecies+k];
            }
        }
        reactionTimer.stop();

        splitTimer.resume();
        // Store the time derivatives for output files
        // TODO: Only do this if we're actually about to write an output file
        dUdtprod = constUprod;
        dTdtprod = constTprod;
        dYdtprod = constYprod;

        // *** Use the time derivatives to calculate the values for the
        //     state variables based on the diagonalized approximation

        // Offset the value of the linear terms at the current solution
        constUprod -= linearUprod*U;
        constTprod -= linearTprod*T;

        for (size_t j=0; j<nPoints; j++) {
            for (size_t k=0; k<nSpec; k++) {
                constYprod(k,j) -= linearYprod(k,j)*Y(k,j);
            }
        }

        // Convection term: Set constant & linear terms
        Array2D YtmpLinear(nSpec, nPoints);
        for (size_t j=0; j<nPoints; j++) {
            for (size_t k=0; k<nSpec; k++) {
                Ytmp(k,j) = constYprod(k,j) + constYdiff(k,j) + constYcross(k,j);
                YtmpLinear(k,j) = linearYprod(k,j) + linearYdiff(k,j) + linearYcross(k,j);
            }
        }
        convectionSystem.setSplitConst(constUprod + constUdiff,
            constTprod + constTdiff + constTcross, Ytmp, true);
        convectionSystem.setSplitLinear(linearUprod + linearUdiff,
            linearTprod + linearTdiff + linearTcross, YtmpLinear);

        // Diffusion terms: Calculate constant & linear terms
        diffusionTerms[kMomentum].splitConst = constUconv + constUprod;
        diffusionTerms[kMomentum].splitLinear = linearUconv + linearUprod;
        diffusionTerms[kEnergy].splitConst = constTconv + constTprod + constTcross;
        diffusionTerms[kEnergy].splitLinear = linearTconv + linearTprod + linearTcross;
        for (size_t k=0; k<nSpec; k++) {
            dvector& splitConst = diffusionTerms[kSpecies+k].splitConst;
            dvector& splitLinear = diffusionTerms[kSpecies+k].splitLinear;
            size_t i = 0;
            for (size_t j = transportStartIndices[k];
                 j <= transportStopIndices[k];
                 j++)
            {
                splitConst[i] = constYconv(k,j) + constYprod(k,j) + constYcross(k,j);
                splitLinear[i] = linearYconv(k,j) + linearYprod(k,j) + linearYcross(k,j);
                i++;
            }
        }
        splitTimer.stop();

        // *** Take one global timestep: diffusion / convection solvers
        int cvode_flag = 0;
        try {
            if (VERY_VERBOSE) {
                cout << "diffusion terms...";
                cout.flush();
            }
            diffusionTimer.start();
            for (size_t k=0; k<nVars; k++) {
                diffusionSolvers[k].integrateToTime(tNext);
            }
            diffusionTimer.stop();

            if (VERY_VERBOSE) {
                cout << "convection term...";
                cout.flush();
            }

            convectionTimer.start();
            convectionSystem.integrateToTime(tNext);
            convectionTimer.stop();
            if (VERY_VERBOSE) {
                cout << "done!" << endl;
            }
        } catch (Cantera::CanteraError) {
            writeStateFile("errorOutput", true);
            cout << "Integration failed at t = " << tNow << std::endl;
        }

        combineTimer.start();
        t = tNext;
        tNow = tNext;
        aPrev = strainfunc.a(t);

        if (cvode_flag == CV_SUCCESS) {
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
        } else {
            cout << format("CVODE Solver failed at time t = %8.6f (dt = %9.3e)") %
                    t % dt << endl;
            writeStateFile("errorOutput",true);
            break;
        }

        // store Interpolation data for V used in the species-split convection
        // solver for the output file
        vInterp.resize(convectionSystem.vInterp->size(), nPoints);
        tvInterp.resize(convectionSystem.vInterp->size());
        size_t i = 0;
        typedef std::pair<double,dvector> dobuledvectorpair;
        foreach (const dobuledvectorpair& item, *convectionSystem.vInterp) {
            for (size_t j=0; j<nPoints; j++) {
                vInterp(i,j) = item.second[j];
                tvInterp[i] = item.first;
            }
            i += 1;
        }

        // This is the exact solution to the linearized problem
        for (size_t j=0; j<nPoints; j++) {
            double K = constUprod[j] + constUdiff[j] + constUconv[j];
            double L = linearUprod[j] + linearUdiff[j] + linearUconv[j];
            double L2; // L2 handles the limit as L->0 for the second term in the solution
            L2 = (abs(L*dt) < 1e-10) ? 1e-10/dt : L;

            Uextrap[j] = U[j]*exp(L*dt) + K/L2*(exp(L2*dt)-1);

            K = constTprod[j] + constTdiff[j] + constTconv[j] + constTcross[j];
            L = linearTprod[j] + linearTdiff[j] + linearTconv[j] + linearTcross[j];
            L2 = (abs(L*dt) < 1e-10) ? 1e-10/dt : L;
            Textrap[j] = T[j]*exp(L*dt) + K/L2*(exp(L2*dt)-1);

            for (size_t k=0; k<nSpec; k++) {
                K = constYprod(k,j) + constYdiff(k,j) + constYconv(k,j) + constYcross(k,j);
                L = linearYconv(k,j) + linearYdiff(k,j) + linearYprod(k,j) + linearYcross(k,j);
                L2 = (abs(L*dt) < 1e-10) ? 1e-10/dt : L;
                Yextrap(k,j) = Y(k,j)*exp(L*dt) + K/L2*(exp(L2*dt)-1);
            }
        }

        // *** Combine the solutions from the split integrators and form
        //     the new state vector
        convectionSystem.unroll_y();
        for (size_t j=0; j<nPoints; j++) {
            sourceTerms[j].unroll_y(sourceSolvers[j].y);
        }

        for (size_t j=0; j<nPoints; j++) {
            U[j] = convectionSystem.U[j] + sourceTerms[j].U + diffusionSolvers[kMomentum].y[j] - 2*Uextrap[j];
            T[j] = convectionSystem.T[j] + sourceTerms[j].T + diffusionSolvers[kEnergy].y[j] - 2*Textrap[j];
            for (size_t k=0; k<nSpec; k++) {
                Y(k,j) = convectionSystem.Y(k,j) + sourceTerms[j].Y[k] - Yextrap(k,j);
            }
        }

        for (size_t k=0; k<nSpec; k++) {
            size_t i = 0;
            const BDFIntegrator& solver = diffusionSolvers[kSpecies+k];
            for (size_t j = transportStartIndices[k];
                 j <= transportStopIndices[k];
                 j++)
            {
                Y(k,j) += solver.y[i] - Yextrap(k,j);
                i++;
            }
        }

        combineTimer.stop();

        setupTimer.resume();
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
        //     Quit now to skip grid adaptation on the last step
        if (nTerminate >= options.terminateStepInterval) {
            nTerminate = 0;
            if (checkTerminationCondition()) {
                break;
            }
        }

        if (t >= tEnd) {
            break;
        }

        // call ckGas::adapBox

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
                dvector dtmp(nPoints);
                for (size_t j=0; j<nPoints; j++) {
                    tmp[j] = Y(k,j);
                    dtmp[j] = dYdt(k,j);
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
                grid.updated = false;
                nIntegrate = 0;
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

                // Assign the new grid to the terms that need it
                for (size_t k=0; k<nVars; k++) {
                    diffusionTerms[k].setGrid(grid);
                }
                convectionSystem.setGrid(grid);

                // Update the mass flux at the left boundary
                rVzero = mathUtils::interp1(x_prev, V_prev, grid.x[0]);

                // Allocate the solvers and arrays for auxiliary variables
                resizeAuxiliary();

                writeStateFile("postAdapt");
            }
            regridTimer.stop();
        }

        setupTimer.resume();
        for (size_t j=0; j<nPoints; j++) {
            // Correct the drift of the total mass fractions and reset
            // any negative mass fractions
            gas.setStateMass(&Y(0,j), T[j]);
            gas.getMassFractions(&Y(0,j));
        }
        setupTimer.stop();

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
        outFile.writeVector("V", convectionSystem.V);
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

        dvector Uprod(nPoints);
        dvector Tprod(nPoints);
        Array2D Ydiff(nSpec, nPoints, 0), Yprod(nSpec, nPoints, 0);

        for (size_t j=0; j<nPoints; j++) {
            Uprod[j] = sourceTerms[j].U;
            Tprod[j] = sourceTerms[j].T;
            for (size_t k=0; k<nSpec; k++) {
                Yprod(k,j) = sourceTerms[j].Y[k];
            }
        }

        for (size_t k=0; k<nSpec; k++) {
            size_t i = 0;
            if (diffusionSolvers[kSpecies+k].y.empty()) {
                continue;
            }
            for (size_t j = transportStartIndices[k];
                 j <= transportStopIndices[k];
                 j++)
            {
                Ydiff(k,j) = diffusionSolvers[kSpecies+k].y[i];
                i++;
            }
        }

        outFile.writeVector("Uconv", convectionSystem.U);
        outFile.writeVector("Udiff", diffusionSolvers[kMomentum].y);
        outFile.writeVector("Uprod", Uprod);
        outFile.writeVector("Tconv", convectionSystem.T);
        outFile.writeVector("Tdiff", diffusionSolvers[kEnergy].y);
        outFile.writeVector("Tprod", Tprod);
        outFile.writeArray2D("Yconv", convectionSystem.Y);
        outFile.writeArray2D("Ydiff", Ydiff);
        outFile.writeArray2D("Yprod", Yprod);
        outFile.writeVector("Wconv", convectionSystem.Wmx);

        outFile.writeVector("Uextrap", Uextrap);
        outFile.writeVector("Textrap", Textrap);
        outFile.writeArray2D("Yextrap", Yextrap);

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

        if (options.transportElimination) {
            dvector jStart, jStop;
            for (size_t k=0; k<nSpec; k++) {
                jStart.push_back(transportStartIndices[k]);
                jStop.push_back(transportStopIndices[k]);
            }
            outFile.writeVector("transportStartIndices", jStart);
            outFile.writeVector("transportStopIndices", jStop);
        }

        outFile.writeVector("tvInterp", tvInterp);
        outFile.writeArray2D("vInterp", vInterp);
    }

    if (options.outputResidualComponents || errorFile) {
        outFile.writeVector("constTdiff", constTdiff);
        outFile.writeVector("constTconv", constTconv);
        outFile.writeVector("constTprod", constTprod);
        outFile.writeVector("constUdiff", constUdiff);
        outFile.writeVector("constUconv", constUconv);
        outFile.writeVector("constUprod", constUprod);
        outFile.writeArray2D("constYdiff", constYdiff);
        outFile.writeArray2D("constYconv", constYconv);
        outFile.writeArray2D("constYprod", constYprod);
        outFile.writeArray2D("constYcross", constYcross);
        outFile.writeVector("constTcross", constTcross);

        outFile.writeVector("dTdtdiff", dTdtdiff);
        outFile.writeVector("dTdtconv", dTdtconv);
        outFile.writeVector("dTdtprod", dTdtprod);
        outFile.writeVector("dTdtcross", dTdtcross);
        outFile.writeVector("dUdtdiff", dUdtdiff);
        outFile.writeVector("dUdtconv", dUdtconv);
        outFile.writeVector("dUdtprod", dUdtprod);
        outFile.writeArray2D("dYdtdiff", dYdtdiff);
        outFile.writeArray2D("dYdtconv", dYdtconv);
        outFile.writeArray2D("dYdtprod", dYdtprod);
        outFile.writeArray2D("dYdtcross", dYdtcross);
        outFile.writeVector("dWdtconv", convectionSystem.dWdt);
        outFile.writeVector("constW", convectionSystem.utwSystem.splitConstW);

        outFile.writeVector("linearTdiff", linearTdiff);
        outFile.writeVector("linearTconv", linearTconv);
        outFile.writeVector("linearTprod", linearTprod);
        outFile.writeVector("linearUdiff", linearUdiff);
        outFile.writeVector("linearUconv", linearUconv);
        outFile.writeVector("linearUprod", linearUprod);
        outFile.writeArray2D("linearYdiff", linearYdiff);
        outFile.writeArray2D("linearYconv", linearYconv);
        outFile.writeArray2D("linearYprod", linearYprod);
        outFile.writeArray2D("linearTcross", linearYcross);
        outFile.writeVector("linearTcross", linearTcross);

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

    constTcross.resize(nPoints,0);
    linearTcross.resize(nPoints,0);
    constYcross.resize(nSpec, nPoints, 0);
    linearYcross.resize(nSpec, nPoints, 0);

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
    transportStartIndices.assign(nSpec, 0);
    transportStopIndices.assign(nSpec, jj);
    nPointsTransport.assign(nSpec, nPoints);
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
    // With transport species elimination enabled, this is done later
    // to handle the varying number of grid points for each species
    if (options.transportElimination) {
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

    convectionSystem.resize(nPoints, nPointsTransport, nSpec);
    convectionSystem.setGrid(grid);

    // Resize the jCorr stabilizer
    jCorrSolver.resize(nPoints, 1, 1);
    jCorrSystem.setGrid(grid);
}

void FlameSolver::updateCrossTerms()
{
    for (size_t j=0; j<jj; j++) {
        jCorr[j] = 0;
        for (size_t k=0; k<nSpec; k++) {
            jFick(k,j) = -0.5*(rhoD(k,j)+rhoD(k,j+1)) * ((Y(k,j+1)-Y(k,j))/hh[j]);
            jSoret(k,j) = -0.5*(Dkt(k,j)/T[j] + Dkt(k,j+1)/T[j+1])
                * (T[j+1]-T[j])/hh[j];
            jCorr[j] -= jFick(k,j) + jSoret(k,j);
        }
    }

    // Add a bit of artificial diffusion to jCorr to improve stability
    jCorrSolver.y = jCorr;
    for (size_t j=0; j<jj; j++) {
        jCorrSystem.B[j] = lambda[j]/(rho[j]*cp[j]); // treat as Le = 1
        jCorrSystem.D[j] = 1.0;
    }

    jCorrSolver.initialize(0, options.diffusionTimestep);
    jCorrSolver.integrateToTime(options.globalTimestep);
    jCorr = jCorrSolver.y;

    for (size_t j=1; j<jj; j++) {
        sumcpj[j] = 0;
        for (size_t k=0; k<nSpec; k++) {
            constYcross(k,j) = -0.5/(r[j]*rho[j]*dlj[j]) *
                (rphalf[j]*(Y(k,j)+Y(k,j+1))*jCorr[j] - rphalf[j-1]*(Y(k,j-1)+Y(k,j))*jCorr[j-1]);
            constYcross(k,j) -= 1/(r[j]*rho[j]*dlj[j]) *
                (rphalf[j]*jSoret(k,j) - rphalf[j-1]*jSoret(k,j-1));
            sumcpj[j] += 0.5*(cpSpec(k,j)+cpSpec(k,j+1))/W[k]*(jFick(k,j)
                    + jSoret(k,j) + 0.5*(Y(k,j)+Y(k,j+1))*jCorr[j]);

            linearYcross(k,j) = -0.5/(r[j]*rho[j]*dlj[j]) *
                    (rphalf[j]*jCorr[j] - rphalf[j-1]*jCorr[j-1]);
        }
        double dTdx = cfm[j]*T[j-1] + cf[j]*T[j] + cfp[j]*T[j+1];
        constTcross[j] = - 0.5*(sumcpj[j]+sumcpj[j-1]) * dTdx / (cp[j]*rho[j]);
        linearTcross[j] = - 0.5*(sumcpj[j]+sumcpj[j-1]) * cf[j] / (cp[j]*rho[j]);
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
    printPerfString("     Split Term Recombination: ", combineTimer);
    cout << endl << " Subcomponents:" << endl;
    printPerfString("               Reaction Rates: ", reactionRatesTimer);
    printPerfString("         Transport Properties: ", transportTimer);
    printPerfString("     Thermodynamic Properties: ", thermoTimer);
    printPerfString("   Source Jacobian Evaluation: ", jacobianTimer);
    cout << endl;
}

void FlameSolver::printPerfString(const std::string& label, const perfTimer& T) const
{
    cout << format("%s %9.3f (%12i)") % label % T.getTime() % T.getCallCount() << endl;
}

void FlameSolver::updateTransportDomain()
{
    if (!options.transportElimination) {
        return;
    }

    // TODO: Use timestep that is based on each component's diffusivity
    diffusionTestSolver.initialize(0, options.diffusionTimestep);

    for (size_t k=0; k<nSpec; k++) {
        // evaluate the full diffusion term for each species
        for (size_t j=0; j<nPoints; j++) {
            diffusionTestSolver.y[j] = Y(k,j);
            diffusionTestTerm.B[j] = 1/rho[j];
            diffusionTestTerm.D[j] = rhoD(k,j);
        }
        const dvector& dYkdt_diff = diffusionTestSolver.get_ydot();
        constYdiff.setRow(k, const_cast<double*>(&dYkdt_diff[0]));
    }

    // Evaluate the full convection term
    // Because the convection solver includes the continuity equation containing
    // drho/dt, we need to include the time derivatives from the other terms when
    // evaluating the derivatives of this term, then subtract them from the output
    convectionSystem.setState(U, T, Y);
    dvector Utmp = constUprod + constUdiff;
    dvector Ttmp = constTprod + constTdiff + constTcross;
    Array2D Ytmp(nSpec, nPoints);
    for (size_t j=0; j<nPoints; j++) {
        for (size_t k=0; k<nSpec; k++) {
            Ytmp(k,j) = constYprod(k,j) + constYdiff(k,j) + constYcross(k,j);
        }
    }
    convectionSystem.setSplitConst(Utmp, Ttmp, Ytmp, false);

    convectionSystem.evaluate();
    constUconv = convectionSystem.dUdt - Utmp;
    constTconv = convectionSystem.dTdt - Ttmp;
    for (size_t j=0; j<nPoints; j++) {
        for (size_t k=0; k<nSpec; k++) {
            constYconv(k,j) = convectionSystem.dYdt(k,j) - Ytmp(k,j);
        }
    }

    for (size_t k=0; k<nSpec; k++) {
        // Find the left boundary for species k
        int jStart = jj;
        for (size_t j=0; j<nPoints; j++) {
            if (abs(constYdiff(k,j)+constYconv(k,j)) > options.adapchem_atol ||
                abs(constYprod(k,j)) > options.adapchem_atol) {
                jStart = j;
                break;
            }
        }
        jStart = max(0, jStart-2);
        transportStartIndices[k] = jStart;

        // Find the right boundary for species k
        size_t jStop = jStart;
        for (int j=jj; j>jStart; j--) {
            if (abs(constYdiff(k,j)+constYconv(k,j)) > options.adapchem_atol ||
                abs(constYprod(k,j)) > options.adapchem_atol) {
                jStop = j;
                break;
            }
        }
        jStop = min(jj, jStop+2);

        transportStopIndices[k] = jStop;
        nPointsTransport[k] = jStop - jStart + 1;
    }

    for (size_t k=0; k<nSpec; k++) {
        // size the Diffusion systems appropriately
        diffusionSolvers[kSpecies+k].resize(nPointsTransport[k], 1, 1);
    }
}
