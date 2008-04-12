#include "flameSolver.h"

void flameSolver::setOptions(const configOptions& theOptions)
{
	options = theOptions;
}

void flameSolver::run(void)
{
	clock_t t1, t2;
	t1 = clock();

	theSys.options = options;

	theSys.copyOptions();

	theSys.gas.initialize();
	// theSys.gas.testFunction();

	// Initial Conditions for ODE
	theSys.setup();
	if (options.haveRestartFile) {
		theSys.loadInitialProfiles();
	} else {
		theSys.generateInitialProfiles();
	}

	double integratorTimestep = 0;
	double t = theSys.tStart;
	double dt;

	int nRegrid = 0;
	int nOutput = 0;
	int nProfile = 0;
	int nFlamePos = 0;
	
	double tOutput = t;
	double tRegrid = t;
	double tProfile = t;
	double tFlamePos = t;

	theSys.grid.updateValues();

	// Get the initial value for rFlameActual:
	theSys.gas.setStateMass(theSys.Y,theSys.T);
	theSys.updateThermoProperties();
	theSys.gas.getReactionRates(theSys.wDot);
	for (int j=0; j<=theSys.nPoints-1; j++) {
		theSys.qDot[j] = 0;
		for (int k=0; k<theSys.nSpec; k++) {
			theSys.qDot[j] -= theSys.wDot(k,j)*theSys.hk(k,j);
		}
	}

	// Flame position (radius) control:
	theSys.tFlamePrev = t;
	theSys.tFlameNext = t + options.rFlameUpdateTimeInterval;
	theSys.rVcenterInitial = theSys.V[0];
	theSys.rVcenterPrev = theSys.rVcenterNext = theSys.rVcenterInitial;
	theSys.writeStateMatFile();
	bool firstIteration = true;

	while (t < theSys.tEnd) {

		theSys.setup();

		// Sundials IDA Solver:
		sundialsIDA theSolver(theSys.N);
		theSolver.reltol = options.idaRelTol;
		theSolver.nRoots = 0;
		theSolver.findRoots = false;

		int N = theSys.nVars;
		// Initial condition:
		theSys.rollY(theSolver.y);
		for (int j=0; j<theSys.nPoints; j++) {
			theSolver.abstol(N*j) = options.idaContinuityAbsTol;
			theSolver.abstol(N*j+1) = options.idaMomentumAbsTol;
			theSolver.abstol(N*j+2) = options.idaEnergyAbsTol;
			for (int k=0; k<theSys.nSpec; k++) {
				theSolver.abstol(N*j+k+3) = options.idaSpeciesAbsTol;
			}
		}

		for (int j=0; j<theSys.N; j++) {
			theSolver.ydot(j) = 0;
		}

		theSys.updateLeftBC();

		theSys.update_rVcenter(t);
		tFlamePos = t + options.rFlameUpdateTimeInterval;
		nFlamePos = 0;

		theSys.updateAlgebraicComponents();
		theSolver.t0 = t;
		int ICflag = -1;
		int ICcount = 0;
		while (ICflag!=0 && ICcount < 5) {
			ICcount++;
			ICflag = theSys.getInitialCondition(t, theSolver.y, theSolver.ydot, theSys.algebraic);
		}

		theSolver.setDAE(&theSys);
		theSolver.calcIC = false;

		theSolver.initialize();
		theSolver.setMaxStepSize(options.maxTimestep);
		
		if (integratorTimestep != 0) {
			theSolver.setInitialStepSize(integratorTimestep);
		}

		int flag;

		while (t < theSys.tEnd) {

			try {
				flag = theSolver.integrateOneStep();
			} catch (Cantera::CanteraError) {
				theSys.writeStateMatFile("errorOutput",true);
				throw;
			}

			firstIteration = false;
			dt = integratorTimestep = theSolver.getStepSize();
			t = theSys.tPrev = theSolver.tInt;

			if (flag == CV_SUCCESS) {
				nOutput++;
				nRegrid++;
				nProfile++;
				nFlamePos++;
				cout << "t = " << t << "  (dt = " << dt << ")" << endl;
			} else {
				cout << "IDA Solver failed at time t = " << t << "  (dt = " << dt << ")" << endl;
				theSys.writeStateMatFile("errorOutput",true);
				integratorTimestep = 0;
				break;
			}
			
			if (t > tOutput || nOutput >= options.outputStepInterval) {
				// Save the time-series data
				timeVector.push_back(t);
				timestepVector.push_back(dt);
				heatReleaseRate.push_back(theSys.getHeatReleaseRate());
				consumptionSpeed.push_back(theSys.getConsumptionSpeed());
				flamePosition.push_back(theSys.getFlamePosition());

				tOutput = t + options.outputTimeInterval;
				nOutput = 0;
			}

			if (t > tProfile || nProfile >= options.profileStepInterval) {
				theSys.writeStateMatFile();

				tProfile = t + options.profileTimeInterval;
				nProfile = 0;
			}

			if (options.flameRadiusControl && 
				(t > tFlamePos || nFlamePos > options.rFlameUpdateStepInterval)) {
				theSys.update_rVcenter(t);

				tFlamePos = t + options.rFlameUpdateTimeInterval;
				nFlamePos = 0;
			}

			if (t > tRegrid || nRegrid >= options.regridStepInterval) {
				tRegrid = t + options.regridTimeInterval;
				nRegrid = 0;

				// Periodic check for terminating the integration
				// (based on steady heat release rate, etc.)
				if (checkTerminationCondition()) {
					theSolver.printStats();
					t2 = clock();
					theSys.writeStateMatFile();
					cout << "Runtime: " << ((double)(t2-t1))/CLOCKS_PER_SEC << " seconds." << endl;
					return;
				}

				// Adapt the grid if necessary

 				for (int j=0; j<theSys.nPoints; j++) {
					double num = min(theSys.mu[j],theSys.lambda[j]/theSys.cp[j]);
					for (int k=0; k<theSys.nSpec; k++) {
						num = min(num,theSys.rhoD(k,j));
					}

					theSys.grid.dampVal[j] = num/abs(theSys.V[j]);
				}
				vector<dvector> currentSolution, currentSolutionDot;
				theSys.rollVectorVector(theSolver.y, currentSolution);
				theSys.rollVectorVector(theSolver.ydot, currentSolutionDot);

				bool adaptFlag = theSys.grid.adapt(currentSolution, currentSolutionDot);
				bool regridFlag = theSys.grid.regrid(currentSolution, currentSolutionDot);

				if (adaptFlag || regridFlag) {
					theSys.nPoints = theSys.grid.jj+1;
					cout << "Grid size: " << theSys.nPoints << " points." << endl;
					theSys.setup();

					theSys.unrollVectorVector(currentSolution);
					theSys.unrollVectorVectorDot(currentSolutionDot);

					// This corrects the drift of the total mass fractions
					theSys.gas.setStateMass(theSys.Y,theSys.T);
					theSys.gas.getMassFractions(theSys.Y);

					break; // exit the inner loop and reinitialize the solver for the new problem size
				}

			}
		}
		
		theSolver.printStats();
		
	}

	t2 = clock();
	theSys.writeStateMatFile();
	cout << "Runtime: " << ((double)(t2-t1))/CLOCKS_PER_SEC << " seconds." << endl;
}

void flameSolver::calculateReactantMixture(void)
{	
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
			cout << "Continuing integration: t < terminationPeriod" << endl;
			return false;
		}
		int j2 = timeVector.size()-1;
		double qMean = mathUtils::mean(heatReleaseRate,j1,j2);
		double hrrError = 0;
		for (int j=j1; j<=j2; j++) {
			hrrError += abs(heatReleaseRate[j]-qMean);
		}
		hrrError /= (j2-j1+1);
		if (hrrError/qMean < options.terminationTolerance) {
			cout << "Terminating integration: Heat release rate deviation = " << hrrError/qMean*100 << "%" << endl;
			return true;
		} else {
			cout << "Continuing integration: Heat release rate deviation = " << hrrError/qMean*100 << "%" << endl;
		}

	}
	return false;
}
