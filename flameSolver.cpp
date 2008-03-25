#include "flameSolver.h"

void flameSolver::setOptions(const configOptions& theOptions)
{
	options = theOptions;
}

void flameSolver::run(void)
{
	clock_t t1, t2;
	t1 = clock();

    flameSys theSys;
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
	int nRegrid = 0;
	int nOutput = 0;
	double tOutput = theSys.tStart;
	double tRegrid = theSys.tStart;
	double t = theSys.tStart;
	double dt;

	theSys.grid.updateValues();

	theSys.writeStateMatFile();

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
		theSys.updateAlgebraicComponents();
		theSolver.t0 = t;
		theSys.getInitialCondition(t, theSolver.y, theSolver.ydot, theSys.algebraic);

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
				theSys.writeErrorFile();
				throw;
			}

			dt = integratorTimestep = theSolver.getStepSize();
			t = theSolver.tInt;

			if (flag == CV_SUCCESS) {
				nOutput++;
				nRegrid++;
				cout << "t = " << t << "  (dt = " << dt << ")" << endl;
			} else {
				cout << "IDA Solver failed at time t = " << t << "  (dt = " << dt << ")" << endl;
				cout << "Writing errorOutput.mat." << endl;
				theSys.writeErrorFile();
				integratorTimestep = 0;
				break;
			}

			if (t > tOutput || nOutput >= options.outputStepInterval) {
				theSys.writeStateMatFile();
				while (t > tOutput) {
					tOutput += options.outputTimeInterval;
				}
				nOutput = 0;
			}

			if (t > tRegrid || nRegrid >= options.regridStepInterval) {
				while (t > tRegrid) {
					tRegrid += options.regridTimeInterval;
				}
				nRegrid = 0;
				// Adapt the grid if necessary

 				for (int j=0; j<theSys.nPoints; j++) {
					theSys.grid.dampVal[j] = abs(theSys.mu[j]/theSys.rhov[j]);
				}
				vector<dvector> currentSolution, currentSolutionDot;

				theSys.rollVectorVector(theSolver.y, currentSolution);
				theSys.rollVectorVector(theSolver.ydot, currentSolutionDot);

				bool adaptFlag = theSys.grid.adapt(currentSolution, currentSolutionDot);
				bool regridFlag = theSys.grid.regrid(currentSolution, currentSolutionDot);

				if (adaptFlag || regridFlag) {
					theSys.nPoints = theSys.grid.jj+1;
					theSys.setup();

					theSys.unrollVectorVector(currentSolution);
					theSys.unrollVectorVectorDot(currentSolutionDot);

					break; // exit the inner loop and reinitialize the solver for the new problem size
				}

			}
		}
		
		theSolver.printStats();
		
	}

	t2 = clock();
	cout << "Runtime: " << ((double)(t2-t1))/CLOCKS_PER_SEC << " seconds." << endl;

}