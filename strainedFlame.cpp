#include "strainedFlame.h"
#include "debugUtils.h"
#include "matlabFile.h"
#include "boost/filesystem.hpp"
#include "sundialsUtils.h"

using namespace mathUtils;

int main(int argc, char** argv)
{

	// Cantera gets angry if it can't find Python
	char* python_cmd = getenv("PYTHON_CMD");
	if (!python_cmd) {
		putenv("PYTHON_CMD=python");
	}

	std::string inputFile;
	if (argc > 1) {
		inputFile = argv[1];
	} else {
		inputFile = "input.txt"; // default input filename
	}

    try {
    	strainedFlame(inputFile);
		//chemistryTest();
		//matlabioTest();
		//miscTest();

    }
	catch (Cantera::CanteraError) {
		Cantera::showErrors(cout);
    }
    
	return 0;
}

void strainedFlame(const std::string& inputFile) 
{
    cout << "**** strainedFlame (1dflame Version 2.0) ****\n" << std::endl;

	clock_t t1, t2;
	t1 = clock();

    strainedFlameSys theSys;
	configOptions& options = theSys.options;

	theSys.readOptionsFile(inputFile);

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

void chemistryTest(void)
{
	using namespace Cantera;
    XML_Node *xc = get_XML_File("gri30.xml");
    XML_Node * const xs = xc->findNameID("phase", "gri30_mix");

	Cantera::XML_Node* foo = NULL;
	delete foo;

	int n = 1;
	clock_t t1, t2;
	t1 = clock();
	Cantera::IdealGasPhase thermoBase;
	Cantera::importPhase(*xs,&thermoBase);
	vector<Cantera::IdealGasPhase> gas(n,thermoBase);
	vector<Cantera::GasKinetics*> kin(n);
	vector<Cantera::Transport*> trans(n);
	for (int i=0; i<n; i++) {
		//gas[i] = new Cantera::IdealGasPhase();
		//Cantera::importPhase(*xs, gas[i]);
		kin[i] = new Cantera::GasKinetics(&gas[i]);
		
		kin[i]->init();
		Cantera::installReactionArrays(*xs,*kin[i],"gri30_mix");
		kin[i]->finalize();

		trans[i] = Cantera::newTransportMgr("Mix",&gas[i],1,0);
	
	}

	t2 = clock();
	cout << "separate: " << t2-t1 << endl;

	t1 = clock();
	gasArray gas2;
	gas2.mechanismFile = "gri30.xml";
	gas2.phaseID = "gri30_mix";
	gas2.initialize();
	gas2.resize(n);
	t2 = clock();
	cout << "gasArray: " << t2-t1 << endl;

	int nSpec = gas[0].nSpecies();
	dvector dkm(nSpec);
	
	dvector y(nSpec);
	gas[0].setState_TPX(300,101325,"O2:1.0, CH4:0.5");
	gas[0].getMassFractions(&y[0]);
	t1 = clock();
	for (int i=0; i<2000; i++) {
		y[1] = 0.005*i;
		gas[0].setMassFractions(&y[0]);
		//trans[0]->getMixDiffCoeffs(&dkm[0]);
	}

	t2 = clock();
	cout << "getMixDiffCoeffs: " << t2-t1 << endl;
//	cout << "mu = " << trans[0]->viscosity() << endl;

	//Cantera::IdealGasPhase gas;
	//Cantera::importPhase(*xs, &gas);
	//Cantera::GasKinetics k;
	//k.addPhase(gas);
	//k.init();
	//Cantera::installReactionArrays(*xs,k,"gri30_mix");
	//k.finalize();

	//gas.setState_TPX(700,101325,"O2:1.0, CH4:0.5");
	//dvector wdot(gas.nSpecies());
	//k.getNetProductionRates(&wdot[0]);
	//cout << wdot << endl;

	//gasArray theGas;
	//theGas.mechanismFile = "gri30.xml";
	//theGas.resize(2);
	//theGas[0].setState_TPX(300,Cantera::OneAtm,"O2:1.0, N2:3.76, CH4:0.5");
	//theGas[1].setState_TPX(500,Cantera::OneAtm,"O2:1.0, N2:3.76, H2:1.0");
	//int n = theGas[0].nSpecies();

	//dvector wdot0(n), wdot1(n);

	//theGas[0].getNetProductionRates(&wdot0[0]);
	//	
	//cout << wdot0 << endl;
	//cout << wdot1 << endl;

	int blargh = 0;
}

void matlabioTest(void)
{
	int n = 2, m=20;
	dvector x(n*m);
	Array2D foo;
	for (int i=0; i<n*m; i++) {
		x[i] = 2*i + 100;
		
	}

	n = 5, m=20;
	foo.resize(n,m);
	for (int i=0; i<n; i++) {
		for (int j=0; j<m; j++) {
			foo(i,j) = i+100*j;
		}
	}
	boost::filesystem::remove_all("output");
	boost::filesystem::create_directory("output");

	matlabFile outFile("output/test.mat");
	outFile.writeVector("hellofriend",x);
	outFile.writeArray2D("hugstiem",foo);
	

	dvector hi = outFile.readVector("hellofriend");
	Array2D bar = outFile.readArray2D("hugstiem");
	cout << hi << endl;	

	outFile.close();
}

void miscTest(void)
{
	Cantera::Array2D a(53,100);
	for (int k=0; k<53; k++) {
		for (int j=0; j<100; j++) {
			a(k,j) = k;
		}
	}

	vector<dvector> v;
	array2DToVectorVector(a,v);

	vectorVectorToArray2D(v,a);

	int blargh = 0;
}
