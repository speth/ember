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
    cout << "**** strainedFlame (1dflame Version 2.0) ****" << std::endl;

	clock_t t1, t2;
	t1 = clock();

    strainedFlameSys theSys;

	theSys.readOptionsFile(inputFile);

	theSys.gas.initialize();
	//theSys.gas.testFunction();

	// output file:
	ofstream outFile;
	std::string outFileName = theSys.options.outputDir+"/out.m";
	outFile.open(outFileName.c_str());

	bool newSolver = true; 
	double nSteps = 0;

	// Initial Conditions for ODE
	theSys.setup();
	if (theSys.options.haveRestartFile) {
		theSys.loadInitialProfiles();
	} else {
		theSys.generateInitialProfiles();
	}
	
	double integratorTimestep = 0;
	int nRegrid = 0;
	int nOutput = 0;
	double tOutput = theSys.tStart;
	double tRegrid = theSys.tStart;

	int i=1;

	double t = theSys.tStart;

	theSys.grid.jj = theSys.nPoints;
	theSys.grid.updateDerivedSizes();

	while (t < theSys.tEnd) {

		theSys.setup();

		// Sundials IDA Solver:
		sundialsIDA theSolver(theSys.N);
		theSolver.reltol = 5e-4;
		theSolver.nRoots = 0;
		theSolver.findRoots = false;
		vector<bool> algebraic(theSys.N);

		int nVars = theSys.nVars;
		// Initial condition:
		theSys.rollY(theSolver.y);
		for (int j=0; j<theSys.nPoints; j++) {
			theSolver.abstol(nVars*j) = 1e-5; // rhov
			theSolver.abstol(nVars*j+1) = 1e-5; // U
			theSolver.abstol(nVars*j+2) = 1e-5; // T
			algebraic[nVars*j] = true;
			algebraic[nVars*j+1] = false;
			algebraic[nVars*j+2] = false;
			for (int k=0; k<theSys.nSpec; k++) {
				theSolver.abstol(nVars*j+k+3) = 1e-1; // Y
				algebraic[nVars*j+3+k] = false;
			}
		}
		algebraic[nVars*theSys.grid.jZero] = false;

		for (int j=0; j<theSys.N; j++) {
			theSolver.ydot(j) = 0;
		}

		theSolver.t0 = t;
		theSys.getInitialCondition(t, theSolver.y, theSolver.ydot, algebraic);

		theSolver.setDAE(&theSys);
		theSolver.calcIC = false;

		theSolver.initialize();
		
		if (integratorTimestep != 0) {
			theSolver.setInitialStepSize(integratorTimestep);
		}

		int flag;

		theSys.writeStateMatFile();
		i++;

		while (t < theSys.tEnd) {
			flag = theSolver.integrateOneStep();

			nOutput++;
			nRegrid++;
			t = theSolver.tInt;

			if (flag == CV_SUCCESS) {
				cout << "t = " << t << endl;
			}

			i++;

			if (t > tOutput || nOutput > theSys.options.outputStepInterval) {
				theSys.writeStateMatFile();
				while (t > tOutput) {
					tOutput += theSys.options.outputTimeInterval;
				}
				nOutput = 0;
			}

			if (t > tRegrid || nRegrid > theSys.options.regridStepInterval) {
				while (t > tRegrid) {
					tRegrid += theSys.options.regridTimeInterval;
				}
				nRegrid = 0;
				// Adapt the grid if necessary

 				for (int j=0; j<theSys.nPoints; j++) {
					theSys.grid.dampVal[j] = abs(theSys.mu/theSys.rhov[j]);
				}
				vector<dvector> currentSolution, currentSolutionDot;

				theSys.rollVectorVector(theSolver.y, currentSolution);
				theSys.rollVectorVector(theSolver.ydot, currentSolutionDot);

				bool adaptFlag = theSys.grid.adapt(currentSolution, currentSolutionDot);
				bool regridFlag = theSys.grid.regrid(currentSolution, currentSolutionDot);
				
				theSys.nPoints = theSys.grid.jj+1;
				theSys.setup();

				theSys.unrollVectorVector(currentSolution);
				theSys.unrollVectorVectorDot(currentSolutionDot);

				if (adaptFlag || regridFlag) {
					break; // exit the inner loop and reinitialize the solver for the new problem size
				}

			}
		}
		integratorTimestep = theSolver.getStepSize();
		theSolver.printStats();
		
	}

	t2 = clock();
	cout << "Runtime: " << ((double)(t2-t1))/CLOCKS_PER_SEC << " seconds." << endl;
	
	// Test of grid adaptation:
	outFile.close();

   	int blargh = 0;
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
