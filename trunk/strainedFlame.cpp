#include "strainedFlame.h"
#include "debugUtils.h"
#include "matlabFile.h"
#include "boost/filesystem.hpp"
#include "sundialsUtils.h"
#include "flameSolver.h"
#include "libconfig.h++"
#include <omp.h>

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
	} catch (libconfig::ParseException err) {
		cout << "Error in config file \"" << inputFile << "\": ";
		cout << err.getError() << " on line " << err.getLine() << endl;
	}

	return 0;
}

void strainedFlame(const std::string& inputFile)
{
    cout << "**** strainedFlame (1dflame Version 2.0) ****\n" << std::endl;

	configOptions mainOptions;
	mainOptions.readOptionsFile(inputFile);

	// OpenMP Configuration
	int nProcs = omp_get_num_procs();
	std::string procStr = (nProcs==1) ? " core." : " cores.";
	cout << "Detected " << nProcs << procStr;
	nProcs = min(mainOptions.numberOfThreads,nProcs);
	procStr = (nProcs==1) ? " core." : " cores.";
	cout << " Running on " << nProcs << procStr << endl;
	omp_set_num_threads(nProcs);

	// Performance Counters

	if (mainOptions.strainRateList.size()!=0) {

		// integral data
		dvector Q, eps, Sc, xFlame;

		for (unsigned int i=0; i<mainOptions.strainRateList.size(); i++) {

			// Instantiate a new solver
			flameSolver theFlameSolver;

			// Set the strain rate and other options
			double a = mainOptions.strainRateList[i];
			std::string restartFile = "prof_eps"+stringify(a,4);
			std::string historyFile = "out_eps"+stringify(a,4);
			if (boost::filesystem::exists(mainOptions.outputDir+"/"+restartFile+".mat") &&
				boost::filesystem::exists(mainOptions.outputDir+"/"+historyFile+".mat"))
			{
				mainOptions.restartFile = mainOptions.outputDir+"/"+restartFile+".mat";
				mainOptions.useRelativeRestartPath = false;
				mainOptions.haveRestartFile = true;
				cout << "Skipping run at strain rate a = " << a << " because the output file \"" <<
					restartFile << "\" already exists." << endl;

				matlabFile oldOutFile(mainOptions.outputDir+"/"+historyFile+".mat");
				dvector oldTime = oldOutFile.readVector("t");
				dvector oldQ = oldOutFile.readVector("Q");
				dvector oldSc = oldOutFile.readVector("Sc");
				dvector oldXflame = oldOutFile.readVector("xFlame");
				double tEnd = *(oldTime.rbegin());

				// Append the integral data from the old run
				int iStart = findFirst(oldTime > tEnd - mainOptions.terminationPeriod);
				int iEnd = oldTime.size()-1;
				eps.push_back(a);
				Q.push_back(mean(oldQ, iStart, iEnd));
				Sc.push_back(mean(oldSc, iStart, iEnd));
				xFlame.push_back(mean(oldXflame, iStart, iEnd));
				oldOutFile.close();
			} else {

				mainOptions.strainRateInitial = a;
				mainOptions.strainRateFinal = a;
				theFlameSolver.setOptions(mainOptions);
				theFlameSolver.calculateReactantMixture();

				// Low Res Run:
				cout << "Begining low-res run at strain rate a = " << a << endl;
	 			theFlameSolver.initialize();
				theFlameSolver.options.idaRelTol = mainOptions.idaRelTolLow;
	 			theFlameSolver.options.terminationTolerance = mainOptions.terminationToleranceLow;
	 			theFlameSolver.options.terminationPeriod = mainOptions.terminationPeriodLow;
	 			theFlameSolver.theSys.tStart = 0;
	 			theFlameSolver.run();
	 			cout << "Completed low-res run at strain rate a = " << a << endl;
	 			double dtLow = theFlameSolver.theSys.tNow - theFlameSolver.theSys.tStart;

				// High Res Run:
				theFlameSolver.theSys.tStart = theFlameSolver.theSys.tNow;
				theFlameSolver.options.idaRelTol = mainOptions.idaRelTol;
				theFlameSolver.options.terminationTolerance = mainOptions.terminationTolerance;
				theFlameSolver.options.terminationPeriod = mainOptions.terminationPeriodHigh + dtLow;
				theFlameSolver.run();
				cout << "Completed high-res run at strain rate a = " << a << endl;

				// copy data needed for the next run
				mainOptions.fileNumberOverride = false;

				theFlameSolver.theSys.writeStateMatFile(restartFile);
				mainOptions.restartFile = mainOptions.outputDir+"/"+restartFile+".mat";
				mainOptions.useRelativeRestartPath = false;
				mainOptions.haveRestartFile = true;

				// Write the time-series data file for this run
				matlabFile outFile(mainOptions.outputDir+"/out_eps"+stringify(a,4)+".mat");
				outFile.writeVector("t",theFlameSolver.timeVector);
				outFile.writeVector("dt",theFlameSolver.timestepVector);
				outFile.writeVector("Q",theFlameSolver.heatReleaseRate);
				outFile.writeVector("Sc",theFlameSolver.consumptionSpeed);
				outFile.writeVector("xFlame",theFlameSolver.flamePosition);
				outFile.close();

				// Append and save the integral data for this run
				int iStart = findFirst(theFlameSolver.timeVector > (theFlameSolver.theSys.tNow - mainOptions.terminationPeriod));
				int iEnd = theFlameSolver.timestepVector.size()-1;
				eps.push_back(a);
				Q.push_back(mean(theFlameSolver.heatReleaseRate, iStart, iEnd));
				Sc.push_back(mean(theFlameSolver.consumptionSpeed, iStart, iEnd));
				xFlame.push_back(mean(theFlameSolver.flamePosition, iStart, iEnd));

			}
			matlabFile dataFile(mainOptions.outputDir+"/integral.mat");
			dataFile.writeVector("a",eps);
			dataFile.writeVector("Q",Q);
			dataFile.writeVector("Sc",Sc);
			dataFile.writeVector("xFlame",xFlame);
			dataFile.close();
		}

	} else {
		flameSolver theFlameSolver;
		theFlameSolver.setOptions(mainOptions);
		theFlameSolver.calculateReactantMixture();
		theFlameSolver.initialize();
		theFlameSolver.run();

		std::string strainString;
		matlabFile outFile(mainOptions.outputDir+"/out.mat");
		outFile.writeVector("t",theFlameSolver.timeVector);
		outFile.writeVector("dt",theFlameSolver.timestepVector);
		outFile.writeVector("Q",theFlameSolver.heatReleaseRate);
		outFile.writeVector("Sc",theFlameSolver.consumptionSpeed);
		outFile.writeVector("xFlame",theFlameSolver.flamePosition);
		outFile.close();
	}
}

void chemistryTest(void)
{
	using namespace Cantera;

	Cantera_CXX::IdealGasMix gas("methane_singlestep.cti");

	gas.setState_TPX(300,101325,"N2:3.76, O2:1.0, CH4:0.5");
	vector<double> wnet(53);
	gas.getNetProductionRates(&wnet[0]);
	cout << wnet << endl;
	return;


//    XML_Node *xc = get_XML_File("gri30.xml");
//    XML_Node * const xs = xc->findNameID("phase", "gri30_mix");
//
//	Cantera::XML_Node* foo = NULL;
//	delete foo;
//
//	int n = 1;
//	clock_t t1, t2;
//	t1 = clock();
//	Cantera::IdealGasPhase thermoBase;
//	Cantera::importPhase(*xs,&thermoBase);
//	vector<Cantera::IdealGasPhase> gas(n,thermoBase);
//	vector<Cantera::GasKinetics*> kin(n);
//	vector<Cantera::Transport*> trans(n);
//	for (int i=0; i<n; i++) {
//		//gas[i] = new Cantera::IdealGasPhase();
//		//Cantera::importPhase(*xs, gas[i]);
//		kin[i] = new Cantera::GasKinetics(&gas[i]);
//
//		kin[i]->init();
//		Cantera::installReactionArrays(*xs,*kin[i],"gri30_mix");
//		kin[i]->finalize();
//
//		trans[i] = Cantera::newTransportMgr("Mix",&gas[i],1,0);
//
//	}
//
//	t2 = clock();
//	cout << "separate: " << t2-t1 << endl;
//
//	t1 = clock();
//	gasArray gas2;
//	gas2.mechanismFile = "gri30.xml";
//	gas2.phaseID = "gri30_mix";
//	gas2.initialize();
//	gas2.resize(n);
//	t2 = clock();
//	cout << "gasArray: " << t2-t1 << endl;
//
//	int nSpec = gas[0].nSpecies();
//	dvector dkm(nSpec);
//
//	dvector y(nSpec);
//	gas[0].setState_TPX(300,101325,"O2:1.0, CH4:0.5");
//	gas[0].getMassFractions(&y[0]);
//	t1 = clock();
//	for (int i=0; i<2000; i++) {
//		y[1] = 0.005*i;
//		gas[0].setMassFractions(&y[0]);
//		//trans[0]->getMixDiffCoeffs(&dkm[0]);
//	}
//
//	t2 = clock();
//	cout << "getMixDiffCoeffs: " << t2-t1 << endl;
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
}

void matlabioTest(void)
{
//	int n = 2, m=20;
//	dvector x(n*m);
//	Array2D foo;
//	for (int i=0; i<n*m; i++) {
//		x[i] = 2*i + 100;
//
//	}
//
//	n = 5, m=20;
//	foo.resize(n,m);
//	for (int i=0; i<n; i++) {
//		for (int j=0; j<m; j++) {
//			foo(i,j) = i+100*j;
//		}
//	}
//	boost::filesystem::remove_all("output");
//	boost::filesystem::create_directory("output");
//
//	matlabFile outFile("output/test.mat");
//	outFile.writeVector("hellofriend",x);
//	outFile.writeArray2D("hugstiem",foo);
//
//
//	dvector hi = outFile.readVector("hellofriend");
//	Array2D bar = outFile.readArray2D("hugstiem");
//	cout << hi << endl;
//
//	outFile.close();
}

void miscTest(void)
{

}
