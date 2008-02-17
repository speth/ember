#include "strainedFlame.h"
#include <libconfig.h++>
#include <iostream>
#include <vector>
#include <cmath>
#include "mathUtils.h"
#include "grid.h"
#include "debugUtils.h"

using namespace mathUtils;
// All Cantera names are in namespace Cantera. You can either
// reference everything as Cantera::<name>, or include the following
// 'using namespace' line.
//using namespace Cantera;

using std::vector; using std::valarray;

int main(int argc, char** argv)
{
	char* python_cmd = getenv("PYTHON_CMD");
	if (!python_cmd) {
		putenv("PYTHON_CMD=python");
	}
    try {
    	strainedFlame();
		//chemistryTest();
    }
	catch (Cantera::CanteraError) {
		Cantera::showErrors(cout);
    }
    
	return 0;
}

void strainedFlame() {

    cout << "**** strainedFlame (1dflame Version 2.0) ****" << std::endl;

    strainedFlameSys theSys;

	theSys.readOptionsFile("input.txt");
	theSys.gas.init();

	// output file:
	ofstream outFile;
	outFile.open("out.m");

	clock_t t1, t2;
	t1 = clock();

	bool newSolver = true; 

	double dt = theSys.tEnd/500;
	double t = 0;
	double dtRegrid = dt*10;
	double tRegrid = dtRegrid;
	int i=1;


	// Initial Conditions for ODE
	theSys.setup();
	theSys.generateInitialProfiles();

	for (int j=0; j<theSys.nPoints; j++) {
		for (int k=0; k<theSys.nSpec; k++) {
			outFile << "Yo( " << k+1 << "," << j+1 << ") = " << theSys.Y(k,j) << ";" << endl;
		}
	}

	theSys.grid.jj = theSys.nPoints;
	theSys.grid.updateDerivedSizes();


	while (t < theSys.tEnd) {

		theSys.setup();

		// Sundials IDA Solver:
		sundialsIDA theSolver(theSys.N);
		theSolver.reltol = 1e-5;
		theSolver.nRoots = 0;
		theSolver.findRoots = false;
		vector<bool> algebraic(theSys.N);

		// Initial condition:
		theSys.rollY(theSolver.y);
		for (int j=0; j<theSys.nPoints; j++) {
			theSolver.abstol(3*j) = 1e-6;
			theSolver.abstol(3*j+1) = 1e-6;
			theSolver.abstol(3*j+2) = 1e-6;
			theSolver.componentId(3*j) = 0.0;
			theSolver.componentId(3*j+1) = 1.0;
			theSolver.componentId(3*j+2) = 1.0;
			algebraic[3*j] = true;
			algebraic[3*j+1] = false;
			algebraic[3*j+2] = false;
		}
		theSolver.componentId(3*theSys.grid.jZero) = 1;
		algebraic[3*theSys.grid.jZero] = false;

		for (int j=0; j<theSys.N; j++) {
			theSolver.ydot(j) = 0;
		}

		theSolver.t0 = t;
		theSys.getInitialCondition(t, theSolver.y, theSolver.ydot, algebraic);

		theSolver.setDAE(&theSys);
		theSolver.calcIC = false;

		if (debugParameters::debugCalcIC) {
			outFile << "xIC{ " << i << "} = [" << theSys.grid.x << "];" << endl;
			outFile << "yIC{ " << i << "} = [" << theSolver.y << "];" << endl;
			outFile << "ydotIC{ " << i << "} = [" << theSolver.ydot << "];" << endl;
		}

		theSolver.initialize();

		int flag;
		theSys.printForMatlab(outFile, theSys.grid.x, i, "x");
		theSys.unrollY(theSolver.y);
		theSys.printForMatlab(outFile, theSys.rhov, i, "rhov");
		theSys.printForMatlab(outFile, theSys.T, i, "T");
		theSys.printForMatlab(outFile, theSys.U, i, "U");
		outFile << "t(" << i << ") = " << t << ";" << endl;
		i++;

		while (t < theSys.tEnd) {
			flag = theSolver.integrateToTime(t+dt);

			if (flag == CV_SUCCESS) {
				t += dt;
			}

			theSys.printForMatlab(outFile, theSys.grid.x, i, "x");
			theSys.unrollY(theSolver.y);
			theSys.printForMatlab(outFile, theSys.rhov, i, "rhov");
			theSys.printForMatlab(outFile, theSys.T, i, "T");
			theSys.printForMatlab(outFile, theSys.U, i, "U");
			outFile << "t(" << i << ") = " << theSolver.tInt << ";" << endl;
			i++;

			if (t > tRegrid) {
				// Adapt the grid if necessary

				tRegrid += dtRegrid;
						theSys.grid.dampVal.resize(theSys.nPoints);

				for (int j=0; j<theSys.nPoints; j++) {
					theSys.grid.dampVal[j] = abs(theSys.mu/theSys.rhov[j]);
				}
				vector<dvector> currentSolution, currentSolutionDot;
				theSys.unrollYdot(theSolver.ydot);
				currentSolution.push_back(theSys.rhov);
				currentSolution.push_back(theSys.U);
				currentSolution.push_back(theSys.T);

				currentSolutionDot.push_back(theSys.drhovdt);
				currentSolutionDot.push_back(theSys.dUdt);
				currentSolutionDot.push_back(theSys.dTdt);

				bool adaptFlag = theSys.grid.adapt(currentSolution, currentSolutionDot);
				bool regridFlag = theSys.grid.regrid(currentSolution, currentSolutionDot);
				

				theSys.rhov = currentSolution[0];
				theSys.U = currentSolution[1];
				theSys.T = currentSolution[2];

				theSys.drhovdt = currentSolutionDot[0];
				theSys.dUdt= currentSolutionDot[1];
				theSys.dTdt = currentSolutionDot[2];
				
				theSys.nPoints = theSys.grid.jj+1;

				if (adaptFlag || regridFlag) {
					break; // exit the inner loop and reinitialize the solver for the new problem size
				}

			}
		}
		theSolver.printStats();
		
	}

	t2 = clock();
	cout << "Runtime: " << ((double)(t2-t1))/CLOCKS_PER_SEC << " seconds." << endl;
	
	// Test of grid adaptation:
	outFile.close();

   	int blargh = 0;
}

void strainedFlameSys::setup(void)
{
	N = 3*nPoints;
	nVars = 3;

	T.resize(nPoints);
	U.resize(nPoints);
	rhov.resize(nPoints);

	dTdt.resize(nPoints,0);
	dUdt.resize(nPoints,0);
	drhovdt.resize(nPoints,0);

	resContinuity.resize(nPoints);
	resMomentum.resize(nPoints);
	resEnergy.resize(nPoints);

	rho.resize(nPoints);
	drhodt.resize(nPoints);

	gas.resize(nPoints);
	nSpec = gas.thermo(0).nSpecies();
	Y.resize(nSpec,nPoints);

	dFdy.resize(N, vector<double> (N));
	dFdydot.resize(N, vector<double> (N));

	if (jacobianIsAllocated) {
		delete bandedJacobian;
	}
	bandedJacobian = new sdBandMatrix(N,nVars+2,nVars+2,2*nVars+4);
	BandZero(bandedJacobian->forSundials());

	pMat.resize(N);
	grid.jj = nPoints-1;
	grid.updateBoundaryIndices();
}

void strainedFlameSys::generateInitialProfiles(void)
{
	grid.x.resize(nPoints);
	int jm = (grid.ju+grid.jb)/2; // midpoint of the profiles.
	int jl = (jm)/3;
	int jr = (5*jm)/3;
	grid.jZero = jm;

	// Reactants
	gas[grid.ju].setState_TPX(Tu,Cantera::OneAtm,reactants);
	rhou = gas[grid.ju].density();

	// Products
	gas[grid.jb].setState_TPX(Tu,Cantera::OneAtm,reactants);
	Cantera::equilibrate(gas[grid.jb],"HP");
	// Tb = gas[grid.jb].temperature();
	
	// Diluent in the center to delay ignition
	gas[jm].setState_TPX(Tu,Cantera::OneAtm,diluent);

	double Tleft = (grid.ju==0) ? Tu : Tb;
	double Tright = (grid.ju==0) ? Tb : Tu;
	T[0] = Tleft; T[grid.jj] = Tright;
	T[jm] = T[grid.ju];

	// Uniform initial grid
	for (int j=0; j<nPoints; j++) {
		grid.x[j] = xLeft + (xRight-xLeft)*((double) j)/((double) nPoints);
	}

	gas[grid.ju].getMassFractions(&Y(0,grid.ju));
	gas[grid.jb].getMassFractions(&Y(0,grid.jb));
	gas[jm].getMassFractions(&Y(0,jm));

	for (int j=1; j<jl; j++) {
		for (int k=0; k<nSpec; k++) {
			Y(k,j) = Y(k,0);
			T[j] = T[0];
		}
	}

	for (int j=jl; j<jm; j++) {
		for (int k=0; k<nSpec; k++) {
			Y(k,j) = Y(k,0) + (Y(k,jm)-Y(k,0))*(grid.x[j]-grid.x[jl])/(grid.x[jm]-grid.x[jl]);
			T[j] = T[0] + (T[jm]-T[0])*(grid.x[j]-grid.x[jl])/(grid.x[jm]-grid.x[jl]);
		}
	}

	for (int j=jm+1; j<jr; j++) {
		for (int k=0; k<nSpec; k++) {
			Y(k,j) = Y(k,jm) + (Y(k,grid.jj)-Y(k,jm))*(grid.x[j]-grid.x[jm])/(grid.x[jr]-grid.x[jm]);
			T[j] = T[jm] + (T[grid.jj]-T[jm])*(grid.x[j]-grid.x[jm])/(grid.x[jr]-grid.x[jm]);
		}
	}

	for (int j=jr; j<nPoints; j++) {
		for (int k=0; k<nSpec; k++) {
			Y(k,j) = Y(k,grid.jj);
			T[j] = T[grid.jj];
		}
	}

	for (int i=0; i<10; i++) {
		mathUtils::smooth(T);
	}

	// Grid and initial profiles of T, U and rhov
	for (int i=0; i<nPoints; i++) {
		grid.x[i] = xLeft + (xRight-xLeft)*((double) i)/((double) nPoints);
		//T[i] = Tleft + (Tright-Tleft)*(grid.x[i]-xLeft)/(xRight-xLeft);
		rho[i] = rhou*Tu/T[i];
		U[i] = sqrt(rho[grid.ju]/rho[i]);
	}

	for (int i=0; i<5; i++) {
		mathUtils::smooth(U);
	}

	rhov[jm] = 0;
	for (int j=jm+1; j<nPoints; j++) {
		rhov[j] = rhov[j-1] - rho[j]*U[j]*strainRate*(grid.x[j]-grid.x[j-1]);
	}

	for (int j=jm-1; j>=0; j--) {
		rhov[j] = rhov[j+1] + rho[j]*U[j]*strainRate*(grid.x[j+1]-grid.x[j]);
	}
}

int strainedFlameSys::f(realtype t, sdVector& y, sdVector& ydot, sdVector& res)
{

	unrollY(y);
	unrollYdot(ydot);

	// Update auxillary data:
	for (int i=0; i<nPoints; i++) {
		rho[i] = rhou*Tu/T[i];
		drhodt[i] = rho[i]/T[i]*dTdt[i];
	}

	int jj = nPoints-1;

	// Left boundary:
	//resContinuity[0] = drhovdt[0];
	resEnergy[0] = dTdt[0];
	resMomentum[0] = dUdt[0];

	// Intermediate points (energy & momentum equations)
	for (int i=1; i<jj; i++) {
		resEnergy[i] = rhov[i]*(T[i-1]*grid.cfm[i] + T[i]*grid.cf[i] + T[i+1]*grid.cfp[i])
			- lambda/cp*(T[i-1]*grid.csm[i] + T[i]*grid.cs[i] + T[i+1]*grid.csp[i])
			+ rho[i]*dTdt[i];

		resMomentum[i] = rhov[i]*(U[i-1]*grid.cfm[i] + U[i]*grid.cf[i] + U[i+1]*grid.cfp[i])
			- mu*(U[i-1]*grid.csm[i] + U[i]*grid.cs[i] + U[i+1]*grid.csp[i])
			+ rho[i]*dUdt[i]
			+ rho[i]*U[i]*U[i]*strainRate
			- rhou*strainRate;
	}

	// Continuity equation
	resContinuity[grid.jZero] = drhovdt[grid.jZero];
	for (int i=grid.jZero+1; i<nPoints; i++) {
		resContinuity[i] = drhodt[i] + rho[i]*U[i]*strainRate
			+ (rhov[i]-rhov[i-1])/grid.hh[i-1];
	}

	for (int i=grid.jZero-1; i>=0; i--) {
		resContinuity[i] = drhodt[i] + rho[i]*U[i]*strainRate
			+ (rhov[i+1]-rhov[i])/grid.hh[i];
	}

	// Right Boundary: fixed values
	resEnergy[jj] = dTdt[jj];
	resMomentum[jj] = dUdt[jj];

	rollResiduals(res);

	//double resNorm = 0;
	//double resC = 0;
	//double resE = 0;
	//double resM = 0;
	//for (int i=0; i<N; i++) {
	//	resNorm += res(i)*res(i);
	//}
	//for (int i=0; i<nPoints; i++) {
	//	resC += resContinuity[i]*resContinuity[i];
	//	resE += resEnergy[i]*resEnergy[i];
	//	resM += resMomentum[i]*resMomentum[i];
	//}
	//resNorm = sqrt(resNorm);
	//resC = sqrt(resC);
	//resE = sqrt(resE);
	//resM = sqrt(resM);
	
	//cout << "resNorm = " << resNorm << endl;
	//cout << "resC = " << resC << endl;
	//cout << "resE = " << resE << endl;
	//cout << "resM = " << resM << endl;

	return 0;
}

int strainedFlameSys::JvProd(realtype t, sdVector& y, sdVector& ydot,
			                 sdVector& res, sdVector& v, sdVector& Jv, realtype c_j)
{
	sdVector resTemp(N);
	sdVector yTemp(N);
	double yNorm = 0;
	double vNorm = 0;
	for (int i=0; i<N; i++) {
		yNorm += y(i)*y(i);
		vNorm += v(i)*v(i);
	}
	yNorm = sqrt(yNorm);
	vNorm = sqrt(vNorm);
	double eps = pow((1+yNorm)*(DBL_EPSILON),0.3333)/vNorm;
	for (int i=0; i<N; i++) {
		yTemp(i) = y(i) + v(i)*eps;
	}
	f(t, yTemp, ydot, resTemp);
	sdVector JvAlt(N);
	for (int i=0; i<N; i++) {
		JvAlt(i) = (resTemp(i) - res(i))/eps;
		Jv(i) = 0;
		for (int j=0; j<N; j++)
			Jv(i) += dFdy[i][j]*v(j) + dFdydot[i][j]*c_j*v(j);
	}

	//for (int i=0; i<N; i++) {
	//	errNorm += pow(JvAlt(i) - Jv(i),2);
	//}
	//errNorm = sqrt(errNorm);
	//cout << "errNorm = " << errNorm << endl;
	return 0;
}

int strainedFlameSys::preconditionerSetup(realtype t, sdVector& y, sdVector& ydot, 
										  sdVector& res, realtype c_j)
{
	double eps = sqrt(DBL_EPSILON)/32;
	sdVector yTemp(N);
	sdVector ydotTemp(N);
	sdVector resTemp(N);
	BandZero(bandedJacobian->forSundials());

	//ofstream outFile;
	//outFile.open("jOut.m");
	//outFile << "y = [" << y << "]; " << endl;
	//outFile << "ydot = [" << ydot << "]; " << endl;
	
	// J = dF/dy
	// Banded, upper & lower bandwidths = nVars. 
	int bw = 2*nVars+1;
	for (int s=0; s<bw; s++)
	{
		for (int i=0; i<N; i++) {
			yTemp(i) = y(i);
		}
		// Form the perturbed y vector
		for (int i=s; i<N; i+=bw) {
			
			if (y(i) != 0) {
				yTemp(i) += yTemp(i)*eps;
			} else {
				yTemp(i) = eps;
			}
		}
		// Evaluate the residuals
		f(t, yTemp, ydot, resTemp);
		
		// Fill in the corresponding Jacobian elements
		int iStart = max(0,s-nVars);
		int counter;
		int col;
		
		if (s >= nVars) {
			counter = 0;
			col = -bw;
			
		} else {
			counter = nVars - s;
			col = 0;
		}
		
		for (int i=iStart; i<N; i++) {
			if (counter++ % bw == 0) {
				col += bw;
			}
			if (col+s>=N) {
				continue;
			}
			(*bandedJacobian)(i,s+col) = (resTemp(i)-res(i))/(yTemp(s+col)-y(s+col));
			//outFile << "J(" << i+1 << "," << s+col+1 << ") = " << (resTemp(i)-res(i))/(yTemp(s+col)-y(s+col)) << ";" << endl;
		}
	}

	// J += c_j*dF/dydot
	// Banded, upper & lower bandwidths = nVars-1. 
	bw = 2*nVars-1;
	for (int s=0; s<bw; s++)
	{
		for (int i=0; i<N; i++) {
			yTemp(i) = y(i);
		}
		// Form the perturbed ydot vector
		for (int i=0; i<N; i++) {
			ydotTemp(i) = ydot(i);
		}
		for (int i=s; i<N; i+=bw) {
			// Form the perturbed y vector
			if (ydot(i) != 0) {
				ydotTemp(i) *= 1+eps;
			} else {
				ydotTemp(i) = eps;
			}
		}

		// Evaluate the residuals
		f(t, y, ydotTemp, resTemp);

		// Fill in the corresponding Jacobian elements
		int iStart = max(0,s-nVars);
		int counter;
		int col;
		if (s >= nVars) {
			counter = 0;
			col = -bw;
			
		} else {
			counter = nVars - s;
			col = 0;
		}
		
		for (int i=iStart; i<N; i++) {
			if (counter++ % bw == 0) {
				col += bw;
			}
			if (col+s>=N) {
				continue;
			}
			(*bandedJacobian)(i,s+col) += c_j*(resTemp(i)-res(i))/(ydotTemp(s+col)-ydot(s+col));
			//outFile << "J2(" << i+1 << "," << s+col+1 << ") = " << (resTemp(i)-res(i))/(ydotTemp(s+col)-ydot(s+col)) << ";" << endl;
		}
	}

	long int iError = BandGBTRF(bandedJacobian->forSundials(),&pMat[0]);

	if (iError!=0) {
		cout << "Error in LU factorization: i = " << iError << endl;
		throw;
	}

//	outFile.close();

 	return 0;
}

int strainedFlameSys::preconditionerSolve(realtype t, sdVector& yIn, sdVector& ydotIn, 
										  sdVector& resIn, sdVector& rhs, 
										  sdVector& outVec, realtype c_j, realtype delta)
{
	vector<double> xVec(N);
	for (int i=0; i<N; i++) {
		xVec[i] = rhs(i);
	}

	BandGBTRS(bandedJacobian->forSundials(),&pMat[0],&xVec[0]);

	for (int i=0; i<N; i++) {
		outVec(i) = xVec[i];
	}

	return 0;
}


strainedFlameSys::strainedFlameSys(void) 
	: jacobianIsAllocated(false)
	, bandedJacobian(NULL)
{
}

strainedFlameSys::~strainedFlameSys(void)
{
	delete bandedJacobian;
}

void strainedFlameSys::unrollYdot(const sdVector& yDot)
{
	for (int i=0; i<nPoints; i++) {
		drhovdt[i] = yDot(3*i);
		dUdt[i] = yDot(3*i+1);
		dTdt[i] = yDot(3*i+2);
	}
}

void strainedFlameSys::rollYdot(sdVector& yDot)
{
	for (int i=0; i<nPoints; i++) {
		yDot(3*i) = drhovdt[i];
		yDot(3*i+1) = dUdt[i];
		yDot(3*i+2) = dTdt[i];
	}
}

void strainedFlameSys::unrollY(const sdVector& y)
{
	for (int i=0; i<nPoints; i++) {
		rhov[i] = y(3*i);
		U[i] = y(3*i+1);
		T[i] = y(3*i+2);
	}
}

void strainedFlameSys::rollY(sdVector& y)
{
	for (int i=0; i<nPoints; i++) {
		y(3*i) = rhov[i];
		y(3*i+1) = U[i];
		y(3*i+2) = T[i];
	}
}

void strainedFlameSys::rollResiduals(sdVector& res)
{
	for (int i=0; i<nPoints; i++) {
		res(3*i) = resContinuity[i];
		res(3*i+1) = resMomentum[i];
		res(3*i+2) = resEnergy[i];
	}
}

void strainedFlameSys::printForMatlab(ofstream& file, vector<double>& v, int index, char* name)
{
	file << name << "{" << index << "} = [";
	for (unsigned int i=0; i<v.size()-1; i++)
	{
		file << v[i] << ", ";
	}
	file << v[v.size()-1] << "];" << endl;
}

void strainedFlameSys::readOptionsFile(const char* filename)
{
	// read input file
	libconfig::Config configFile;
	
	configFile.readFile(filename);
	configFile.setAutoConvert(true);

	bool cfgFlag = 
		configFile.lookupValue("chemistry.mechanismFile",gas.mechanismFile) &&
		configFile.lookupValue("chemistry.phaseID",gas.phaseID) &&
		configFile.lookupValue("chemistry.reactants",reactants) &&
		configFile.lookupValue("chemistry.diluent",diluent) &&
		configFile.lookupValue("grid.nPoints",nPoints) &&
		configFile.lookupValue("grid.xLeft",xLeft) &&
		configFile.lookupValue("grid.xRight",xRight) &&
		configFile.lookupValue("flow.Tu",Tu) &&
		configFile.lookupValue("flow.Tb",Tb) &&
		configFile.lookupValue("flow.rhou",rhou) &&
		configFile.lookupValue("flow.strainRate",strainRate) &&
		configFile.lookupValue("fluid.mu",mu) &&
		configFile.lookupValue("fluid.lambda",lambda) &&
		configFile.lookupValue("fluid.cp",cp) &&
		configFile.lookupValue("tStart",tStart) &&
		configFile.lookupValue("tEnd",tEnd) &&
		configFile.lookupValue("adaptation.vtol",grid.vtol) &&
		configFile.lookupValue("adaptation.dvtol",grid.dvtol) &&
		configFile.lookupValue("adaptation.rmTol",grid.rmTol) &&
		configFile.lookupValue("adaptation.dampConst",grid.dampConst) &&
		configFile.lookupValue("adaptation.gridMin",grid.gridMin) &&
		configFile.lookupValue("adaptation.gridMax",grid.gridMax) &&
		configFile.lookupValue("adaptation.uniformityTol",grid.uniformityTol) &&
		configFile.lookupValue("general.fixedBurnedValFlag",grid.fixedBurnedValFlag) &&
		configFile.lookupValue("general.fixedLeftLocation",grid.fixedLeftLoc) &&
		configFile.lookupValue("general.unburnedLeft",grid.unburnedLeft) &&
		configFile.lookupValue("regridding.boundaryTol",grid.boundaryTol) &&
		configFile.lookupValue("regridding.boundaryTolRm",grid.boundaryTolRm);

	if (!cfgFlag) {
		cout << "Failed to read required settings from input.txt" << endl;
		throw;
	}
	
	grid.alpha = 0;

	grid.kContinuity = 0;
	grid.kMomentum = 1;
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
	gas2.init();
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

void strainedFlameSys::getInitialCondition(double t, sdVector& y, sdVector& ydot, vector<bool>& algebraic)
{
	sdBandMatrix ICmatrix(N,nVars+2,nVars+2,2*nVars+4);
	double eps = 1; // eps might be a bad name for this.
	sdVector yTemp(N);
	sdVector ydotTemp(N);
	sdVector resTemp(N);
	sdVector res(N);

	//ofstream outFile;
	//outFile.open("IC.m");

	BandZero(ICmatrix.forSundials());
	f(t, y, ydot, res);
	
	// ICmatrix = mix of dF/dy and dF/dydot
	// Banded, upper & lower bandwidths = nVars. 
	int bw = 2*nVars+1;
	for (int s=0; s<bw; s++)
	{
		for (int i=0; i<N; i++) {
			yTemp(i) = y(i);
			ydotTemp(i) = ydot(i);
		}
		// Form the perturbed y vector
		for (int i=s; i<N; i+=bw) {
			if (algebraic[i]) {
				if (abs(y(i)) > DBL_EPSILON) {
					yTemp(i) += yTemp(i)*eps;
				} else {
					yTemp(i) = eps;
				}
			} else {
				if (abs(ydot(i)) > DBL_EPSILON) {
					ydotTemp(i) += ydotTemp(i)*eps;
				} else {
					ydotTemp(i) = eps;
				}
			}
		}
		// Evaluate the residuals
		f(t, yTemp, ydotTemp, resTemp);
		
		// Fill in the corresponding Jacobian elements
		int iStart = max(0,s-nVars);
		int counter;
		int col;

		if (s >= nVars) {
			counter = 0;
			col = -bw;
			
		} else {
			counter = nVars - s;
			col = 0;
		}
		
		for (int i=iStart; i<N; i++) {
			if (counter++ % bw == 0) {
				col += bw;
			}
			if (col+s>=N) {
				continue;
			}
			if (algebraic[s+col]) {
				ICmatrix(i,s+col) = (resTemp(i)-res(i))/(yTemp(s+col)-y(s+col));
			} else {
				ICmatrix(i,s+col) = (resTemp(i)-res(i))/(ydotTemp(s+col)-ydot(s+col));
			}
			//outFile << "J(" << i+1 << "," << s+col+1 << ") = " << ICmatrix(i,s+col) << ";" << endl;
		}
	}


	// Solve the linear system to get the IC:

	// LU factorization
	vector<long int> ICpermMat(N);
	long int iError = BandGBTRF(ICmatrix.forSundials(),&ICpermMat[0]);
	if (iError!=0) {
		cout << "iError = " << iError << endl;
		throw;
	}

	// Generate the RHS
	for (int i=0; i<N; i++) {
		if (algebraic[i]) {
			y(i) = 0;
		} else {
			ydot(i) = 0;
		}
	}

	f(t, y, ydot, resTemp);
	for (int i=0; i<N; i++) {
		resTemp(i) = -resTemp(i);
	}

	// Backsubstitute to find the initial conditions:
	BandGBTRS(ICmatrix.forSundials(),&ICpermMat[0],&resTemp(0));

	for (int i=0; i<N; i++) {
		if (algebraic[i]) {
			y(i) = resTemp(i);
		} else {
			ydot(i) = resTemp(i);
		}
	}

	f(t, y, ydot, resTemp);
}