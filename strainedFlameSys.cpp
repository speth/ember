#include "strainedFlameSys.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include "debugUtils.h"
#include "mathUtils.h"
#include "grid.h"
#include "matlabFile.h"
#include "boost/filesystem.hpp"

using namespace mathUtils;

int strainedFlameSys::f(realtype t, sdVector& y, sdVector& ydot, sdVector& res)
{
	tNow = t;
	unrollY(y);
	unrollYdot(ydot);
	double a = strainRate(t);

	// Update auxillary data:
	for (int j=0; j<nPoints; j++) {
		rho[j] = rhou*Tu/T[j];
		drhodt[j] = -rho[j]/T[j]*dTdt[j];
	}

	int jj = nPoints-1;

	// Left boundary:
	//resContinuity[0] = drhovdt[0];
	resEnergy[0] = dTdt[0];
	resMomentum[0] = dUdt[0];
	for (int k=0; k<nSpec; k++) {
		resSpecies(k,0) = dYdt(k,0);
	}
	double diffCoeff = lambda/cp; // Le = 1;

	//// Upwinded spatial derivatives:
	//for (int j=1; j<jj; j++) {
	//	if (rhov[j] > 0) {
	//		dUdx[j] = (U[j]-U[j-1])/grid.hh[j-1];
	//		dTdx[j] = (T[j]-T[j-1])/grid.hh[j-1];
	//		for (int k=0; k<nSpec; k++) {
	//			dYdx(k,j) = (Y(k,j)-Y(k,j-1))/grid.hh[j-1];
	//		}
	//	} else {
	//		dUdx[j] = (U[j+1]-U[j])/grid.hh[j];
	//		dTdx[j] = (T[j+1]-T[j])/grid.hh[j];
	//		for (int k=0; k<nSpec; k++) {
	//			dYdx(k,j) = (Y(k,j+1)-Y(k,j))/grid.hh[j];
	//		}
	//	}
	//}


	// Intermediate points (energy & momentum equations)
	for (int j=1; j<jj; j++) {
		dTdx[j] = (T[j-1]*grid.cfm[j] + T[j]*grid.cf[j] + T[j+1]*grid.cfp[j]);
		dUdx[j] = (U[j-1]*grid.cfm[j] + U[j]*grid.cf[j] + U[j+1]*grid.cfp[j]);

		resEnergy[j] = rhov[j]*dTdx[j]
			- lambda/cp*(T[j-1]*grid.csm[j] + T[j]*grid.cs[j] + T[j+1]*grid.csp[j])
			+ rho[j]*dTdt[j];

		resMomentum[j] = rhov[j]*dUdx[j]
			- mu*(U[j-1]*grid.csm[j] + U[j]*grid.cs[j] + U[j+1]*grid.csp[j])
			+ rho[j]*dUdt[j]
			+ rho[j]*U[j]*U[j]*a
			- rhou*a;

		for (int k=0; k<nSpec; k++) {
			resSpecies(k,j) = rhov[j]*dYdx(k,j)
				- diffCoeff*(Y(k,j-1)*grid.csm[j] + Y(k,j)*grid.cs[j] + Y(k,j+1)*grid.csp[j])
				+ rho[j]*dYdt(k,j);
		}
	}

	// Continuity equation
	resContinuity[grid.jZero] = drhovdt[grid.jZero];
	for (int j=grid.jZero+1; j<nPoints; j++) {
		resContinuity[j] = drhodt[j] + rho[j]*U[j]*a
			+ (rhov[j]-rhov[j-1])/grid.hh[j-1];
	}

	for (int j=grid.jZero-1; j>=0; j--) {
		resContinuity[j] = drhodt[j] + rho[j]*U[j]*a
			+ (rhov[j+1]-rhov[j])/grid.hh[j];
	}

	// Right Boundary: fixed values
	resEnergy[jj] = dTdt[jj];
	resMomentum[jj] = dUdt[jj];
	for (int k=0; k<nSpec; k++) {
		resSpecies(k,jj) = dYdt(k,jj);
	}

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

int strainedFlameSys::preconditionerSetup(realtype t, sdVector& y, sdVector& ydot, 
										  sdVector& res, realtype c_j)
{
	double eps = sqrt(DBL_EPSILON)*10;
	sdVector yTemp(N);
	sdVector ydotTemp(N);
	sdVector resTemp(N);
	BandZero(bandedJacobian->forSundials());

	ofstream outFile;
	if (debugParameters::debugJacobian) {
		outFile.open("jOut.m");
		outFile << "J = sparse( " << N << "," << N << ");" << endl;
		outFile << "J2 = sparse( " << N << "," << N << ");" << endl;
	}

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
			
			if (abs(y(i)) > eps*eps) {
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
			if (debugParameters::debugJacobian && abs(resTemp(i)-res(i)) > DBL_EPSILON) {
				outFile << "J(" << i+1 << "," << s+col+1 << ") = " << (resTemp(i)-res(i))/(yTemp(s+col)-y(s+col)) << ";" << endl;
			}
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
			if (debugParameters::debugJacobian && abs(resTemp(i)-res(i)) > DBL_EPSILON) {
				outFile << "J2(" << i+1 << "," << s+col+1 << ") = " << (resTemp(i)-res(i))/(ydotTemp(s+col)-ydot(s+col)) << ";" << endl;
			}

		}
	}

	long int iError = BandGBTRF(bandedJacobian->forSundials(),&pMat[0]);

	if (iError!=0) {
		cout << "Error in LU factorization: i = " << iError << endl;
		throw;
	}

	if (debugParameters::debugJacobian) {
		outFile.close();
	}

 	return 0;
}

int strainedFlameSys::preconditionerSolve(realtype t, sdVector& yIn, sdVector& ydotIn, 
										  sdVector& resIn, sdVector& rhs, 
										  sdVector& outVec, realtype c_j, realtype delta)
{
	dvector xVec(N);
	for (int i=0; i<N; i++) {
		xVec[i] = rhs(i);
	}

	BandGBTRS(bandedJacobian->forSundials(),&pMat[0],&xVec[0]);

	for (int i=0; i<N; i++) {
		outVec(i) = xVec[i];
	}

	return 0;
}

void strainedFlameSys::setup(void)
{
	int nPointsOld = T.size();
	if (nPoints == nPointsOld) {
		return;
	}

	gas.resize(nPoints);
	nSpec = gas.thermo(0).nSpecies();

	nVars = 3+nSpec;
	N = nVars*nPoints;

	rhov.resize(nPoints);
	U.resize(nPoints);
	T.resize(nPoints);
	Y.resize(nSpec,nPoints);

	drhovdt.resize(nPoints,0);
	dUdt.resize(nPoints,0);
	dTdt.resize(nPoints,0);
	dYdt.resize(nSpec,nPoints);

	dUdx.resize(nPoints,0);
	dTdx.resize(nPoints,0);
	dYdx.resize(nSpec,nPoints,0);

	resContinuity.resize(nPoints);
	resMomentum.resize(nPoints);
	resEnergy.resize(nPoints);
	resSpecies.resize(nSpec,nPoints);

	rho.resize(nPoints);
	drhodt.resize(nPoints);

	delete bandedJacobian;
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

	mu = gas.trans(grid.ju).viscosity();
	lambda = gas.trans(grid.ju).thermalConductivity();
	cp = gas.thermo(grid.ju).cp_mass();

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

	dvector yTemp(nPoints);
	for (int k=0; k<nSpec; k++) {
		for (int j=0; j<nPoints; j++) {
			yTemp[j] = Y(k,j);
		}
		
		for (int i=0; i<10; i++) {
			mathUtils::smooth(yTemp);
		}

		for (int j=0; j<nPoints; j++) {
			Y(k,j) = yTemp[j];
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
		rhov[j] = rhov[j-1] - rho[j]*U[j]*strainRate(tStart)*(grid.x[j]-grid.x[j-1]);
	}

	for (int j=jm-1; j>=0; j--) {
		rhov[j] = rhov[j+1] + rho[j]*U[j]*strainRate(tStart)*(grid.x[j+1]-grid.x[j]);
	}
}

void strainedFlameSys::loadInitialProfiles(void)
{
	std::string inputFilename = options.inputDir + "/" + options.restartFile;
	matlabFile infile(inputFilename);
	grid.x = infile.readVector("x");
	
	nPoints = grid.x.size();
	setup();

	U = infile.readVector("U");
	rhov = infile.readVector("V");
	T = infile.readVector("T");
	Y = infile.readArray2D("Y");

	cout << "Read profiles from restart file: U = " << U << endl;
	cout << "rhov = " << rhov << endl;
	cout << "T = " << T << endl;
	
	tStart = infile.readScalar("t");

	infile.close();

	if (options.overrideTu) {
		T[grid.ju] = Tu;
	} else {
		Tu = T[grid.ju];
	}

	if (options.overrideReactants) {
		gas.thermo(grid.ju).setMoleFractionsByName(reactants);
		dvector yu(nSpec);
		gas.thermo(grid.ju).getMassFractions(&yu[0]);
		for (int k=0; k<nSpec; k++) {
			Y(k,grid.ju) = yu[k];
		}
	} else {
		gas.thermo(grid.ju).setMoleFractions(&Y(0,grid.ju));
	}
	gas.thermo(grid.ju).setState_TP(T[grid.ju],gas.pressure);

	mu = gas.trans(grid.ju).viscosity()*5;


	lambda = gas.trans(grid.ju).thermalConductivity()*5;
	cp = gas.thermo(grid.ju).cp_mass();
	rhou = gas.thermo(grid.ju).density();

	cout << "mu = " << mu << endl;
	cout << "lambda = " << lambda << endl;
	cout << "cp = " << cp << endl;
	cout << "rhou = " << rhou << endl;

	grid.update_jZero(rhov);

}

strainedFlameSys::strainedFlameSys(void) 
	: bandedJacobian(NULL)
	, outputFileNumber(0)
{
}

strainedFlameSys::~strainedFlameSys(void)
{
	delete bandedJacobian;
}

void strainedFlameSys::unrollYdot(const sdVector& yDot)
{
	for (int j=0; j<nPoints; j++) {
		drhovdt[j] = yDot(nVars*j);
		dUdt[j] = yDot(nVars*j+1);
		dTdt[j] = yDot(nVars*j+2);
		for (int k=0; k<nSpec; k++) {
			dYdt(k,j) = yDot(nVars*j+3+k);
		}
	}
}

void strainedFlameSys::rollYdot(sdVector& yDot)
{
	for (int j=0; j<nPoints; j++) {
		yDot(nVars*j) = drhovdt[j];
		yDot(nVars*j+1) = dUdt[j];
		yDot(nVars*j+2) = dTdt[j];
		for (int k=0; k<nSpec; k++) {
			yDot(nVars*j+3+k) = dYdt(k,j);
		}
	}
}

void strainedFlameSys::unrollY(const sdVector& y)
{
	for (int j=0; j<nPoints; j++) {
		rhov[j] = y(nVars*j);
		U[j] = y(nVars*j+1);
		T[j] = y(nVars*j+2);
		for (int k=0; k<nSpec; k++) {
			Y(k,j) = y(nVars*j+3+k);
		}
	}
}

void strainedFlameSys::rollY(sdVector& y)
{
	for (int j=0; j<nPoints; j++) {
		y(nVars*j) = rhov[j];
		y(nVars*j+1) = U[j];
		y(nVars*j+2) = T[j];
		for (int k=0; k<nSpec; k++) {
			y(nVars*j+k+3) = Y(k,j);
		}
	}
}

void strainedFlameSys::rollResiduals(sdVector& res)
{
	for (int j=0; j<nPoints; j++) {
		res(nVars*j) = resContinuity[j];
		res(nVars*j+1) = resMomentum[j];
		res(nVars*j+2) = resEnergy[j];
		for (int k=0; k<nSpec; k++) {
			res(nVars*j+k+3) = resSpecies(k,j);
		}
	}
}

void strainedFlameSys::rollVectorVector(const sdVector& y, vector<dvector>& v)
{
	v.resize(nVars);
	for (int i=0; i<nVars; i++) {
		v[i].resize(nPoints);
		for (int j=0; j<nPoints; j++) {
			v[i][j] = y(i+nVars*j);
		}
	}
}

void strainedFlameSys::unrollVectorVector(const vector<dvector>& v)
{
	for (int j=0; j<nPoints; j++) {
		rhov[j] = v[grid.kContinuity][j];
		U[j] = v[grid.kMomentum][j];
		T[j] = v[grid.kEnergy][j];
	}

	for (int k=0; k<nSpec; k++) {
		for (int j=0; j<nPoints; j++) {
			Y(k,j) = v[grid.kSpecies+k][j];
		}
	}
}

void strainedFlameSys::unrollVectorVectorDot(const vector<dvector>& v)
{
	for (int j=0; j<nPoints; j++) {
		drhovdt[j] = v[grid.kContinuity][j];
		dUdt[j] = v[grid.kMomentum][j];
		dTdt[j] = v[grid.kEnergy][j];
	}

	for (int k=0; k<nSpec; k++) {
		for (int j=0; j<nPoints; j++) {
			dYdt(k,j) = v[grid.kSpecies+k][j];
		}
	}
}

void strainedFlameSys::printForMatlab(ofstream& file, dvector& v, int index, char* name)
{
	file << name << "{" << index << "} = [";
	for (unsigned int i=0; i<v.size()-1; i++)
	{
		file << v[i] << ", ";
	}
	file << v[v.size()-1] << "];" << endl;
}


void strainedFlameSys::getInitialCondition(double t, sdVector& y, sdVector& ydot, std::vector<bool>& algebraic)
{
	sdBandMatrix ICmatrix(N,nVars+2,nVars+2,2*nVars+4);
	double eps = 1; // eps might be a bad name for this.
	sdVector yTemp(N);
	sdVector ydotTemp(N);
	sdVector resTemp(N);
	sdVector res(N);

	ofstream outFile;
	if (debugParameters::debugCalcIC) {
		outFile.open("IC.m");
		outFile << "J = sparse( " << N << "," << N << ");" << endl;
	}
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
			if (debugParameters::debugCalcIC) {
				outFile << "J(" << i+1 << "," << s+col+1 << ") = " << ICmatrix(i,s+col) << ";" << endl;
			}
		}
	}

	if (debugParameters::debugCalcIC) {
		outFile.close();
	}
	// Solve the linear system to get the IC:

	// LU factorization
	std::vector<long int> ICpermMat(N);
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

void strainedFlameSys::writeStateMatFile(void)
{
	// Determine the name of the output file (outXXXXXX.mat)
	std::ostringstream fileName(ostringstream::out);
	fileName << options.outputDir << "/out";
	fileName.flags(ios_base::right);
	fileName.fill('0');
	fileName.width(6);
	fileName << outputFileNumber++ << ".mat";

	cout << "Writing output file: " << fileName.str() << endl;
	// Erase the existing file and create a new one
	if (boost::filesystem::exists(fileName.str())) {
		boost::filesystem::remove(fileName.str());
	}
	matlabFile outFile(fileName.str());

	// Write the system data to the output file:
	outFile.writeScalar("t", tNow);
	outFile.writeVector("x", grid.x);
	outFile.writeVector("T", T);
	outFile.writeVector("U", U);
	outFile.writeVector("V", rhov);
	outFile.writeArray2D("Y", Y);
	outFile.writeVector("rho", rho);

	outFile.close();

}

double strainedFlameSys::strainRate(const double t)
{	
	// Strain rate is at initial value until time strainRateT0,
	// then increases linearly to the final value at time strainRateT0+strainRateDt
	return (t <= strainRateT0) ? strainRateInitial 
		:  (t >= strainRateT0+strainRateDt) ? strainRateFinal
		: strainRateInitial + (strainRateFinal-strainRateInitial)*(t-tStart)/strainRateDt;
}

double strainedFlameSys::dStrainRateDt(const double t)
{
	return (t <= strainRateT0) ? 0
		:  (t >= strainRateT0+strainRateDt) ? 0
		: (strainRateFinal-strainRateInitial)/strainRateDt;
}