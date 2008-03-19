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
	
	int jj = nPoints-1;

	// Update the thermodynamic state, and evaluate the 
	// thermodynamic, transport and kinetic parameters

	gas.setState(Y,T);
	gas.getMoleFractions(X);

	if (!inJacobianUpdate && !inGetIC) {
		updateTransportProperties();
	}

	if (!inGetIC) {
		updateThermoProperties();
		gas.getReactionRates(wDot);
		for (int j=0; j<=jj; j++) {
			qDot[j] = 0;
			for (int k=0; k<nSpec; k++) {
				qDot[j] += wDot(k,j)*hk(k,j);
			}
		}
	}
	
	// Update auxillary data:
	for (int j=0; j<=jj; j++) {
		double sum = 0;
		for (int k=0; k<nSpec; k++) {
			sum += dYdt(k,j)/W[k];
		}
		drhodt[j] = -rho[j]*(dTdt[j]/T[j] + Wmx[j]*sum);
	}

	double a = strainRate(t);
	double daDt = dStrainRateDt(t);

	// Calculate diffusion mass fluxes, heat flux, enthalpy flux
	for (int j=0; j<jj; j++) {
		sumcpj[j] = 0;
		for (int k=0; k<nSpec; k++) {
			jFick(k,j) = -W[k]*0.5
				* (rho[j]*Dkm(k,j)/Wmx[j] + rho[j+1]*Dkm(k,j+1)/Wmx[j+1])
				* (X(k,j+1)-X(k,j))/grid.hh[j];
			jSoret(k,j) = 0.5*(Dkt(k,j)/T[j] + Dkt(k,j+1)/T[j+1])
				* (T[j+1]-T[j])/grid.hh[j];
			sumcpj[j] += cpSpec(k,j)*(jFick(k,j) + jSoret(k,j));
		}
		qFourier[j] = -0.5*(lambda[j]+lambda[j+1])*(T[j+1]-T[j])/grid.hh[j];
	}

	// Left Boundary values for U, T, Y
	if (grid.leftBoundaryConfig == grid.lbFixedValNonCenter) {
		// Fixed values for T, Y and U
		momentumUnst[0] = dUdt[0];
		energyUnst[0] = dTdt[0];
		energyDiff[0] = energyConv[0] = energyProd[0] = 0;
		momentumDiff[0] = momentumConv[0] = momentumProd[0] = 0;

		for (int k=0; k<nSpec; k++) {
			speciesUnst(k,0) = dYdt(k,0);
			speciesDiff(k,0) = speciesConv(k,0) = speciesProd(k,0) = 0;
		}

	} else if (grid.leftBoundaryConfig == grid.lbZeroGradNonCenter) {
		// Zero gradient condition, not at centerline
		energyDiff[0] = (T[1]-T[0])/grid.hh[0];
		momentumDiff[0] = (U[1]-U[0])/grid.hh[0];
		energyProd[0] = energyConv[0] = energyProd[0] = 0;
		momentumProd[0] = momentumConv[0] = momentumProd[0] = 0;

		for (int k=0; k<nSpec; k++) {
			speciesDiff(k,0) = (Y(k,1)-Y(k,0))/grid.hh[0];
			speciesUnst(k,0) = speciesConv(k,0) = speciesProd(k,0) = 0;
		}

	} else if (grid.leftBoundaryConfig == grid.lbZeroGradCenter) {
		// Zero-gradient condition at domain centerline
		energyUnst[0] = rho[0]*dTdt[0];
		energyDiff[0] = 2*(grid.alpha+1)/(cp[0]*grid.hh[0])*qFourier[0];
		energyProd[0] = qDot[0]/cp[0];
		energyConv[0] = 0;

		momentumUnst[0] = rho[0]*dUdt[0];
		momentumDiff[0] = (grid.alpha+1)*(mu[1]+mu[0])
			/(grid.hh[0]*grid.hh[0])*(U[1]-U[0]);
		momentumProd[0] = daDt*(rho[0]*U[0]-rhou)/a + (rho[0]*U[0]*U[0]-rhou)/a;
		momentumConv[0] = 0;

		for (int k=0; k<nSpec; k++) {
			speciesUnst(k,0) = rho[0]*dYdt(k,0);
			speciesDiff(k,0) = 2*(grid.alpha+1)*(jFick(k,0)+jSoret(k,0))/grid.hh[0];
			speciesProd(k,0) = -wDot(k,0)*W[k];
			speciesConv(k,0) = 0;
		}

	} else if (grid.leftBoundaryConfig == grid.lbControlVolume) {
		double centerVol = pow(grid.x[1],grid.alpha+1)/(grid.alpha+1);
		double centerArea = pow(grid.x[1],grid.alpha+1);

		energyUnst[0] = centerVol*rho[0]*dUdt[0];
		energyDiff[0] = centerArea*qFourier[0]/cp[0];
		energyProd[0] = centerVol*qDot[0]/cp[0];
		energyConv[0] = -rrhov[0]*(Tleft-T[0]);

		momentumUnst[0] = centerVol*rho[0]*dUdt[0];
		momentumDiff[0] = -centerArea*0.5*(mu[0]+mu[1])*(U[1]-U[0])/grid.hh[0];
		momentumProd[0] = centerVol*daDt*(rho[0]*U[0]-rhou)/a + (rho[0]*U[0]*U[0]-rhou)*a;
		momentumConv[0] = -rrhov[0]*(Uleft-U[0]);

		for (int k=0; k<nSpec; k++) {
			speciesUnst(k,0) = centerVol*rho[0]*dYdt(k,0);
			speciesDiff(k,0) = centerArea*(jFick(k,0) + jSoret(k,0));
			speciesProd(k,0) = -centerVol*wDot(k,0)*W[k];
			speciesConv(k,0) = -rrhov[0]*(Yleft[k]-Y(k,0));
		}
	}

	// Intermediate points for U, T, Y
	for (int j=1; j<jj; j++) {
		// First derivative for convective terms: centered difference
		//dTdx[j] = (T[j-1]*grid.cfm[j] + T[j]*grid.cf[j] + T[j+1]*grid.cfp[j]);
		//dUdx[j] = (U[j-1]*grid.cfm[j] + U[j]*grid.cf[j] + U[j+1]*grid.cfp[j]);
		//for (int k=0; k<nSpec; k++) {
		//	dYdx(k,j) = (Y(k,j-1)*grid.cfm[j] + Y(k,j)*grid.cf[j] + Y(k,j+1)*grid.cfp[j]);
		//}

		//// Upwinded convective derivatives:
		if (rhov[j] > 0) {
			dUdx[j] = (U[j]-U[j-1])/grid.hh[j-1];
			dTdx[j] = (T[j]-T[j-1])/grid.hh[j-1];
			for (int k=0; k<nSpec; k++) {
				dYdx(k,j) = (Y(k,j)-Y(k,j-1))/grid.hh[j-1];
			}
		} else {
			dUdx[j] = (U[j+1]-U[j])/grid.hh[j];
			dTdx[j] = (T[j+1]-T[j])/grid.hh[j];
			for (int k=0; k<nSpec; k++) {
				dYdx(k,j) = (Y(k,j+1)-Y(k,j))/grid.hh[j];
			}
		}

		// Momentum Equation
		momentumConv[j] = rhov[j]*dUdx[j];
		momentumDiff[j] = -0.5*( grid.rphalf[j]*(mu[j]+mu[j+1])*(U[j+1]-U[j])/grid.hh[j] -
			grid.rphalf[j-1]*(mu[j-1]+mu[j])*(U[j]-U[j-1])/grid.hh[j-1])/(grid.dlj[j]*grid.r[j]);
		momentumProd[j] = daDt*(rho[j]*U[j]-rhou)/a + (rho[j]*U[j]*U[j]-rhou)*a;
		momentumUnst[j] = rho[j]*dUdt[j];

		// Energy Equation
		energyConv[j] = rhov[j]*dTdx[j];
		energyDiff[j] = sumcpj[j]*(T[j-1]*grid.cfm[j] + T[j]*grid.cf[j] + T[j+1]*grid.cfp[j])/cp[j] +
			(grid.rphalf[j]*qFourier[j] - grid.rphalf[j-1]*qFourier[j-1])/(grid.dlj[j]*cp[j]*grid.r[j]);
		energyProd[j] = qDot[j]/cp[j];
		energyUnst[j] = rho[j]*dTdt[j];

		// Species Equations
		for (int k=0; k<nSpec; k++) {
			speciesUnst(k,j) = rho[j]*dYdt(k,j);
			speciesConv(k,j) = rhov[j]*dYdx(k,j);
			speciesDiff(k,j) = (grid.rphalf[j]*(jFick(k,j)+jSoret(k,j))
				- grid.rphalf[j-1]*(jFick(k,j-1)+jSoret(k,j-1)))/(grid.dlj[j]*grid.r[j]);
			speciesProd(k,j) = -wDot(k,j)*W[k];
		}
	}

	// Right boundary values for U, T, Y
	if (grid.unburnedLeft && !grid.fixedBurnedVal) {
		// zero gradient condition
		energyDiff[jj] = (T[jj]-T[jj-1])/grid.hh[jj-1];
		momentumDiff[jj] = (U[jj]-U[jj-1])/grid.hh[jj-1];
		energyProd[jj] = energyConv[jj] = energyUnst[jj] = 0;
		momentumProd[jj] = momentumConv[jj] = momentumUnst[jj] = 0;

		for (int k=0; k<nSpec; k++) {
			speciesDiff(k,jj) = (Y(k,jj)-Y(k,jj-1))/grid.hh[jj-1];
			speciesProd(k,jj) = speciesConv(k,jj) = speciesUnst(k,jj) = 0;
		}
	} else {
		// Fixed values
		energyUnst[jj] = dTdt[jj];
		momentumUnst[jj] = dUdt[jj];
		energyDiff[jj] = energyProd[jj] = energyConv[jj] = 0;
		momentumDiff[jj] = momentumProd[jj] = momentumConv[jj] = 0;

		for (int k=0; k<nSpec; k++) {
			speciesUnst(k,jj) = dYdt(k,jj);
			speciesProd(k,jj) = speciesConv(k,jj) = speciesDiff(k,jj) = 0;
		}
	}

	// Continuity Equation
	continuityUnst[grid.jZero] = rhov[grid.jZero];
	continuityRhov[grid.jZero] = continuityStrain[grid.jZero] = 0;

	for (int j=grid.jZero+1; j<nPoints; j++) {
		continuityRhov[j] = (rrhov[j]-rrhov[j-1])/(grid.hh[j-1]*grid.rphalf[j-1]);
		continuityUnst[j] = drhodt[j];
		continuityStrain[j] = rho[j]*U[j]*a;
	}

	for (int j=grid.jZero-1; j>=0; j--) {
		continuityRhov[j] = (rrhov[j+1]-rrhov[j])/(grid.hh[j]*grid.rphalf[j]);
		continuityUnst[j] = drhodt[j];
		continuityStrain[j] = rho[j]*U[j]*a;
	}

	rollResiduals(res);

	return 0;
}

int strainedFlameSys::preconditionerSetup(realtype t, sdVector& y, sdVector& ydot, 
										  sdVector& res, realtype c_j)
{
	inJacobianUpdate = true;
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

	inJacobianUpdate = false;
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
		return; // nothing to do
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

	energyUnst.resize(nPoints,0); energyDiff.resize(nPoints,0);
	energyConv.resize(nPoints,0); energyProd.resize(nPoints,0);
	momentumUnst.resize(nPoints,0); momentumDiff.resize(nPoints,0);
	momentumConv.resize(nPoints,0); momentumProd.resize(nPoints,0);
	speciesUnst.resize(nSpec,nPoints,0); speciesDiff.resize(nSpec,nPoints,0);
	speciesConv.resize(nSpec,nPoints,0); speciesProd.resize(nSpec,nPoints,0);
	continuityUnst.resize(nPoints,0); continuityRhov.resize(nPoints,0);
	continuityStrain.resize(nPoints,0);

	rrhov.resize(nPoints);
	rho.resize(nPoints);
	drhodt.resize(nPoints);
	Wmx.resize(nPoints);
	W.resize(nSpec);

	mu.resize(nPoints);
	lambda.resize(nPoints);
	cp.resize(nPoints);
	cpSpec.resize(nSpec,nPoints);
	qFourier.resize(nPoints);

	X.resize(nSpec,nPoints);
	dXdx.resize(nSpec,nPoints);
	Dkm.resize(nSpec,nPoints);
	Dkt.resize(nSpec,nPoints);
	sumcpj.resize(nPoints);
	
	jFick.resize(nSpec,nPoints);
	jSoret.resize(nSpec,nPoints);
	wDot.resize(nSpec,nPoints);
	hk.resize(nSpec,nPoints);
	qDot.resize(nPoints);

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
	Tb = gas[grid.jb].temperature();
	rhob = gas[grid.jb].density();
	
	// Diluent in the center to delay ignition
	gas[jm].setState_TPX(Tu,Cantera::OneAtm,diluent);

	double Tleft = (grid.ju==0) ? Tu : Tb;
	double Tright = (grid.ju==0) ? Tb : Tu;

	double rhoLeft = (grid.ju==0) ? rhou : rhob;
	double rhoRight = (grid.ju==0) ? rhob : rhou;

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
		gas[i].setState_TPY(T[i],gas.pressure,&Y(0,i));
		rho[i] = gas[i].density();
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

	gas.setState(Y,T);
	gas.getMolecularWeights(W);
	
	grid.update_jZero(rhov);
	grid.x -= grid.x[grid.jZero];
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
	}

	gas.setState(Y,T);
	gas.getMolecularWeights(W);
	rhou = gas.thermo(grid.ju).density();

	grid.update_jZero(rhov);

}

strainedFlameSys::strainedFlameSys(void) 
	: bandedJacobian(NULL)
	, outputFileNumber(0)
	, inJacobianUpdate(false)
	, inGetIC(false)
	, grid(options)
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

	if (grid.alpha == 0) {
		for (int j=0; j<nPoints; j++) {
			rrhov[j] = rhov[j];
		} 
	} else {
		for (int j=0; j<nPoints; j++) {
			rrhov[j] = grid.x[j]*rhov[j];
		}
	}

	// In the case where the centerline is part of the domain,
	// rhov[0] actually stores rrhov[0] (since rhov = Inf there)
	if (grid.x[0] == 0) {
		rrhov[0] = rhov[0];
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
		res(nVars*j) = continuityRhov[j] + continuityStrain[j] + continuityUnst[j];
		res(nVars*j+1) = momentumConv[j] + momentumDiff[j] + momentumProd[j] + momentumUnst[j];
		res(nVars*j+2) = energyConv[j] + energyDiff[j] + energyProd[j] + energyUnst[j];
		for (int k=0; k<nSpec; k++) {
			res(nVars*j+k+3) = speciesConv(k,j) + speciesDiff(k,j) + speciesProd(k,j) + speciesUnst(k,j);
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

void strainedFlameSys::updateTransportProperties(void)
{
	gas.getViscosity(mu);
	gas.getThermalConductivity(lambda);
	gas.getDiffusionCoefficients(Dkm);
	gas.getThermalDiffusionCoefficients(Dkt);
}

void strainedFlameSys::updateThermoProperties(void)
{
	gas.getDensity(rho);
	gas.getSpecificHeatCapacity(cp);
	gas.getSpecificHeatCapacities(cpSpec);
	gas.getEnthalpies(hk);
	gas.getMixtureMolecularWeight(Wmx);
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
	
	inGetIC = true; // Turns off property evaluations in f
	
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
	double resnorm = 0;
	for (int i=0; i<N; i++) {
		resnorm += resTemp(i)*resTemp(i);
	}
	resnorm = sqrt(resnorm);
	cout << "Residual norm after IC calculation: " << resnorm << endl;

	inGetIC = false;

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
	outFile.writeVector("q",qDot);
	
	outFile.writeVector("dUdt", dUdt);
	outFile.writeVector("dTdt", dTdt);
	outFile.writeVector("dVdt", drhovdt);
	outFile.writeArray2D("dYdt", dYdt);
	outFile.writeArray2D("wdot", wDot);

	outFile.close();

}

void strainedFlameSys::writeErrorFile(void)
{
	std::string fileName = "errorOutput.mat";
	cout << "Writing error output file: " << fileName << endl;
	// Erase the existing file and create a new one
	if (boost::filesystem::exists(fileName)) {
		boost::filesystem::remove(fileName);
	}
	matlabFile outFile(options.outputDir+"/"+fileName);

	// Write the system data to the output file:
	outFile.writeScalar("t", tNow);
	outFile.writeVector("x", grid.x);
	outFile.writeVector("T", T);
	outFile.writeVector("U", U);
	outFile.writeVector("V", rhov);
	outFile.writeArray2D("Y", Y);

	outFile.writeVector("dUdt", dUdt);
	outFile.writeVector("dTdt", dTdt);
	outFile.writeVector("dVdt", drhovdt);
	outFile.writeArray2D("dYdt", dYdt);
	
	outFile.writeVector("rho", rho);
	outFile.writeVector("drhodt", drhodt);
	outFile.writeArray2D("dkm",Dkm);
	outFile.writeVector("lambda",lambda);
	outFile.writeVector("cp",cp);
	outFile.writeVector("mu",mu);

	outFile.writeScalar("a",strainRate(tNow));

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

void strainedFlameSys::updateAlgebraicComponents(void)
{
	int jj = nPoints-1;
	algebraic.resize(N);

	algebraic[0] = true; // continuity
	if (grid.leftBoundaryConfig == grid.lbFixedValNonCenter) {
		algebraic[1] = false; // momentum
		algebraic[2] = false; // energy
		for (int k=0; k<nSpec; k++) {
			algebraic[3+k] = false; // species
		}
	} else if (grid.leftBoundaryConfig == grid.lbZeroGradNonCenter) {
		algebraic[1] = true; // momentum
		algebraic[2] = true; // energy
		for (int k=0; k<nSpec; k++) {
			algebraic[3+k] = true; // species
		}
	} else if (grid.leftBoundaryConfig == grid.lbZeroGradCenter) {
		algebraic[1] = false; // momentum
		algebraic[2] = false; // energy
		for (int k=0; k<nSpec; k++) {
			algebraic[3+k] = false; // species
		}
	} else if (grid.leftBoundaryConfig == grid.lbControlVolume) {
		algebraic[1] = false; // momentum
		algebraic[2] = false; // energy
		for (int k=0; k<nSpec; k++) {
			algebraic[3+k] = false; // species
		}
	}

	for (int j=1; j<jj; j++) {
		algebraic[nVars*j] = true; // continuity
		algebraic[nVars*j+1] = false; // momentum
		algebraic[nVars*j+2] = false; // energy
		for (int k=0; k<nSpec; k++) {
			algebraic[nVars*j+3+k] = false; //species
		}
	}

	algebraic[nVars*jj] = true; // continuity
	if (grid.unburnedLeft && !grid.fixedBurnedVal) {
		algebraic[nVars*jj+1] = true; // momentum;
		algebraic[nVars*jj+2] = true; // energy
		for (int k=0; k<nSpec; k++) {
			algebraic[nVars*jj+3+k] = true; // species
		}
	} else {
		algebraic[nVars*jj+1] = false; // momentum;
		algebraic[nVars*jj+2] = false; // energy
		for (int k=0; k<nSpec; k++) {
			algebraic[nVars*jj+3+k] = false; // species
		}
	}

}

void strainedFlameSys::updateLeftBC(void)
{
	// Boundary condition for left edge of domain
	if ( (grid.ju == 0 && grid.x[0] != 0) || (grid.jb == 0 && grid.fixedBurnedVal) ) {
		grid.leftBoundaryConfig = grid.lbFixedValNonCenter;
	} else if (grid.x[0] != 0) {
		grid.leftBoundaryConfig = grid.lbZeroGradNonCenter;
	} else if (rrhov[0] <= 0) {
		grid.leftBoundaryConfig = grid.lbZeroGradCenter;
	} else {
		grid.leftBoundaryConfig = grid.lbControlVolume;
	}
}