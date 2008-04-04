#include "flameSys.h"
#include <sstream>
#include <vector>
#include <cmath>
#include "debugUtils.h"
#include "mathUtils.h"
#include "matlabFile.h"
#include "boost/filesystem.hpp"

using namespace mathUtils;

int flameSys::f(realtype t, sdVector& y, sdVector& ydot, sdVector& res)
{
	tNow = t;
	unrollY(y);
	unrollYdot(ydot);
	
	int jj = nPoints-1;

	// Update the thermodynamic state, and evaluate the 
	// thermodynamic, transport and kinetic parameters

	gas.setStateMass(Y,T);
	gas.getMoleFractions(X);
	gas.getMassFractions(Y);

	if (!inGetIC) {
		updateTransportProperties();
		updateThermoProperties();
		gas.getReactionRates(wDot);

		for (int j=0; j<=jj; j++) {
			qDot[j] = 0;
			for (int k=0; k<nSpec; k++) {
				qDot[j] -= wDot(k,j)*hk(k,j);
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
	double dadt = dStrainRatedt(t);

	// Calculate diffusion mass fluxes, heat flux, enthalpy flux
	for (int j=0; j<jj; j++) {
		sumcpj[j] = 0;
		for (int k=0; k<nSpec; k++) {
			jFick(k,j) = -W[k]*0.5
				* (rhoD(k,j)/Wmx[j] + rhoD(k,j+1)/Wmx[j+1])
				* (X(k,j+1)-X(k,j))/grid.hh[j];
			jSoret(k,j) = 0.5*(Dkt(k,j)/T[j] + Dkt(k,j+1)/T[j+1])
				* (T[j+1]-T[j])/grid.hh[j];
			sumcpj[j] += cpSpec(k,j)/W[k]*(jFick(k,j) + jSoret(k,j));
		}
		qFourier[j] = -0.5*(lambda[j]+lambda[j+1])*(T[j+1]-T[j])/grid.hh[j];
	}

	// Left Boundary values for U, T, Y
	if (grid.leftBoundaryConfig == grid.lbFixedVal) {
		// Fixed values for T, Y
		
		energyUnst[0] = dTdt[0];
		energyDiff[0] = energyConv[0] = energyProd[0] = 0;
		
		// Momentum equation is always zero gradient on the burned side
		if (grid.unburnedLeft) {
			momentumUnst[0] = dUdt[0];
			momentumDiff[0] = momentumConv[0] = momentumProd[0] = 0;
		} else {
			momentumDiff[0] = (U[1]-U[0])/grid.hh[0];
			momentumProd[0] = momentumConv[0] = momentumProd[0] = 0;
		}

		for (int k=0; k<nSpec; k++) {
			speciesUnst(k,0) = dYdt(k,0);
			speciesDiff(k,0) = speciesConv(k,0) = speciesProd(k,0) = 0;
		}

	} else if (grid.leftBoundaryConfig == grid.lbZeroGradNonCenter) {
		// Zero gradient condition for T, Y, U, not at centerline
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
		energyProd[0] = -qDot[0]/cp[0];
		energyConv[0] = 0;

		momentumUnst[0] = rho[0]*dUdt[0];
		momentumDiff[0] = -(grid.alpha+1)*(mu[1]+mu[0])
			/(grid.hh[0]*grid.hh[0])*(U[1]-U[0]);
		momentumProd[0] = dadt*(rho[0]*U[0]-rhou)/a + (rho[0]*U[0]*U[0]-rhou)*a;
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

		energyUnst[0] = centerVol*rho[0]*dTdt[0];
		energyDiff[0] = centerArea*qFourier[0]/cp[0];
		energyProd[0] = -centerVol*qDot[0]/cp[0];
		energyConv[0] = -rV[0]*(Tleft-T[0]);

		momentumUnst[0] = centerVol*rho[0]*dUdt[0];
		momentumDiff[0] = -centerArea*0.5*(mu[0]+mu[1])*(U[1]-U[0])/grid.hh[0];
		momentumProd[0] = centerVol*dadt*(rho[0]*U[0]-rhou)/a + (rho[0]*U[0]*U[0]-rhou)*a;
		momentumConv[0] = -rV[0]*(Uleft-U[0]);

		for (int k=0; k<nSpec; k++) {
			speciesUnst(k,0) = centerVol*rho[0]*dYdt(k,0);
			speciesDiff(k,0) = centerArea*(jFick(k,0) + jSoret(k,0));
			speciesProd(k,0) = -centerVol*wDot(k,0)*W[k];
			speciesConv(k,0) = -rV[0]*(Yleft[k]-Y(k,0));
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

		dTdxCen[j] = T[j-1]*grid.cfm[j] + T[j]*grid.cf[j] + T[j+1]*grid.cfp[j];
		// Upwinded convective derivatives:
		if (V[j] > 0) {
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
		momentumUnst[j] = rho[j]*dUdt[j];
		momentumConv[j] = V[j]*dUdx[j];
		momentumDiff[j] = -0.5*( grid.rphalf[j]*(mu[j]+mu[j+1])*(U[j+1]-U[j])/grid.hh[j] -
			grid.rphalf[j-1]*(mu[j-1]+mu[j])*(U[j]-U[j-1])/grid.hh[j-1])/(grid.dlj[j]*grid.r[j]);
		momentumProd[j] = dadt*(rho[j]*U[j]-rhou)/a + (rho[j]*U[j]*U[j]-rhou)*a;
		
		// Energy Equation
		energyUnst[j] = rho[j]*dTdt[j];
		energyConv[j] = V[j]*dTdx[j];
		energyDiff[j] = sumcpj[j]*dTdxCen[j]/cp[j] +
			(grid.rphalf[j]*qFourier[j] - grid.rphalf[j-1]*qFourier[j-1])/(grid.dlj[j]*cp[j]*grid.r[j]);
		energyProd[j] = -qDot[j]/cp[j];

		// Species Equations
		for (int k=0; k<nSpec; k++) {
			speciesUnst(k,j) = rho[j]*dYdt(k,j);
			speciesConv(k,j) = V[j]*dYdx(k,j);
			speciesDiff(k,j) = (grid.rphalf[j]*(jFick(k,j)+jSoret(k,j))
				- grid.rphalf[j-1]*(jFick(k,j-1)+jSoret(k,j-1)))/(grid.dlj[j]*grid.r[j]);
			speciesProd(k,j) = -wDot(k,j)*W[k];
		}
	}

	// Right boundary values for T, Y
	if (grid.unburnedLeft && !grid.fixedBurnedVal) {
		// zero gradient condition
		energyDiff[jj] = (T[jj]-T[jj-1])/grid.hh[jj-1];
		energyProd[jj] = energyConv[jj] = energyUnst[jj] = 0;
		
		for (int k=0; k<nSpec; k++) {
			speciesDiff(k,jj) = (Y(k,jj)-Y(k,jj-1))/grid.hh[jj-1];
			speciesProd(k,jj) = speciesConv(k,jj) = speciesUnst(k,jj) = 0;
		}
	} else {
		// Fixed values
		energyUnst[jj] = dTdt[jj];
		energyDiff[jj] = energyProd[jj] = energyConv[jj] = 0;
		
		for (int k=0; k<nSpec; k++) {
			speciesUnst(k,jj) = dYdt(k,jj);
			speciesProd(k,jj) = speciesConv(k,jj) = speciesDiff(k,jj) = 0;
		}
	}

	// Right boundary values for U
	if (grid.unburnedLeft) {
		momentumDiff[jj] = (U[jj]-U[jj-1])/grid.hh[jj-1];
		momentumProd[jj] = momentumConv[jj] = momentumUnst[jj] = 0;
	} else {
		momentumUnst[jj] = dUdt[jj];
		momentumDiff[jj] = momentumProd[jj] = momentumConv[jj] = 0;
	}

	// Continuity Equation

	if (options.flameRadiusControl) {
		rVcenter = rVcenterPrev + (rVcenterNext-rVcenterPrev)/(tFlameNext-tFlamePrev)*(t-tFlamePrev);
		continuityUnst[0] = rV[0] - rVcenter;
	} else if (options.stagnationRadiusControl) {
		continuityUnst[0] = rV[0] - pow(options.rStag,grid.alpha+1)/((grid.alpha+1)*rhoLeft*a);
	} else {
		continuityUnst[0] = drhovdt[0];
	}
	
	continuityRhov[0] = continuityStrain[0] = 0;

	for (int j=1; j<=jj; j++) {
		continuityRhov[j] = (rV[j]-rV[j-1])/(grid.hh[j-1]*grid.rphalf[j-1]);
		continuityUnst[j] = drhodt[j];
		continuityStrain[j] = rho[j]*U[j]*a;
	}

	rollResiduals(res);
	return 0;
}

/*
int flameSys::preconditionerSetup(realtype t, sdVector& y, sdVector& ydot, 
										  sdVector& res, realtype c_j)
{
	inJacobianUpdate = true;
	// The constant "10" here has been empirically determined
	// to give the best performance for typical test cases.
	double eps = sqrt(DBL_EPSILON)*10;

	BandZero(bandedJacobian->forSundials());

	sdVector& yTemp = *yTempJac;
	sdVector& ydotTemp = *ydotTempJac;
	sdVector& resTemp = *resTempJac;

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

	if (debugParameters::debugJacobian) {
		outFile.close();
	}

	long int iError = BandGBTRF(bandedJacobian->forSundials(),&pMat[0]);

	if (iError!=0) {
		cout << "Error in LU factorization: i = " << iError << endl;
		debugParameters::debugJacobian = true;
		return 1;
	} else {
		debugParameters::debugJacobian = false;
	}

	inJacobianUpdate = false;
 	return 0;
}
*/


int flameSys::preconditionerSetup(realtype t, sdVector& y, sdVector& ydot, 
										  sdVector& res, realtype c_j)
{
	int j;
	int jj = nPoints-1;
	// The constant "800" here has been empirically determined to give the 
	// best performance for typical test cases. This value can have 
	// a substantial impact on the convergence rate of the solver.
	double eps = sqrt(DBL_EPSILON)*800;

	BandZero(bandedJacobian->forSundials());

	double a = strainRate(t);
	double dadt = dStrainRatedt(t);
	
	dvector TplusdT(nPoints), UplusdU(nPoints), VplusdV(nPoints);
	
	for (j=0; j<nPoints; j++) {
		TplusdT[j] = T[j]*(1+eps);
		UplusdU[j] = (abs(U[j]) > eps) ? U[j]*(1+eps) : eps;
		VplusdV[j] = (abs(V[j]) > eps) ? V[j]*(1+eps) : eps;
	}
	
	dvector dwhdT(nPoints,0);
	Array2D dwdT(nSpec,nPoints);
	Array2D wDot2(nSpec,nPoints);

	gas.setStateMass(Y,TplusdT);
	gas.getReactionRates(wDot2);

	for (j=0; j<nPoints; j++) {
		for (int k=0; k<nSpec; k++) {
			dwdT(k,j) = (wDot2(k,j)-wDot(k,j))/(TplusdT[j]-T[j]);
			dwhdT[j] += hk(k,j)*dwdT(k,j) + cpSpec(k,j)*wDot(k,j);
		}
	}

	vector<Array2D> dwdY(nPoints, Array2D(nSpec,nSpec));
	Array2D hdwdY(nSpec,nPoints,0);
	Array2D YplusdY(nSpec,nPoints);
	for (int k=0; k<nSpec; k++) {
		YplusdY = Y;
		for (j=0; j<nPoints; j++) {
			YplusdY(k,j) = (abs(Y(k,j)) > eps) ? Y(k,j)*(1+eps) : eps;
		}
		gas.setStateMass(YplusdY,T);
		gas.getReactionRates(wDot2);
		for (j=0; j<nPoints; j++) {
			for (int i=0; i<nSpec; i++) {
				dwdY[j](i,k) = (wDot2(i,j)-wDot(i,j))/(YplusdY(k,j)-Y(k,j));
				hdwdY(k,j) += hk(k,j)*dwdY[j](i,k);
			}
		}
	}

	// J = dF/dy + c_j*dF/dydot
	j=0;
	// Left Boundary values for U, T, Y
	if (grid.leftBoundaryConfig == grid.lbFixedVal) {
		// Fixed values for T, Y
		jacB(j,kEnergy,kEnergy) = c_j;
		
		// Momentum equation is always zero gradient on the burned side
		if (grid.unburnedLeft) {
			jacB(j,kMomentum,kMomentum) = c_j;
		} else {
			jacB(j,kMomentum,kMomentum) = -1/grid.hh[0];
			jacC(j,kMomentum,kMomentum) = 1/grid.hh[0];
		}

		for (int k=0; k<nSpec; k++) {
			jacB(j,kSpecies+k,kSpecies+k) = c_j;
		}

	} else if (grid.leftBoundaryConfig == grid.lbZeroGradNonCenter) {
		// Zero gradient condition for T, Y, U, not at centerline
		jacB(j,kEnergy,kEnergy) = -1/grid.hh[0];
		jacC(j,kEnergy,kEnergy) = 1/grid.hh[0];

		jacB(j,kMomentum,kMomentum) = -1/grid.hh[0];
		jacC(j,kMomentum,kMomentum) = 1/grid.hh[0];

		for (int k=0; k<nSpec; k++) {
			jacB(j,kSpecies+k,kSpecies+k) = -1/grid.hh[0];
			jacC(j,kSpecies+k,kSpecies+k) = 1/grid.hh[0];
		}

	} else if (grid.leftBoundaryConfig == grid.lbZeroGradCenter) {
		// Zero-gradient condition at domain centerline

		// dEnergy/dT_j
		jacB(j,kEnergy,kEnergy) = c_j*rho[j] - rho[j]*dTdt[j]/T[j] +
			(grid.alpha+1)*(lambda[j+1]+lambda[j])/(grid.hh[j]*grid.hh[j]*cp[j]) -
			dwhdT[j]/cp[j];
		// dEnergy/dT_(j+1)
		jacC(j,kEnergy,kEnergy) = - (grid.alpha+1)*(lambda[j+1]+lambda[j])/(grid.hh[j]*grid.hh[j]*cp[j]);

		// dEnergy/dY
		for (int k=0; k<nSpec; k++) {
			jacB(j,kEnergy,kSpecies+k) = -rho[j]*Wmx[j]/W[k]*dTdt[j] -
				hdwdY(k,j)/cp[j];
		}

		// dMomentum/dT
		jacB(j,kMomentum,kEnergy) = -rho[j]/T[j] *
			(dUdt[j] + U[j]*dadt/a + a*U[j]*U[j]);

		// dMomentum/dY
		for (int k=0; k<nSpec; k++) {
			jacB(j,kMomentum,kSpecies+k) = -rho[j]*Wmx[j]/W[k] *
				(dUdt[j] + U[j]*dadt/a + a*U[j]*U[j]);
		}

		// dMomentum/dU_j
		jacB(j,kMomentum,kMomentum) = c_j*rho[j] + dadt/a*rho[j] + 2*rho[j]*U[j]*a 
			+ (grid.alpha+1)*(mu[1]+mu[0])/(grid.hh[0]*grid.hh[0]);

		// dMomentum/dU_(j+1)
		jacC(j,kMomentum,kMomentum) = - (grid.alpha+1)*(mu[1]+mu[0])/(grid.hh[0]*grid.hh[0]);

		for (int k=0; k<nSpec; k++) {
			//dSpecies/dT
			jacB(j,kSpecies+k,kEnergy) = - rho[j]/T[j]*dYdt(k,j)
				- dwdT(k,j)*W[k];

			//dSpecies/dY_j
			for (int i=0; i<nSpec; i++) {
				jacB(j,kSpecies+k,kSpecies+i) = -dwdY[j](k,i)*W[k] - rho[j]*Wmx[j]/W[i]*dYdt(k,j);
			}
			jacB(j,kSpecies+k,kSpecies+k) += c_j*rho[j] +
				(grid.alpha+1)*(rhoD(k,j)+rhoD(k,j+1))*Wmx[j]/Wmx[j+1]/(grid.hh[j]*grid.hh[j]);
			jacC(j,kSpecies+k,kSpecies+k) = - (grid.alpha+1)*(rhoD(k,j)+rhoD(k,j+1))*Wmx[j]/Wmx[j+1]/(grid.hh[j]*grid.hh[j]);
		}

	} else if (grid.leftBoundaryConfig == grid.lbControlVolume) {
		double centerVol = pow(grid.x[1],grid.alpha+1)/(grid.alpha+1);
		double centerArea = pow(grid.x[1],grid.alpha+1);

		// dEnergy/dT_j
		jacB(j,kEnergy,kEnergy) = c_j*centerVol*rho[j] - centerVol*rho[j]/T[j] + 
			centerArea*(grid.alpha+1)*(lambda[j]+lambda[j+1])/(grid.hh[j]*cp[j]) +
			centerVol*dwhdT[j]/cp[j] + rV[0];

		// dEnergy/dT_(j+1)
		jacC(j,kEnergy,kEnergy) = -centerArea*(grid.alpha+1)*(lambda[j]+lambda[j+1])/(grid.hh[j]*cp[j]);

		// dEnergy/dY
		for (int k=0; k<nSpec; k++) {
			jacB(j,kEnergy,kSpecies+k) = hdwdY(k,j)/cp[j] - rho[j]*Wmx[j]/W[k]*dTdt[j];
		}

		// dMomentum/dT
		jacB(j,kMomentum,kEnergy) = -rho[j]/T[j] *
			(centerVol*dUdt[j] - U[j]*dadt/a + U[j]*U[j]*a);

		// dMomentum/dY
		for (int k=0; k<nSpec; k++) {
			jacB(j,kMomentum,kSpecies+k) = -rho[j]*Wmx[j]/W[k] *
				(centerVol*dUdt[j] - U[j]*dadt/a + U[j]*U[j]*a);
		}

		// dMomentum/dU_j
		jacB(j,kMomentum,kMomentum) = centerVol*rho[j]*c_j
			+ centerVol*rho[j]*dadt/a + centerVol*2*rho[j]*U[j]*a
			+ centerArea*(grid.alpha+1)*(mu[j]+mu[j+1])/grid.hh[j] + rV[0];

		// dMomentum/dU_j+1
		jacC(j,kMomentum,kMomentum) = -centerArea*(grid.alpha+1)*(mu[j]+mu[j+1])/grid.hh[j];

		for (int k=0; k<nSpec; k++) {
			// dSpecies/dT
			jacB(j,kSpecies+k,kEnergy) = -centerVol*dwdT(k,j)*W[k] -
				centerVol*rho[j]/T[j]*dYdt(k,j);

			// dSpecies/dY_j
			for (int i=0; i<nSpec; i++) {
				jacB(j,kSpecies+k,kSpecies+i) = -centerVol*dwdY[j](k,i)*W[k] - 
					centerVol*rho[j]*Wmx[j]/W[i]*dYdt(k,j);
			}
			jacB(j,kSpecies+k,kSpecies+k) += centerVol*rho[j]*c_j
				+ centerArea*(grid.alpha+1)*(rhoD(k,j)+rhoD(k,j+1))*Wmx[j]/Wmx[j+1]/grid.hh[j] + rVcenter;

			// dSpecies/dY_(j+1)
			jacC(j,kSpecies+k,kSpecies+k) = -centerArea*(grid.alpha+1)*(rhoD(k,j)+rhoD(k,j+1))*Wmx[j]/Wmx[j+1]/grid.hh[j];
		}
	}

	// Intermediate points for U, T, Y
	for (j=1; j<jj; j++) {

		// Upwinded convective derivatives:
		if (V[j] > 0) {
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
		dTdxCen[j] = T[j-1]*grid.cfm[j] + T[j]*grid.cf[j] + T[j+1]*grid.cfp[j];

		
		for (int k=0; k<nSpec; k++) {
			// dSpecies/dY
			for (int i=0; i<nSpec; i++) {
				jacB(j,kSpecies+k,kSpecies+i) = -dwdY[j](k,i)*W[k] -
					rho[j]*Wmx[j]/W[i]*dYdt(k,j);
			}

			// dSpecies/dT
			jacB(j,kSpecies+k,kEnergy) = -dwdT(k,j)*W[k] - rho[j]/T[j]*dYdt(k,j);

			// dSpecies/dV
			jacB(j,kSpecies+k,kContinuity) = dYdx(k,j);

			// dSpeciesk/dYk_j
			jacB(j,kSpecies+k,kSpecies+k) += rho[j]*c_j  +
				0.5*grid.rphalf[j]*(rhoD(k,j)+rhoD(k,j+1))/(grid.hh[j]*grid.dlj[j]*grid.r[j]) +
				0.5*grid.rphalf[j-1]*(rhoD(k,j-1)+rhoD(k,j))/(grid.hh[j-1]*grid.dlj[j]*grid.r[j]);

			// dSpeciesk/dYk_j-1
			jacA(j,kSpecies+k,kSpecies+k) = -0.5*grid.rphalf[j-1]*(rhoD(k,j-1)+rhoD(k,j))/(grid.hh[j-1]*grid.dlj[j]*grid.r[j]);

			// dSpeciesk/dYk_j+1
			jacC(j,kSpecies+k,kSpecies+k) = -0.5*grid.rphalf[j]*(rhoD(k,j)+rhoD(k,j+1))/(grid.hh[j]*grid.dlj[j]*grid.r[j]);

			if (V[j] > 0) {
				jacB(j,kSpecies+k,kSpecies+k) += V[j]/grid.hh[j-1];
				jacA(j,kSpecies+k,kSpecies+k) -= V[j]/grid.hh[j-1];
			} else {
				jacB(j,kSpecies+k,kSpecies+k) -= V[j]/grid.hh[j];
				jacC(j,kSpecies+k,kSpecies+k) += V[j]/grid.hh[j];
			}
		}

		// dEnergy/dY
		for (int k=0; k<nSpec; k++) {
			jacB(j,kEnergy,kSpecies+k) = hdwdY(k,j)/cp[j] - rho[j]*Wmx[j]/W[k]*dTdt[j];

			// Enthalpy flux term
			if (V[j] > 0) {
				jacB(j,kEnergy,kSpecies+k) += dTdxCen[j]*(rhoD(k,j)/Wmx[j]+rhoD(k,j-1)/Wmx[j-1])/cp[j]/grid.hh[j-1];
				jacA(j,kEnergy,kSpecies+k) -= dTdxCen[j]*(rhoD(k,j)/Wmx[j]+rhoD(k,j-1)/Wmx[j-1])/cp[j]/grid.hh[j-1];
			} else {
				jacB(j,kEnergy,kSpecies+k) -= dTdxCen[j]*(rhoD(k,j)/Wmx[j]+rhoD(k,j+1)/Wmx[j+1])/cp[j]/grid.hh[j];
				jacC(j,kEnergy,kSpecies+k) += dTdxCen[j]*(rhoD(k,j)/Wmx[j]+rhoD(k,j+1)/Wmx[j+1])/cp[j]/grid.hh[j];
			}
		}

		// dEnergy/dT_j
		jacB(j,kEnergy,kEnergy) = rho[j]*c_j - rho[j]/T[j]*dTdt[j] + 
			(0.5*grid.rphalf[j]*(lambda[j]+lambda[j+1])/grid.hh[j] +
			0.5*grid.rphalf[j-1]*(lambda[j-1]+lambda[j])/grid.hh[j-1])/(grid.dlj[j]*cp[j]*grid.r[j]) +
			dwhdT[j]/cp[j] + sumcpj[j]*grid.cf[j]/cp[j];

		// dEnergy/dT_j-1
		jacA(j,kEnergy,kEnergy) = -0.5*grid.rphalf[j-1]*(lambda[j-1]+lambda[j])/(grid.hh[j-1]*grid.dlj[j]*cp[j]*grid.r[j]) +
			sumcpj[j]*grid.cfm[j]/cp[j];

		// dEnergy/dT_j+1
		jacC(j,kEnergy,kEnergy) = -0.5*grid.rphalf[j]*(lambda[j]+lambda[j+1])/(grid.hh[j]*grid.dlj[j]*cp[j]*grid.r[j]) +
			sumcpj[j]*grid.cfp[j]/cp[j];

		if (V[j] > 0) {
			jacB(j,kEnergy,kEnergy) += V[j]/grid.hh[j-1];
			jacA(j,kEnergy,kEnergy) -= V[j]/grid.hh[j-1];
		} else {
			jacB(j,kEnergy,kEnergy) -= V[j]/grid.hh[j];
			jacC(j,kEnergy,kEnergy) += V[j]/grid.hh[j];
		}
		
		// dEnergy/dV
		jacB(j,kEnergy,kContinuity) = dTdx[j];

		// dMomentum/dY
		for (int k=0; k<nSpec; k++) {
			jacB(j,kMomentum,kSpecies+k) = -rho[j]*Wmx[j]/W[k] *
				(dUdt[j] + a*U[j]*U[j] + U[j]*dadt/a);
		}

		// dMomentum/dT
		jacB(j,kMomentum,kEnergy) = -rho[j]/T[j] * 
			(dUdt[j] + a*U[j]*U[j] + U[j]*dadt/a);

		// dMomentum/dU_j
		jacB(j,kMomentum,kMomentum) = rho[j]*c_j +
			0.5*(grid.rphalf[j]*(mu[j]+mu[j+1])/grid.hh[j] +
			grid.rphalf[j]*(mu[j-1]+mu[j])/grid.hh[j-1])/(grid.dlj[j]*grid.r[j]) +
			2*a*U[j]*rho[j] + rho[j]*dadt/a;
		
		// dMomentum/dU_j-1
		jacA(j,kMomentum,kMomentum) = -0.5*grid.rphalf[j]*(mu[j-1]+mu[j])/(grid.hh[j-1]*grid.dlj[j]*grid.r[j]);

		// dMomentum/dU_j+1
		jacC(j,kMomentum,kMomentum) = -0.5*grid.rphalf[j]*(mu[j]+mu[j+1])/(grid.hh[j]*grid.dlj[j]*grid.r[j]);

		if (V[j] > 0) {
			jacB(j,kMomentum,kMomentum) += V[j]/grid.hh[j-1];
			jacA(j,kMomentum,kMomentum) -= V[j]/grid.hh[j-1];
		} else {
			jacB(j,kMomentum,kMomentum) -= V[j]/grid.hh[j];
			jacC(j,kMomentum,kMomentum) += V[j]/grid.hh[j];
		}

		// dMomentum/dV
		jacB(j,kMomentum,kContinuity) = dUdx[j];
	}

	j = jj;
	// Right boundary values for T, Y
	if (grid.unburnedLeft && !grid.fixedBurnedVal) {
		// zero gradient condition
		jacA(j,kEnergy,kEnergy) = -1/grid.hh[j-1];
		jacB(j,kEnergy,kEnergy) = 1/grid.hh[j-1];
		for (int k=0; k<nSpec; k++) {
			jacA(j,kSpecies+k,kSpecies+k) = -1/grid.hh[j-1];
			jacB(j,kSpecies+k,kSpecies+k) = 1/grid.hh[j-1];
		}
	} else {
		// Fixed values
		jacB(j,kEnergy,kEnergy) = c_j;
		for (int k=0; k<nSpec; k++) {
			jacB(j,kSpecies+k,kSpecies+k) = c_j;
		}
	}

	// Right boundary values for U
	if (grid.unburnedLeft) {
		// zero gradient
		jacA(j,kMomentum,kMomentum) = -1/grid.hh[j-1];
		jacB(j,kMomentum,kMomentum) = 1/grid.hh[j-1];
	} else {
		// fixed value
		jacB(j,kMomentum,kMomentum) = c_j;
	}

	// Continuity Equation
	if (options.stagnationRadiusControl || options.flameRadiusControl) {
		jacB(0,kContinuity,kContinuity) = 1;
	} else {
		jacB(0,kContinuity,kContinuity) = c_j;
	}

	for (j=1; j<=jj; j++) {
		for (int k=0; k<nSpec; k++) {
			// dContinuity/dY
			jacB(j,kContinuity,kSpecies+k) = -rho[j]*Wmx[j]/W[k]*(c_j + U[j]*a);
		}

		// dContinuity/dT
		jacB(j,kContinuity,kEnergy) = -rho[j]/T[j]*(c_j + U[j]*a);

		// dContinuity/dU
		jacB(j,kContinuity,kMomentum) = rho[j]*a;

		// dContinuity/dV
		jacB(j,kContinuity,kContinuity) = grid.r[j]/grid.hh[j-1]/grid.rphalf[j-1];

		// dContinuity/dV_(j-1)
		jacA(j,kContinuity,kContinuity) = - grid.r[j]/grid.hh[j-1]/grid.rphalf[j-1];
	}

	long int iError = BandGBTRF(bandedJacobian->forSundials(),&pMat[0]);

	if (iError!=0) {
		cout << "Error in LU factorization: i = " << iError << endl;
		debugParameters::debugJacobian = true;
		return 1;
	} else {
		debugParameters::debugJacobian = false;
	}

	return 0;
}


int flameSys::preconditionerSolve(realtype t, sdVector& yIn, sdVector& ydotIn, 
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

void flameSys::setup(void)
{
	int nPointsOld = T.size();

	if (nPoints == nPointsOld) {
		return; // nothing to do
	}

	gas.resize(nPoints);
	nSpec = gas.thermo(0).nSpecies();

	nVars = 3+nSpec;
	N = nVars*nPoints;

	V.resize(nPoints);
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
	dTdxCen.resize(nPoints,0);

	energyUnst.resize(nPoints,0); energyDiff.resize(nPoints,0);
	energyConv.resize(nPoints,0); energyProd.resize(nPoints,0);
	momentumUnst.resize(nPoints,0); momentumDiff.resize(nPoints,0);
	momentumConv.resize(nPoints,0); momentumProd.resize(nPoints,0);
	speciesUnst.resize(nSpec,nPoints,0); speciesDiff.resize(nSpec,nPoints,0);
	speciesConv.resize(nSpec,nPoints,0); speciesProd.resize(nSpec,nPoints,0);
	continuityUnst.resize(nPoints,0); continuityRhov.resize(nPoints,0);
	continuityStrain.resize(nPoints,0);

	rV.resize(nPoints);
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
	rhoD.resize(nSpec,nPoints);
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

void flameSys::copyOptions(void)
{
	strainRateInitial = options.strainRateInitial;
	strainRateFinal = options.strainRateFinal;
	strainRateT0 = options.strainRateT0;
	strainRateDt = options.strainRateDt;

	tStart = options.tStart;
	tEnd = options.tEnd;

	gas.mechanismFile = options.gasMechanismFile;
	gas.phaseID = options.gasPhaseID;
	gas.pressure = options.pressure;

	nPoints = options.nPoints;
	grid.unburnedLeft = options.unburnedLeft;

	grid.vtol = options.vtol;
	grid.dvtol = options.dvtol;
	grid.absvtol = options.absvtol;
	grid.rmTol = options.rmTol;
	grid.uniformityTol = options.uniformityTol;
	grid.gridMin = options.gridMin;
	grid.gridMax = options.gridMax;
	grid.dampConst = options.dampConst;
	grid.fixedBurnedVal = options.fixedBurnedVal;
	grid.unburnedLeft = options.unburnedLeft;
	grid.boundaryTol = options.boundaryTol;
	grid.boundaryTolRm = options.boundaryTolRm;
	grid.addPointCount = options.addPointCount;

	grid.alpha = options.gridAlpha;
	grid.kMomentum = options.kMomentum;
	grid.kContinuity = options.kContinuity;
	grid.kEnergy = options.kEnergy;
	grid.kSpecies = options.kSpecies;

	kMomentum = options.kMomentum;
	kContinuity = options.kContinuity;
	kEnergy = options.kEnergy;
	kSpecies = options.kSpecies;
}

void flameSys::generateInitialProfiles(void)
{
	grid.x.resize(nPoints);
	int jm = (grid.ju+grid.jb)/2; // midpoint of the profiles.
	int jl = jm - 4;
	int jr = jm + 4;
	//int jl = (jm)/3;
	//int jr = (5*jm)/3;
	grid.jZero = jm;
	Yb.resize(nSpec); Yu.resize(nSpec);

	// Reactants
	Tu = options.Tu;
	gas[grid.ju].setState_TPX(Tu,Cantera::OneAtm,&options.reactants[0]);
	rhou = gas[grid.ju].density();

	// Products
	Tb = options.Tb;
	gas[grid.jb].setState_TPX(Tu,Cantera::OneAtm,&options.reactants[0]);
	Cantera::equilibrate(gas[grid.jb],"HP");
	Tb = gas[grid.jb].temperature();
	rhob = gas[grid.jb].density();
	gas[grid.jb].getMassFractions(&Yb[0]);
	
	// Diluent in the center to delay ignition
	gas[jm].setState_TPY(Tu,Cantera::OneAtm,options.oxidizer);

	if (options.unburnedLeft) {
		Yleft = Yu;
		Uleft = 1;
		rhoRight = rhob;
		Tright = Tb;
		Yright = Yb;
		Uright = Ub;
	} else {
		rhoLeft = rhob;
		Tleft = Tb;
		Yleft = Yb;
		Uleft = Ub;
		rhoRight = rhou;
		Tright = Tu;
		Yright = Yu;
		Uright = 1;
	}

	T[0] = Tleft; T[grid.jj] = Tright;
	T[jm] = T[grid.ju];



	xLeft = (options.twinFlame || options.curvedFlame) ?
		max(options.xLeft,0.0) : options.xLeft;
	xRight = options.xRight;

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

	for (int i=0; i<5; i++) {
		mathUtils::smooth(T);
	}

	// Grid and initial profiles of T, U and V
	for (int i=0; i<nPoints; i++) {
		gas[i].setState_TPY(T[i],gas.pressure,&Y(0,i));
		rho[i] = gas[i].density();
		U[i] = sqrt(rhou/rho[i]);
	}

	for (int i=0; i<2; i++) {
		mathUtils::smooth(U);
	}

	if (options.fixedLeftLoc)
	{
		jm = 0;
	}

	V[jm] = 0;
	for (int j=jm+1; j<nPoints; j++) {
		V[j] = V[j-1] - rho[j]*U[j]*strainRate(tStart)*(grid.x[j]-grid.x[j-1]);
	}

	for (int j=jm-1; j>=0; j--) {
		V[j] = V[j+1] + rho[j]*U[j]*strainRate(tStart)*(grid.x[j+1]-grid.x[j]);
	}

	gas.setStateMass(Y,T);
	gas.getMoleFractions(X);
	gas.getMolecularWeights(W);
	
	grid.update_jZero(V);
	grid.x -= grid.x[grid.jZero];
}

void flameSys::loadInitialProfiles(void)
{

	std::string inputFilename;
	if (options.useRelativeRestartPath) {
		inputFilename = options.inputDir + "/" + options.restartFile;
	} else {
		inputFilename = options.restartFile;
	}

	matlabFile infile(inputFilename);
	grid.x = infile.readVector("x");
	
	nPoints = grid.x.size();
	setup();

	U = infile.readVector("U");
	V = infile.readVector("V");
	T = infile.readVector("T");
	Y = infile.readArray2D("Y");
	tStart = infile.readScalar("t");
	if (!options.fileNumberOverride) {
		options.outputFileNumber = (int) infile.readScalar("fileNumber");
	}

	infile.close();

	if (options.overrideTu) {
		T[grid.ju] = options.Tu;
	}
	Tu = T[grid.ju];

	if (options.overrideReactants) {
		gas.thermo(grid.ju).setState_TPX(Tu,gas.pressure,&options.reactants[0]);
		gas.thermo(grid.ju).getMassFractions(&Y(0,grid.ju));
		gas.thermo(grid.jb).setState_TPX(Tu,gas.pressure,&options.reactants[0]);
		Cantera::equilibrate(gas.thermo(grid.jb),"HP");
		gas.thermo(grid.jb).getMassFractions(&Y(0,grid.jb));
		T[grid.jb] = gas.thermo(grid.jb).temperature();
	}
	Tb = T[grid.jb];

	gas.setStateMass(Y,T);
	gas.getMolecularWeights(W);
	rhou = gas.thermo(grid.ju).density();
	rhob = gas.thermo(grid.jb).density();
	Yu.resize(nSpec); Yb.resize(nSpec);
	gas[grid.ju].getMassFractions(&Yu[0]);
	gas[grid.jb].getMassFractions(&Yb[0]);
	Ub = sqrt(rhou/rhob);

	if (options.unburnedLeft) {
		rhoLeft = rhou;
		Tleft = Tu;
		Yleft = Yu;
		Uleft = 1;
		rhoRight = rhob;
		Tright = Tb;
		Yright = Yb;
		Uright = Ub;

	} else {
		rhoLeft = rhob;
		Tleft = Tb;
		Yleft = Yb;
		Uleft = Ub;
		rhoRight = rhou;
		Tright = Tu;
		Yright = Yu;
		Uright = 1;
	}

	grid.update_jZero(V);
}

flameSys::flameSys(void) 
	: grid(options)
	, bandedJacobian(NULL)
	, flamePosIntegralError(0)
	, inGetIC(false)
{
}

flameSys::~flameSys(void)
{
	delete bandedJacobian;
}

void flameSys::unrollYdot(const sdVector& yDot)
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

void flameSys::rollYdot(sdVector& yDot)
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

void flameSys::unrollY(const sdVector& y)
{
	for (int j=0; j<nPoints; j++) {
		V[j] = y(nVars*j);
		U[j] = y(nVars*j+1);
		T[j] = y(nVars*j+2);
		for (int k=0; k<nSpec; k++) {
			Y(k,j) = y(nVars*j+3+k);
		}
	}

	if (grid.alpha == 0) {
		for (int j=0; j<nPoints; j++) {
			rV[j] = V[j];
		} 
	} else {
		for (int j=0; j<nPoints; j++) {
			rV[j] = grid.x[j]*V[j];
		}
	}

	// In the case where the centerline is part of the domain,
	// V[0] actually stores rV[0] (since V = Inf there)
	if (grid.x[0] == 0) {
		rV[0] = V[0];
	}
}

void flameSys::rollY(sdVector& y)
{
	for (int j=0; j<nPoints; j++) {
		y(nVars*j) = V[j];
		y(nVars*j+1) = U[j];
		y(nVars*j+2) = T[j];
		for (int k=0; k<nSpec; k++) {
			y(nVars*j+k+3) = Y(k,j);
		}
	}
}

void flameSys::rollResiduals(sdVector& res)
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

void flameSys::rollVectorVector(const sdVector& y, vector<dvector>& v)
{
	v.resize(nVars);
	for (int i=0; i<nVars; i++) {
		v[i].resize(nPoints);
		for (int j=0; j<nPoints; j++) {
			v[i][j] = y(i+nVars*j);
		}
	}
}

void flameSys::unrollVectorVector(const vector<dvector>& v)
{
	for (int j=0; j<nPoints; j++) {
		V[j] = v[grid.kContinuity][j];
		U[j] = v[grid.kMomentum][j];
		T[j] = v[grid.kEnergy][j];
	}

	for (int k=0; k<nSpec; k++) {
		for (int j=0; j<nPoints; j++) {
			Y(k,j) = v[grid.kSpecies+k][j];
		}
	}
}

void flameSys::unrollVectorVectorDot(const vector<dvector>& v)
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

void flameSys::updateTransportProperties(void)
{
	gas.getViscosity(mu);
	gas.getThermalConductivity(lambda);
	gas.getWeightedDiffusionCoefficients(rhoD);
	gas.getThermalDiffusionCoefficients(Dkt);
}

void flameSys::updateThermoProperties(void)
{
	gas.getDensity(rho);
	gas.getSpecificHeatCapacity(cp);
	gas.getSpecificHeatCapacities(cpSpec);
	gas.getEnthalpies(hk);
	gas.getMixtureMolecularWeight(Wmx);
}

void flameSys::printForMatlab(ofstream& file, dvector& v, int index, char* name)
{
	file << name << "{" << index << "} = [";
	for (unsigned int i=0; i<v.size()-1; i++)
	{
		file << v[i] << ", ";
	}
	file << v[v.size()-1] << "];" << endl;
}

int flameSys::getInitialCondition(double t, sdVector& y, sdVector& ydot, std::vector<bool>& algebraic)
{
	
	sdBandMatrix ICmatrix(N,nVars+2,nVars+2,2*nVars+4);
	double eps = 0.1; // eps might be a bad name for this.
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
		debugParameters::debugCalcIC = true;
		return 1;
	} else {
		debugParameters::debugCalcIC = false;
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
	if (resnorm > 1) {

		matlabFile outfile("ICmiss.mat");
		dvector resVec(N);
		dvector yVec(N);
		dvector ydotVec(N);
		for (int i=0; i<N; i++) {
			resVec[i] = resTemp(i);
			yVec[i] = y(i);
			ydotVec[i] = ydot(i);
		}
		outfile.writeVector("y",yVec);
		outfile.writeVector("ydot",ydotVec);
		outfile.writeVector("res",resVec);
		outfile.writeVector("eUnst",energyUnst);
		outfile.writeVector("eDiff",energyDiff);
		outfile.writeVector("eConv",energyConv);
		outfile.writeVector("eProd",energyProd);

		outfile.close();
		inGetIC = false;
		return 1;
	}


	inGetIC = false;
	return 0;
}

void flameSys::writeStateMatFile(const std::string fileNameStr, bool errorFile)
{
	std::ostringstream fileName(ostringstream::out);
	bool incrementFileNumber = false;

	if (fileNameStr.length() == 0) {
	// Determine the name of the output file (outXXXXXX.mat)
		incrementFileNumber = true;
		fileName << options.outputDir << "/prof";
		fileName.flags(ios_base::right);
		fileName.fill('0');
		fileName.width(6);
		fileName << options.outputFileNumber << ".mat";
	} else {
		fileName << options.outputDir << "/" << fileNameStr << ".mat";
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
	matlabFile outFile(fileName.str());

	// Write the state data to the output file:
	outFile.writeScalar("t", tNow);
	outFile.writeVector("x", grid.x);
	outFile.writeVector("T", T);
	outFile.writeVector("U", U);
	outFile.writeVector("V", V);
	outFile.writeArray2D("Y", Y);	
	outFile.writeScalar("a",strainRate(tNow));
	outFile.writeScalar("dadt",dStrainRatedt(tNow));
	outFile.writeScalar("fileNumber", options.outputFileNumber);

	if (options.outputHeatReleaseRate || errorFile) {
		outFile.writeVector("q",qDot);
	}

	if (options.outputTimeDerivatives || errorFile) {
		outFile.writeVector("dUdt", dUdt);
		outFile.writeVector("dTdt", dTdt);
		outFile.writeVector("dVdt", drhovdt);
		outFile.writeArray2D("dYdt", dYdt);
		outFile.writeVector("drhodt",drhodt);
	}

	if (options.outputAuxiliaryVariables || errorFile) {
		outFile.writeVector("rho", rho);
		outFile.writeArray2D("wdot", wDot);
		outFile.writeArray2D("rhoD",rhoD);
		outFile.writeVector("lambda",lambda);
		outFile.writeVector("cp",cp);
		outFile.writeVector("mu",mu);
	}

	outFile.close();
	if (incrementFileNumber) {
		options.outputFileNumber++;
	}
}

double flameSys::strainRate(const double t)
{	
	// Strain rate is at initial value until time strainRateT0,
	// then increases linearly to the final value at time strainRateT0+strainRateDt
	return (t <= strainRateT0) ? strainRateInitial 
		:  (t >= strainRateT0+strainRateDt) ? strainRateFinal
		: strainRateInitial + (strainRateFinal-strainRateInitial)*(t-tStart)/strainRateDt;
}

double flameSys::dStrainRatedt(const double t)
{
	return (t <= strainRateT0) ? 0
		:  (t >= strainRateT0+strainRateDt) ? 0
		: (strainRateFinal-strainRateInitial)/strainRateDt;
}

void flameSys::updateAlgebraicComponents(void)
{
	int jj = nPoints-1;
	algebraic.resize(N);
	
	algebraic[0] = options.stagnationRadiusControl 
		|| options.flameRadiusControl; // continuity

	if (grid.leftBoundaryConfig == grid.lbFixedVal) {
		algebraic[1] = !grid.unburnedLeft; // momentum
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

	algebraic[nVars*jj+1] = grid.unburnedLeft;

}

void flameSys::updateLeftBC(void)
{
	// Boundary condition for left edge of domain
	if ( (grid.ju == 0 && grid.x[0] != 0) || (grid.jb == 0 && grid.fixedBurnedVal) ) {
		grid.leftBoundaryConfig = grid.lbFixedVal;
	} else if (grid.x[0] != 0) {
		grid.leftBoundaryConfig = grid.lbZeroGradNonCenter;
	} else if (rV[0] <= 0) {
		grid.leftBoundaryConfig = grid.lbZeroGradCenter;
	} else {
		grid.leftBoundaryConfig = grid.lbControlVolume;
	}
}

double& flameSys::jacA(const int j, const int k1, const int k2)
{
	return (*bandedJacobian)(j*nVars+k1, (j-1)*nVars+k2);
}

double& flameSys::jacB(const int j, const int k1, const int k2)
{
	return (*bandedJacobian)(j*nVars+k1, j*nVars+k2);
}

double& flameSys::jacC(const int j, const int k1, const int k2)
{
	return (*bandedJacobian)(j*nVars+k1, (j+1)*nVars+k2);
}


double flameSys::getHeatReleaseRate(void)
{
	return mathUtils::integrate(grid.x, qDot);
}

double flameSys::getConsumptionSpeed(void)
{
	double QoverCp = mathUtils::integrate(grid.x,qDot/cp);
	double rhouDeltaT = rho[grid.ju]*(T[grid.jb]-T[grid.ju]);
	return QoverCp/rhouDeltaT;
}

double flameSys::getFlamePosition(void)
{
	return mathUtils::trapz(grid.x,grid.x*qDot)/mathUtils::trapz(grid.x,qDot);
}

double flameSys::targetFlamePosition(double t)
{
	return (t <= options.rFlameT0) ? options.rFlameInitial 
		:  (t >= options.rFlameT0+options.rFlameDt) ? options.rFlameFinal
		: options.rFlameInitial + (options.rFlameFinal-options.rFlameInitial)*(t-options.rFlameT0)/options.rFlameDt;
}

void flameSys::update_rVcenter(const double t)
{
	rFlameActual = getFlamePosition();
	rFlameTarget = targetFlamePosition(t);
	flamePosIntegralError += (rFlameTarget-rFlameActual)*(t-tFlamePrev);
	double flamePosDerivativeError = (t > tFlamePrev) ? -(rFlameActual-rFlamePrev)/(t-tFlamePrev) : 0;

	rVcenterPrev = rVcenterPrev + (rVcenterNext-rVcenterPrev)/(tFlameNext-tFlamePrev)*(t-tFlamePrev);

	rFlamePrev = rFlameActual;
	tFlamePrev = t;
	tFlameNext = t + options.rFlameUpdateTimeInterval;
	double a = strainRate(t);

	rVcenterNext = rVcenterInitial*a/options.strainRateInitial + 
		rhou*a*options.rFlameProportionalGain*( (rFlameTarget-rFlameActual) +
		flamePosDerivativeError*options.rFlameDerivativeGain +
		flamePosIntegralError*options.rFlameIntegralGain);
}