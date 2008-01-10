#include "strainedFlame.h"
#include <iostream>
#include <vector>

#include <cmath>

typedef vector<double> dvector;

// All Cantera names are in namespace Cantera. You can either
// reference everything as Cantera::<name>, or include the following
// 'using namespace' line.
//using namespace Cantera;

int main(int argc, char** argv) {
		
	char* python_cmd = getenv("PYTHON_CMD");
	if (!python_cmd) {
		putenv("PYTHON_CMD=python");
	}
    try {
    	strainedFlame();
    }
	catch (Cantera::CanteraError) {
		Cantera::showErrors(cout);
    }
    
	return 0;
}

void strainedFlame() {

    cout << "**** strainedFlame (1dflame Version 2.0) ****" << std::endl;

	// output file:
	//ofstream outFile;
	//outFile.open("out.m");
	clock_t t1, t2;
	t1 = clock();

    strainedFlameSys theSys;
	theSys.nPoints = 20;
	theSys.xLeft = -0.05;
	theSys.xRight = 0.05;
	theSys.Tleft = 300;
	theSys.Tright = 300;
	theSys.strainRate = 100;
	theSys.mu = 1.8e-5;
	theSys.lambda = 0.0243;
	theSys.rhoLeft = 1.20;
	theSys.tStart = 0;
	theSys.tEnd = 10;
	// Initial Conditions for ODE
	theSys.setup();

	// Sundials IDA Solver:
	sundialsIDA theSolver(theSys.N);
	theSolver.reltol = 1e-5;
	theSolver.nRoots = 0;
	theSolver.findRoots = false;

	// Initial condition:
	theSys.rollY(theSolver.y);
	for (int i=0; i<theSys.nPoints; i++) {
		theSolver.abstol(3*i) = 1e-6;
		theSolver.abstol(3*i+1) = 1e-6;
		theSolver.abstol(3*i+2) = 1e-6;
	}
	for (int i=0; i<theSys.N; i++) {
		theSolver.ydot(i) = 0;
	}
	theSolver.t0 = theSys.tStart;

	theSolver.setDAE(&theSys);
	theSolver.initialize();
	

	// In loop, call CVode and print results
	//   Break out of loop when NOUT preset output times have been reached.

	int iout = 0;  
	double dt = theSys.tEnd/500;
	double t = dt;
	int i=0;
	
	int flag;
	double R = Cantera::GasConstant;
	
	while(t <= theSys.tEnd) {
		
		flag = theSolver.integrateToTime(t);
		//cout << theSolver.tInt << ": " << theSolver.y << endl;


		if (flag == CV_SUCCESS) {
			i++;
			t += dt;
		}
	}
	theSolver.printStats();
	t2 = clock();
	cout << "Runtime: " << ((double)(t2-t1))/CLOCKS_PER_SEC << " seconds." << endl;
	
	//outFile.close();
	int blargh = 0;
}

void strainedFlameSys::setup(void)
{
	N = 3*nPoints;

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

	x.resize(nPoints);
	dx.resize(nPoints);

	// Grid and initial profiles of T, U and rhov
	for (int i=0; i<nPoints; i++) {
		x[i] = xLeft + (xRight-xLeft)*((double) i)/((double) nPoints);
		T[i] = Tleft + (Tright-Tleft)*(x[i]-xLeft)/(xRight-xLeft);
		rho[i] = Tleft/T[i];
		U[i] = T[i]/Tleft;
	}
	rhov[0] = -strainRate*rho[0]*x[0];
	for (int i=1; i<nPoints; i++) {
		rhov[i] = rhov[i-1] - rho[i]*U[i]*strainRate*(x[i]-x[i-1]);
		dx[i-1] = x[i]-x[i-1];
	}

}

int strainedFlameSys::f(realtype t, sdVector& y, sdVector& ydot, sdVector& res)
{
	unrollY(y);
	unrollYdot(ydot);

	// Update auxillary data:
	for (int i=0; i<nPoints; i++) {
		rho[i] = rho[0]*T[0]/T[i];
		drhodt[i] = rho[i]/T[i]*dTdt[i];
	}

	int jj = nPoints-1;

	// Left boundary:
	resContinuity[0] = drhovdt[0];
	resEnergy[0] = dTdt[0];
	resMomentum[0] = dUdt[0];

	// Intermediate points (energy & momentum equations)
	for (int i=1; i<jj; i++) {
		resEnergy[i] = rhov[i]*(T[i+1]-T[i-1])/(dx[i-1]+dx[i])
			- lambda*(T[i+1]-2*T[i]+T[i-1])/(dx[i-1]*dx[i])
			+ rho[i]*dTdt[i];

		resMomentum[i] = rhov[i]*(U[i+1]-U[i-1])/(dx[i-1]+dx[i])
			- mu*(U[i+1]-2*U[i]+U[i-1])/(dx[i-1]*dx[i])
			+ rho[i]*dUdt[i]
			+ rho[i]*U[i]*U[i]*strainRate
			- rho[1]*strainRate;
	}

	// Intermediate points (continuity equation)
	for (int i=1; i<nPoints; i++) {
		resContinuity[i] = drhodt[i] + rho[i]*U[i]*strainRate
			+ (rhov[i]-rhov[i-1])/dx[i-1];
	}

	// Right Boundary: fixed values
	resEnergy[jj] = dTdt[jj];
	resMomentum[jj] = dUdt[jj];

	rollResiduals(res);
	return 0;
}

// A general, numerical method for finding the jacobian
int strainedFlameSys::Jac(realtype t, sdVector& y, sdVector& ydot,
						  sdVector& res, realtype c_j, sdMatrix& J)
{
	double eps = sqrt(DBL_EPSILON);
//		sdVector ydot(N);
		sdVector yTemp(N);
		sdVector ydotTemp(N);
		sdVector resTemp(N);

		// dF/dy
		for (int j=0; j<N; j++) {

			for (int i=0; i<N; i++) {
				yTemp(i) = y(i);
			}
			double dy;
			if (yTemp(j) != 0) {
				dy = y(j)*(eps);
			} else {
				dy = eps;
			}
			yTemp(j) += dy;
			f(t,yTemp,ydot,resTemp);

			for (int i=0; i<N; i++) {
				J(i,j) = (resTemp(i)-res(i))/dy;
			}
		}

		// cj*dF/dydot
		for (int j=0; j<N; j++) {

			for (int i=0; i<N; i++) {
				ydotTemp(i) = ydot(i);
			}
			double dy;
			if (ydotTemp(j) != 0) {
				dy = ydot(j)*(eps);
			} else {
				dy = eps;
			}
			ydotTemp(j) += dy;
			f(t,y,ydotTemp,resTemp);

			for (int i=0; i<N; i++) {
				J(i,j) += c_j*(resTemp(i)-res(i))/dy;
			}
		}
	return 0;
}

strainedFlameSys::strainedFlameSys(void) 
{
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
