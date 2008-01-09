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

    strainedFlameODE theODE;
	theODE.nPoints = 20;
	theODE.xLeft = -0.05;
	theODE.xRight = 0.05;
	theODE.Tleft = 300;
	theODE.Tright = 400;
	theODE.strainRate = 100;
	theODE.mu = 1.8e-5;
	theODE.lambda = 0.0243;
	theODE.rhoLeft = 1.20;
	theODE.tStart = 0;
	theODE.tEnd = 10;
	// Initial Conditions for ODE
	theODE.setup();

	// Sundials CVODE Solver:
	sundialsSolver theSolver(3*theODE.N);
	theSolver.reltol = 1e-5;
	theSolver.nRoots = 0;
	theSolver.findRoots = false;
	theSolver.linearMultistepMethod = CV_BDF;
	theSolver.nonlinearSolverMethod = CV_NEWTON;

	// Initial condition:
	theODE.rollY(theSolver.y0);
	for (int i=0; i<theODE.nPoints; i++) {
		theSolver.abstol(3*i) = 1e-6;
		theSolver.abstol(3*i+1) = 1e-6;
		theSolver.abstol(3*i+2) = 1e-6;
	}
	theSolver.t0 = theODE.tStart;

	theSolver.setODE(&theODE);
	theSolver.initialize();

	// In loop, call CVode and print results
	//   Break out of loop when NOUT preset output times have been reached.

	int iout = 0;  
	double dt = theODE.tEnd/500;
	double t = dt;
	int i=0;
	
	int flag;
	double R = Cantera::GasConstant;
	
	while(t <= theODE.tEnd) {
		
		flag = theSolver.integrateToTime(t);
		//std::cout << theSolver.tInt << ": " << theSolver.y << std::endl;


		if (flag == CV_SUCCESS) {
			i++;
			t += dt;
		}
	}
	theSolver.printStats();
	t2 = clock();
	std::cout << "Runtime: " << ((double)(t2-t1))/CLOCKS_PER_SEC << " seconds." << std::endl;
	
	//outFile.close();
	int blargh = 0;
}

void strainedFlameODE::setup(void)
{
	N = 3*nPoints;
	T.resize(nPoints);
	U.resize(nPoints);
	rhov.resize(nPoints);
	
	dTdt.resize(nPoints,0);
	dUdt.resize(nPoints,0);
	drhovdt.resize(nPoints,0);

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

int strainedFlameODE::f(realtype t, sdVector& y, sdVector& ydot)
{
	unrollY(y);

	// Update auxillary data:
	for (int i=0; i<nPoints; i++) {
		rho[i] = rho[0]*T[0]/T[i];
	}
	int jj = N-1;

	// Left boundary:
	drhovdt[0] = 0;
	dTdt[0] = 0;
	dUdt[0] = 0;

	// Intermediate points (energy & momentum equations)
	for (int i=1; i<jj; i++) {
		dTdt[i] = -rhov[i]*(T[i+1]-T[i-1])/(dx[i-1]+dx[i])
			+ lambda*(T[i+1]-2*T[i]+T[i-1])/(dx[i-1]*dx[i]);

		dUdt[i] = -rhov[i]*(U[i+1]-U[i-1])/(dx[i-1]+dx[i])
			+ mu*(U[i+1]-2*U[i]+U[i-1])/(dx[i-1]*dx[i]);
	}

	// Intermediate points (continuity equation)
	for (int i=1; i<nPoints; i++) {

	}


	// Right Boundary:
	dTdt[jj] = 0;
	dUdt[jj] = (U[jj]-U[jj-1])/dx[jj-1];
	
	rollYdot(ydot);
	return 0;
}

int strainedFlameODE::Jac(realtype t, sdVector& y, sdVector& fy, sdMatrix& J)
{
	double eps = 10*DBL_EPSILON;
		sdVector ydot(N);
		sdVector yTemp(N);
		sdVector ydotTemp(N);
		ydot(0) = 2;
		int foo = f(t,y,ydot);

		for (int j=0; j<N; j++) {

			for (int i=0; i<N; i++) {
				yTemp(i) = y(i);
			}
			double dy;
			if (yTemp(j) != 0) {
				dy = y(j)*(eps);
			} else {
				dy = eps*eps;
			}
			yTemp(j) += dy;
			f(t,yTemp,ydotTemp);

			for (int i=0; i<N; i++) {
				J(i,j) = (ydotTemp(i)-ydot(i))/dy;
			}
		}
	return 0;
}

strainedFlameODE::strainedFlameODE(void) 
{
}

void strainedFlameODE::rollYdot(sdVector& yDot)
{
	for (int i=0; i<nPoints; i++) {
		yDot(3*i) = drhovdt[i];
		yDot(3*i+1) = dUdt[i];
		yDot(3*i+2) = dTdt[i];
	}
}

void strainedFlameODE::unrollY(const sdVector& y)
{
	for (int i=0; i<nPoints; i++) {
		rhov[i] = y(3*i);
		U[i] = y(3*i+1);
		T[i] = y(3*i+2);
	}
}

void strainedFlameODE::rollY(sdVector& y)
{
	for (int i=0; i<nPoints; i++) {
		y(3*i) = rhov[i];
		y(3*i+1) = U[i];
		y(3*i+2) = T[i];
	}
}
