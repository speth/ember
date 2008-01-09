#pragma once
#define _USE_MATH_DEFINES

#include "sundialsUtils.h"

#include <cantera/Cantera.h>
#include <cantera/IdealGasMix.h>    // defines class IdealGasMix
#include <cantera/equilibrium.h>    // chemical equilibrium
#include <cantera/transport.h>      // transport properties

void strainedFlame(void);

class strainedFlameODE : public sdODE 
{
public:
	strainedFlameODE(void);
	
	// the functions for solving the ODE
	int f(realtype t, sdVector& y, sdVector& ydot);
	int Jac(realtype t, sdVector& y, sdVector& fy, sdMatrix& J);
	
	// Chemistry:
	//Cantera::IdealGasMix gas;

	// Problem definition
	double Tleft, Tright;
	double xLeft, xRight;
	int nPoints;
	double strainRate;
	double tStart;
	double tEnd;
	double mu; // viscosity
	double lambda; // thermal conductivity
	double rhoLeft;
	void setup();

	// Utility functions
	void unrollY(const sdVector& y);
	void rollYdot(sdVector& yDot);
	void rollY(sdVector& y);

	// these should be read-only:
	int N; // total problem size;

	// State variables:
	vector<double> U;
	vector<double> rhov;
	vector<double> T;

	// Derivatives of state variables:
	vector<double> dUdt;
	vector<double> drhovdt;
	vector<double> dTdt;

	// extra variables:
	vector<double> rho;

	// the grid:
	vector<double> x;
	vector<double> dx;

private:


};

