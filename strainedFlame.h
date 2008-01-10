#pragma once
#define _USE_MATH_DEFINES

#include "sundialsUtils.h"

#include <cantera/Cantera.h>
#include <cantera/IdealGasMix.h>    // defines class IdealGasMix
#include <cantera/equilibrium.h>    // chemical equilibrium
#include <cantera/transport.h>      // transport properties

void strainedFlame(void);

class strainedFlameSys : public sdDAE 
{
public:
	strainedFlameSys(void);
	
	// the functions for solving the ODE
	int f(realtype t, sdVector& y, sdVector& ydot, sdVector& res);
	int Jac(realtype t, sdVector& y, sdVector& ydot, sdVector& res,
		    realtype c_j, sdMatrix& J);
	
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
	double cp; // specific heat capacity
	double rhoLeft;
	void setup();

	// Utility functions
	void unrollY(const sdVector& y);
	void unrollYdot(const sdVector& yDot);
	void rollY(sdVector& y);
	void rollYdot(sdVector& yDot);
	void rollResiduals(sdVector& res);
	
	void printForMatlab(ofstream& file, vector<double>& v, int index, char* name);

	// these should be read-only:
	int N; // total problem size;

	// State variables:
	vector<double> rhov;
	vector<double> U;
	vector<double> T;

	// Auxillary variables:
	vector<double> rho;
	vector<double> drhodt;

	// the grid:
	vector<double> x;
	vector<double> dx;

private:
	// Derivatives of state variables:
	vector<double> drhovdt;
	vector<double> dUdt;
	vector<double> dTdt;

	// Residuals of governing equations:
	vector<double> resContinuity;
	vector<double> resMomentum;
	vector<double> resEnergy;
};

