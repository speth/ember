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
	~strainedFlameSys(void);
	
	// the functions for solving the ODE
	int f(realtype t, sdVector& y, sdVector& ydot, sdVector& res);

	int JvProd(realtype t, sdVector& yIn, sdVector& ydotIn, sdVector& resIn, 
			   sdVector& vIn, sdVector& JvIn, realtype c_j);

	int preconditionerSetup(realtype t, sdVector& yIn, sdVector& ydotIn, 
		                    sdVector& resIn, realtype c_j);

	int preconditionerSolve(realtype t, sdVector& yIn, sdVector& ydotIn, sdVector& resIn,
						    sdVector& rhs, sdVector& outVec, realtype c_j, realtype delta);
	
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

	// Jacobian data
	vector< vector<double> > dFdy;
	vector< vector<double> > dFdydot;

	int nVars; // Number of solution variables at each point
	int jacBW; // Bandwidth of the Jacobian (number of filled blocks per row, 
			   // dependent on the order of the finite difference stencil)
	int jacBWdot; // Bandwidth of dF/dydot component of Jacobian

	sdMatrix* jacMatrix;
	sdMatrix* jacMatrix2;
	sdMatrix* jacMatrixErr;
	vector<long int> pMat;

	double cjSave;
};

