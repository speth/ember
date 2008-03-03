#pragma once

#include <iostream>
#include "sundialsUtils.h"
#include "chemistry.h"
#include "grid.h"
#include "readConfig.h"

using Cantera::Array2D;
using std::string;

class strainedFlameSys : public sdDAE 
{
public:
	strainedFlameSys(void);
	~strainedFlameSys(void);
	
	// the functions for solving the ODE
	int f(realtype t, sdVector& y, sdVector& ydot, sdVector& res);

	//int JvProd(realtype t, sdVector& yIn, sdVector& ydotIn, sdVector& resIn, 
	//		   sdVector& vIn, sdVector& JvIn, realtype c_j);

	int preconditionerSetup(realtype t, sdVector& yIn, sdVector& ydotIn, 
		                    sdVector& resIn, realtype c_j);

	int preconditionerSolve(realtype t, sdVector& yIn, sdVector& ydotIn, sdVector& resIn,
						    sdVector& rhs, sdVector& outVec, realtype c_j, realtype delta);
	
	// Finds a consistent solution for the DAE to begin the integration
	void getInitialCondition(double t, sdVector& y, sdVector& ydot, vector<bool>& algebraic);

	// Problem definition
	double Tu, Tb;
	double xLeft, xRight;
	int nPoints;

	double strainRateInitial;
	double strainRateFinal;
	double strainRateDt;
	double strainRateT0;

	double tStart;
	double tEnd;
	double tNow;

	double mu; // viscosity
	double lambda; // thermal conductivity
	double cp; // specific heat capacity

	double rhou;
	std::string reactants;
	std::string diluent;

	void setup(void);
	void readOptionsFile(std::string filename);
	
	void generateInitialProfiles(void);
	void loadInitialProfiles(void);
	
	// Utility functions
	void unrollY(const sdVector& y);
	void unrollYdot(const sdVector& yDot);
	void rollY(sdVector& y);
	void rollYdot(sdVector& yDot);
	void rollResiduals(sdVector& res);
	
	void printForMatlab(ofstream& file, vector<double>& v, int index, char* name);
	void writeStateMatFile(void);

	// these should be read-only:
	int N; // total problem size;
	int nVars; // Number of solution variables at each point
	int nSpec; // Number of chemical species

	// State variables:
	vector<double> rhov; // mass flux normal to flame per unit area (rho*v) 
	vector<double> U; // normalized tangential velocity (u/u_inf)
	vector<double> T; // temperature
	Array2D Y; // species mass fractions; Y(k,j)
	
	// Derivatives of state variables:
	vector<double> drhovdt;
	vector<double> dUdt;
	vector<double> dTdt;
	Array2D dYdt;

	// Auxillary variables:
	vector<double> rho; // density [kg/m^3]
	vector<double> drhodt;

	// the grid:
	oneDimGrid grid;

	// Cantera data
	gasArray gas;

	// Miscellaneous configuration options
	configOptions options;

private:
	// Residuals of governing equations:
	vector<double> resContinuity;
	vector<double> resMomentum;
	vector<double> resEnergy;
	Array2D resSpecies;
	
	// Jacobian data
	int jacBW; // Bandwidth of the Jacobian (number of filled blocks per row, 
			   // dependent on the order of the finite difference stencil)
	int jacBWdot; // Bandwidth of dF/dydot component of Jacobian
	sdBandMatrix* bandedJacobian;
	bool jacobianIsAllocated;
	vector<long int> pMat;

	int outputFileNumber; // number of output files written

	double strainRate(const double t);
	double dStrainRateDt(const double t);
};
