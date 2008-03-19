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
	std::string reactants;
	std::string diluent;

	double xLeft, xRight;
	int nPoints;

	double tStart;
	double tEnd;
	double tNow;

	// Boundary values
	double rhou, rhob, rhoLeft, rhoRight;
	double Tu, Tb, Tleft, Tright;
	double Ub, Uleft, Uright; // Uu = 1 by definition
	dvector Yu, Yb, Yleft, Yright;

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

	// Utility functions for adaptation & regridding
	void rollVectorVector(const sdVector& y, vector<dvector>& v);
	void unrollVectorVector(const vector<dvector>& v);
	void unrollVectorVectorDot(const vector<dvector>& v);

	void updateTransportProperties(void);
	void updateThermoProperties(void);
	void updateLeftBC(void);

	void printForMatlab(ofstream& file, vector<double>& v, int index, char* name);
	void writeStateMatFile(void);
	void writeErrorFile(void);

	// these should be read-only:
	int N; // total problem size;
	int nVars; // Number of solution variables at each point
	int nSpec; // Number of chemical species

	// State variables:
	vector<double> rhov; // mass flux normal to flame per unit area (rho*v) 
	vector<double> U; // normalized tangential velocity (u/u_inf)
	vector<double> T; // temperature
	Array2D Y; // species mass fractions, Y(k,j)
	
	// Time derivatives of state variables:
	vector<double> drhovdt;
	vector<double> dUdt;
	vector<double> dTdt;
	Array2D dYdt;

	// Spatial derivatives of state variables:
	dvector dUdx;
	dvector dTdx;
	Array2D dYdx;

	// Auxillary variables:
	dvector rrhov;
	dvector rho; // density [kg/m^3]
	dvector drhodt;
	Array2D X; // species mole fractions, X(k,j)
	Array2D dXdx;
	dvector W; // species molecular weights
	dvector Wmx; // mixture molecular weight
	dvector sumcpj; // for enthalpy flux term
	
	dvector mu; // viscosity

	dvector lambda; // thermal conductivity
	dvector cp; // specific heat capacity (average)
	Array2D cpSpec; // species specific heat capacity
	dvector qFourier; // heat flux

	// Diffusion coefficients
	Array2D Dkm; // mixture-averaged diffusion coefficients
	Array2D Dkt; // thermal diffusion coefficients

	// Diffusion mass fluxes
	Array2D jFick; // Normal diffusion (Fick's Law)
	Array2D jSoret; // Soret Effect
	
	Array2D wDot; // species net production rates
	Array2D hk; // species enthalpies
	dvector qDot; // heat release rate per unit volume

	// the grid:
	oneDimGrid grid;

	// Strain rate parameters
	double strainRateInitial, strainRateFinal;
	double strainRateDt;
	double strainRateT0;
	double rStag; // nominal stagnation point radius

	 // Algebraic components of state, for IC calculation
	vector<bool> algebraic;
	void updateAlgebraicComponents(void);

	// Cantera data
	gasArray gas;

	// Miscellaneous options
	configOptions options;

private:
	// Subdivided governing equation components
	dvector energyUnst, energyDiff, energyConv, energyProd;
	dvector momentumUnst, momentumDiff, momentumConv, momentumProd;
	Array2D speciesUnst, speciesDiff, speciesConv, speciesProd;
	dvector continuityUnst, continuityRhov, continuityStrain;
	
	// Jacobian data
	int jacBW; // Bandwidth of the Jacobian (number of filled blocks per row, 
			   // dependent on the order of the finite difference stencil)
	int jacBWdot; // Bandwidth of dF/dydot component of Jacobian
	sdBandMatrix* bandedJacobian;
	vector<long int> pMat;
	bool inJacobianUpdate;
	bool inGetIC;

	int outputFileNumber; // number of output files written

	double strainRate(const double t);
	double dStrainRateDt(const double t);
};
