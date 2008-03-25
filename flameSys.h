#pragma once

#include <iostream>
#include "sundialsUtils.h"
#include "chemistry.h"
#include "grid.h"
#include "readConfig.h"

using Cantera::Array2D;
using std::string;

class flameSys : public sdDAE 
{
public:
	flameSys(void);
	~flameSys(void);
	
	// the functions for solving the ODE
	int f(realtype t, sdVector& y, sdVector& ydot, sdVector& res);

	//int JvProd(realtype t, sdVector& yIn, sdVector& ydotIn, sdVector& resIn, 
	//		   sdVector& vIn, sdVector& JvIn, realtype c_j);

	int preconditionerSetup(realtype t, sdVector& yIn, sdVector& ydotIn, 
		                    sdVector& resIn, realtype c_j);

	int preconditionerSolve(realtype t, sdVector& yIn, sdVector& ydotIn, sdVector& resIn,
						    sdVector& rhs, sdVector& outVec, realtype c_j, realtype delta);
	
	// Finds a consistent solution for the DAE to begin the integration
	int getInitialCondition(double t, sdVector& y, sdVector& ydot, vector<bool>& algebraic);

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
	
	void generateInitialProfiles(void);
	void loadInitialProfiles(void);
	void copyOptions(void);
	
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

	void printForMatlab(ofstream& file, dvector& v, int index, char* name);
	void writeStateMatFile(void);
	void writeErrorFile(void);

	// these should be read-only:
	int N; // total problem size;
	int nVars; // Number of solution variables at each point
	int nSpec; // Number of chemical species

	// State variables:
	dvector V; // mass flux normal to flame per unit area (rho*v) 
	dvector U; // normalized tangential velocity (u/u_inf)
	dvector T; // temperature
	Array2D Y; // species mass fractions, Y(k,j)
	
	// Time derivatives of state variables:
	dvector drhovdt;
	dvector dUdt;
	dvector dTdt;
	Array2D dYdt;

	// Spatial derivatives of state variables:
	dvector dUdx; // upwinded
	dvector dTdx; // upwinded
	Array2D dYdx; // upwinded
	dvector dTdxCen; // centered difference

	// Auxillary variables:
	dvector rrhov; // mass flux [kg/m^2*s or kg/m*rad*s]
	double rrhovC; // mass flux at centerline [kg/m^2 or kg/m*rad*s]
	dvector rho; // density [kg/m^3]
	dvector drhodt;
	Array2D X; // species mole fractions, X(k,j)
	Array2D dXdx;
	dvector W; // species molecular weights [kg/kmol]
	dvector Wmx; // mixture molecular weight [kg/kmol]
	dvector sumcpj; // for enthalpy flux term
	
	
	dvector mu; // viscosity [Pa*s]

	dvector lambda; // thermal conductivity [W/m*K]
	dvector cp; // specific heat capacity (average) [J/kg*K]
	Array2D cpSpec; // species specific heat capacity [J/mol*K]
	dvector qFourier; // heat flux [W/m^2]

	// Diffusion coefficients
	Array2D rhoD; // density-weighted, mixture-averaged diffusion coefficients [kg/m*s]
	Array2D Dkt; // thermal diffusion coefficients [m^2/s]

	// Diffusion mass fluxes
	Array2D jFick; // Normal diffusion (Fick's Law)
	Array2D jSoret; // Soret Effect diffusion
	
	Array2D wDot; // species net production rates [kmol/m^3*s]
	Array2D hk; // species enthalpies [J/kmol]
	dvector qDot; // heat release rate per unit volume [W/m^3]

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

	bool forceTransportUpdate;
	int transportUpdateCounter;

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
	bool inGetIC;
	sdVector* yTempJac;
	sdVector* ydotTempJac;
	sdVector* resTempJac;

	double& jacA(const int j, const int k1, const int k2);
	double& jacB(const int j, const int k1, const int k2);
	double& jacC(const int j, const int k1, const int k2);

	int outputFileNumber; // number of output files written

	double strainRate(const double t);
	double dStrainRateDt(const double t);

	int kMomentum, 	kContinuity, kEnergy, kSpecies;
};
