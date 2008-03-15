#pragma once

#include "mathUtils.h"
#include "readConfig.h"

using std::vector;

class oneDimGrid
{
public:
	oneDimGrid(configOptions& theOptions);

	dvector x; // The grid points
	dvector dampVal; // ratio of convective to diffusive coefficients (e.g. nu/v)
	
	int alpha;	// domain curvature parameter. 0: planar, 1: cylindrical
	int ju, jb; 	// indices of burned / unburned boundaries
	int kMomentum; // index into solutionState of the Momentum equation
	int kContinuity; // index into solutionState of the continuity equation
	int kEnergy; // index into solutionState of the Energy equation
	int kSpecies; // index into solutionState of the first species equation

	// Parameters for controlling internal grid points:
	double vtol; // relative solution variable tolerance for point insertion
	double dvtol; // derivative solution variable tolerance for point insertion
	double absvtol; // absolute tolerance (ignore components with range smaller than this)
	double rmTol; // relative grid point removal tolerance
	double uniformityTol; // maximum ratio of adjacent grid point separation distances
	double gridMin, gridMax; // minimum and maximum grid point separation
	double dampConst; // relative allowable numerical damping

	// Parameters for controlling exterior grid points:
	bool fixedBurnedVal;
	bool unburnedLeft;
	double boundaryTol;
	double boundaryTolRm;
	int addPointCount; // number of points to add when regridding

	// Derived mesh parameters (calculated by updateDerivedSizes)
	dvector cfm, cf, cfp; // first centered difference
	dvector csm, cs, csp; // second centered difference
	dvector dlj;
	dvector hh;
	dvector rphalf;
	int jj; // number of grid points

	int jZero; // index of the stagnation point

	// Update the grid based on the solutionState, and adjust it to fit 
	// on the new grid. Each element of the solutionState is a vector 
	// containing the value of a solution variable at each grid point.
	// Return value is true if the grid has been modified
	// adapt inserts and removes points within the problem domain
	// regrid inserts and removes points at the ends of the problem domain
	bool adapt(vector<dvector>& y, vector<dvector>& ydot);
	bool regrid(vector<dvector>& y, vector<dvector>& ydot);
	void updateDerivedSizes(void);
	void update_jZero(dvector& rhov);

	void updateBoundaryIndices(void);

private:

	configOptions& options;
	int nVars; // number of solution variables at each grid point
	
	void removePoint(int jRemove, vector<dvector>& y, vector<dvector>& ydot);
	void addPoint(int jInsert, vector<dvector>& y, vector<dvector>& ydot);
	
};