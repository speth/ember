#pragma once

#include "mathUtils.h"

using std::vector;

class oneDimGrid
{
public:
	dvector x; // The grid points
	dvector dampVal; // ratio of convective to diffusive coefficients (e.g. nu/v)
	
	// domain curvature parameter. 0: planar, 1: cylindrical
	int alpha;

	double vtol; // relative solution variable tolerance
	double dvtol; // derivative solution variable tolerance
	double absvtol; // absolute tolerance (ignore components with range smaller than this)
	double rmTol; // relative removal tolerance

	double uniformityTol;
	double gridMin, gridMax;
	double dampConst;

	// Derived mesh parameters (calculated by updateDerivedSizes)
	dvector cfm;
	dvector cf;
	dvector cfp;
	dvector dlj;
	dvector hh;
	dvector rphalf;

	double jZero; // index adjacent to the stagnation point

	// Update the grid based on the solutionState, and adjust it to fit 
	// on the new grid. Each element of the solutionState is a vector 
	// containing the value of a solution variable at each grid point.
	// Return value is true if the grid has been modified
	// adapt inserts and removes points within the problem domain
	// regrid inserts and removes points at the ends of the problem domain
	bool adapt(vector<dvector>& solutionState);
	bool regrid(vector<dvector>& solutionState);

private:	
	int jj; // number of grid points
	int nVars; // number of solution variables at each grid point

	void updateDerivedSizes(void);
	void removePoint(int jRemove, vector<dvector>& y);
	void addPoint(int jInsert, vector<dvector>& y);
	void update_jZero(void);

};