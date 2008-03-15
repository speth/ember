#include "grid.h"
#include "debugUtils.h"

using std::cout;
using std::endl;
using namespace mathUtils;

oneDimGrid::oneDimGrid(configOptions& theOptions)
	: options(theOptions)
{
}

void oneDimGrid::updateDerivedSizes()
{
	hh.resize(jj);
	cfm.resize(jj);
	cf.resize(jj);
	cfp.resize(jj);
	csm.resize(jj);
	cs.resize(jj);
	csp.resize(jj);
	dlj.resize(jj);
	rphalf.resize(jj);
	dampVal.resize(jj+1);

	for (unsigned int j=0; j<x.size()-1; j++) {
		hh[j] = x[j+1]-x[j];
		rphalf[j] =  pow(0.5*(x[j]+x[j+1]),alpha);
	}
	hh[jj-1] = hh[jj-2];

	for (unsigned int j=1; j<x.size()-1; j++) {
		cfp[j] = hh[j-1]/(hh[j]*(hh[j]+hh[j-1]));
		cf[j] = (hh[j]-hh[j-1])/(hh[j]*hh[j-1]);
		cfm[j] = -hh[j]/(hh[j-1]*(hh[j]+hh[j-1]));

		csp[j] = 2/(hh[j]*(hh[j]+hh[j-1]));
		cs[j] = -2/(hh[j]*hh[j-1]);
		csm[j] = 2/(hh[j-1]*(hh[j]+hh[j-1]));

		dlj[j]= 0.5 * (x[j+1]-x[j-1]);
	}
}

bool oneDimGrid::adapt(vector<dvector>& y, vector<dvector>& ydot)
{
	// This function takes the unadapted solution vector y, analyzes it,
	// and returns an updated solution vector. It tries to remove
	// unnescessary grid points located in the regions of small gradients,
	// insert new grid points in the regions of large gradients, and,
	// at the same time maintain relative uniformity of the grid.

	// Adaptation algorithm:
	//
	// The insertion of the grid points is performed first. For each
	// component of the solution vector:
	// 1. Find its range and the range of its derivative.
	// 2. Apply four criteria and and find where insertions are 
	//    needed. the criteria for a component f(j)
	//    a. |f[j+1]-f[j]| < vtol*range(f)
	//    b. |dfdy[j+1]-dfdy[j] | < dvtol*range(dfdy)
	//    c. 1/uniformityTol < hh[j]/hh[j-1] < uniformityTol
	// 3. If any of these criteria is not satisfied, a grid point j
	//    is inserted.
	// Next, the unnecessary grid points are removed, and the algorithm
	// is applied in reverse. If the criteria
	//    a. |f[j]-f[j-1]| > rmTol*vtol*range(f)
	//    b. |dfdy[j]-dfdy[j-1]| > rmTol*dvtol*range(dfdy)
	//    c. hh[j]+hh[j-1] < uniformityTol*hh[j-2]
	//    d. hh[j]+hh[j-1] < uniformityTol*hh[j+1]
	//  are satisfied for all components at a point, it is removed.

	nVars = y.size();
	jj = y[0].size()-1;

	bool gridUpdated = false; // flag if adaption has occured

	// *** Grid point insertion algorithm

	int j = 0;
	vector<double> dv(jj+1); // dv/dx

	while (j < jj) {
		updateDerivedSizes();
		dv.resize(jj+1);
		bool insert = false;

		// Consider tolerances for each variable v in the solution y
		for (int k=0; k<nVars; k++) {

			vector<double>& v = y[k];
			for (int i=1; i<jj; i++) {
				dv[i] = cfp[i]*v[i+1] + cf[i]*v[i] + cfm[i]*v[i-1];
			}

			double vRange = mathUtils::range(v);
			double dvRange = mathUtils::range(dv,1,jj-1);

			if (vRange < absvtol) {
				continue; // Ignore minor species
			}

			// Apply grid point addition criteria:

			// resolution of v
			if (abs(v[j+1]-v[j]) > vtol*vRange) {
				insert = true;
				if (debugParameters::debugAdapt) {
					cout << "Adapt: v resolution wants grid point j = " << j << ", k = " << k;
					cout << " |v(j+1)-v(j)|/vrange = " << abs(v[j+1]-v[j])/vRange << " > " << vtol << endl;
				}
			}

			// resolution of dv
			if (j!=0 && j!=jj-1 && abs(dv[j+1]-dv[j]) > dvtol*dvRange) {
				insert = true;
				if (debugParameters::debugAdapt) {
					cout << "Adapt: dv resolution wants grid point j = " << j << ", k = " << k;
					cout << " |dv(j+1)-dv(j)|/vrange = " << abs(dv[j+1]-dv[j])/dvRange << " > " << dvtol << endl;
				}
			}
		}
		
		// Left uniformity
		if (j!=0 && hh[j]/hh[j-1] > uniformityTol) {
			insert = true;
			if (debugParameters::debugAdapt) { 
				cout << "Adapt: left uniformity wants grid point j = " << j;
				cout << " hh(j)/hh(j-1) = " << hh[j]/hh[j-1] << ' > ' << uniformityTol << endl;
			}
		}

		// Right uniformity
		if (j!=jj-1 && hh[j]/hh[j+1] > uniformityTol) {
			insert = true;
			if (debugParameters::debugAdapt) {
				cout << "Adapt: right uniformity wants grid point j = " << j;
				cout << " hh(j)/hh(j+1) = " << hh[j]/hh[j+1] << " > " << uniformityTol << endl;
			}
		}

		// Maximum grid size
		if (hh[j] > gridMax) {
			insert = true;
			if (debugParameters::debugAdapt) {
				cout << "Adapt: Maximum grid size criterion wants a grid point j = " << j;
				cout << " hh(j) = " << hh[j] << " > " << gridMax << endl;
			}
		}

		// Damping of high-frequency numerical error
		if (hh[j] > dampConst*dampVal[j]) {
			insert = true;
			if (debugParameters::debugAdapt) {
				cout << "Adapt: damping criterion wants a grid point j = " << j;
				cout << " hh(j) = " << hh[j] << " > " << dampConst*dampVal[j] << endl;
			}
		}
		
		// Enforce minimum grid size
		if (hh[j] < 2*gridMin) {
			insert = false;
			if (debugParameters::debugAdapt) {
				cout << "Adapt: grid point addition cancelled by minimum grid size j = " << j;
				cout << " hh(j) = " << hh[j] << " < " << 2*gridMin << endl;
			}
		}

		if (insert) {
			// Insert a new point
			cout << "Adapt: inserting grid point j = " << j << endl;
			addPoint(j, y, ydot);
			gridUpdated = true;
			jj++;
			j+=2;
		} else {
			// No insertion; step to the next point.
			j++;
		}
	}

	// *** Grid point removal algorithm

	j = 1;
	while (j < jj) {
		updateDerivedSizes();

		// Assume removal, then look for a conditon which prevents removal
		bool remove = true;

		// Consider tolerances each variable v in the solution y
		for (int k=0; k<nVars; k++) {

			vector<double>& v = y[k];
			for (int i=1; i<jj; i++) {
				dv[i] = cfp[i]*v[i+1] + cf[i]*v[i] + cfm[i]*v[i-1];
			}

			double vRange = mathUtils::range(v);
			double dvRange = mathUtils::range(dv,1,jj-1);
			if (vRange < absvtol) {
				continue; // ignore minor componenents
			}

			// Apply grid point removal criteria:

			// resolution of v
			if (abs(v[j+1]-v[j-1]) > rmTol*vtol*vRange) {
				if (debugParameters::debugAdapt) {
					cout << "Adapt: no removal - v res. j = " << j << ", k = " << k;
					cout << " |v(j+1)-v(j-1)|/vtrange = " << abs(v[j+1]-v[j-1])/vRange << " > " << vtol*rmTol << endl;
				}
				remove = false;
			}

			// resolution of dv
			if (j!=2 && j!=jj-1 && abs(dv[j+1]-dv[j-1]) > rmTol*dvtol*dvRange) {
				if (debugParameters::debugAdapt) {
					cout << "Adapt: no removal - dv res. j = " << j << ", k = " << k;
					cout << " |dv(j+1)-dv(j-1)|/dvrange = " << abs(dv[j+1]-dv[j-1])/dvRange << " > " << dvtol*rmTol << endl;
				}
				remove = false;
			}
		}

		// Enforce maximum grid size
		if (hh[j]+hh[j-1] > gridMax) {
			if (debugParameters::debugAdapt) {
				cout << "Adapt: no removal - maximum grid size. j = " << j;
				cout << " hh(j)+hh(j-1) = " << hh[j]+hh[j-1] << " > " << gridMax << endl;
			}
			remove = false;
		}

		// Damping of high-frequency numerical error
		if (hh[j]+hh[j-1] >= rmTol*dampConst*dampVal[j]) {
			if (debugParameters::debugAdapt) {
				cout << "Adapt: no removal - damping criterion. j = " << j;
				cout << "hh(j) = " << hh[j] << " > " << dampConst*dampVal[j] << endl;
			}
			remove = false;
		}

		// Enforce left uniformity
		if (j>=2 && hh[j]+hh[j-1] > uniformityTol*hh[j-2]) {
			if (debugParameters::debugAdapt) {
				cout << "Adapt: no removal - left uniformity. j = " << j;
				cout << "(hh(j)+hh(j-1))/hh(j-2) = " << (hh[j]+hh[j-1])/hh[j-2] << " > " << uniformityTol << endl;
			}
			remove = false;
		}

		// Enforce right uniformity
		if (j<=jj-2 && hh[j]+hh[j-1] > uniformityTol*hh[j+1]) {
			if (debugParameters::debugAdapt) {
				cout << "Adapt: no removal - right uniformity. j = " << j;
				cout << "(hh(j)+hh(j-1))/hh(j+1) = " << (hh[j]+hh[j-1])/hh[j+1] << " > " << uniformityTol << endl;
			}
			remove = false;
		}

		if (j==jZero) {
			if (debugParameters::debugAdapt) {
				cout << "Adapt: no removal - not allowed to remove jZero" << endl;
			}
			remove = false;
		}


		if (remove) {
			cout << "Adapt: removing grid node j = " << j << endl;
			removePoint(j, y, ydot);
            jj--;
		    gridUpdated = true;
		} else {
			j++;
		}
	}

	if (gridUpdated) {
		updateDerivedSizes();
		updateBoundaryIndices();
		update_jZero(y[kContinuity]);
	}

	return gridUpdated;

}

void oneDimGrid::addPoint(int jInsert, vector<dvector>& y, vector<dvector>& ydot)
{
  	dvector::iterator iter;
	double xInsert = 0.5*(x[jInsert+1]+x[jInsert]);
	
	dampVal.insert(dampVal.begin()+jInsert+1, mathUtils::splines(x,dampVal, xInsert));

	vector<dvector>::iterator i;
	double yNew, ydotNew;

	for (i=y.begin(); i!=y.end(); i++) {
		yNew = mathUtils::splines(x,*i, xInsert);
		i->insert(i->begin()+jInsert+1, yNew);
	}

	for (i=ydot.begin(); i!=ydot.end(); i++) {
		ydotNew = mathUtils::splines(x,*i, xInsert);
		i->insert(i->begin()+jInsert+1, ydotNew);
	}

	x.insert(x.begin()+jInsert+1, xInsert);
}

void oneDimGrid::removePoint(int jRemove, vector<dvector>& y, vector<dvector>& ydot)
{
	x.erase(x.begin() + jRemove);
	dampVal.erase(dampVal.begin() + jRemove);

	vector<dvector>::iterator i;
	for (i=y.begin(); i!=y.end(); i++) {
		i->erase(i->begin() + jRemove);
	}

	for (i=ydot.begin(); i!=ydot.end(); i++) {
		i->erase(i->begin() + jRemove);
	}
}

bool oneDimGrid::regrid(vector<dvector>& y, vector<dvector>& ydot)
{

	bool gridUpdated = false;

	nVars = y.size();
	jj = y[0].size()-1;

	// *** Criteria for addition to right (j==jj) ***

	// Pick the comparison point for the grid flatness criterion,
	// depending on whether a zero-gradient condition is being enforced

	// djMom is for momentum equation, which always uses a
	// zero-gradient condition on the burned side:
	int djMom = (jb==jj) ? 2 : 1;

	// All other variables use fixed or zero gradient conditions
	// depending on fixedBurnedVal
	int djOther = (jb==jj && !fixedBurnedVal) ? 2 : 1;

	bool addRight = false; // Assume no addition

	// check flatness of temperature, velocity and species profiles at the boundary
	for (unsigned int k=0; k<y.size(); k++) {
		if (k == kContinuity) {
			continue;
		}

		int dj = (k == kMomentum) ? djMom : djOther;
		if (abs(y[k][jj]-y[k][jj-dj])/range(y[k]) > boundaryTol && range(y[k]) > absvtol ) {
			if (!addRight && debugParameters::debugRegrid) {
				cout << "Regrid: Flatness of component " << k << " requires right addition: ";
				cout << abs(y[k][jj]-y[k][jj-dj])/mathUtils::range(y[k]) << " > " << boundaryTol << endl;
			}
			addRight = true;
			break;
		}
	}

	if (!addRight && debugParameters::debugRegrid) {
		cout << "Regrid: no addition on right" << endl;
	}

	// *** Criteria for addition to the left (j=1) ***

	djOther = (jb==1 && !fixedBurnedVal) ? 2 : 1;
	djMom = (jb==1) ? 2 : 1;

	bool addLeft = false;
	for (unsigned int k=0; k<y.size(); k++) {
		if (k==kContinuity) {
			continue;
		}

		int dj = (k==kMomentum) ? djMom : djOther;
		if (abs(y[k][dj]-y[k][0])/range(y[k]) > boundaryTol && range(y[k]) > absvtol) {
			if (!addLeft && debugParameters::debugRegrid) {
				cout << "Regrid: Flatness of component " << k << " requires left addition: ";
				cout << abs(y[k][dj]-y[k][0])/range(y[k]) << " > " << boundaryTol << endl;
			}
			addLeft = true;
			break;
		}
	}

	if (options.fixedLeftLoc && x[0]!=0.0) {
		if (!addLeft && debugParameters::debugRegrid) {
			cout << "Regrid: Adding point to force left boundary toward x=0" << endl;
		}
		addLeft = 1;
	}

	if (abs(x[0]) <= gridMin) {
		if (addLeft && debugParameters::debugRegrid) {
			cout << "Regrid: Cannot add point to the left because x[0] == 0" << endl;
		}
		addLeft = false;
	}

	if (!addLeft && debugParameters::debugRegrid) {
		cout << "Regrid: No addition to left." << endl;
	}

	// *** Criteria for removal from the right (j==jj) ***
	
	// Comparison point for flatness criteria, depending on 
	// zero-gradient or fixed value boundary condition
	djMom = (jb==jj) ? 3 : 2;
	djOther = (jb==jj && !fixedBurnedVal) ? 3 : 2;

	bool removeRight = true; // assume removal
	for (unsigned int k=0; k<y.size(); k++) {
		if (k==kContinuity) {
			continue; // no flatness criterion for continuity equation
		}
		int dj = (k==kMomentum) ? djMom : djOther;
		if (abs(y[k][jj]-y[k][jj-dj])/range(y[k]) > boundaryTolRm && range(y[k]) > absvtol) {
			if (removeRight && debugParameters::debugRegrid) {
				cout << "Regrid: Right removal prevented by component " << k << " flatness: ";
				cout << abs(y[k][jj]-y[k][jj-dj])/range(y[k]) << " > " << boundaryTolRm << endl;
			}
			removeRight = false;
			break;
		}
	}

	if (removeRight && debugParameters::debugRegrid) {
		cout << "Regrid: Removal of rightmost point allowed." << endl;
	}

	// *** Criteria for removal from the left (j==1) ***
	djMom = (jb==1) ? 3 : 2;
	djOther = (jb==1 && !fixedBurnedVal) ? 3 : 2;
	
	// Don't remove points if the location of the left boundary is fixed
	bool removeLeft = true; // assume removal
	if (options.fixedLeftLoc) {
		if (removeLeft && debugParameters::debugRegrid) {
			cout << "Regrid: left removal prevented by fixed left boundary" << endl;
		}
		removeLeft = false;
	}

	for (unsigned int k=0; k<y.size(); k++) {
		int dj = (k==kMomentum) ? djMom : djOther;
		if (k==kContinuity) {
			continue;
		}

		if (abs(y[k][dj]-y[k][0])/range(y[k]) > boundaryTolRm && range(y[k]) > absvtol) {
			if (removeLeft && debugParameters::debugRegrid) {
				cout << "Regrid: Left removal prevented component " << k << " flatness requirement: ";
				cout << abs(y[k][dj]-y[k][0])/range(y[k]) << " > " << boundaryTolRm << endl;
			}
			removeLeft = false;
		}
	}

	if (removeLeft && debugParameters::debugRegrid) {
		cout << "Regrid: Removal of leftmost point allowed." << endl;
	}

	// *** Perform the indicated additions and removals
	
	if (addRight) {
		gridUpdated = true;
		cout << "Regrid: Adding point to right side." << endl;

		for (int i=0; i<addPointCount; i++) {
			x.push_back(x[jj] + (x[jj]-x[jj-1]));
			for (unsigned int k=0; k<y.size(); k++) {
				if (k==kContinuity) {
					// linear extrapolation for rhov
					y[k].push_back(y[k][jj] + (y[k][jj]-y[k][jj-1])*(x[jj+1]-x[jj])/(x[jj]-x[jj-1]));
					ydot[k].push_back(ydot[k][jj] + (ydot[k][jj]-ydot[k][jj-1])*(x[jj+1]-x[jj])/(x[jj]-x[jj-1]));
				} else {
					// keep constant boundary value
					y[k].push_back(y[k][jj]);
					ydot[k].push_back(ydot[k][jj]);
				}
			}
			jj++;
		}

	} else if (removeRight) {
		gridUpdated = true;
		cout << "Regrid: Removing point from right side." << endl;
		// Remove point from the right
		removePoint(jj,y,ydot);
		//for (unsigned int k=0; k<y.size(); k++) {
		//	if (k != kContinuity) {
		//		y[k][jj-1] = y[k][jj];
		//		ydot[k][jj-1] = ydot[k][jj];
		//	}
		//	y[k].erase(y[k].end()-1);
		//	ydot[k].erase(ydot[k].end()-1);

		//}
		//x.erase(x.end()-1);
		jj--;
	}

	if (addLeft) {
		gridUpdated = true;
		cout << "Regrid: Adding point to left side." << endl;

		// Add point to the left.
		for (int i=0; i<addPointCount; i++) {
			x.insert(x.begin(),x[0]-(x[1]-x[0]));
			for (unsigned int k=0; k<y.size(); k++) {
				if (k==kContinuity) {
					// linear extrapolation
					y[k].insert(y[k].begin(),y[k][0] + (y[k][1]-y[k][0])*(x[0]-x[1])/(x[2]-x[1]));
					ydot[k].insert(ydot[k].begin(),ydot[k][0] + (ydot[k][1]-ydot[k][0])*(x[1]-x[0])/(x[2]-x[1]));
				} else {
					// keep constant boundary value
					y[k].insert(y[k].begin(),y[k][0]);
					ydot[k].insert(ydot[k].begin(),ydot[k][0]);
				}
			}
			jj++;
		}

	} else if (removeLeft) {
		gridUpdated = true;
		cout << "Regrid: Removing point from left side." << endl;

		//for (unsigned int k=0; k<y.size(); k++) {
		//	if (k != kContinuity) {
		//		y[k][1] = y[k][0];
		//		ydot[k][1] = ydot[k][0];
		//	}
		//	y[k].erase(y[k].begin());
		//	ydot[k].erase(ydot[k].begin());
		//}
		//x.erase(x.begin());
		removePoint(0,y,ydot);
		jj--;
	}

	if (gridUpdated) {
		updateDerivedSizes();
		updateBoundaryIndices();
		update_jZero(y[kContinuity]);
	}

	return gridUpdated;

}

void oneDimGrid::update_jZero(dvector& rhov)
{
	jZero = mathUtils::minloc(mathUtils::abs(rhov));
	rhov[jZero] = 0;
}

void oneDimGrid::updateBoundaryIndices(void) {
	if (unburnedLeft) {
		ju = 0;
		jb = jj;
	} else {
		jb = 0;
		ju = jj;
	}
}