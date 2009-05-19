#include "grid.h"
#include "debugUtils.h"

using std::cout;
using std::endl;
using namespace mathUtils;
using std::max; using std::min;

oneDimGrid::oneDimGrid(configOptions& theOptions)
    : options(theOptions)
{
}

void oneDimGrid::updateValues()
{
    // Derived mesh sizes
    hh.resize(jj);
    cfm.resize(jj);
    cf.resize(jj);
    cfp.resize(jj);
    csm.resize(jj);
    cs.resize(jj);
    csp.resize(jj);
    dlj.resize(jj);
    rphalf.resize(jj);
    r.resize(jj+1);
    dampVal.resize(jj+1);

    for (unsigned int j=0; j<x.size()-1; j++) {
        hh[j] = x[j+1]-x[j];
        rphalf[j] =  pow(0.5*(x[j]+x[j+1]),alpha);
        r[j] = pow(x[j],alpha);
    }
    r[jj] = pow(x[jj],alpha);

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
    // unnecessary grid points located in the regions of small gradients,
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
    //    b. |dfdy[j+1]-dfdy[j]| < dvtol*range(dfdy)
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

    // Used for informational purposes only
    std::vector<int> insertionIndicies;
    std::vector<int> removalIndices;

    bool gridUpdated = false; // flag if adaption has occurred

    // *** Grid point insertion algorithm

    int j = 0;
    dvector dv(jj+1); // dv/dx

    while (j < jj) {
        updateValues();
        dv.resize(jj+1);
        bool insert = false;

        // Consider tolerances for each variable v in the solution y
        for (int k=0; k<nVars; k++) {

            if (k==kContinuity) {
                continue;
            }
            dvector& v = y[k];
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
                    cout << "Adapt: dv resolution (global) wants grid point j = " << j << ", k = " << k;
                    cout << " |dv(j+1)-dv(j)|/vrange = " << abs(dv[j+1]-dv[j])/dvRange << " > " << dvtol << endl;
                }
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

        // Maximum grid size
        if (hh[j] > gridMax) {
            insert = true;
            if (debugParameters::debugAdapt) {
                cout << "Adapt: Maximum grid size criterion wants a grid point j = " << j;
                cout << " hh(j) = " << hh[j] << " > " << gridMax << endl;
            }
        }

        // Left uniformity
        if (j!=0 && hh[j]/hh[j-1] > uniformityTol) {
            insert = true;
            if (debugParameters::debugAdapt) {
                cout << "Adapt: left uniformity wants grid point j = " << j;
                cout << " hh(j)/hh(j-1) = " << hh[j]/hh[j-1] << " > " << uniformityTol << endl;
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

        // Special minimum grid size for flames pinned at x=0
        if (x[0]==0.0 && x[j] < centerGridMin) {
            insert = false;
            if (debugParameters::debugAdapt) {
                cout << "Adapt: grid point addition canceled by minimum center grid size j = " << j;
                cout << " hh(j) = " << hh[j] << " < " << 2*centerGridMin << endl;
            }
        }

        // Enforce minimum grid size
        if (insert && hh[j] < 2*gridMin) {
            insert = false;
            if (debugParameters::debugAdapt) {
                cout << "Adapt: grid point addition canceled by minimum grid size j = " << j;
                cout << " hh(j) = " << hh[j] << " < " << 2*gridMin << endl;
            }
        }

        if (insert) {
            // Insert a new point
            insertionIndicies.push_back(j);
            addPoint(j, y, ydot);
            gridUpdated = true;
            jj++;
            j+=2;
        } else {
            // No insertion; step to the next point.
            j++;
        }
    }

    if (insertionIndicies.size() != 0) {
        cout << "Grid points inserted at j = ";
        for (unsigned int i=0; i<insertionIndicies.size()-1; i++) {
            cout << insertionIndicies[i] << ", ";
        }
        cout << insertionIndicies[insertionIndicies.size()-1] << endl;
    }

    // *** Grid point removal algorithm

    j = 1;
    while (j < jj) {
        updateValues();

        // Assume removal, then look for a condition which prevents removal
        bool remove = true;

        // Consider tolerances each variable v in the solution y
        for (int k=0; k<nVars; k++) {

            if (k==kContinuity) {
                continue;
            }

            dvector& v = y[k];
            for (int i=1; i<jj; i++) {
                dv[i] = cfp[i]*v[i+1] + cf[i]*v[i] + cfm[i]*v[i-1];
            }

            double vRange = mathUtils::range(v);
            double dvRange = mathUtils::range(dv,1,jj-1);

            if (vRange < absvtol) {
                continue; // Ignore minor species
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

        // Damping of high-frequency numerical error
        if (hh[j]+hh[j-1] >= rmTol*dampConst*dampVal[j]) {
            if (debugParameters::debugAdapt) {
                cout << "Adapt: no removal - damping criterion. j = " << j;
                cout << " hh(j)+hh(j-1) = " << hh[j]+hh[j-1] << " > " << dampConst*dampVal[j] << endl;
            }
            remove = false;
        }

        // Enforce maximum grid size
        if (hh[j]+hh[j-1] > gridMax) {
            if (debugParameters::debugAdapt) {
                cout << "Adapt: no removal - maximum grid size. j = " << j;
                cout << " hh(j)+hh(j-1) = " << hh[j]+hh[j-1] << " > " << gridMax << endl;
            }
            remove = false;
        }

        // Enforce left uniformity
        if (j>=2 && hh[j]+hh[j-1] > uniformityTol*hh[j-2]) {
            if (debugParameters::debugAdapt) {
                cout << "Adapt: no removal - left uniformity. j = " << j;
                cout << " (hh(j)+hh(j-1))/hh(j-2) = " << (hh[j]+hh[j-1])/hh[j-2] << " > " << uniformityTol << endl;
            }
            remove = false;
        }

        // Enforce right uniformity
        if (j<=jj-2 && hh[j]+hh[j-1] > uniformityTol*hh[j+1]) {
            if (debugParameters::debugAdapt) {
                cout << "Adapt: no removal - right uniformity. j = " << j;
                cout << " (hh(j)+hh(j-1))/hh(j+1) = " << (hh[j]+hh[j-1])/hh[j+1] << " > " << uniformityTol << endl;
            }
            remove = false;
        }

        if (j==jZero) {
            if (debugParameters::debugAdapt) {
                cout << "Adapt: no removal - not allowed to remove jZero" << endl;
            }
            remove = false;
        }

        // Special fixd grid for flames pinned at x=0
        if (x[0]==0.0 && x[j] < centerGridMin) {
            if (debugParameters::debugAdapt) {
                cout << "Adapt: no removal - fixed grid near r = 0." << endl; 
            }
            remove = false;
        }

        if (remove) {
            removalIndices.push_back(j);
            removePoint(j, y, ydot);
            jj--;
            gridUpdated = true;
        } else {
            j++;
        }
    }

    if (removalIndices.size() != 0) {
        cout << "Grid points removed at j = ";
        for (unsigned int i=0; i<removalIndices.size()-1; i++) {
            cout << removalIndices[i] << ", ";
        }
        cout << removalIndices[removalIndices.size()-1] << endl;
    }

    if (gridUpdated) {
        updateValues();
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
    vector<dvector>::iterator i;
    for (i=y.begin(); i!=y.end(); i++) {
        i->erase(i->begin() + jRemove);
    }
    for (i=ydot.begin(); i!=ydot.end(); i++) {
        i->erase(i->begin() + jRemove);
    }
}

bool oneDimGrid::addRight(void)
{
    vector<dvector>& y = *yIn;
    vector<dvector>& ydot = *ydotIn;

// *** Criteria for addition to right (j==jj) ***

    // Pick the comparison point for the grid flatness criterion,
    // depending on whether a zero-gradient condition is being enforced

    // djMom is for momentum equation, which always uses a
    // zero-gradient condition on the burned side:
    int djMom = (jb==jj) ? 2 : 1;

    // All other variables use fixed or zero gradient conditions
    // depending on fixedBurnedVal
    int djOther = (jb==jj && !fixedBurnedVal) ? 2 : 1;

    bool pointAdded = false; // Assume no addition

    // check flatness of temperature, velocity and species profiles at the boundary
    for (int k=0; k<nVars; k++) {
        if (k == kContinuity || k == kQdot) {
            continue;
        }

        int dj = (k == kMomentum) ? djMom : djOther;
        double ymax = maxval(y[k]);
        if (abs(y[k][jj]-y[k][jj-dj])/ymax > boundaryTol && ymax > absvtol ) {
            if (!pointAdded && debugParameters::debugRegrid) {
                cout << "Regrid: Flatness of component " << k << " requires right addition: ";
                cout << abs(y[k][jj]-y[k][jj-dj])/ymax << " > " << boundaryTol << endl;
            }
            pointAdded = true;
            break;
        }
    }

    if (!pointAdded && debugParameters::debugRegrid) {
        cout << "Regrid: no addition on right" << endl;
    }

    if (pointAdded) {
        cout << "Regrid: Adding points to right side." << endl;

        for (int i=0; i<addPointCount; i++) {
            x.push_back(x[jj] + (x[jj]-x[jj-1]));
            for (int k=0; k<nVars; k++) {
                if (k==kContinuity) {
                    // linear extrapolation for V
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
        updateBoundaryIndices();
    }

    return pointAdded;

}

bool oneDimGrid::addLeft(void)
{
    vector<dvector>& y = *yIn;
    vector<dvector>& ydot = *ydotIn;

    // *** Criteria for addition to the left (j=1) ***

    int djOther = (jb==1 && !fixedBurnedVal) ? 2 : 1;
    int djMom = (jb==1) ? 2 : 1;

    bool pointAdded = false;
    if (!options.fixedLeftLoc) {
        for (int k=0; k<nVars; k++) {
            if (k==kContinuity || k==kQdot) {
                continue;
            }

            double ymax = maxval(y[k]);
            int dj = (k==kMomentum) ? djMom : djOther;
            if (abs(y[k][dj]-y[k][0])/ymax > boundaryTol && ymax > absvtol) {
                if (!pointAdded && debugParameters::debugRegrid) {
                    cout << "Regrid: Flatness of component " << k << " requires left addition: ";
                    cout << abs(y[k][dj]-y[k][0])/ymax << " > " << boundaryTol << endl;
                }
                pointAdded = true;
                break;
            }
        }
    }

    if (options.fixedLeftLoc && x[0]!=0.0) {
        if (!pointAdded && debugParameters::debugRegrid) {
            cout << "Regrid: Adding point to force left boundary toward x=0" << endl;
        }
        pointAdded = true;
    }

    if (!pointAdded && debugParameters::debugRegrid) {
        cout << "Regrid: No addition to left." << endl;
    }

    if (pointAdded) {
        cout << "Regrid: Adding points to left side." << endl;
        // Add point to the left.
        for (int i=0; i<addPointCount; i++) {
            double xLeft = x[0]-hh[0];
            if (options.twinFlame || options.curvedFlame) {
                xLeft = std::max(xLeft,0.0e0);
            }
            x.insert(x.begin(),xLeft);

            for (int k=0; k<nVars; k++) {
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
        updateBoundaryIndices();
    }

    return pointAdded;
}

bool oneDimGrid::removeRight(void)
{
    vector<dvector>& y = *yIn;
    vector<dvector>& ydot = *ydotIn;

    // *** Criteria for removal from the right (j==jj) ***

    // Comparison point for flatness criteria, depending on
    // zero-gradient or fixed value boundary condition
    int djMom = (jb==jj) ? 3 : 2;
    int djOther = (jb==jj && !fixedBurnedVal) ? 3 : 2;

    bool pointRemoved = true; // assume removal
    for (int k=0; k<nVars; k++) {
        if (k==kContinuity || k==kQdot) {
            continue; // no flatness criterion for continuity equation
        }
        int dj = (k==kMomentum) ? djMom : djOther;
        double ymax = maxval(y[k]);
        if (abs(y[k][jj]-y[k][jj-dj])/ymax > boundaryTolRm && ymax > absvtol) {
            if (pointRemoved && debugParameters::debugRegrid) {
                cout << "Regrid: Right removal prevented by component " << k << " flatness: ";
                cout << abs(y[k][jj]-y[k][jj-dj])/ymax << " > " << boundaryTolRm << endl;
            }
            pointRemoved = false;
            break;
        }
    }

    if (pointRemoved) {
        removePoint(jj,y,ydot);
        jj--;
        updateBoundaryIndices();
    }

    return pointRemoved;
}

bool oneDimGrid::removeLeft(void)
{
    vector<dvector>& y = *yIn;
    vector<dvector>& ydot = *ydotIn;

    // *** Criteria for removal from the left (j==1) ***
    int djMom = (jb==1) ? 3 : 2;
    int djOther = (jb==1 && !fixedBurnedVal) ? 3 : 2;

    // Don't remove points if the location of the left boundary is fixed
    bool pointRemoved = true; // assume removal
    if (options.fixedLeftLoc) {
        if (pointRemoved && debugParameters::debugRegrid) {
            cout << "Regrid: left removal prevented by fixed left boundary" << endl;
        }
        pointRemoved = false;
    }

    for (int k=0; k<nVars; k++) {
        int dj = (k==kMomentum) ? djMom : djOther;
        if (k==kContinuity || k==kQdot) {
            continue;
        }

        double ymax = maxval(y[k]);
        if (abs(y[k][dj]-y[k][0])/ymax > boundaryTolRm && ymax > absvtol) {
            if (pointRemoved && debugParameters::debugRegrid) {
                cout << "Regrid: Left removal prevented component " << k << " flatness requirement: ";
                cout << abs(y[k][dj]-y[k][0])/ymax << " > " << boundaryTolRm << endl;
            }
            pointRemoved = false;
        }
    }

    if (pointRemoved) {
        removePoint(0,y,ydot);
        jj--;
        updateBoundaryIndices();
    }

    return pointRemoved;
}


bool oneDimGrid::regrid(vector<dvector>& y, vector<dvector>& ydot)
{
    nVars = y.size();
    kQdot = nVars-1;

    jj = y[0].size()-1;

    yIn = &y;
    ydotIn = &ydot;

    bool rightAddition = addRight();
    bool leftAddition = addLeft();

    bool rightRemoval = false;
    bool leftRemoval = false;

    int rightRemovalCount = 0;

    if (!rightAddition) {
        bool continueRemoval = true;
        while (continueRemoval) {
            continueRemoval = removeRight();
            if (continueRemoval) {
                rightRemoval = true;
                rightRemovalCount++;
            }
        }
    }
    if (rightRemovalCount > 0) {
        std::string suffix = (rightRemovalCount == 1) ? "" : "s";
        cout << "Removed " << rightRemovalCount << " point" << suffix << " from the right side. ";
    }

    int leftRemovalCount = 0;
    if (!leftAddition) {
        bool continueRemoval = true;
        while (continueRemoval) {
            continueRemoval = removeLeft();
            if (continueRemoval) {
                leftRemoval = true;
                leftRemovalCount++;
            }
        }
    }

    if (leftRemovalCount > 0) {
        std::string suffix = (leftRemovalCount == 1) ? "" : "s";
        cout << "Removed " << leftRemovalCount << " point" << suffix << " from the left side.";
    }
    if (leftRemovalCount > 0 || rightRemovalCount > 0) {
        cout << endl;
    }

    bool gridUpdated = (leftAddition || rightAddition || leftRemoval || rightRemoval);

    if (gridUpdated) {
        updateValues();
        updateBoundaryIndices();
        update_jZero(y[kContinuity]);
    }

    return gridUpdated;
}

void oneDimGrid::update_jZero(dvector& V)
{
    jZero = mathUtils::minloc(mathUtils::abs(V));
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
