#include "grid.h"
#include "debugUtils.h"

#include <assert.h>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace mathUtils;
using std::max; using std::min;

oneDimGrid::oneDimGrid()
    : updated(true)
    , leftBC(BoundaryCondition::FixedValue)
    , rightBC(BoundaryCondition::FixedValue)
{
}

void oneDimGrid::setOptions(const configOptions& options)
{
    vtol_in = options.vtol;
    dvtol_in = options.dvtol;
    absvtol = options.absvtol;
    rmTol = options.rmTol;
    uniformityTol = options.uniformityTol;
    gridMin = options.gridMin;
    gridMax = options.gridMax;
    dampConst = options.dampConst;
    centerGridMin = options.centerGridMin;

    fixedBurnedVal = options.fixedBurnedVal;
    unburnedLeft = options.unburnedLeft;
    fixedLeftLoc = options.fixedLeftLoc;
    twinFlame = options.twinFlame;
    curvedFlame = options.curvedFlame;

    boundaryTol = options.boundaryTol;
    boundaryTolRm = options.boundaryTolRm;
    addPointCount = options.addPointCount;
    alpha = options.gridAlpha;
}

void oneDimGrid::updateValues()
{
    // Derived mesh sizes
    hh.resize(jj);
    cfm.resize(jj);
    cf.resize(jj);
    cfp.resize(jj);

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

        dlj[j]= 0.5 * (x[j+1]-x[j-1]);
    }
}

void oneDimGrid::setSize(const size_t N)
{
    nPoints = N;
    jj = N-1;
}

void oneDimGrid::adapt(vector<dvector>& y)
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
    setSize(y[0].size());

    vtol.resize(nVars);
    dvtol.resize(nVars);
    for (size_t k=0; k<nVars; k++) {
        vtol[k] = vtol_in;
        dvtol[k] = dvtol_in;
    }

    // Used for informational purposes only
    std::vector<int> insertionIndicies;
    std::vector<int> removalIndices;

    // *** Grid point insertion algorithm

    size_t j = 0;
    dvector dv(jj+1); // dv/dx

    while (j < jj) {
        updateValues();
        dv.resize(jj+1);
        bool insert = false;

        // Consider tolerances for each variable v in the solution y
        for (size_t k=0; k<nVars; k++) {

            dvector& v = y[k];
            for (size_t i=1; i<jj; i++) {
                dv[i] = cfp[i]*v[i+1] + cf[i]*v[i] + cfm[i]*v[i-1];
            }

            double vRange = mathUtils::range(v);
            double dvRange = mathUtils::range(dv,1,jj-1);

            if (vRange < absvtol) {
                continue; // Ignore minor species
            }

            // Apply grid point addition criteria:

            // resolution of v
            if (abs(v[j+1]-v[j]) > vtol[k]*vRange) {
                insert = true;
                if (debugParameters::debugAdapt) {
                    logFile.write(format("Adapt: v resolution wants grid point"
                            " j = %i, k = %i; |v(j+1)-v(j)|/vrange = %g > %g") %
                            j % k % (abs(v[j+1]-v[j])/vRange) % vtol[k]);
                }
            }

            // resolution of dv
            if (j!=0 && j!=jj-1 && abs(dv[j+1]-dv[j]) > dvtol[k]*dvRange) {
                insert = true;
                if (debugParameters::debugAdapt) {
                    logFile.write(format("Adapt: dv resolution (global) wants grid point"
                            " j = %i, k = %i; |dv(j+1)-dv(j)|/vrange = %g > % g") %
                            j % k % (abs(dv[j+1]-dv[j])/dvRange) % dvtol[k]);
                }
            }
        }

        // Damping of high-frequency numerical error
        if (hh[j] > dampConst*dampVal[j]) {
            insert = true;
            if (debugParameters::debugAdapt) {
                logFile.write(format("Adapt: damping criterion wants a grid point"
                        " j = %i; hh[j] = %g > %g") %
                        j % hh[j] % (dampConst*dampVal[j]));
            }
        }

        // Maximum grid size
        if (hh[j] > gridMax) {
            insert = true;
            if (debugParameters::debugAdapt) {
                logFile.write(format("Adapt: Maximum grid size criterion wants a grid point"
                        " j = %i; hh[j] = %g < %g") %
                        j % hh[j] % gridMax);
            }
        }

        // Left uniformity
        if (j!=0 && hh[j]/hh[j-1] > uniformityTol) {
            insert = true;
            if (debugParameters::debugAdapt) {
                logFile.write(format("Adapt: left uniformity wants grid point"
                        " j = %i; hh(j)/hh(j-1) = %g > %g") %
                        j % (hh[j]/hh[j-1]) % uniformityTol);
            }
        }

        // Right uniformity
        if (j!=jj-1 && hh[j]/hh[j+1] > uniformityTol) {
            insert = true;
            if (debugParameters::debugAdapt) {
                logFile.write(format("Adapt: right uniformity wants grid point"
                        " j = %i; hh(j)/hh(j+1) = %g > %g") %
                        j % (hh[j]/hh[j+1]) % uniformityTol);
            }
        }

        // Special minimum grid size for flames pinned at x=0
        if (j == 0 && (leftBC == BoundaryCondition::ControlVolume ||
                       leftBC == BoundaryCondition::WallFlux))
        {
            double xLeftMin = min(centerGridMin, 0.02*x[jj]);
            if (hh[j] < 2*xLeftMin) {
                insert = false;
                if (debugParameters::debugAdapt) {
                    logFile.write(format("Adapt: grid point addition canceled"
                            " by minimum center grid size j = %i; hh[j] = %g < %g") %
                            j % hh[j] % (2*xLeftMin));
                }
            }
        }

        // Enforce minimum grid size
        if (insert && hh[j] < 2*gridMin) {
            insert = false;
            if (debugParameters::debugAdapt) {
                logFile.write(format("Adapt: grid point addition canceled"
                        " by minimum grid size j = %i; hh[j] = %g < %g") %
                        j % hh[j] % (2*gridMin));
            }
        }

        if (insert) {
            // Insert a new point
            insertionIndicies.push_back(j);
            addPoint(j, y);
            updated = true;
            setSize(nPoints+1);
            j+=2;
        } else {
            // No insertion; step to the next point.
            j++;
        }
    }

    if (insertionIndicies.size() != 0) {
        logFile.write("Grid points inserted at j = ", false);
        for (unsigned int i=0; i<insertionIndicies.size()-1; i++) {
            logFile.write(format("%i, ") % insertionIndicies[i], false);
        }
        logFile.write(insertionIndicies[insertionIndicies.size()-1]);
    }

    // *** Grid point removal algorithm

    j = 1;
    while (j < jj) {
        updateValues();

        // Assume removal, then look for a condition which prevents removal
        bool remove = true;

        // Consider tolerances each variable v in the solution y
        for (size_t k=0; k<nVars; k++) {

            dvector& v = y[k];
            for (size_t i=1; i<jj; i++) {
                dv[i] = cfp[i]*v[i+1] + cf[i]*v[i] + cfm[i]*v[i-1];
            }

            double vRange = mathUtils::range(v);
            double dvRange = mathUtils::range(dv,1,jj-1);

            if (vRange < absvtol) {
                continue; // Ignore minor species
            }

            // Apply grid point removal criteria:

            // resolution of v
            if (abs(v[j+1]-v[j-1]) > rmTol*vtol[k]*vRange) {
                if (debugParameters::debugAdapt) {
                    logFile.write(format("Adapt: no removal - v res. j = %i, k = %i;"
                            " |v[j+1]-v[j-1]|/vtrange = %g > %g") %
                            j % k % (abs(v[j+1]-v[j-1])/vRange) % (vtol[k]*rmTol));
                }
                remove = false;
            }

            // resolution of dv
            if (j!=2 && j!=jj-1 && abs(dv[j+1]-dv[j-1]) > rmTol*dvtol[k]*dvRange) {
                if (debugParameters::debugAdapt) {
                    logFile.write(format("Adapt: no removal - dv res. j = %i, k = %i;"
                            " |dv(j+1)-dv(j-1)|/dvrange = %g > %g") %
                            j % k % (abs(dv[j+1]-dv[j-1])/dvRange) % (dvtol[k]*rmTol));
                }
                remove = false;
            }
        }

        // Damping of high-frequency numerical error
        if (hh[j]+hh[j-1] >= rmTol*dampConst*dampVal[j]) {
            if (debugParameters::debugAdapt) {
                logFile.write(format("Adapt: no removal - damping criterion. j = %i;"
                        " hh(j)+hh(j-1) = %g > %g") %
                        j % (hh[j]+hh[j-1]) % (dampConst*dampVal[j]));
            }
            remove = false;
        }

        // Enforce maximum grid size
        if (hh[j]+hh[j-1] > gridMax) {
            if (debugParameters::debugAdapt) {
                logFile.write(format("Adapt: no removal - maximum grid size. j = %i;"
                        " hh(j)+hh(j-1) = %g > %g") %
                        j % (hh[j]+hh[j-1]) % gridMax);
            }
            remove = false;
        }

        // Enforce left uniformity
        if (j>=2 && hh[j]+hh[j-1] > uniformityTol*hh[j-2]) {
            if (debugParameters::debugAdapt) {
                logFile.write(format("Adapt: no removal - left uniformity. j = %i;"
                        " (hh(j)+hh(j-1))/hh(j-2) = %g > %g") %
                        j % ((hh[j]+hh[j-1])/hh[j-2]) % uniformityTol);
            }
            remove = false;
        }

        // Enforce right uniformity
        if (j<=jj-2 && hh[j]+hh[j-1] > uniformityTol*hh[j+1]) {
            if (debugParameters::debugAdapt) {
                logFile.write(format("Adapt: no removal - right uniformity. j = %i;"
                        " (hh(j)+hh(j-1))/hh(j+1) = %g > %g") %
                        j % ((hh[j]+hh[j-1])/hh[j+1]) % uniformityTol);
            }
            remove = false;
        }

        // Special fixed grid for flames pinned at x=0
        if (j == 1 && (leftBC == BoundaryCondition::ControlVolume ||
                       leftBC == BoundaryCondition::WallFlux))
        {
            if (debugParameters::debugAdapt) {
                logFile.write("Adapt: no removal - fixed grid near r = 0.");
            }
            remove = false;
        }

        if (remove) {
            removalIndices.push_back(j);
            removePoint(j, y);
            setSize(nPoints-1);
            updated = true;
        } else {
            j++;
        }
    }

    if (removalIndices.size() != 0) {
        logFile.write("Grid points removed at j = ", false);
        for (unsigned int i=0; i<removalIndices.size()-1; i++) {
            logFile.write(format("%i, ") % removalIndices[i], false);
        }
        logFile.write(removalIndices[removalIndices.size()-1]);
    }

    if (updated) {
        updateValues();
        updateBoundaryIndices();
    }
}

void oneDimGrid::addPoint(int jInsert, vector<dvector>& y)
{
    assert(x.size() == dampVal.size());

    double xInsert = 0.5*(x[jInsert+1]+x[jInsert]);

    dampVal.insert(dampVal.begin()+jInsert+1, mathUtils::splines(x,dampVal, xInsert));

    foreach(dvector& row, y) {
        double yNew = mathUtils::splines(x, row, xInsert);
        row.insert(row.begin()+jInsert+1, yNew);
    }

    x.insert(x.begin()+jInsert+1, xInsert);
}

void oneDimGrid::removePoint(int jRemove, vector<dvector>& y)
{
    x.erase(x.begin() + jRemove);
    foreach(dvector& row, y) {
        row.erase(row.begin() + jRemove);
    }
}

bool oneDimGrid::addRight(vector<dvector>& y)
{
    // *** Criteria for addition to right (j==jj) ***

    // Pick the comparison point for the grid flatness criterion,
    // depending on whether a zero-gradient condition is being enforced

    // djMom is for momentum equation, which always uses a special zero-gradient condition
    int djMom = 1;

    // All other variables use fixed or zero gradient conditions
    // depending on fixedBurnedVal
    size_t djOther = (jb==jj && !fixedBurnedVal) ? 2 : 1;

    bool pointAdded = false; // Assume no addition

    // check flatness of temperature, velocity and species profiles at the boundary
    for (size_t k=0; k<nVars; k++) {
        if (!rightComponents[k]) {
            continue;
        }
        int dj = (k == kMomentum) ? djMom : djOther;
        double ymax = maxval(y[k]);
        if (abs(y[k][jj]-y[k][jj-dj])/ymax > boundaryTol && ymax > absvtol ) {
            if (!pointAdded && debugParameters::debugRegrid) {
                logFile.write(format("Regrid: Flatness of component %i requires"
                        " right addition: %g > %g") %
                        k % (abs(y[k][jj]-y[k][jj-dj])/ymax) % boundaryTol);
            }
            pointAdded = true;
            break;
        }
    }

    if (!pointAdded && debugParameters::debugRegrid) {
        logFile.write("Regrid: no addition on right");
    }

    if (pointAdded) {
        logFile.write("Regrid: Adding points to right side.");

        for (size_t i=0; i<addPointCount; i++) {
            x.push_back(x[jj] + (x[jj]-x[jj-1]));
            for (size_t k=0; k<nVars; k++) {
                // keep constant boundary value
                y[k].push_back(y[k][jj]);
            }
            setSize(nPoints+1);
        }
        updateBoundaryIndices();
    }

    return pointAdded;

}

bool oneDimGrid::addLeft(vector<dvector>& y)
{
    // *** Criteria for addition to the left (j==0) ***

    int djOther = (jb==1 && !fixedBurnedVal) ? 2 : 1;
    int djMom = 1;

    bool pointAdded = false;
    if (!fixedLeftLoc) {
        for (size_t k=0; k<nVars; k++) {
            if (!leftComponents[k]) {
                continue;
            }
            double ymax = maxval(y[k]);
            int dj = (k==kMomentum) ? djMom : djOther;
            if (abs(y[k][dj]-y[k][0])/ymax > boundaryTol && ymax > absvtol) {
                if (!pointAdded && debugParameters::debugRegrid) {
                    logFile.write(format("Regrid: Flatness of component %i"
                            " requires left addition: %g > %g") %
                            k % (abs(y[k][dj]-y[k][0])/ymax) % boundaryTol);
                }
                pointAdded = true;
                break;
            }
        }
    }

    if (fixedLeftLoc &&
        leftBC != BoundaryCondition::ControlVolume &&
        leftBC != BoundaryCondition::WallFlux)
    {
        if (!pointAdded && debugParameters::debugRegrid) {
            logFile.write("Regrid: Adding point to force left boundary toward x = 0");
        }
        pointAdded = true;
    }

    if (!pointAdded && debugParameters::debugRegrid) {
        logFile.write("Regrid: No addition to left.");
    }

    if (pointAdded) {
        logFile.write("Regrid: Adding points to left side.");
        // Add point to the left.
        for (size_t i=0; i<addPointCount; i++) {
            double xLeft = x[0]-hh[0];
            if (twinFlame || curvedFlame) {
                if (x[0] == 0) {
                    break;
                }

                double xLeftMin = std::min(centerGridMin, 0.005*x[jj]);
                if (xLeft < 0.0) {
                    xLeft = 0.0;
                } else if (xLeft < xLeftMin) {
                    xLeft = xLeftMin;
                }
            }
            if (debugParameters::debugRegrid) {
                logFile.write(format("adding point at %i") % xLeft);
            }
            x.insert(x.begin(),xLeft);

            for (size_t k=0; k<nVars; k++) {
                // keep constant boundary value
                y[k].insert(y[k].begin(),y[k][0]);
            }
            setSize(nPoints+1);
        }
        updateBoundaryIndices();
    }

    return pointAdded;
}

bool oneDimGrid::removeRight(vector<dvector>& y)
{
    // *** Criteria for removal from the right (j==jj) ***

    // Comparison point for flatness criteria, depending on
    // zero-gradient or fixed value boundary condition
    int djMom = 2;
    size_t djOther = (jb==jj && !fixedBurnedVal) ? 3 : 2;

    bool pointRemoved = true; // assume removal
    for (size_t k=0; k<nVars; k++) {
        if (!rightComponents[k]) {
            continue;
        }
        int dj = (k==kMomentum) ? djMom : djOther;
        double ymax = maxval(y[k]);
        if (abs(y[k][jj]-y[k][jj-dj])/ymax > boundaryTolRm && ymax > absvtol) {
            if (pointRemoved && debugParameters::debugRegrid) {
                logFile.write(format("Regrid: Right removal prevented by component %i"
                        " flatness: %g > %g") %
                        k % (abs(y[k][jj]-y[k][jj-dj])/ymax) % boundaryTolRm);
            }
            pointRemoved = false;
            break;
        }
    }

    if (jj < 3) {
        pointRemoved = false;
    }

    if (pointRemoved) {
        removePoint(jj,y);
        setSize(nPoints-1);
        updateBoundaryIndices();
    }

    return pointRemoved;
}

bool oneDimGrid::removeLeft(vector<dvector>& y)
{
    // *** Criteria for removal from the left (j==0) ***
    int djMom = 2;
    int djOther = (jb==1 && !fixedBurnedVal) ? 3 : 2;

    // Don't remove points if the location of the left boundary is fixed
    bool pointRemoved = true; // assume removal
    if (fixedLeftLoc) {
        if (pointRemoved && debugParameters::debugRegrid) {
            logFile.write("Regrid: left removal prevented by fixed left boundary");
        }
        pointRemoved = false;
    }

    for (size_t k=0; k<nVars; k++) {
        if (!leftComponents[k]) {
            continue;
        }
        size_t dj = (k==kMomentum) ? djMom : djOther;
        double ymax = maxval(y[k]);
        if (abs(y[k][dj]-y[k][0])/ymax > boundaryTolRm && ymax > absvtol) {
            if (pointRemoved && debugParameters::debugRegrid) {
                logFile.write(format("Regrid: Left removal prevented component %i"
                        " flatness requirement: %g > %g") %
                        k % (abs(y[k][dj]-y[k][0])/ymax) % boundaryTolRm);
            }
            pointRemoved = false;
        }
    }

    if (jj < 3) {
        pointRemoved = false;
    }

    if (pointRemoved) {
        removePoint(0,y);
        setSize(nPoints-1);
        updateBoundaryIndices();
    }

    return pointRemoved;
}


void oneDimGrid::regrid(vector<dvector>& y)
{
    nVars = y.size();

    setSize(y[0].size());

    // By default, all components count for regridding
    leftComponents.resize(nVars, true);
    rightComponents.resize(nVars, true);


    bool rightAddition = addRight(y);
    bool leftAddition = addLeft(y);

    bool rightRemoval = false;
    bool leftRemoval = false;

    int rightRemovalCount = 0;

    if (!rightAddition) {
        bool continueRemoval = true;
        while (continueRemoval) {
            continueRemoval = removeRight(y);
            if (continueRemoval) {
                rightRemoval = true;
                rightRemovalCount++;
            }
        }
    }
    if (rightRemovalCount > 0) {
        std::string suffix = (rightRemovalCount == 1) ? "" : "s";
        logFile.write(format("Removed %i point%s from the right side.") %
                rightRemovalCount % suffix);
    }

    int leftRemovalCount = 0;
    if (!leftAddition) {
        bool continueRemoval = true;
        while (continueRemoval) {
            continueRemoval = removeLeft(y);
            if (continueRemoval) {
                leftRemoval = true;
                leftRemovalCount++;
            }
        }
    }

    if (leftRemovalCount > 0) {
        std::string suffix = (leftRemovalCount == 1) ? "" : "s";
        logFile.write(format("Removed %i point%s from the left side.") %
                leftRemovalCount % suffix);
    }

    updated = (updated || leftAddition || rightAddition || leftRemoval || rightRemoval);

    if (updated) {
        updateValues();
        updateBoundaryIndices();
    }
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

GridBased::GridBased()
    : x(grid.x)
    , r(grid.r)
    , rphalf(grid.rphalf)
    , hh(grid.hh)
    , dlj(grid.dlj)
    , cfm(grid.cfm)
    , cf(grid.cf)
    , cfp(grid.cfp)
    , alpha(grid.alpha)
    , nPoints(grid.nPoints)
    , jj(grid.jj)
{
}

void GridBased::setGrid(const oneDimGrid& other)
{
    grid = other;
}
